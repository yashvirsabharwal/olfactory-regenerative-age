"""Leakage-safe stacked ORA age models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .age_model import (
    _performance_row,
    add_oraa,
    donor_cv_folds,
    fit_model_predictions,
    fit_preprocessor,
    summarize_repeated_performance,
    transform_preprocessor,
)
from .tuning import _training_frame


DEFAULT_BASE_MODELS = ["ridge", "random_forest", "xgboost", "catboost"]


@dataclass
class StackingResult:
    performance: pd.DataFrame
    performance_summary: pd.DataFrame
    predictions: pd.DataFrame
    meta_weights: pd.DataFrame


def run_stacked_ora(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    model_config: dict[str, Any] | None = None,
    *,
    base_models: list[str] | None = None,
    repeats: int = 5,
    inner_folds: int = 3,
    meta_alphas: list[float] | None = None,
) -> StackingResult:
    """Run donor-level nested stacking with inner OOF base predictions."""

    model_config = dict(model_config or {})
    if base_models is None:
        base_models = DEFAULT_BASE_MODELS
    base_models = [str(model) for model in base_models]
    if not base_models:
        raise ValueError("At least one base model is required for stacking.")
    repeats = max(1, int(repeats))
    inner_folds = max(2, int(inner_folds))
    train, feature_cols = _training_frame(features, manifest, model_config)
    y = train["age"].astype(float).to_numpy()
    base_seed = int(model_config.get("random_seed", 42))
    performance_rows = []
    prediction_tables = []
    weight_rows = []

    for repeat in range(repeats):
        repeat_config = dict(model_config)
        repeat_config["random_seed"] = base_seed + repeat
        outer_folds = donor_cv_folds(train, repeat_config)
        pred = np.full(train.shape[0], np.nan, dtype=float)
        for outer_fold, (outer_train_idx, outer_test_idx) in enumerate(outer_folds):
            meta_train = _inner_oof_base_predictions(
                train=train,
                feature_cols=feature_cols,
                y=y,
                outer_train_idx=outer_train_idx,
                base_models=base_models,
                model_config=repeat_config,
                inner_folds=inner_folds,
                random_seed=base_seed + repeat * 1000 + outer_fold,
            )
            meta_model = _fit_meta_model(meta_train, y[outer_train_idx], meta_alphas=meta_alphas)
            meta_test = _outer_base_predictions(
                train=train,
                feature_cols=feature_cols,
                y=y,
                outer_train_idx=outer_train_idx,
                outer_test_idx=outer_test_idx,
                base_models=base_models,
                model_config=repeat_config,
            )
            pred[outer_test_idx] = _predict_meta(meta_model, meta_test)
            for model, weight in zip(base_models, meta_model["coef"], strict=True):
                weight_rows.append(
                    {
                        "repeat": repeat,
                        "outer_fold": outer_fold,
                        "base_model": model,
                        "weight": float(weight),
                        "meta_alpha": meta_model["alpha"],
                    }
                )
            weight_rows.append(
                {
                    "repeat": repeat,
                    "outer_fold": outer_fold,
                    "base_model": "intercept",
                    "weight": float(meta_model["intercept"]),
                    "meta_alpha": meta_model["alpha"],
                }
            )

        row = _performance_row("stacked_ensemble", y, pred)
        row["repeat"] = repeat
        performance_rows.append(row)
        prediction = pd.DataFrame(
            {
                "donor_id": train["donor_id"].to_numpy(),
                "model": "stacked_ensemble",
                "chronological_age": y,
                "ora": pred,
            }
        )
        prediction = add_oraa(prediction, train)
        prediction.insert(0, "repeat", repeat)
        prediction_tables.append(prediction)

    performance = pd.DataFrame(performance_rows)
    predictions = pd.concat(prediction_tables, ignore_index=True)
    meta_weights = pd.DataFrame(weight_rows)
    return StackingResult(
        performance=performance,
        performance_summary=summarize_repeated_performance(performance),
        predictions=predictions,
        meta_weights=meta_weights,
    )


def _inner_oof_base_predictions(
    *,
    train: pd.DataFrame,
    feature_cols: list[str],
    y: np.ndarray,
    outer_train_idx: np.ndarray,
    base_models: list[str],
    model_config: dict[str, Any],
    inner_folds: int,
    random_seed: int,
) -> np.ndarray:
    inner_data = train.iloc[outer_train_idx].reset_index(drop=True)
    inner_y = y[outer_train_idx]
    inner_config = dict(model_config)
    inner_config["outer_cv_folds"] = inner_folds
    inner_config["random_seed"] = random_seed
    folds = donor_cv_folds(inner_data, inner_config)
    output = np.full((inner_data.shape[0], len(base_models)), np.nan, dtype=float)
    for model_idx, model_name in enumerate(base_models):
        for inner_train_idx, inner_valid_idx in folds:
            prep = fit_preprocessor(inner_data.iloc[inner_train_idx][feature_cols])
            x_train = transform_preprocessor(inner_data.iloc[inner_train_idx][feature_cols], prep)
            x_valid = transform_preprocessor(inner_data.iloc[inner_valid_idx][feature_cols], prep)
            fold_pred, _ = fit_model_predictions(model_name, x_train, inner_y[inner_train_idx], x_valid, inner_config)
            output[inner_valid_idx, model_idx] = fold_pred
    if np.isnan(output).any():
        raise ValueError("Inner OOF base predictions contain missing values.")
    return output


def _outer_base_predictions(
    *,
    train: pd.DataFrame,
    feature_cols: list[str],
    y: np.ndarray,
    outer_train_idx: np.ndarray,
    outer_test_idx: np.ndarray,
    base_models: list[str],
    model_config: dict[str, Any],
) -> np.ndarray:
    output = np.full((outer_test_idx.size, len(base_models)), np.nan, dtype=float)
    for model_idx, model_name in enumerate(base_models):
        prep = fit_preprocessor(train.iloc[outer_train_idx][feature_cols])
        x_train = transform_preprocessor(train.iloc[outer_train_idx][feature_cols], prep)
        x_test = transform_preprocessor(train.iloc[outer_test_idx][feature_cols], prep)
        pred, _ = fit_model_predictions(model_name, x_train, y[outer_train_idx], x_test, model_config)
        output[:, model_idx] = pred
    return output


def _fit_meta_model(x: np.ndarray, y: np.ndarray, *, meta_alphas: list[float] | None = None) -> dict[str, object]:
    alphas = meta_alphas or [0.1, 1.0, 10.0, 100.0]
    try:
        from sklearn.linear_model import RidgeCV  # type: ignore

        model = RidgeCV(alphas=alphas)
        model.fit(x, y)
        return {
            "coef": np.asarray(model.coef_, dtype=float),
            "intercept": float(model.intercept_),
            "alpha": float(model.alpha_),
        }
    except ModuleNotFoundError:
        design = np.column_stack([np.ones(x.shape[0]), x])
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)
        return {
            "coef": np.asarray(coef[1:], dtype=float),
            "intercept": float(coef[0]),
            "alpha": np.nan,
        }


def _predict_meta(meta_model: dict[str, object], x: np.ndarray) -> np.ndarray:
    coef = np.asarray(meta_model["coef"], dtype=float)
    intercept = float(meta_model["intercept"])
    return intercept + x @ coef
