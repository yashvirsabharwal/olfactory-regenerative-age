"""Leakage-safe nested tuning for ORA age models."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from itertools import product
import json
from typing import Any

import numpy as np
import pandas as pd

from .age_model import (
    _boolean_series,
    _combine_backend_info,
    _performance_row,
    add_oraa,
    biological_feature_columns,
    donor_cv_folds,
    fit_model_predictions,
    fit_preprocessor,
    model_names_from_config,
    summarize_repeated_performance,
    transform_preprocessor,
)


DEFAULT_SEARCH_SPACES: dict[str, list[dict[str, object]]] = {
    "xgboost": [
        {"max_depth": 1, "min_child_weight": 5.0, "reg_lambda": 10.0, "subsample": 0.8, "colsample_bytree": 0.8},
        {"max_depth": 2, "min_child_weight": 3.0, "reg_lambda": 5.0, "subsample": 0.8, "colsample_bytree": 0.8},
        {"max_depth": 2, "min_child_weight": 5.0, "reg_lambda": 10.0, "subsample": 0.7, "colsample_bytree": 0.8},
        {"max_depth": 3, "min_child_weight": 5.0, "reg_lambda": 10.0, "subsample": 0.8, "colsample_bytree": 0.7},
        {"n_estimators": 200, "learning_rate": 0.05, "max_depth": 2, "min_child_weight": 5.0, "reg_lambda": 10.0},
        {"n_estimators": 500, "learning_rate": 0.02, "max_depth": 2, "min_child_weight": 3.0, "reg_lambda": 5.0},
    ],
    "catboost": [
        {"depth": 2, "l2_leaf_reg": 5.0, "learning_rate": 0.03, "iterations": 300},
        {"depth": 2, "l2_leaf_reg": 10.0, "learning_rate": 0.03, "iterations": 300},
        {"depth": 3, "l2_leaf_reg": 10.0, "learning_rate": 0.03, "iterations": 300},
        {"depth": 3, "l2_leaf_reg": 20.0, "learning_rate": 0.03, "iterations": 300},
        {"depth": 3, "l2_leaf_reg": 10.0, "learning_rate": 0.02, "iterations": 500},
        {"depth": 4, "l2_leaf_reg": 20.0, "learning_rate": 0.02, "iterations": 500},
    ],
    "random_forest": [
        {"n_estimators": 500, "max_depth": 4},
        {"n_estimators": 500, "max_depth": 5},
        {"n_estimators": 800, "max_depth": 5},
        {"n_estimators": 500, "max_depth": 6},
        {"n_estimators": 800, "max_depth": 6},
    ],
}


@dataclass
class NestedTuningResult:
    performance: pd.DataFrame
    performance_summary: pd.DataFrame
    predictions: pd.DataFrame
    tuning_trace: pd.DataFrame
    selected_params: pd.DataFrame


def run_nested_tuning(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    model_config: dict[str, Any] | None = None,
    tuning_config: dict[str, Any] | None = None,
    *,
    repeats: int = 1,
    inner_folds: int = 3,
    max_candidates: int | None = None,
) -> NestedTuningResult:
    """Tune model parameters inside each outer donor fold and evaluate held-out donors."""

    model_config = dict(model_config or {})
    tuning_config = tuning_config or {}
    repeats = max(1, int(repeats))
    inner_folds = max(2, int(inner_folds))
    train, feature_cols = _training_frame(features, manifest, model_config)
    y = train["age"].astype(float).to_numpy()
    model_names = model_names_from_config(model_config)
    base_seed = int(model_config.get("random_seed", 42))
    performance_rows = []
    prediction_tables = []
    trace_rows = []
    selected_rows = []
    for repeat in range(repeats):
        repeat_config = dict(model_config)
        repeat_config["random_seed"] = base_seed + repeat
        outer_folds = donor_cv_folds(train, repeat_config)
        for model_name in model_names:
            candidates = candidate_params_for_model(model_name, model_config, tuning_config)
            if max_candidates is not None:
                candidates = candidates[: max(1, int(max_candidates))]
            pred = np.full(train.shape[0], np.nan, dtype=float)
            fold_backends = []
            for outer_fold, (outer_train_idx, outer_test_idx) in enumerate(outer_folds):
                selected, candidate_trace = tune_outer_fold(
                    train=train,
                    feature_cols=feature_cols,
                    y=y,
                    outer_train_idx=outer_train_idx,
                    model_name=model_name,
                    model_config=repeat_config,
                    candidates=candidates,
                    inner_folds=inner_folds,
                    random_seed=base_seed + repeat * 1000 + outer_fold,
                )
                for row in candidate_trace:
                    row.update({"repeat": repeat, "model": model_name, "outer_fold": outer_fold})
                    trace_rows.append(row)
                selected_rows.append(
                    {
                        "repeat": repeat,
                        "model": model_name,
                        "outer_fold": outer_fold,
                        "candidate_id": selected["candidate_id"],
                        "inner_mae_mean": selected["inner_mae_mean"],
                        "inner_rmse_mean": selected["inner_rmse_mean"],
                        "inner_spearman_r_mean": selected["inner_spearman_r_mean"],
                        "params_json": selected["params_json"],
                        **{f"param__{key}": value for key, value in selected["params"].items()},
                    }
                )
                selected_config = _config_with_params(repeat_config, model_name, selected["params"])
                prep = fit_preprocessor(train.iloc[outer_train_idx][feature_cols])
                x_train = transform_preprocessor(train.iloc[outer_train_idx][feature_cols], prep)
                x_test = transform_preprocessor(train.iloc[outer_test_idx][feature_cols], prep)
                fold_pred, _, backend_info = fit_model_predictions(
                    model_name,
                    x_train,
                    y[outer_train_idx],
                    x_test,
                    selected_config,
                )
                pred[outer_test_idx] = fold_pred
                fold_backends.append(backend_info)

            row = _performance_row(model_name, y, pred, _combine_backend_info(model_name, fold_backends))
            row["repeat"] = repeat
            performance_rows.append(row)
            prediction = pd.DataFrame(
                {
                    "donor_id": train["donor_id"].to_numpy(),
                    "model": model_name,
                    "chronological_age": y,
                    "ora": pred,
                }
            )
            prediction = add_oraa(prediction, train)
            prediction.insert(0, "repeat", repeat)
            prediction_tables.append(prediction)

    performance = pd.DataFrame(performance_rows)
    predictions = pd.concat(prediction_tables, ignore_index=True)
    tuning_trace = pd.DataFrame(trace_rows)
    selected_params = pd.DataFrame(selected_rows)
    return NestedTuningResult(
        performance=performance,
        performance_summary=summarize_repeated_performance(performance),
        predictions=predictions,
        tuning_trace=tuning_trace,
        selected_params=selected_params,
    )


def candidate_params_for_model(
    model_name: str,
    model_config: dict[str, Any] | None = None,
    tuning_config: dict[str, Any] | None = None,
) -> list[dict[str, object]]:
    """Return de-duplicated candidate parameter dictionaries for one model."""

    model_config = model_config or {}
    tuning_config = tuning_config or {}
    candidates: list[dict[str, object]] = []
    base_params = model_config.get("models", {}).get(model_name, {})
    if isinstance(base_params, dict):
        candidates.append({key: value for key, value in base_params.items() if key != "enabled"})
    spaces = tuning_config.get("search_spaces", {})
    if model_name in spaces:
        candidates.extend(_expand_candidates(spaces[model_name]))
    else:
        candidates.extend(DEFAULT_SEARCH_SPACES.get(model_name, [{}]))
    return _deduplicate_candidates(candidates)


def tune_outer_fold(
    *,
    train: pd.DataFrame,
    feature_cols: list[str],
    y: np.ndarray,
    outer_train_idx: np.ndarray,
    model_name: str,
    model_config: dict[str, Any],
    candidates: list[dict[str, object]],
    inner_folds: int,
    random_seed: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Evaluate candidate params on inner donor folds for one outer training split."""

    inner_data = train.iloc[outer_train_idx].reset_index(drop=True)
    inner_y = y[outer_train_idx]
    inner_config = dict(model_config)
    inner_config["outer_cv_folds"] = inner_folds
    inner_config["random_seed"] = random_seed
    folds = donor_cv_folds(inner_data, inner_config)
    rows = []
    for candidate_id, params in enumerate(candidates):
        candidate_config = _config_with_params(model_config, model_name, params)
        metrics = []
        for inner_train_idx, inner_valid_idx in folds:
            prep = fit_preprocessor(inner_data.iloc[inner_train_idx][feature_cols])
            x_train = transform_preprocessor(inner_data.iloc[inner_train_idx][feature_cols], prep)
            x_valid = transform_preprocessor(inner_data.iloc[inner_valid_idx][feature_cols], prep)
            pred, _, _ = fit_model_predictions(
                model_name,
                x_train,
                inner_y[inner_train_idx],
                x_valid,
                candidate_config,
            )
            metrics.append(_metric_row(inner_y[inner_valid_idx], pred))
        metrics_frame = pd.DataFrame(metrics)
        row = {
            "candidate_id": candidate_id,
            "params": params,
            "params_json": _params_json(params),
            "inner_folds": len(metrics),
            "inner_mae_mean": float(metrics_frame["mae"].mean()),
            "inner_mae_sd": float(metrics_frame["mae"].std(ddof=1)) if metrics_frame.shape[0] > 1 else 0.0,
            "inner_rmse_mean": float(metrics_frame["rmse"].mean()),
            "inner_spearman_r_mean": float(metrics_frame["spearman_r"].mean()),
        }
        rows.append(row)
    selected = sorted(rows, key=lambda row: (row["inner_mae_mean"], row["inner_rmse_mean"], row["candidate_id"]))[0]
    return selected, rows


def _training_frame(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    model_config: dict[str, Any],
) -> tuple[pd.DataFrame, list[str]]:
    donor_meta = (
        manifest.sort_values(["donor_id", "sample_id"] if "sample_id" in manifest else ["donor_id"])
        .drop_duplicates("donor_id")
        .set_index("donor_id")
    )
    data = features.set_index("donor_id").join(donor_meta, how="inner", rsuffix="_meta").reset_index()
    train_mask = _boolean_series(data["usable_for_ora_training"]) & data["age"].notna()
    train = data[train_mask].copy().reset_index(drop=True)
    feature_cols = biological_feature_columns(train, model_config)
    if not feature_cols:
        raise ValueError("No biological numeric feature columns available for nested tuning.")
    max_missing = float(model_config.get("missingness_max_fraction", 0.30))
    missing_fraction = train[feature_cols].isna().mean()
    feature_cols = [col for col in feature_cols if missing_fraction[col] <= max_missing]
    if not feature_cols:
        raise ValueError("All feature columns exceeded the missingness threshold.")
    return train, feature_cols


def _config_with_params(model_config: dict[str, Any], model_name: str, params: dict[str, object]) -> dict[str, Any]:
    config = deepcopy(model_config)
    config.setdefault("models", {})
    config["models"].setdefault(model_name, {})
    config["models"][model_name].update(params)
    return config


def _expand_candidates(spec: object) -> list[dict[str, object]]:
    if isinstance(spec, list):
        return [dict(item) for item in spec if isinstance(item, dict)]
    if isinstance(spec, dict):
        keys = list(spec.keys())
        values = [value if isinstance(value, list) else [value] for value in spec.values()]
        return [dict(zip(keys, combo, strict=True)) for combo in product(*values)]
    raise ValueError("Tuning search space must be a list of parameter maps or a parameter grid map.")


def _deduplicate_candidates(candidates: list[dict[str, object]]) -> list[dict[str, object]]:
    seen = set()
    output = []
    for params in candidates:
        key = _params_json(params)
        if key in seen:
            continue
        seen.add(key)
        output.append(dict(params))
    return output


def _metric_row(y_true: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    err = pred - y_true
    return {
        "mae": float(np.mean(np.abs(err))),
        "rmse": float(np.sqrt(np.mean(err**2))),
        "spearman_r": _spearman(y_true, pred),
    }


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2:
        return np.nan
    rx = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    ry = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    if float(np.std(rx)) == 0.0 or float(np.std(ry)) == 0.0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def _params_json(params: dict[str, object]) -> str:
    return json.dumps(params, sort_keys=True, separators=(",", ":"))
