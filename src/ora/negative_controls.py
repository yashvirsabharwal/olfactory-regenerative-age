"""Negative-control and technical-confound analyses for ORA."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .age_model import biological_feature_columns, donor_cv_folds, fit_preprocessor, transform_preprocessor


DEFAULT_TECHNICAL_COVARIATES = ("log10_total_cells", "sex", "chemistry", "collection_method", "site")


@dataclass
class NegativeControlResult:
    """Container for M1.3 negative-control outputs."""

    performance: pd.DataFrame
    baseline_comparison: pd.DataFrame
    covariate_explainability: pd.DataFrame


def run_negative_controls(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    scores: pd.DataFrame,
    model_config: dict[str, Any] | None = None,
    *,
    n_shuffles: int = 50,
    ridge_alpha: float = 10.0,
) -> NegativeControlResult:
    """Run donor-level negative controls and technical baseline checks."""

    model_config = dict(model_config or {})
    data = _training_frame(features, manifest)
    biological = _biological_matrix(data, model_config)
    technical = _technical_matrix(data, DEFAULT_TECHNICAL_COVARIATES)
    yield_only = _technical_matrix(data, ("log10_total_cells",))
    y = pd.to_numeric(data["age"], errors="coerce").to_numpy(dtype=float)

    performance_rows = [
        _cv_metric_row(
            data,
            biological,
            y,
            model_config,
            control="biological_ridge_cv",
            feature_set="biological",
            target_strategy="chronological_age",
            ridge_alpha=ridge_alpha,
        ),
        _cv_metric_row(
            data,
            technical,
            y,
            model_config,
            control="technical_only_ridge_cv",
            feature_set="technical_covariates",
            target_strategy="chronological_age",
            ridge_alpha=ridge_alpha,
        ),
        _cv_metric_row(
            data,
            yield_only,
            y,
            model_config,
            control="yield_only_ridge_cv",
            feature_set="log10_total_cells",
            target_strategy="chronological_age",
            ridge_alpha=ridge_alpha,
        ),
        _cv_metric_row(
            data,
            pd.DataFrame(index=data.index),
            y,
            model_config,
            control="null_mean_cv",
            feature_set="none",
            target_strategy="chronological_age",
            ridge_alpha=ridge_alpha,
            null_model=True,
        ),
        _not_applicable_row(
            control="disease_label_negative_control",
            feature_set="disease_group",
            target_strategy="healthy_training_only",
            n=int(data.shape[0]),
            detail="Primary ORA training contains healthy age-known donors only, so disease-label prediction is not applicable.",
        ),
    ]
    performance_rows.extend(
        _age_shuffle_rows(
            data,
            biological,
            y,
            model_config,
            n_shuffles=n_shuffles,
            ridge_alpha=ridge_alpha,
        )
    )
    performance = pd.DataFrame(performance_rows)
    comparison = _baseline_comparison(performance)
    explainability = covariate_explainability(scores, manifest)
    return NegativeControlResult(
        performance=performance,
        baseline_comparison=comparison,
        covariate_explainability=explainability,
    )


def covariate_explainability(scores: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    """Summarize how much technical covariates explain ORA predictions."""

    if scores.empty:
        return pd.DataFrame()
    meta = _manifest_meta(manifest)
    frame = scores.copy()
    frame["donor_id"] = frame["donor_id"].astype(str)
    frame = frame.merge(meta.reset_index(), on="donor_id", how="left", suffixes=("", "_manifest"))
    frame["chronological_age"] = pd.to_numeric(
        frame.get("chronological_age", frame.get("age")),
        errors="coerce",
    )
    frame["ora"] = pd.to_numeric(frame.get("ora"), errors="coerce")
    frame["oraa"] = pd.to_numeric(frame.get("oraa"), errors="coerce")
    frame["log10_total_cells"] = np.log10(pd.to_numeric(frame.get("total_cells"), errors="coerce") + 1.0)

    rows: list[dict[str, Any]] = []
    covariate_sets: list[tuple[str, tuple[str, ...]]] = [
        ("total_cell_yield", ("log10_total_cells",)),
        ("chemistry", ("chemistry",)),
        ("collection_method", ("collection_method",)),
        ("site", ("site",)),
        ("all_technical", DEFAULT_TECHNICAL_COVARIATES),
    ]
    for model, model_frame in frame.groupby("model", observed=True, sort=False):
        if str(model) == "null_model":
            continue
        for covariate_name, covariates in covariate_sets:
            rows.append(
                _explainability_row(
                    model_frame,
                    model=str(model),
                    outcome="ora",
                    covariate_name=covariate_name,
                    covariates=covariates,
                    age_adjusted=True,
                )
            )
            rows.append(
                _explainability_row(
                    model_frame,
                    model=str(model),
                    outcome="oraa",
                    covariate_name=covariate_name,
                    covariates=covariates,
                    age_adjusted=False,
                )
            )
    return pd.DataFrame(rows).sort_values(["model", "outcome", "covariate_set"]).reset_index(drop=True)


def _training_frame(features: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    meta = _manifest_meta(manifest)
    data = features.copy()
    data["donor_id"] = data["donor_id"].astype(str)
    joined = data.set_index("donor_id").join(meta, how="inner", rsuffix="_manifest").reset_index()
    joined["age"] = pd.to_numeric(joined["age"], errors="coerce")
    if "usable_for_ora_training" in joined.columns:
        mask = _as_bool(joined["usable_for_ora_training"]) & joined["age"].notna()
    else:
        mask = joined["age"].notna()
    train = joined.loc[mask].copy().reset_index(drop=True)
    if train.shape[0] < 5:
        raise ValueError("At least five ORA training donors are required for negative controls.")
    train["log10_total_cells"] = np.log10(pd.to_numeric(train.get("total_cells"), errors="coerce") + 1.0)
    return train


def _manifest_meta(manifest: pd.DataFrame) -> pd.DataFrame:
    required = {"donor_id", "age"}
    missing = sorted(required.difference(manifest.columns))
    if missing:
        raise KeyError(f"Manifest missing required columns: {missing}")
    sort_cols = [col for col in ["donor_id", "sample_id"] if col in manifest.columns]
    meta = manifest.sort_values(sort_cols).drop_duplicates("donor_id").copy()
    meta["donor_id"] = meta["donor_id"].astype(str)
    return meta.set_index("donor_id")


def _biological_matrix(data: pd.DataFrame, model_config: dict[str, Any]) -> pd.DataFrame:
    cols = biological_feature_columns(data, model_config)
    max_missing = float(model_config.get("missingness_max_fraction", 0.30))
    if cols:
        missing_fraction = data[cols].isna().mean()
        cols = [col for col in cols if missing_fraction[col] <= max_missing]
    return data[cols].copy() if cols else pd.DataFrame(index=data.index)


def _technical_matrix(data: pd.DataFrame, covariates: tuple[str, ...]) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for covariate in covariates:
        if covariate not in data.columns:
            continue
        series = data[covariate]
        if pd.api.types.is_numeric_dtype(series):
            numeric = pd.to_numeric(series, errors="coerce")
            if numeric.notna().sum() and numeric.nunique(dropna=True) > 1:
                parts.append(pd.DataFrame({covariate: numeric}, index=data.index))
            continue
        clean = _clean_category(series)
        if clean.nunique(dropna=True) <= 1:
            continue
        dummies = pd.get_dummies(clean, prefix=covariate, drop_first=True, dtype=float)
        dummies.index = data.index
        parts.append(dummies)
    if not parts:
        return pd.DataFrame(index=data.index)
    return pd.concat(parts, axis=1)


def _cv_metric_row(
    data: pd.DataFrame,
    x_frame: pd.DataFrame,
    y: np.ndarray,
    model_config: dict[str, Any],
    *,
    control: str,
    feature_set: str,
    target_strategy: str,
    ridge_alpha: float,
    repeat: int = 0,
    null_model: bool = False,
) -> dict[str, Any]:
    if not null_model and x_frame.shape[1] == 0:
        return _empty_metric_row(
            control=control,
            feature_set=feature_set,
            target_strategy=target_strategy,
            repeat=repeat,
            n=int(data.shape[0]),
            status="no_features",
        )
    pred = np.full(data.shape[0], np.nan, dtype=float)
    folds = donor_cv_folds(data[["donor_id", "age"]].copy(), model_config)
    for train_idx, test_idx in folds:
        y_train = y[train_idx]
        if null_model:
            pred[test_idx] = float(np.mean(y_train))
            continue
        prep = fit_preprocessor(x_frame.iloc[train_idx])
        x_train = transform_preprocessor(x_frame.iloc[train_idx], prep)
        x_test = transform_preprocessor(x_frame.iloc[test_idx], prep)
        coef, intercept = _ridge_closed_form(x_train, y_train, alpha=ridge_alpha)
        pred[test_idx] = x_test @ coef + intercept
    row = _metrics(y, pred)
    row.update(
        {
            "control": control,
            "feature_set": feature_set,
            "target_strategy": target_strategy,
            "repeat": repeat,
            "n_features": int(x_frame.shape[1]),
            "ridge_alpha": float(ridge_alpha),
            "status": "ok",
            "spearman_vs_chronological_age": _spearman(pd.to_numeric(data["age"], errors="coerce"), pred),
            "detail": "",
        }
    )
    return row


def _age_shuffle_rows(
    data: pd.DataFrame,
    biological: pd.DataFrame,
    y: np.ndarray,
    model_config: dict[str, Any],
    *,
    n_shuffles: int,
    ridge_alpha: float,
) -> list[dict[str, Any]]:
    rows = []
    base_seed = int(model_config.get("random_seed", 42))
    strata = _shuffle_strata(data)
    for repeat in range(n_shuffles):
        rng = np.random.default_rng(base_seed + 1000 + repeat)
        shuffled = pd.Series(y.copy(), index=data.index, dtype=float)
        for _, idx in strata.groupby(strata, dropna=False).groups.items():
            index = np.array(list(idx), dtype=int)
            if index.size > 1:
                shuffled.iloc[index] = rng.permutation(shuffled.iloc[index].to_numpy(dtype=float))
        row = _cv_metric_row(
            data,
            biological,
            shuffled.to_numpy(dtype=float),
            model_config,
            control="age_shuffle_within_technical_strata",
            feature_set="biological",
            target_strategy="age_shuffled_within_sex_chemistry_collection",
            ridge_alpha=ridge_alpha,
            repeat=repeat,
        )
        rows.append(row)
    return rows


def _baseline_comparison(performance: pd.DataFrame) -> pd.DataFrame:
    ok = performance[performance["status"].eq("ok")].copy()
    observed = ok[ok["control"].eq("biological_ridge_cv")]
    biological_mae = float(observed["mae"].iloc[0]) if not observed.empty else np.nan
    rows: list[dict[str, Any]] = []
    for control in ["technical_only_ridge_cv", "yield_only_ridge_cv", "null_mean_cv"]:
        frame = ok[ok["control"].eq(control)]
        if frame.empty:
            continue
        mae = float(frame["mae"].iloc[0])
        rows.append(
            {
                "comparison": control,
                "n": int(frame["n"].iloc[0]),
                "mae": mae,
                "rmse": float(frame["rmse"].iloc[0]),
                "r2": float(frame["r2"].iloc[0]),
                "spearman_r": float(frame["spearman_r"].iloc[0]),
                "biological_ridge_mae": biological_mae,
                "delta_mae_vs_biological_ridge": mae - biological_mae,
                "interpretation_vs_biological_ridge": _baseline_interpretation(mae, biological_mae),
            }
        )
    shuffle = ok[ok["control"].eq("age_shuffle_within_technical_strata")]
    if not shuffle.empty:
        shuffle_mae = pd.to_numeric(shuffle["mae"], errors="coerce")
        empirical_p = float((shuffle_mae <= biological_mae).mean()) if np.isfinite(biological_mae) else np.nan
        rows.append(
            {
                "comparison": "age_shuffle_within_technical_strata",
                "n": int(shuffle["n"].median()),
                "mae": float(shuffle_mae.mean()),
                "rmse": float(pd.to_numeric(shuffle["rmse"], errors="coerce").mean()),
                "r2": float(pd.to_numeric(shuffle["r2"], errors="coerce").mean()),
                "spearman_r": float(pd.to_numeric(shuffle["spearman_r"], errors="coerce").mean()),
                "biological_ridge_mae": biological_mae,
                "delta_mae_vs_biological_ridge": float(shuffle_mae.mean() - biological_mae),
                "shuffle_mae_ci_low": float(shuffle_mae.quantile(0.025)),
                "shuffle_mae_ci_high": float(shuffle_mae.quantile(0.975)),
                "shuffle_empirical_p_mae_le_observed": empirical_p,
                "interpretation_vs_biological_ridge": _baseline_interpretation(float(shuffle_mae.mean()), biological_mae),
            }
        )
    return pd.DataFrame(rows)


def _explainability_row(
    frame: pd.DataFrame,
    *,
    model: str,
    outcome: str,
    covariate_name: str,
    covariates: tuple[str, ...],
    age_adjusted: bool,
) -> dict[str, Any]:
    y = pd.to_numeric(frame.get(outcome), errors="coerce")
    base_covariates: tuple[str, ...] = ("chronological_age",) if age_adjusted else ()
    base_design = _design_for_covariates(frame, base_covariates)
    full_design = _design_for_covariates(frame, base_covariates + covariates)
    if full_design.shape[1] <= base_design.shape[1]:
        return {
            "model": model,
            "outcome": outcome,
            "covariate_set": covariate_name,
            "age_adjusted": bool(age_adjusted),
            "n": int(y.notna().sum()),
            "base_r2": np.nan,
            "full_r2": np.nan,
            "incremental_r2": np.nan,
            "status": "invariant_or_missing_covariate",
        }
    valid = y.notna() & np.isfinite(base_design).all(axis=1) & np.isfinite(full_design).all(axis=1)
    if int(valid.sum()) <= full_design.shape[1] + 1:
        return {
            "model": model,
            "outcome": outcome,
            "covariate_set": covariate_name,
            "age_adjusted": bool(age_adjusted),
            "n": int(valid.sum()),
            "base_r2": np.nan,
            "full_r2": np.nan,
            "incremental_r2": np.nan,
            "status": "underdetermined",
        }
    yv = y[valid].to_numpy(dtype=float)
    base_r2 = _ols_r2(yv, base_design[valid])
    full_r2 = _ols_r2(yv, full_design[valid])
    return {
        "model": model,
        "outcome": outcome,
        "covariate_set": covariate_name,
        "age_adjusted": bool(age_adjusted),
        "n": int(valid.sum()),
        "base_r2": base_r2,
        "full_r2": full_r2,
        "incremental_r2": full_r2 - base_r2,
        "status": "ok",
    }


def _design_for_covariates(frame: pd.DataFrame, covariates: tuple[str, ...]) -> np.ndarray:
    parts = [np.ones(frame.shape[0], dtype=float)]
    for covariate in covariates:
        if covariate not in frame.columns:
            continue
        series = frame[covariate]
        if pd.api.types.is_numeric_dtype(series):
            numeric = pd.to_numeric(series, errors="coerce")
            if numeric.notna().sum() and numeric.nunique(dropna=True) > 1:
                parts.append(numeric.fillna(numeric.median()).to_numpy(dtype=float))
            continue
        clean = _clean_category(series)
        if clean.nunique(dropna=True) <= 1:
            continue
        dummies = pd.get_dummies(clean, drop_first=True, dtype=float)
        for col in dummies.columns:
            parts.append(dummies[col].to_numpy(dtype=float))
    return np.vstack(parts).T


def _metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, Any]:
    valid = np.isfinite(y) & np.isfinite(pred)
    yv = y[valid]
    pv = pred[valid]
    resid = yv - pv
    ss_tot = float(np.sum((yv - np.mean(yv)) ** 2))
    ss_res = float(np.sum(resid**2))
    return {
        "n": int(yv.size),
        "mae": float(np.mean(np.abs(resid))) if yv.size else np.nan,
        "rmse": float(np.sqrt(np.mean(resid**2))) if yv.size else np.nan,
        "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan,
        "spearman_r": _spearman(yv, pv) if yv.size > 1 else np.nan,
    }


def _ridge_closed_form(x: np.ndarray, y: np.ndarray, alpha: float) -> tuple[np.ndarray, float]:
    y_mean = float(np.mean(y))
    coef = np.linalg.pinv(x.T @ x + alpha * np.eye(x.shape[1])) @ x.T @ (y - y_mean)
    return coef, y_mean


def _ols_r2(y: np.ndarray, x: np.ndarray) -> float:
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    pred = x @ coef
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan


def _spearman(a: Any, b: Any) -> float:
    ar = pd.Series(a).rank(method="average").to_numpy(dtype=float)
    br = pd.Series(b).rank(method="average").to_numpy(dtype=float)
    valid = np.isfinite(ar) & np.isfinite(br)
    if valid.sum() < 2 or np.std(ar[valid]) == 0 or np.std(br[valid]) == 0:
        return np.nan
    return float(np.corrcoef(ar[valid], br[valid])[0, 1])


def _shuffle_strata(data: pd.DataFrame) -> pd.Series:
    parts = []
    for covariate in ["sex", "chemistry", "collection_method"]:
        if covariate in data.columns:
            parts.append(_clean_category(data[covariate]))
    if not parts:
        return pd.Series("all", index=data.index)
    strata = parts[0].copy()
    for part in parts[1:]:
        strata = strata + "|" + part
    return strata


def _baseline_interpretation(control_mae: float, biological_mae: float) -> str:
    if not np.isfinite(control_mae) or not np.isfinite(biological_mae):
        return "not_evaluable"
    if control_mae <= biological_mae:
        return "control_matches_or_beats_biological_ridge"
    if control_mae - biological_mae < 1.0:
        return "control_close_to_biological_ridge"
    return "control_worse_than_biological_ridge"


def _empty_metric_row(
    *,
    control: str,
    feature_set: str,
    target_strategy: str,
    repeat: int,
    n: int,
    status: str,
) -> dict[str, Any]:
    return {
        "control": control,
        "feature_set": feature_set,
        "target_strategy": target_strategy,
        "repeat": repeat,
        "n": n,
        "mae": np.nan,
        "rmse": np.nan,
        "r2": np.nan,
        "spearman_r": np.nan,
        "spearman_vs_chronological_age": np.nan,
        "n_features": 0,
        "ridge_alpha": np.nan,
        "status": status,
        "detail": "",
    }


def _not_applicable_row(*, control: str, feature_set: str, target_strategy: str, n: int, detail: str) -> dict[str, Any]:
    row = _empty_metric_row(
        control=control,
        feature_set=feature_set,
        target_strategy=target_strategy,
        repeat=0,
        n=n,
        status="not_applicable",
    )
    row["detail"] = detail
    return row


def _clean_category(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("unknown").str.strip().replace("", "unknown")


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).astype(float).ne(0)
    return series.fillna("").astype(str).str.strip().str.lower().isin({"1", "true", "t", "yes", "y"})
