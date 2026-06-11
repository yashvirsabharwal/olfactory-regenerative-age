"""Lightweight statistics for ORA MVP tables."""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
import pandas as pd


def run_age_associations(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    feature_columns: Iterable[str] | None = None,
    covariates: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Fit one linear model per feature against age and optional covariates."""

    donor_meta = (
        manifest.sort_values(["donor_id", "sample_id"])
        .drop_duplicates("donor_id")
        .set_index("donor_id")
    )
    data = features.set_index("donor_id").join(donor_meta, how="inner", rsuffix="_meta")
    data = data[data["age"].notna()]
    covariates = list(covariates or [])
    if feature_columns is None:
        feature_columns = [
            col
            for col in features.columns
            if col != "donor_id" and pd.api.types.is_numeric_dtype(features[col])
        ]

    rows = []
    for feature in feature_columns:
        cols = ["age", feature, *covariates]
        if feature not in data.columns:
            rows.append(_empty_result(feature, 0, "missing_feature"))
            continue
        frame = data[[col for col in cols if col in data.columns]].copy()
        frame["age"] = pd.to_numeric(frame["age"], errors="coerce")
        frame[feature] = pd.to_numeric(frame[feature], errors="coerce")
        frame = frame.replace([np.inf, -np.inf], np.nan).dropna(subset=["age", feature])
        if frame.shape[0] < 5:
            rows.append(_empty_result(feature, frame.shape[0], "too_few_samples"))
            continue
        y = pd.to_numeric(frame[feature], errors="coerce").to_numpy(dtype=float)
        x_parts = [np.ones(frame.shape[0]), pd.to_numeric(frame["age"], errors="coerce").to_numpy(dtype=float)]
        design_names = ["intercept", "age"]
        for covariate in covariates:
            if covariate not in frame.columns:
                continue
            series = frame[covariate]
            if pd.api.types.is_numeric_dtype(series):
                numeric = pd.to_numeric(series, errors="coerce")
                if not numeric.notna().any():
                    continue
                x_parts.append(numeric.fillna(numeric.median()).to_numpy(dtype=float))
                design_names.append(covariate)
            else:
                categorical = series.astype("string").fillna("unknown").str.strip().replace("", "unknown")
                if categorical.nunique(dropna=True) <= 1:
                    continue
                dummies = pd.get_dummies(categorical, prefix=covariate, drop_first=True, dtype=float)
                for dummy_col in dummies.columns:
                    x_parts.append(dummies[dummy_col].to_numpy(dtype=float))
                    design_names.append(dummy_col)
        x = np.vstack(x_parts).T
        valid = np.isfinite(x).all(axis=1) & np.isfinite(y)
        x = x[valid]
        y = y[valid]
        if x.shape[0] <= x.shape[1] + 1:
            rows.append(_empty_result(feature, x.shape[0], "underdetermined"))
            continue
        result = _ols(y, x)
        age_idx = design_names.index("age")
        beta = result["coef"][age_idx] * 10.0
        se = result["se"][age_idx] * 10.0
        t_value = result["coef"][age_idx] / result["se"][age_idx] if result["se"][age_idx] > 0 else np.nan
        rows.append(
            {
                "feature": feature,
                "n": int(x.shape[0]),
                "beta_per_10_years": beta,
                "standard_error": se,
                "t_value": t_value,
                "p_value": _two_sided_normal_p(t_value),
                "direction": "positive" if beta > 0 else "negative" if beta < 0 else "zero",
                "status": "ok",
            }
        )
    output = pd.DataFrame(rows)
    output["fdr"] = bh_fdr(output["p_value"].to_numpy(dtype=float))
    return output.sort_values(["fdr", "p_value", "feature"], na_position="last").reset_index(drop=True)


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    order = np.argsort(p[valid])
    ranked = p[valid][order]
    n = ranked.size
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    valid_indices = np.flatnonzero(valid)
    out[valid_indices[order]] = adjusted
    return out


def _ols(y: np.ndarray, x: np.ndarray) -> dict[str, np.ndarray]:
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ coef
    df = max(x.shape[0] - x.shape[1], 1)
    sigma2 = float((resid @ resid) / df)
    xtx_inv = np.linalg.pinv(x.T @ x)
    se = np.sqrt(np.clip(np.diag(xtx_inv) * sigma2, 0, np.inf))
    return {"coef": coef, "se": se}


def _two_sided_normal_p(t_value: float) -> float:
    if not np.isfinite(t_value):
        return np.nan
    return math.erfc(abs(float(t_value)) / math.sqrt(2.0))


def _empty_result(feature: str, n: int, status: str) -> dict[str, object]:
    return {
        "feature": feature,
        "n": n,
        "beta_per_10_years": np.nan,
        "standard_error": np.nan,
        "t_value": np.nan,
        "p_value": np.nan,
        "direction": "",
        "status": status,
    }
