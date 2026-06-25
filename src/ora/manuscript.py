"""Manuscript-readiness summary tables."""

from __future__ import annotations

import pandas as pd


def build_model_card(
    *,
    feature_set_comparison: pd.DataFrame | None = None,
    calibration: pd.DataFrame | None = None,
    permutation: pd.DataFrame | None = None,
    nested_tuning: pd.DataFrame | None = None,
    stacking: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build a compact ORA model-card table from existing benchmark outputs."""

    rows = []
    comparison = feature_set_comparison.copy() if feature_set_comparison is not None else pd.DataFrame()
    if comparison.empty:
        return pd.DataFrame(columns=_model_card_columns())
    for _, row in comparison.iterrows():
        model = str(row.get("model", ""))
        feature_set = str(row.get("feature_set", ""))
        rows.append(
            {
                "model": model,
                "feature_set": feature_set,
                "role": _model_role(model, feature_set, row),
                "training_cohort": "healthy donors with valid age only",
                "excluded_from_training": "AD/PD donors; donors missing age",
                "cv_design": "donor-level repeated outer CV",
                "n": row.get("n", pd.NA),
                "repeats": row.get("repeats", pd.NA),
                "backend": row.get("backend", pd.NA),
                "backend_package": row.get("backend_package", pd.NA),
                "backend_version": row.get("backend_version", pd.NA),
                "fallback_used": row.get("fallback_used", pd.NA),
                "fallback_reason": row.get("fallback_reason", pd.NA),
                "mae_mean": row.get("mae_mean", pd.NA),
                "mae_ci_low": row.get("mae_ci_low", pd.NA),
                "mae_ci_high": row.get("mae_ci_high", pd.NA),
                "spearman_r_mean": row.get("spearman_r_mean", pd.NA),
                "calibration_slope": _lookup(calibration, model, "calibration_slope_ora_on_age"),
                "permutation_p_mae": _lookup(permutation, model, "empirical_p_mae"),
                "nested_tuning_mae_mean": _lookup(nested_tuning, model, "mae_mean"),
                "stacking_mae_mean": _lookup(stacking, model, "mae_mean") if model == "stacked_ensemble" else pd.NA,
                "limitations": _limitations(model, feature_set),
            }
        )
    return pd.DataFrame(rows, columns=_model_card_columns()).sort_values(["role", "mae_mean", "model"]).reset_index(drop=True)


def _model_role(model: str, feature_set: str, row: pd.Series) -> str:
    if model == "null_model":
        return "negative_control"
    if bool(row.get("is_best_overall", False)):
        return "preferred_benchmark"
    if feature_set == "composition":
        return "interpretable_baseline"
    return "secondary_benchmark"


def _lookup(frame: pd.DataFrame | None, model: str, column: str) -> object:
    if frame is None or frame.empty or "model" not in frame or column not in frame:
        return pd.NA
    rows = frame[frame["model"].astype(str).eq(model)]
    return rows.iloc[0][column] if not rows.empty else pd.NA


def _limitations(model: str, feature_set: str) -> str:
    notes = [
        "under-dispersed age predictions",
        "donor n=187 limits model complexity",
        "chemistry and collection method require sensitivity interpretation",
    ]
    if feature_set == "composition_plus_modules":
        notes.append("module gain is modest with overlapping intervals")
    if model in {"xgboost", "lightgbm", "catboost", "boosted_ensemble"}:
        notes.append("native booster result is benchmark-only unless interpretability is secondary")
    if model == "null_model":
        notes = ["negative-control baseline only"]
    return "; ".join(notes)


def _model_card_columns() -> list[str]:
    return [
        "model",
        "feature_set",
        "role",
        "training_cohort",
        "excluded_from_training",
        "cv_design",
        "n",
        "repeats",
        "backend",
        "backend_package",
        "backend_version",
        "fallback_used",
        "fallback_reason",
        "mae_mean",
        "mae_ci_low",
        "mae_ci_high",
        "spearman_r_mean",
        "calibration_slope",
        "permutation_p_mae",
        "nested_tuning_mae_mean",
        "stacking_mae_mean",
        "limitations",
    ]
