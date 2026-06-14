"""ORA calibration and residual diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .age_model import _spearman


@dataclass
class ORADiagnosticsResult:
    calibration: pd.DataFrame
    age_bin_errors: pd.DataFrame
    residual_diagnostics: pd.DataFrame
    calibrated_scores: pd.DataFrame


def summarize_ora_diagnostics(
    scores: pd.DataFrame,
    model_config: dict[str, Any] | None = None,
    manifest: pd.DataFrame | None = None,
) -> ORADiagnosticsResult:
    """Summarize out-of-fold ORA calibration and residual diagnostics."""

    model_config = model_config or {}
    frame = _prepare_scores(scores, manifest)
    frame["error"] = frame["ora"] - frame["chronological_age"]
    frame["abs_error"] = frame["error"].abs()
    frame["squared_error"] = frame["error"] ** 2
    frame["age_bin"] = _assign_age_bins(frame["chronological_age"], model_config)
    frame["total_cells_bin"] = _quantile_bins(frame.get("total_cells"))

    calibrated = _add_calibrated_ora(frame)
    return ORADiagnosticsResult(
        calibration=_calibration_summary(calibrated),
        age_bin_errors=_age_bin_summary(calibrated),
        residual_diagnostics=_residual_summary(calibrated),
        calibrated_scores=calibrated.drop(columns=["squared_error"], errors="ignore"),
    )


def _prepare_scores(scores: pd.DataFrame, manifest: pd.DataFrame | None) -> pd.DataFrame:
    required = {"donor_id", "model", "chronological_age", "ora"}
    missing = required.difference(scores.columns)
    if missing:
        raise ValueError(f"ORA scores are missing required columns: {sorted(missing)}")
    frame = scores.copy()
    if manifest is not None and not manifest.empty and "donor_id" in manifest:
        donor_meta = manifest.sort_values(["donor_id", "sample_id"] if "sample_id" in manifest else ["donor_id"])
        donor_meta = donor_meta.drop_duplicates("donor_id")
        extra_cols = ["donor_id"] + [
            col
            for col in ["race_ethnicity", "disease_group", "total_cells", "lineage_cells", "mature_neurons"]
            if col in donor_meta.columns and col not in frame.columns
        ]
        if len(extra_cols) > 1:
            frame = frame.merge(donor_meta[extra_cols], on="donor_id", how="left")
    frame["chronological_age"] = pd.to_numeric(frame["chronological_age"], errors="coerce")
    frame["ora"] = pd.to_numeric(frame["ora"], errors="coerce")
    frame["oraa"] = pd.to_numeric(frame.get("oraa"), errors="coerce")
    for col in ["total_cells", "lineage_cells", "mature_neurons"]:
        if col in frame:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame[frame["chronological_age"].notna() & frame["ora"].notna()].copy()


def _add_calibrated_ora(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    output["calibrated_ora"] = np.nan
    output["calibrated_error"] = np.nan
    for model, sub in output.groupby("model", sort=False):
        idx = sub.index
        pred = sub["ora"].to_numpy(dtype=float)
        age = sub["chronological_age"].to_numpy(dtype=float)
        valid = np.isfinite(pred) & np.isfinite(age)
        if valid.sum() < 2 or float(np.nanstd(pred[valid])) == 0.0:
            calibrated = pred.copy()
        else:
            slope, intercept = np.polyfit(pred[valid], age[valid], 1)
            calibrated = intercept + slope * pred
        output.loc[idx, "calibrated_ora"] = calibrated
        output.loc[idx, "calibrated_error"] = calibrated - age
    output["calibrated_abs_error"] = output["calibrated_error"].abs()
    output["calibrated_squared_error"] = output["calibrated_error"] ** 2
    return output


def _calibration_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model, sub in frame.groupby("model", sort=False):
        age = sub["chronological_age"].to_numpy(dtype=float)
        pred = sub["ora"].to_numpy(dtype=float)
        valid = np.isfinite(age) & np.isfinite(pred)
        if valid.sum() >= 2 and float(np.nanstd(age[valid])) > 0.0:
            raw_slope, raw_intercept = np.polyfit(age[valid], pred[valid], 1)
        else:
            raw_slope, raw_intercept = np.nan, np.nan
        rows.append(
            {
                "model": model,
                "n": int(valid.sum()),
                "calibration_slope_ora_on_age": raw_slope,
                "calibration_intercept_ora_on_age": raw_intercept,
                "mean_error": float(np.nanmean(sub["error"])),
                "mae": float(np.nanmean(sub["abs_error"])),
                "rmse": float(np.sqrt(np.nanmean(sub["squared_error"]))),
                "spearman_r": _spearman(age[valid], pred[valid]) if valid.sum() > 1 else np.nan,
                "calibrated_mean_error": float(np.nanmean(sub["calibrated_error"])),
                "calibrated_mae": float(np.nanmean(sub["calibrated_abs_error"])),
                "calibrated_rmse": float(np.sqrt(np.nanmean(sub["calibrated_squared_error"]))),
            }
        )
    return pd.DataFrame(rows)


def _age_bin_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (model, age_bin), sub in frame.groupby(["model", "age_bin"], observed=True, sort=False):
        rows.append(_summary_row(sub, model=model, group="age_bin", level=age_bin))
    return pd.DataFrame(rows)


def _residual_summary(frame: pd.DataFrame) -> pd.DataFrame:
    covariates = ["sex", "chemistry", "collection_method", "site", "race_ethnicity", "age_bin", "total_cells_bin"]
    rows = []
    for covariate in covariates:
        if covariate not in frame.columns:
            continue
        values = frame[covariate].fillna("unknown").astype(str)
        if values.nunique(dropna=False) == 0:
            continue
        work = frame.assign(_level=values)
        for (model, level), sub in work.groupby(["model", "_level"], observed=True, sort=False):
            rows.append(_summary_row(sub, model=model, group=covariate, level=level))
    return pd.DataFrame(rows)


def _summary_row(sub: pd.DataFrame, *, model: str, group: str, level: str) -> dict[str, object]:
    return {
        "model": model,
        "group": group,
        "level": level,
        "n": int(sub.shape[0]),
        "mean_age": float(sub["chronological_age"].mean()),
        "mean_ora": float(sub["ora"].mean()),
        "mean_error": float(sub["error"].mean()),
        "mae": float(sub["abs_error"].mean()),
        "rmse": float(np.sqrt(sub["squared_error"].mean())),
        "mean_oraa": float(sub["oraa"].mean()) if "oraa" in sub else np.nan,
        "calibrated_mean_error": float(sub["calibrated_error"].mean()),
        "calibrated_mae": float(sub["calibrated_abs_error"].mean()),
    }


def _assign_age_bins(age: pd.Series, model_config: dict[str, Any]) -> pd.Series:
    bins = model_config.get("age_bins", {})
    output = pd.Series("unknown", index=age.index, dtype=object)
    for label, bounds in bins.items():
        if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
            continue
        low, high = bounds
        output[(age >= low) & (age <= high)] = str(label)
    return output


def _quantile_bins(values: pd.Series | None) -> pd.Series:
    if values is None:
        return pd.Series(dtype=object)
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().sum() < 4 or numeric.nunique(dropna=True) < 2:
        return pd.Series("unknown", index=numeric.index, dtype=object)
    try:
        bins = pd.qcut(numeric, q=4, duplicates="drop")
    except ValueError:
        return pd.Series("unknown", index=numeric.index, dtype=object)
    return bins.astype(str).fillna("unknown")
