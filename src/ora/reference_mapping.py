"""Helpers for reference-mapped external single-cell labels."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def mapped_label_donor_features(
    obs: pd.DataFrame,
    *,
    label_column: str,
    confidence_column: str = "",
    group_columns: tuple[str, ...] = ("dataset_id", "sample_id", "donor_id", "age", "disease_group"),
) -> pd.DataFrame:
    """Summarize mapped cell labels as ORA-compatible donor/sample features."""

    if obs.empty or label_column not in obs:
        return pd.DataFrame()
    available_groups = tuple(column for column in group_columns if column in obs)
    rows: list[dict[str, Any]] = []
    for keys, group in obs.groupby(list(available_groups), dropna=False, observed=True):
        key_values = keys if isinstance(keys, tuple) else (keys,)
        row = dict(zip(available_groups, key_values, strict=False))
        labels = group[label_column].astype(str)
        counts = labels.value_counts()
        total = float(counts.sum())
        fractions = counts / total if total else counts.astype(float)
        row["n_cells"] = int(total)
        if confidence_column and confidence_column in group:
            row["mean_label_confidence"] = float(pd.to_numeric(group[confidence_column], errors="coerce").mean())
        label_names = sorted(fractions.index)
        geometric = _geometric_fraction_mean(fractions, label_names)
        for label in label_names:
            feature_label = _feature_safe_label(label)
            fraction = float(fractions.get(label, 0.0))
            row[f"prop__{feature_label}"] = fraction
            row[f"clr__{feature_label}"] = float(np.log((fraction + 1e-6) / geometric)) if geometric > 0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values(list(available_groups)).reset_index(drop=True)


def mapping_qc_by_sample(
    obs: pd.DataFrame,
    *,
    label_column: str,
    confidence_column: str,
    entropy_column: str,
    group_columns: tuple[str, ...] = ("dataset_id", "sample_id", "donor_id", "age", "disease_group"),
) -> pd.DataFrame:
    """Summarize per-sample reference mapping confidence and label diversity."""

    if obs.empty:
        return pd.DataFrame(columns=[*group_columns, "n_cells", "n_labels", "mean_confidence", "median_confidence", "mean_entropy"])
    available_groups = tuple(column for column in group_columns if column in obs)
    rows: list[dict[str, Any]] = []
    for keys, group in obs.groupby(list(available_groups), dropna=False, observed=True):
        key_values = keys if isinstance(keys, tuple) else (keys,)
        row = dict(zip(available_groups, key_values, strict=False))
        row["n_cells"] = int(group.shape[0])
        row["n_labels"] = int(group[label_column].astype(str).nunique()) if label_column in group else 0
        row["mean_confidence"] = _numeric_mean(group, confidence_column)
        row["median_confidence"] = _numeric_median(group, confidence_column)
        row["mean_entropy"] = _numeric_mean(group, entropy_column)
        row["status"] = _mapping_qc_status(row["mean_confidence"], row["mean_entropy"])
        rows.append(row)
    return pd.DataFrame(rows).sort_values(list(available_groups)).reset_index(drop=True)


def normalized_entropy(probabilities: pd.DataFrame | np.ndarray) -> np.ndarray:
    """Return row-wise entropy scaled to [0, 1]."""

    probs = probabilities.to_numpy(dtype=float) if isinstance(probabilities, pd.DataFrame) else np.asarray(probabilities, dtype=float)
    if probs.ndim != 2 or probs.shape[1] <= 1:
        return np.zeros(probs.shape[0] if probs.ndim else 0, dtype=float)
    clipped = np.clip(probs, 1e-12, 1.0)
    entropy = -(clipped * np.log(clipped)).sum(axis=1)
    return entropy / np.log(probs.shape[1])


def _feature_safe_label(label: str) -> str:
    return str(label).strip().lower().replace(" ", "_").replace("-", "_").replace("/", "_")


def _geometric_fraction_mean(fractions: pd.Series, labels: list[str]) -> float:
    if not labels:
        return np.nan
    return float(np.exp(np.mean(np.log([float(fractions.get(label, 0.0)) + 1e-6 for label in labels]))))


def _numeric_mean(frame: pd.DataFrame, column: str) -> float:
    if column not in frame:
        return np.nan
    return float(pd.to_numeric(frame[column], errors="coerce").mean())


def _numeric_median(frame: pd.DataFrame, column: str) -> float:
    if column not in frame:
        return np.nan
    return float(pd.to_numeric(frame[column], errors="coerce").median())


def _mapping_qc_status(mean_confidence: float, mean_entropy: float) -> str:
    if np.isnan(mean_confidence):
        return "missing_confidence"
    if mean_confidence >= 0.70 and (np.isnan(mean_entropy) or mean_entropy <= 0.55):
        return "high_confidence"
    if mean_confidence >= 0.45:
        return "moderate_confidence"
    return "low_confidence"
