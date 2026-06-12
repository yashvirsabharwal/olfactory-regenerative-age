"""Feature matrix construction helpers."""

from __future__ import annotations

import pandas as pd


COMPOSITION_PREFIXES = ("prop__", "clr__", "ratio__")
MODULE_PREFIXES = ("module_score__",)


def build_ora_feature_matrix(
    cell_features: pd.DataFrame,
    module_features: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build a donor-level ORA feature matrix from composition and optional modules."""

    _require_donor_id(cell_features, "cell feature table")
    keep = ["donor_id"] + [
        col for col in cell_features.columns if col.startswith(COMPOSITION_PREFIXES)
    ]
    matrix = cell_features[keep].copy()
    if module_features is None:
        return matrix

    _require_donor_id(module_features, "module feature table")
    module_cols = [
        col
        for col in module_features.columns
        if col.startswith(MODULE_PREFIXES) and pd.api.types.is_numeric_dtype(module_features[col])
    ]
    if not module_cols:
        raise ValueError("Module feature table has no numeric module_score__ columns.")
    matrix = matrix.merge(module_features[["donor_id", *module_cols]], on="donor_id", how="left")
    return matrix


def feature_kind_counts(matrix: pd.DataFrame) -> dict[str, int]:
    """Count major feature families in a donor-level matrix."""

    return {
        "composition": sum(col.startswith(COMPOSITION_PREFIXES) for col in matrix.columns),
        "module": sum(col.startswith(MODULE_PREFIXES) for col in matrix.columns),
    }


def _require_donor_id(frame: pd.DataFrame, context: str) -> None:
    if "donor_id" not in frame.columns:
        raise KeyError(f"{context} must include donor_id.")
