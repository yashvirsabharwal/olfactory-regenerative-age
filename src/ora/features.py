"""Feature matrix construction helpers."""

from __future__ import annotations

import pandas as pd


COMPOSITION_PREFIXES = ("prop__", "clr__", "ratio__")
MODULE_PREFIXES = ("module_score__",)
SCVI_GLOBAL_PREFIXES = ("scvi_global_",)
SCVI_STATE_PREFIXES = ("scvi_state_",)

FEATURE_FAMILY_COLUMNS = [
    "model",
    "feature_family",
    "features",
    "selected_features",
    "mean_abs_importance",
    "total_abs_importance",
    "max_selection_fraction",
    "family_rank_within_model",
]


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
        "scvi_global": sum(col.startswith(SCVI_GLOBAL_PREFIXES) for col in matrix.columns),
        "scvi_cell_state": sum(col.startswith(SCVI_STATE_PREFIXES) for col in matrix.columns),
    }


def merge_donor_feature_matrices(tables: list[pd.DataFrame]) -> pd.DataFrame:
    """Merge donor-level feature matrices, rejecting duplicate feature columns."""

    if not tables:
        raise ValueError("At least one donor feature table is required.")
    merged: pd.DataFrame | None = None
    seen: set[str] = set()
    for idx, table in enumerate(tables):
        _require_donor_id(table, f"donor feature table {idx + 1}")
        columns = [col for col in table.columns if col != "donor_id"]
        duplicates = sorted(seen.intersection(columns))
        if duplicates:
            raise ValueError(f"Duplicate donor feature columns: {', '.join(duplicates)}")
        seen.update(columns)
        deduped = table.drop_duplicates("donor_id").copy()
        if merged is None:
            merged = deduped
        else:
            merged = merged.merge(deduped, on="donor_id", how="outer")
    assert merged is not None
    return merged.sort_values("donor_id").reset_index(drop=True)


def summarize_feature_family_stability(feature_stability: pd.DataFrame) -> pd.DataFrame:
    """Summarize repeated-CV feature stability by broad biological feature family."""

    if feature_stability.empty:
        return pd.DataFrame(columns=FEATURE_FAMILY_COLUMNS)
    required = {"model", "feature", "abs_mean_importance", "selection_fraction"}
    missing = sorted(required.difference(feature_stability.columns))
    if missing:
        raise ValueError(f"Feature stability table is missing required columns: {', '.join(missing)}")

    frame = feature_stability.copy()
    frame["feature_family"] = frame["feature"].map(feature_family)
    frame["abs_mean_importance"] = pd.to_numeric(frame["abs_mean_importance"], errors="coerce").fillna(0.0)
    frame["selection_fraction"] = pd.to_numeric(frame["selection_fraction"], errors="coerce").fillna(0.0)
    summary = (
        frame.groupby(["model", "feature_family"], observed=True)
        .agg(
            features=("feature", "nunique"),
            selected_features=("selection_fraction", lambda s: int((pd.to_numeric(s, errors="coerce") > 0).sum())),
            mean_abs_importance=("abs_mean_importance", "mean"),
            total_abs_importance=("abs_mean_importance", "sum"),
            max_selection_fraction=("selection_fraction", "max"),
        )
        .reset_index()
    )
    summary["family_rank_within_model"] = (
        summary.groupby("model", observed=True)["total_abs_importance"]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    return summary.sort_values(
        ["model", "family_rank_within_model", "total_abs_importance", "feature_family"],
        ascending=[True, True, False, True],
    )[FEATURE_FAMILY_COLUMNS].reset_index(drop=True)


def feature_family(feature: object) -> str:
    """Classify a donor-level feature name into a broad family."""

    text = str(feature)
    if text.startswith(COMPOSITION_PREFIXES):
        return "composition"
    if text.startswith(MODULE_PREFIXES):
        return "module"
    if text.startswith(SCVI_GLOBAL_PREFIXES):
        return "scvi_global"
    if text.startswith(SCVI_STATE_PREFIXES):
        return "scvi_cell_state"
    return "other"


def _require_donor_id(frame: pd.DataFrame, context: str) -> None:
    if "donor_id" not in frame.columns:
        raise KeyError(f"{context} must include donor_id.")
