"""Annotation helpers for Milo-style neighborhood differential-abundance tables."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd

from ora.interpretation import classify_feature_theme


DEFAULT_FDR_THRESHOLD = 0.10


def annotate_neighborhood_table(
    neighborhoods: pd.DataFrame,
    *,
    run_name: str,
    fdr_threshold: float = DEFAULT_FDR_THRESHOLD,
) -> pd.DataFrame:
    """Add manuscript-facing theme and claim-gate annotations to one DA table."""

    required = {"top_fine_celltype", "top_coarse_celltype", "age_coef", "age_fdr", "status"}
    missing = required.difference(neighborhoods.columns)
    if missing:
        raise ValueError(f"Neighborhood table is missing required columns: {sorted(missing)}")

    frame = neighborhoods.copy()
    frame.insert(0, "run", run_name)
    frame["age_coef"] = pd.to_numeric(frame["age_coef"], errors="coerce")
    frame["age_fdr"] = pd.to_numeric(frame["age_fdr"], errors="coerce")
    frame["age_direction"] = np.where(frame["age_coef"].gt(0), "positive", np.where(frame["age_coef"].lt(0), "negative", "zero"))
    frame["is_tested"] = frame["status"].astype(str).eq("tested")
    frame["is_age_associated_fdr_0_10"] = frame["is_tested"] & frame["age_fdr"].lt(fdr_threshold)
    frame["fine_theme"] = frame["top_fine_celltype"].astype(str).map(lambda label: classify_feature_theme(f"prop__{label}"))
    frame["coarse_theme"] = frame["top_coarse_celltype"].astype(str).map(lambda label: classify_feature_theme(f"prop__{label}"))
    frame["claim_gate"] = frame.apply(_claim_gate, axis=1)
    return frame


def summarize_neighborhood_themes(annotated: pd.DataFrame) -> pd.DataFrame:
    """Summarize significant neighborhoods by run, theme, direction, and top labels."""

    columns = [
        "run",
        "fine_theme",
        "age_direction",
        "n_neighborhoods",
        "n_significant",
        "median_age_coef",
        "min_age_fdr",
        "top_fine_celltypes",
        "top_coarse_celltypes",
        "claim_gate",
    ]
    if annotated.empty:
        return pd.DataFrame(columns=columns)

    sig = annotated[annotated["is_age_associated_fdr_0_10"]].copy()
    if sig.empty:
        return pd.DataFrame(columns=columns)

    rows = []
    for (run, theme, direction), group in sig.groupby(["run", "fine_theme", "age_direction"], observed=True):
        rows.append(
            {
                "run": run,
                "fine_theme": theme,
                "age_direction": direction,
                "n_neighborhoods": int(group.shape[0]),
                "n_significant": int(group["is_age_associated_fdr_0_10"].sum()),
                "median_age_coef": float(group["age_coef"].median()),
                "min_age_fdr": float(group["age_fdr"].min()),
                "top_fine_celltypes": _top_labels(group["top_fine_celltype"]),
                "top_coarse_celltypes": _top_labels(group["top_coarse_celltype"]),
                "claim_gate": _theme_claim_gate(group),
            }
        )
    result = pd.DataFrame(rows, columns=columns)
    return result.sort_values(["run", "n_significant", "min_age_fdr"], ascending=[True, False, True]).reset_index(drop=True)


def build_neighborhood_annotation(
    tables: Mapping[str, pd.DataFrame],
    *,
    fdr_threshold: float = DEFAULT_FDR_THRESHOLD,
    top_n: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build annotated top-neighborhood and theme-summary tables from multiple runs."""

    annotated_frames = [
        annotate_neighborhood_table(table, run_name=run_name, fdr_threshold=fdr_threshold)
        for run_name, table in tables.items()
    ]
    if not annotated_frames:
        return pd.DataFrame(), pd.DataFrame()
    annotated = pd.concat(annotated_frames, ignore_index=True)
    top = (
        annotated.sort_values(["run", "is_age_associated_fdr_0_10", "age_fdr"], ascending=[True, False, True])
        .groupby("run", observed=True, group_keys=False)
        .head(top_n)
        .reset_index(drop=True)
    )
    summary = summarize_neighborhood_themes(annotated)
    return top, summary


def _claim_gate(row: pd.Series) -> str:
    if not bool(row.get("is_age_associated_fdr_0_10", False)):
        return "not_significant"
    run = str(row.get("run", ""))
    theme = str(row.get("fine_theme", ""))
    if "matched" in run and "neuronal lineage" in theme:
        return "matched_regenerative_neuronal_support"
    if "matched" in run and "immune" in theme:
        return "matched_immune_support_interpret_with_balance"
    if "matched" in run:
        return "matched_exploratory"
    if "secretory" in run:
        return "all_donor_secretory_exploratory"
    return "all_donor_hypothesis_map"


def _theme_claim_gate(group: pd.DataFrame) -> str:
    gates = group["claim_gate"].astype(str).value_counts()
    if gates.empty:
        return "not_significant"
    return str(gates.index[0])


def _top_labels(values: pd.Series, *, limit: int = 6) -> str:
    counts = values.astype(str).value_counts().head(limit)
    return ";".join(f"{label}:{count}" for label, count in counts.items())
