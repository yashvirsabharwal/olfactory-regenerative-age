"""Publication-facing comparison helpers for scVI validation runs."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

PRIMARY_MODEL = "full_4m_reduced"
CORE_MIXING_METRICS = (
    "neighbor_mixing_entropy__flex_version",
    "neighbor_mixing_entropy__device_guided",
    "neighbor_mixing_entropy__condition",
    "neighbor_mixing_entropy__sex",
)
MARKER_CHECK_PREFIX = "marker_continuity__"
EXPECTED_MARKER_LABELS = {
    "basal": ("hbc", "basal", "suprabasal"),
    "progenitor": ("gbc", "inp", "progenitor", "neurod", "neurog", "ascl"),
    "immature_osn": ("iosn", "inp", "immature", "osn"),
    "mature_osn": ("mosn", "mature", "osn"),
    "sustentacular": ("sus", "sustentacular", "gland", "secretory"),
    "immune": ("dc", "immune", "mono", "mac", "tcell", "bcell", "ptprc"),
}


def compare_scvi_validation_runs(validation_paths: dict[str, str | Path]) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """Summarize validation runs into model and marker claim-gate outputs."""

    validations = {model: _read_validation(path) for model, path in validation_paths.items()}
    model_summary = pd.DataFrame([_summarize_model(model, table) for model, table in validations.items()])
    marker_rows = []
    markers = sorted(
        {
            check.replace(MARKER_CHECK_PREFIX, "", 1)
            for table in validations.values()
            for check in table["check"].astype(str)
            if check.startswith(MARKER_CHECK_PREFIX)
        }
    )
    for marker in markers:
        marker_rows.append(_summarize_marker(marker, validations))
    marker_summary = pd.DataFrame(marker_rows)
    markdown = render_scvi_comparison_note(model_summary, marker_summary)
    return model_summary, marker_summary, markdown


def render_scvi_comparison_note(model_summary: pd.DataFrame, marker_summary: pd.DataFrame) -> str:
    """Render a concise Markdown claim-gate note from comparison tables."""

    primary = model_summary[model_summary["model"].eq(PRIMARY_MODEL)]
    primary_line = "Primary full 4M model was not found in the comparison."
    if not primary.empty:
        row = primary.iloc[0]
        primary_line = (
            f"The primary full 4M reduced scVI model represents {int(row['cells']):,} cells, "
            f"has {int(row['latent_dimensions'])} latent dimensions, fine-label purity "
            f"{row['fine_label_purity']:.3f}, and coarse-label purity {row['coarse_label_purity']:.3f}."
        )

    caveats = marker_summary[marker_summary["claim_gate"].ne("supported")]["marker"].tolist()
    caveat_line = (
        "No marker-specific caveats were detected."
        if not caveats
        else "Marker-specific caveats remain for: " + ", ".join(caveats) + "."
    )
    supported = marker_summary[marker_summary["claim_gate"].eq("supported")]["marker"].tolist()
    support_line = (
        "No marker panels are fully supported across the comparison."
        if not supported
        else "Supported marker panels across the current comparison: " + ", ".join(supported) + "."
    )
    return "\n".join(
        [
            "# scVI Embedding Comparison Claim Gate",
            "",
            "Updated: 2026-06-22",
            "",
            "## Verdict",
            "",
            primary_line,
            "",
            "The full 4M reduced model remains the manuscript-primary latent substrate. The 250k seed models and the 100k lineage model should be used as sensitivity anchors, not as competing primary analyses.",
            "",
            support_line,
            caveat_line,
            "",
            "## Manuscript Rule",
            "",
            "- Use full 4M reduced scVI for the main latent/neighborhood methods.",
            "- Use 250k seed and lineage-focused runs to qualify marker-continuity and lineage-specific interpretation.",
            "- Keep progenitor and immune latent-mechanism claims guarded unless follow-up marker/program or compartment-specific analyses support them.",
            "- Describe the Early iOSN neighborhood as an exact-neighborhood finding supported by edgeR, age-bin directionality, and curated program scoring, not as a fully replicated official-MiloR or all-embedding claim.",
            "",
        ]
    )


def _read_validation(path_like: str | Path) -> pd.DataFrame:
    path = Path(path_like)
    if not path.exists():
        return pd.DataFrame(columns=["check", "status", "detail", "recommendation"])
    return pd.read_csv(path, sep="\t")


def _summarize_model(model: str, table: pd.DataFrame) -> dict[str, object]:
    cells, genes = _shape_from_check(table, "pilot_h5ad")
    _, dims = _shape_from_check(table, "embedding_dimensions")
    marker_statuses = table[table["check"].astype(str).str.startswith(MARKER_CHECK_PREFIX)]["status"].astype(str)
    limited_markers = int(marker_statuses.eq("limited").sum())
    failed_markers = int(marker_statuses.isin(["missing", "failed"]).sum())
    role = _model_role(model)
    min_core_mixing = _min_metric_value(table, CORE_MIXING_METRICS, "normalized_entropy")
    claim_gate = _model_claim_gate(model, dims, failed_markers, min_core_mixing)
    return {
        "model": model,
        "role": role,
        "cells": cells,
        "genes": genes,
        "latent_dimensions": dims,
        "fine_label_purity": _metric_value(table, "neighbor_label_purity__fine_celltype", "mean_same_label"),
        "coarse_label_purity": _metric_value(table, "neighbor_label_purity__coarse_celltype", "mean_same_label"),
        "min_core_mixing_entropy": min_core_mixing,
        "ok_marker_panels": int(marker_statuses.eq("ok").sum()),
        "limited_marker_panels": limited_markers,
        "failed_marker_panels": failed_markers,
        "claim_gate": claim_gate,
    }


def _summarize_marker(marker: str, validations: dict[str, pd.DataFrame]) -> dict[str, object]:
    row: dict[str, object] = {"marker": marker}
    primary_label = ""
    primary_status = "missing"
    primary_enrichment = np.nan
    labels = []
    statuses = []
    for model, table in validations.items():
        check = f"{MARKER_CHECK_PREFIX}{marker}"
        detail = _detail_for_check(table, check)
        status = _status_for_check(table, check)
        label = _extract_label(detail)
        enrichment = _extract_float(detail, "top_decile_enrichment")
        row[f"{model}__status"] = status
        row[f"{model}__top_label"] = label
        row[f"{model}__enrichment"] = enrichment
        if model == PRIMARY_MODEL:
            primary_label = label
            primary_status = status
            primary_enrichment = enrichment
        if label:
            labels.append(label)
        statuses.append(status)
    row["primary_status"] = primary_status
    row["primary_top_label"] = primary_label
    row["primary_enrichment"] = primary_enrichment
    row["label_concordance"] = _label_concordance(marker, labels)
    row["claim_gate"] = _marker_claim_gate(marker, primary_status, primary_label, labels, statuses)
    row["interpretation"] = _marker_interpretation(marker, str(row["claim_gate"]))
    return row


def _model_role(model: str) -> str:
    if model == PRIMARY_MODEL:
        return "primary_latent_substrate"
    if "250k" in model:
        return "scaled_seed_sensitivity"
    if "lineage" in model:
        return "lineage_sensitivity"
    return "pilot_or_other"


def _model_claim_gate(model: str, dims: float, failed_markers: int, min_core_mixing: float) -> str:
    if not np.isfinite(dims) or dims < 10 or failed_markers:
        return "blocked"
    if model == PRIMARY_MODEL:
        if np.isfinite(min_core_mixing) and min_core_mixing < 0.4:
            return "primary_with_technical_caveat"
        return "primary"
    return "sensitivity"


def _marker_claim_gate(marker: str, primary_status: str, primary_label: str, labels: list[str], statuses: list[str]) -> str:
    if primary_status != "ok":
        return "guarded"
    if not _label_concordance(marker, labels):
        return "guarded"
    if marker == "immature_osn" and not _label_matches(primary_label, ("iosn", "immature", "inp")):
        return "guarded"
    if statuses.count("ok") > 1:
        return "supported"
    return "guarded"


def _marker_interpretation(marker: str, claim_gate: str) -> str:
    if claim_gate == "supported":
        return "Marker continuity supports manuscript-level latent interpretation across primary and sensitivity runs."
    if marker == "immature_osn":
        return "Continuity is biologically related but not label-identical across runs; keep Early iOSN wording narrow."
    if marker in {"progenitor", "immune"}:
        return "Continuity is limited in the full 4M model; do not promote this as a latent-mechanism claim without extra support."
    return "Use as supportive context rather than a standalone latent claim."


def _label_concordance(marker: str, labels: list[str]) -> bool:
    expected = EXPECTED_MARKER_LABELS.get(marker, ())
    if not labels or not expected:
        return False
    return all(_label_matches(label, expected) for label in labels)


def _label_matches(label: str, expected: tuple[str, ...]) -> bool:
    normalized = label.lower().replace("_", " ")
    return any(term in normalized for term in expected)


def _shape_from_check(table: pd.DataFrame, check: str) -> tuple[float, float]:
    detail = _detail_for_check(table, check)
    h5ad_match = re.search(r"(\d+)\s+cells\s+x\s+(\d+)", detail)
    if h5ad_match:
        return float(h5ad_match.group(1)), float(h5ad_match.group(2))
    embedding_match = re.search(r"\((\d+),\s*(\d+)\)", detail)
    if embedding_match:
        return float(embedding_match.group(1)), float(embedding_match.group(2))
    return np.nan, np.nan


def _metric_value(table: pd.DataFrame, check: str, key: str) -> float:
    return _extract_float(_detail_for_check(table, check), key)


def _min_metric_value(table: pd.DataFrame, checks: tuple[str, ...], key: str) -> float:
    values = [_metric_value(table, check, key) for check in checks]
    values = [value for value in values if np.isfinite(value)]
    return min(values) if values else np.nan


def _detail_for_check(table: pd.DataFrame, check: str) -> str:
    if table.empty or "check" not in table:
        return ""
    rows = table[table["check"].astype(str).eq(check)]
    return str(rows["detail"].iloc[0]) if not rows.empty else ""


def _status_for_check(table: pd.DataFrame, check: str) -> str:
    if table.empty or "check" not in table:
        return "missing"
    rows = table[table["check"].astype(str).eq(check)]
    return str(rows["status"].iloc[0]) if not rows.empty else "missing"


def _extract_float(detail: str, key: str) -> float:
    match = re.search(rf"{re.escape(key)}=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", detail)
    return float(match.group(1)) if match else np.nan


def _extract_label(detail: str) -> str:
    match = re.search(r"top_label=([^;]+)", detail)
    return match.group(1) if match else ""
