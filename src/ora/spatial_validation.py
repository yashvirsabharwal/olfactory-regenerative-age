"""Spatial and histology validation design helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .utils import ensure_parent


SPATIAL_TERMS = ("spatial", "geomx", "merfish", "visium", "xenium", "cosmx")


def build_spatial_candidate_matrix(external_config: dict[str, Any]) -> pd.DataFrame:
    """Build a public-data triage table for spatial/histology validation."""

    rows = [_direct_spatial_sentinel()]
    for candidate in external_config.get("public_data_exhaustion", {}).get("candidates", []):
        if not isinstance(candidate, dict):
            continue
        assay = str(candidate.get("assay", ""))
        tissue = str(candidate.get("tissue", ""))
        species = str(candidate.get("species", ""))
        text = f"{assay} {tissue}".lower()
        if not any(term in text for term in SPATIAL_TERMS):
            continue
        rows.append(
            {
                "dataset_id": _dataset_id(candidate.get("accession_or_dataset", "")),
                "accession_or_dataset": candidate.get("accession_or_dataset", ""),
                "source_url": candidate.get("source_url", ""),
                "tissue": tissue,
                "assay": assay,
                "species": species,
                "donor_or_sample_count": candidate.get("donor_or_sample_count", ""),
                "age_availability": candidate.get("age_availability", ""),
                "counts_availability": candidate.get("counts_availability", ""),
                "labels_availability": candidate.get("labels_availability", ""),
                "spatial_validation_role": _spatial_role(candidate),
                "usable_for_primary_spatial_validation": _primary_usable(candidate),
                "recommended_action": _recommended_action(candidate),
                "limitations": _limitations(candidate),
                "notes": candidate.get("notes", ""),
            }
        )
    columns = [
        "dataset_id",
        "accession_or_dataset",
        "source_url",
        "tissue",
        "assay",
        "species",
        "donor_or_sample_count",
        "age_availability",
        "counts_availability",
        "labels_availability",
        "spatial_validation_role",
        "usable_for_primary_spatial_validation",
        "recommended_action",
        "limitations",
        "notes",
    ]
    return pd.DataFrame(rows, columns=columns)


def build_spatial_search_log(external_config: dict[str, Any]) -> pd.DataFrame:
    """Return the spatial-search subset of public-data search records plus this refresh."""

    rows = []
    for row in external_config.get("public_data_exhaustion", {}).get("search_log", []):
        query = str(row.get("query_or_filter", ""))
        resource = str(row.get("database_or_resource", ""))
        if any(term in query.lower() for term in ["spatial", "geomx", "visium", "xenium", "histology"]):
            rows.append(
                {
                    "search_date": row.get("search_date", ""),
                    "database_or_resource": resource,
                    "query_or_filter": query,
                    "result_summary": "No direct adult human olfactory spatial aging dataset found in project registry.",
                }
            )
    rows.extend(
        [
            {
                "search_date": "2026-06-25",
                "database_or_resource": "NCBI GEO / web search",
                "query_or_filter": "human olfactory epithelium spatial transcriptomics; human nasal GeoMx spatial; GSE235714; GSE292993",
                "result_summary": "Confirmed nasal CRS GeoMx context dataset GSE235714 and airway/lung spatial context dataset GSE292993; no direct adult human olfactory aging spatial dataset found.",
            },
            {
                "search_date": "2026-06-25",
                "database_or_resource": "Local ORA registry",
                "query_or_filter": "public_data_exhaustion spatial/MERFISH/GeoMx candidates",
                "result_summary": "GSE303809 is fetal/developmental olfactory MERFISH context only; adult olfactory spatial validation requires new data or targeted histology.",
            },
        ]
    )
    return pd.DataFrame(rows)


def build_spatial_marker_panel(marker_config: dict[str, Any]) -> pd.DataFrame:
    """Flatten the configured marker panel into a table."""

    rows = []
    for panel in marker_config.get("panels", []):
        markers = [str(marker) for marker in panel.get("markers", []) if str(marker)]
        rows.append(
            {
                "panel_id": panel.get("panel_id", ""),
                "theme": panel.get("theme", ""),
                "compartment": panel.get("compartment", ""),
                "priority": panel.get("priority", ""),
                "expected_age_direction": panel.get("expected_age_direction", ""),
                "markers": ",".join(markers),
                "n_markers": len(markers),
                "readout": panel.get("readout", ""),
                "rationale": panel.get("rationale", ""),
            }
        )
    return pd.DataFrame(rows)


def build_spatial_readout_plan(marker_config: dict[str, Any]) -> pd.DataFrame:
    """Flatten configured spatial readouts into a table."""

    return pd.DataFrame(
        [
            {
                "readout_id": row.get("readout_id", ""),
                "metric": row.get("metric", ""),
                "model": row.get("model", ""),
                "unit_of_inference": "donor",
                "claim_boundary": "Orthogonal localization/support only; not lineage flux or causal mechanism.",
            }
            for row in marker_config.get("readouts", [])
        ]
    )


def render_spatial_validation_plan(
    *,
    candidate_matrix: pd.DataFrame,
    marker_panel: pd.DataFrame,
    readout_plan: pd.DataFrame,
    search_log: pd.DataFrame,
) -> str:
    """Render the spatial/histology validation plan."""

    direct = candidate_matrix[
        candidate_matrix["dataset_id"].eq("direct_adult_human_olfactory_spatial_not_found")
    ]
    direct_status = "not_found" if not direct.empty else "unknown"
    lines = [
        "# Spatial And Histology Validation Plan",
        "",
        "Updated: 2026-06-25",
        "",
        "## Decision",
        "",
        f"Direct public adult human olfactory epithelial spatial aging dataset status: `{direct_status}`.",
        "",
        "The best public spatial resources found are context datasets, not primary ORA validation datasets. GSE235714 is human nasal CRS/healthy-control GeoMx spatial transcriptomics; GSE292993 is human lung airway/parenchymal/vessel spatial context; GSE303809 is fetal olfactory/head-section MERFISH and is developmental only. Therefore the next strong validation step is a targeted adult olfactory histology/spatial experiment.",
        "",
        "## Public Spatial Candidate Triage",
        "",
        "| Dataset | Tissue | Assay | Role | Primary usable | Limitation |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in candidate_matrix.iterrows():
        label = str(row.get("accession_or_dataset", ""))
        url = str(row.get("source_url", ""))
        dataset = f"[{label}]({url})" if url and not label.startswith("No direct") else label
        lines.append(
            "| "
            f"{dataset} | {row.get('tissue', '')} | {row.get('assay', '')} | "
            f"{row.get('spatial_validation_role', '')} | {row.get('usable_for_primary_spatial_validation', '')} | "
            f"{row.get('limitations', '')} |"
        )
    lines.extend(
        [
            "",
            "## Experimental Design",
            "",
            "- Cohort: at least 10 donors per age bin if feasible: 18-40, 41-65, and 66+ years.",
            "- Sampling: olfactory cleft/olfactory epithelium with histologic confirmation; record exact biopsy site, CRS/allergy/COVID/smoking history, sex, race/ethnicity, PMI or processing delay, and section quality.",
            "- Exclusions or strata: active CRS/nasal polyps, acute infection, tumor-adjacent tissue, severe technical artifact, and ambiguous respiratory-only sections.",
            "- Preferred assay: Xenium or CosMx for single-cell spatial resolution with the expanded marker panel. Practical fallback: RNAscope/IF plus a minimal marker set on adjacent sections.",
            "- Unit of inference: donor. Regions, fields, cells, ROIs, or sections are nested observations and should not be treated as independent donors.",
            "",
            "## Marker Panel",
            "",
            "| Panel | Theme | Priority | Expected age direction | Markers | Readout |",
            "| --- | --- | ---: | --- | --- | --- |",
        ]
    )
    for _, row in marker_panel.sort_values(["priority", "panel_id"]).iterrows():
        lines.append(
            "| "
            f"{row.get('panel_id', '')} | {row.get('theme', '')} | {row.get('priority', '')} | "
            f"{row.get('expected_age_direction', '')} | {row.get('markers', '')} | {row.get('readout', '')} |"
        )
    lines.extend(
        [
            "",
            "## Quantitative Readouts",
            "",
            "| Readout | Metric | Model |",
            "| --- | --- | --- |",
        ]
    )
    for _, row in readout_plan.iterrows():
        lines.append(f"| {row.get('readout_id', '')} | {row.get('metric', '')} | {row.get('model', '')} |")
    lines.extend(
        [
            "",
            "## Expected Support Patterns",
            "",
            "- Strong support: older donors show reduced immature/mature OSN spatial density or transduction intensity, altered HBC/INP localization, increased respiratory-metaplasia or stress regions, and immune/stress neighborhoods that align with ORA feature families.",
            "- Partial support: only some ORA themes localize spatially, especially neuronal decline and respiratory metaplasia, while immune/LR signals remain context-dependent.",
            "- Negative result: ORA-associated features do not localize to coherent epithelial or immune regions after donor-level modeling and histologic QC.",
            "",
            "## Search Log",
            "",
            "| Date | Resource | Query/filter | Result |",
            "| --- | --- | --- | --- |",
        ]
    )
    for _, row in search_log.iterrows():
        lines.append(
            "| "
            f"{row.get('search_date', '')} | {row.get('database_or_resource', '')} | "
            f"{row.get('query_or_filter', '')} | {row.get('result_summary', '')} |"
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "This plan can validate localization and compartment-level consistency of ORA biology. It cannot by itself prove lineage flux, regeneration rate, or causality; perturbation or longitudinal injury-repair data would be needed for those claims.",
            "",
        ]
    )
    return "\n".join(lines)


def write_spatial_validation_outputs(
    *,
    candidate_matrix: pd.DataFrame,
    marker_panel: pd.DataFrame,
    readout_plan: pd.DataFrame,
    search_log: pd.DataFrame,
    plan_md: str,
    candidate_out: str | Path,
    marker_out: str | Path,
    readout_out: str | Path,
    search_log_out: str | Path,
    plan_out: str | Path,
) -> None:
    """Write spatial validation design artifacts."""

    candidate_matrix.to_csv(ensure_parent(candidate_out), sep="\t", index=False)
    marker_panel.to_csv(ensure_parent(marker_out), sep="\t", index=False)
    readout_plan.to_csv(ensure_parent(readout_out), sep="\t", index=False)
    search_log.to_csv(ensure_parent(search_log_out), sep="\t", index=False)
    ensure_parent(plan_out).write_text(plan_md, encoding="utf-8")


def _direct_spatial_sentinel() -> dict[str, object]:
    return {
        "dataset_id": "direct_adult_human_olfactory_spatial_not_found",
        "accession_or_dataset": "No direct adult human olfactory epithelial spatial aging dataset found",
        "source_url": "",
        "tissue": "adult human olfactory epithelium",
        "assay": "spatial transcriptomics or histology",
        "species": "human",
        "donor_or_sample_count": "not found",
        "age_availability": "not found",
        "counts_availability": "not found",
        "labels_availability": "not found",
        "spatial_validation_role": "blocking_gap",
        "usable_for_primary_spatial_validation": False,
        "recommended_action": "run targeted adult olfactory spatial/histology validation",
        "limitations": "public data gap",
        "notes": "Context datasets exist, but no public adult olfactory spatial aging dataset was found.",
    }


def _spatial_role(candidate: dict[str, Any]) -> str:
    tissue = str(candidate.get("tissue", "")).lower()
    species = str(candidate.get("species", "")).lower()
    assay = str(candidate.get("assay", "")).lower()
    if "fetal" in species or "fetal" in tissue or "pcw" in str(candidate.get("age_availability", "")).lower():
        return "developmental_spatial_context_only"
    if "olfactory" in tissue:
        return "olfactory_spatial_context"
    if "nasal" in tissue or "sinus" in tissue or "ethmoid" in tissue:
        return "nasal_spatial_context"
    if "lung" in tissue or "airway" in tissue or "vessel" in tissue:
        return "airway_lung_spatial_context"
    if "geomx" in assay:
        return "spatial_context"
    return "spatial_context"


def _primary_usable(candidate: dict[str, Any]) -> bool:
    role = _spatial_role(candidate)
    species = str(candidate.get("species", "")).lower()
    age = str(candidate.get("age_availability", "")).lower()
    return role == "olfactory_spatial_context" and "human" in species and "fetal" not in species and "adult" in age and "aging" in age


def _recommended_action(candidate: dict[str, Any]) -> str:
    role = _spatial_role(candidate)
    if role == "nasal_spatial_context":
        return "Use only as nasal inflammation/metaplasia context; do not treat as adult olfactory aging validation."
    if role == "airway_lung_spatial_context":
        return "Use only for airway/lung spatial specificity context."
    if role == "developmental_spatial_context_only":
        return "Use only for developmental marker localization; exclude from adult aging validation."
    if role == "olfactory_spatial_context":
        return "Audit labels and ages before any primary use."
    return "Context only."


def _limitations(candidate: dict[str, Any]) -> str:
    role = _spatial_role(candidate)
    limits = []
    if role == "developmental_spatial_context_only":
        limits.append("developmental/fetal only")
    if role == "nasal_spatial_context":
        limits.append("nasal/CRS context, not olfactory aging")
    if role == "airway_lung_spatial_context":
        limits.append("airway/lung context, not olfactory tissue")
    if "not obvious" in str(candidate.get("age_availability", "")).lower():
        limits.append("age metadata not obvious")
    if "fastq raw data files are unavailable" in str(candidate.get("notes", "")).lower():
        limits.append("raw FASTQ unavailable")
    if not limits:
        limits.append("manual audit required")
    return "; ".join(limits)


def _dataset_id(accession_or_dataset: object) -> str:
    text = str(accession_or_dataset).strip()
    if not text:
        return ""
    return (
        text.lower()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("+", "plus")
        .replace(".", "")
        .replace("-", "_")
    )
