"""Perturbation/organoid/ALI validation planning helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .utils import ensure_parent


def build_perturbation_candidate_matrix(config: dict[str, Any]) -> pd.DataFrame:
    """Flatten configured perturbation candidates into a ranked matrix."""

    rows = []
    for candidate in config.get("candidates", []):
        if not isinstance(candidate, dict):
            continue
        rows.append(
            {
                "accession": candidate.get("accession", ""),
                "source_url": candidate.get("source_url", ""),
                "title": candidate.get("title", ""),
                "tissue_model": candidate.get("tissue_model", ""),
                "assay": candidate.get("assay", ""),
                "perturbations": candidate.get("perturbations", ""),
                "mechanism_match": candidate.get("mechanism_match", ""),
                "data_access": candidate.get("data_access", ""),
                "adapter_decision": candidate.get("adapter_decision", ""),
                "priority": int(candidate.get("priority", 99)),
                "directness": candidate.get("directness", ""),
                "usable_for_direct_ora_mechanism": _direct_usable(candidate),
                "limitations": candidate.get("limitations", ""),
                "recommended_next_step": candidate.get("recommended_next_step", ""),
            }
        )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["priority", "accession"]).reset_index(drop=True)


def build_perturbation_search_log(config: dict[str, Any]) -> pd.DataFrame:
    """Flatten perturbation search records."""

    return pd.DataFrame(
        [
            {
                "search_date": row.get("search_date", ""),
                "database_or_resource": row.get("database_or_resource", ""),
                "query_or_filter": row.get("query_or_filter", ""),
                "result_summary": row.get("result_summary", ""),
            }
            for row in config.get("search_log", [])
        ]
    )


def build_minimum_experiment_table(config: dict[str, Any]) -> pd.DataFrame:
    """Flatten minimum experiment designs."""

    return pd.DataFrame(
        [
            {
                "experiment_id": row.get("experiment_id", ""),
                "model": row.get("model", ""),
                "perturbations": row.get("perturbations", ""),
                "timepoints": row.get("timepoints", ""),
                "readout": row.get("readout", ""),
                "target_n": row.get("target_n", ""),
                "claim_boundary": "Causal support only after donor-level perturbation response is measured; not implied by cross-sectional ORA alone.",
            }
            for row in config.get("minimum_experiment", [])
        ]
    )


def render_perturbation_validation_plan(
    *,
    config: dict[str, Any],
    candidates: pd.DataFrame,
    search_log: pd.DataFrame,
    minimum_experiment: pd.DataFrame,
) -> str:
    """Render a perturbation/organoid/ALI validation plan."""

    high = candidates[candidates["priority"].le(1)] if not candidates.empty else pd.DataFrame()
    direct = candidates[candidates["usable_for_direct_ora_mechanism"].astype(bool)] if not candidates.empty else pd.DataFrame()
    conclusion = config.get("metadata", {}).get("conclusion", "")
    lines = [
        "# Perturbation, Organoid, And ALI Validation Plan",
        "",
        "Updated: 2026-06-25",
        "",
        "## Decision",
        "",
        conclusion or "No direct perturbation decision recorded.",
        "",
        f"High-priority public candidates: {high.shape[0]}. Direct adult olfactory aging perturbation datasets: {direct.shape[0]}.",
        "",
        "## Public Candidate Triage",
        "",
        "| Dataset | Model | Assay | Perturbation | Priority | Decision | Limitation |",
        "| --- | --- | --- | --- | ---: | --- | --- |",
    ]
    for _, row in candidates.iterrows():
        accession = str(row.get("accession", ""))
        url = str(row.get("source_url", ""))
        dataset = f"[{accession}]({url})" if url else accession
        lines.append(
            "| "
            f"{dataset} | {row.get('tissue_model', '')} | {row.get('assay', '')} | "
            f"{row.get('perturbations', '')} | {row.get('priority', '')} | "
            f"{row.get('adapter_decision', '')} | {row.get('limitations', '')} |"
        )
    lines.extend(
        [
            "",
            "## Adapter Order",
            "",
            "1. Audit GSE309325 first because it is the only human organoid candidate explicitly containing olfactory and respiratory epithelial compartments.",
            "2. Build a GSE299529 adapter next if Seurat RDS conversion is practical; it is the best cytokine-driven nasal ALI single-cell perturbation candidate.",
            "3. Use GSE286616 or GSE309353 for bulk IFN/NF-kB/viral-injury module direction only if single-cell adapters are blocked.",
            "4. Keep GSE175541, GSE324335, and GSE271245 as lower-priority context only.",
            "",
            "## Minimum New Experiment",
            "",
            "| Experiment | Model | Perturbations | Timepoints | Readout | Target n |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for _, row in minimum_experiment.iterrows():
        lines.append(
            "| "
            f"{row.get('experiment_id', '')} | {row.get('model', '')} | {row.get('perturbations', '')} | "
            f"{row.get('timepoints', '')} | {row.get('readout', '')} | {row.get('target_n', '')} |"
        )
    lines.extend(
        [
            "",
            "## ORA Scoring Readout",
            "",
            "- Score ORA regeneration modules, respiratory metaplasia, IFN/TNF/IL17, senescence/SASP, oxidative stress, EGFR/AREG, Wnt, and Notch programs per condition.",
            "- For single-cell data, summarize by donor, condition, timepoint, and harmonized epithelial state before testing perturbation effects.",
            "- For bulk data, use module-level contrasts only; do not infer cell-state composition without single-cell or histology support.",
            "- Treat donor or independent organoid line as the unit of inference; cells, ROIs, or technical replicates are nested observations.",
            "",
            "## Claim Boundary",
            "",
            "Public perturbation candidates can strengthen mechanism plausibility, especially inflammation and epithelial remodeling, but they cannot prove that ORA age associations are causal unless donor-level olfactory epithelial perturbation responses reproduce the ORA directions.",
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
    lines.append("")
    return "\n".join(lines)


def write_perturbation_validation_outputs(
    *,
    candidates: pd.DataFrame,
    search_log: pd.DataFrame,
    minimum_experiment: pd.DataFrame,
    plan_md: str,
    candidate_out: str | Path,
    search_log_out: str | Path,
    experiment_out: str | Path,
    plan_out: str | Path,
) -> None:
    """Write perturbation validation design outputs."""

    candidates.to_csv(ensure_parent(candidate_out), sep="\t", index=False)
    search_log.to_csv(ensure_parent(search_log_out), sep="\t", index=False)
    minimum_experiment.to_csv(ensure_parent(experiment_out), sep="\t", index=False)
    ensure_parent(plan_out).write_text(plan_md, encoding="utf-8")


def _direct_usable(candidate: dict[str, Any]) -> bool:
    directness = str(candidate.get("directness", "")).lower()
    perturbations = str(candidate.get("perturbations", "")).lower()
    limitations = str(candidate.get("limitations", "")).lower()
    return (
        "olfactory" in directness
        and "aging" in perturbations
        and "infection" not in limitations
        and "not aging" not in limitations
    )
