"""Organoid perturbation module scoring adapters."""

from __future__ import annotations

import gzip
import re
import tarfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .modules import GeneSet, parse_gene_sets
from .utils import ensure_parent


def score_gse309325_organoid_modules(
    archive_path: str | Path,
    gene_set_config: dict[str, Any],
    *,
    chunksize: int = 500,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score ORA/regeneration modules on GSE309325 gene-by-cell CSV samples."""

    gene_sets = parse_gene_sets(gene_set_config)
    if not gene_sets:
        raise ValueError("No gene sets were configured.")
    archive = Path(archive_path)
    if not archive.exists():
        raise FileNotFoundError(f"Missing GSE309325 archive: {archive}")
    qc_rows = []
    score_rows = []
    coverage_rows = []
    with tarfile.open(archive) as tar:
        members = [member for member in tar.getmembers() if member.isfile() and member.name.endswith(".csv.gz")]
        for member in sorted(members, key=lambda item: item.name):
            sample = parse_gse309325_sample_name(member.name)
            extracted = tar.extractfile(member)
            if extracted is None:
                continue
            with gzip.GzipFile(fileobj=extracted) as handle:
                gene_means, n_cells, n_genes = _selected_gene_means(handle, gene_sets, chunksize=chunksize)
            qc_rows.append(
                {
                    **sample,
                    "n_cells": n_cells,
                    "n_genes": n_genes,
                    "n_selected_genes_present": len(gene_means),
                    "archive_member": member.name,
                }
            )
            for gene_set in gene_sets:
                present = [gene for gene in gene_set.genes if gene.upper() in gene_means]
                missing = [gene for gene in gene_set.genes if gene.upper() not in gene_means]
                module_score = float(np.mean([gene_means[gene.upper()] for gene in present])) if present else np.nan
                score_rows.append(
                    {
                        **sample,
                        "module": gene_set.name,
                        "n_requested": len(gene_set.genes),
                        "n_present": len(present),
                        "coverage_fraction": len(present) / len(gene_set.genes) if gene_set.genes else np.nan,
                        "mean_log1p_expression": module_score,
                    }
                )
                coverage_rows.append(
                    {
                        "sample_id": sample["sample_id"],
                        "module": gene_set.name,
                        "n_requested": len(gene_set.genes),
                        "n_present": len(present),
                        "coverage_fraction": len(present) / len(gene_set.genes) if gene_set.genes else np.nan,
                        "present_genes": ",".join(present),
                        "missing_genes": ",".join(missing),
                    }
                )
    qc = pd.DataFrame(qc_rows)
    scores = pd.DataFrame(score_rows)
    coverage = pd.DataFrame(coverage_rows)
    contrasts = gse309325_organoid_module_contrasts(scores)
    return qc, scores, coverage, contrasts


def parse_gse309325_sample_name(member_name: str) -> dict[str, object]:
    """Parse GSE309325 sample metadata from GEO supplementary member names."""

    name = Path(member_name).name.replace(".csv.gz", "")
    sample_id = name.split("_", 1)[0]
    label = name.split("_", 1)[1] if "_" in name else name
    if "SARSCoV2" in label:
        match = re.search(r"day(\d+)", label)
        day = int(match.group(1)) if match else pd.NA
        condition = "sars_cov_2"
        replicate = pd.NA
    else:
        match = re.search(r"_(\d+)$", label)
        day = 0
        condition = "mock"
        replicate = int(match.group(1)) if match else pd.NA
    return {
        "sample_id": sample_id,
        "sample_label": label,
        "condition": condition,
        "timepoint_day": day,
        "replicate": replicate,
    }


def gse309325_organoid_module_contrasts(scores: pd.DataFrame) -> pd.DataFrame:
    """Compare each infected timepoint to the mock baseline mean."""

    if scores.empty:
        return pd.DataFrame()
    mock = scores[scores["condition"].eq("mock")].copy()
    baseline = (
        mock.groupby("module")["mean_log1p_expression"]
        .agg(mock_mean="mean", mock_sd="std", n_mock="count")
        .reset_index()
    )
    infected = scores[~scores["condition"].eq("mock")].copy()
    merged = infected.merge(baseline, on="module", how="left")
    merged["delta_vs_mock"] = merged["mean_log1p_expression"] - merged["mock_mean"]
    merged["mock_sd"] = merged["mock_sd"].replace(0, np.nan)
    merged["z_vs_mock"] = merged["delta_vs_mock"] / merged["mock_sd"]
    merged["direction_vs_mock"] = np.where(
        merged["delta_vs_mock"] > 0,
        "increased",
        np.where(merged["delta_vs_mock"] < 0, "decreased", "flat"),
    )
    merged["evidence_status"] = "single_timepoint_no_replicate_p_value"
    return merged[
        [
            "sample_id",
            "sample_label",
            "condition",
            "timepoint_day",
            "module",
            "mean_log1p_expression",
            "mock_mean",
            "mock_sd",
            "n_mock",
            "delta_vs_mock",
            "z_vs_mock",
            "direction_vs_mock",
            "evidence_status",
        ]
    ].sort_values(["module", "timepoint_day"]).reset_index(drop=True)


def render_gse309325_organoid_status(
    *,
    qc: pd.DataFrame,
    scores: pd.DataFrame,
    coverage: pd.DataFrame,
    contrasts: pd.DataFrame,
) -> str:
    """Render a concise adapter status report."""

    top = contrasts.assign(abs_delta=contrasts["delta_vs_mock"].abs()).sort_values(
        ["abs_delta", "module"],
        ascending=[False, True],
    ).head(12)
    lines = [
        "# GSE309325 Organoid Perturbation Adapter Status",
        "",
        "Updated: 2026-06-25",
        "",
        "## Scope",
        "",
        "GSE309325 is a human nasal organoid scRNA-seq perturbation dataset containing mock samples and SARS-CoV-2 timepoints. It is olfactory-relevant because the organoid model includes olfactory and nasal respiratory epithelial compartments, but it is infection context, not aging or adult donor validation.",
        "",
        "## Adapter Outputs",
        "",
        f"- Samples scored: {qc.shape[0]}",
        f"- Modules scored per sample: {scores['module'].nunique() if not scores.empty else 0}",
        f"- Perturbation contrasts: {contrasts.shape[0]}",
        "",
        "## Sample QC",
        "",
        "| Sample | Condition | Day | Cells | Selected genes present |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for _, row in qc.iterrows():
        lines.append(
            f"| {row.get('sample_id', '')} | {row.get('condition', '')} | {row.get('timepoint_day', '')} | "
            f"{int(row.get('n_cells', 0))} | {int(row.get('n_selected_genes_present', 0))} |"
        )
    lines.extend(
        [
            "",
            "## Largest Module Shifts Versus Mock",
            "",
            "| Sample | Day | Module | Direction | Delta vs mock | z vs mock |",
            "| --- | ---: | --- | --- | ---: | ---: |",
        ]
    )
    for _, row in top.iterrows():
        lines.append(
            "| "
            f"{row.get('sample_id', '')} | {row.get('timepoint_day', '')} | {row.get('module', '')} | "
            f"{row.get('direction_vs_mock', '')} | {row.get('delta_vs_mock', np.nan):.4g} | "
            f"{row.get('z_vs_mock', np.nan):.4g} |"
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "Use this as perturbation-context evidence only. Infected timepoints are single samples without matched biological replicates, so effect sizes are descriptive module shifts rather than inferential p-values.",
            "",
        ]
    )
    return "\n".join(lines)


def write_gse309325_organoid_outputs(
    *,
    qc: pd.DataFrame,
    scores: pd.DataFrame,
    coverage: pd.DataFrame,
    contrasts: pd.DataFrame,
    report_md: str,
    qc_out: str | Path,
    scores_out: str | Path,
    coverage_out: str | Path,
    contrasts_out: str | Path,
    report_out: str | Path,
) -> None:
    """Write GSE309325 adapter outputs."""

    qc.to_csv(ensure_parent(qc_out), sep="\t", index=False)
    scores.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    coverage.to_csv(ensure_parent(coverage_out), sep="\t", index=False)
    contrasts.to_csv(ensure_parent(contrasts_out), sep="\t", index=False)
    ensure_parent(report_out).write_text(report_md, encoding="utf-8")


def _selected_gene_means(handle: Any, gene_sets: list[GeneSet], *, chunksize: int) -> tuple[dict[str, float], int, int]:
    selected = {gene.upper() for gene_set in gene_sets for gene in gene_set.genes}
    gene_means: dict[str, float] = {}
    n_cells = 0
    n_genes = 0
    for chunk in pd.read_csv(handle, index_col=0, chunksize=chunksize):
        if n_cells == 0:
            n_cells = int(chunk.shape[1])
        n_genes += int(chunk.shape[0])
        upper_index = pd.Index(str(item).upper() for item in chunk.index)
        mask = upper_index.isin(selected)
        if not mask.any():
            continue
        selected_chunk = chunk.loc[mask].copy()
        selected_chunk.index = upper_index[mask]
        log_values = np.log1p(selected_chunk.apply(pd.to_numeric, errors="coerce").fillna(0.0))
        for gene, value in log_values.mean(axis=1).items():
            gene_means[str(gene)] = float(value)
    return gene_means, n_cells, n_genes
