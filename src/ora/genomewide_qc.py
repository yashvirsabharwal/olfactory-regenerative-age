"""QC summaries for genome-wide pseudobulk exports."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class GenomewidePseudobulkQC:
    summary: pd.DataFrame
    gene_qc: pd.DataFrame
    group_qc: pd.DataFrame
    disease_summary: pd.DataFrame
    cell_state_summary: pd.DataFrame


def summarize_genomewide_pseudobulk(
    counts_path: str | Path,
    metadata_path: str | Path,
    genes_path: str | Path | None = None,
    *,
    chunksize: int = 500,
) -> GenomewidePseudobulkQC:
    """Summarize a gene-by-pseudobulk count matrix without loading all rows at once."""

    metadata = pd.read_csv(metadata_path, sep="\t")
    genes = pd.read_csv(genes_path, sep="\t") if genes_path else pd.DataFrame()
    with _open_counts(counts_path) as reader:
        header = reader.readline().rstrip("\n").split("\t")
    id_columns = header[:3]
    group_ids = header[3:]
    expected_ids = metadata["pseudobulk_id"].astype(str).tolist() if "pseudobulk_id" in metadata else []
    ids_match = group_ids == expected_ids
    gene_rows = []
    group_totals = np.zeros(len(group_ids), dtype=np.float64)
    group_detected_genes = np.zeros(len(group_ids), dtype=np.int64)
    total_counts = 0.0

    for chunk in pd.read_csv(counts_path, sep="\t", chunksize=chunksize):
        count_values = chunk[group_ids].to_numpy(dtype=np.float64, copy=False)
        gene_total = count_values.sum(axis=1)
        detected_groups = (count_values > 0).sum(axis=1)
        mean_count = count_values.mean(axis=1)
        variance_log1p = np.log1p(count_values).var(axis=1)
        total_counts += float(gene_total.sum())
        group_totals += count_values.sum(axis=0)
        group_detected_genes += (count_values > 0).sum(axis=0)
        frame = chunk[id_columns].copy()
        frame["total_count"] = gene_total
        frame["detected_groups"] = detected_groups
        frame["detected_group_fraction"] = detected_groups / max(len(group_ids), 1)
        frame["mean_count"] = mean_count
        frame["variance_log1p"] = variance_log1p
        gene_rows.append(frame)

    gene_qc = pd.concat(gene_rows, ignore_index=True) if gene_rows else pd.DataFrame()
    if not genes.empty:
        gene_qc = gene_qc.merge(genes, on=[col for col in ["gene_id", "gene_symbol", "gene_index"] if col in gene_qc and col in genes], how="left")
    group_qc = metadata.copy()
    group_qc["matrix_total_count"] = group_totals
    group_qc["detected_genes"] = group_detected_genes
    group_qc["detected_gene_fraction"] = group_detected_genes / max(gene_qc.shape[0], 1)
    if "sum_n_counts" in group_qc:
        denom = pd.to_numeric(group_qc["sum_n_counts"], errors="coerce").replace(0, np.nan)
        group_qc["matrix_to_metadata_count_ratio"] = group_qc["matrix_total_count"] / denom

    disease_summary = _summarize_groups(group_qc, ["disease_group"])
    cell_state_summary = _summarize_groups(group_qc, ["disease_group", "fine_cell_type"])
    summary = pd.DataFrame(
        [
            {
                "n_genes": int(gene_qc.shape[0]),
                "n_groups": int(len(group_ids)),
                "metadata_rows": int(metadata.shape[0]),
                "matrix_columns_match_metadata": bool(ids_match),
                "matrix_total_counts": total_counts,
                "metadata_total_counts": float(pd.to_numeric(metadata.get("sum_n_counts"), errors="coerce").sum())
                if "sum_n_counts" in metadata
                else np.nan,
                "median_group_detected_genes": float(np.median(group_detected_genes)) if len(group_detected_genes) else np.nan,
                "median_gene_detected_group_fraction": float(gene_qc["detected_group_fraction"].median()) if not gene_qc.empty else np.nan,
            }
        ]
    )
    return GenomewidePseudobulkQC(
        summary=summary,
        gene_qc=gene_qc.sort_values(["variance_log1p", "total_count"], ascending=[False, False]).reset_index(drop=True),
        group_qc=group_qc,
        disease_summary=disease_summary,
        cell_state_summary=cell_state_summary,
    )


def _summarize_groups(frame: pd.DataFrame, groupby: list[str]) -> pd.DataFrame:
    if frame.empty or any(col not in frame for col in groupby):
        return pd.DataFrame()
    return (
        frame.groupby(groupby, observed=True, dropna=False)
        .agg(
            groups=("pseudobulk_id", "size"),
            donors=("donor_id", "nunique"),
            cells=("n_cells", "sum"),
            matrix_total_count=("matrix_total_count", "sum"),
            median_detected_genes=("detected_genes", "median"),
        )
        .reset_index()
        .sort_values(["groups", "cells"], ascending=[False, False])
    )


def _open_counts(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        import gzip

        return gzip.open(path, "rt")
    return path.open("rt", encoding="utf-8")
