#!/usr/bin/env python3
"""Summarize genome-wide pseudobulk export QC and top variable genes."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.genomewide_qc import summarize_genomewide_pseudobulk
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--counts", default=None)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--genes", default=None)
    parser.add_argument("--chunksize", type=int, default=500)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--gene-qc-out", default=None)
    parser.add_argument("--group-qc-out", default=None)
    parser.add_argument("--disease-summary-out", default=None)
    parser.add_argument("--cell-state-summary-out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    counts = args.counts or outputs.get("pseudobulk_genomewide_counts_tsv", "data/processed/pseudobulk_genomewide_counts.tsv.gz")
    metadata = args.metadata or outputs.get("pseudobulk_genomewide_metadata_tsv", "data/processed/pseudobulk_genomewide_metadata.tsv")
    genes = args.genes or outputs.get("pseudobulk_genomewide_genes_tsv", "data/processed/pseudobulk_genomewide_genes.tsv")
    summary_out = args.summary_out or outputs.get("pseudobulk_genomewide_qc_summary_tsv", "results/tables/pseudobulk_genomewide_qc_summary.tsv")
    gene_qc_out = args.gene_qc_out or outputs.get("pseudobulk_genomewide_gene_qc_tsv", "results/tables/pseudobulk_genomewide_gene_qc.tsv")
    group_qc_out = args.group_qc_out or outputs.get("pseudobulk_genomewide_group_qc_tsv", "results/tables/pseudobulk_genomewide_group_qc.tsv")
    disease_summary_out = args.disease_summary_out or outputs.get(
        "pseudobulk_genomewide_disease_summary_tsv",
        "results/tables/pseudobulk_genomewide_disease_summary.tsv",
    )
    cell_state_summary_out = args.cell_state_summary_out or outputs.get(
        "pseudobulk_genomewide_cell_state_summary_tsv",
        "results/tables/pseudobulk_genomewide_cell_state_summary.tsv",
    )

    result = summarize_genomewide_pseudobulk(counts, metadata, genes, chunksize=args.chunksize)
    result.summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.gene_qc.to_csv(ensure_parent(gene_qc_out), sep="\t", index=False)
    result.group_qc.to_csv(ensure_parent(group_qc_out), sep="\t", index=False)
    result.disease_summary.to_csv(ensure_parent(disease_summary_out), sep="\t", index=False)
    result.cell_state_summary.to_csv(ensure_parent(cell_state_summary_out), sep="\t", index=False)
    row = result.summary.iloc[0].to_dict()
    print(
        "Wrote genome-wide pseudobulk QC: "
        f"{int(row['n_genes'])} genes x {int(row['n_groups'])} groups; "
        f"metadata match={bool(row['matrix_columns_match_metadata'])}"
    )
    print(f"Wrote gene QC: {gene_qc_out}")
    print(f"Wrote group QC: {group_qc_out}")


if __name__ == "__main__":
    main()
