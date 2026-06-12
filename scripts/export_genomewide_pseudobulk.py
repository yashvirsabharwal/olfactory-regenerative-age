#!/usr/bin/env python3
"""Export genome-wide pseudobulk counts for edgeR/limma/DESeq2 workflows."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.pseudobulk import DEFAULT_GROUPBY, export_genomewide_pseudobulk_h5ad
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--groupby", nargs="+", default=list(DEFAULT_GROUPBY))
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--gene-chunk-size", type=int, default=500)
    parser.add_argument("--min-cells-per-group", type=int, default=10)
    parser.add_argument("--min-donors-per-cell-state", type=int, default=3)
    parser.add_argument("--apply-qc", action="store_true")
    parser.add_argument("--counts-out", default=None)
    parser.add_argument("--metadata-out", default=None)
    parser.add_argument("--genes-out", default=None)
    parser.add_argument("--summary-out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    h5ad = args.h5ad or config.get("source", {}).get("h5ad_path", "data/raw/gateway.h5ad")
    counts_out = args.counts_out or outputs.get("pseudobulk_genomewide_counts_tsv", "data/processed/pseudobulk_genomewide_counts.tsv.gz")
    metadata_out = args.metadata_out or outputs.get("pseudobulk_genomewide_metadata_tsv", "data/processed/pseudobulk_genomewide_metadata.tsv")
    genes_out = args.genes_out or outputs.get("pseudobulk_genomewide_genes_tsv", "data/processed/pseudobulk_genomewide_genes.tsv")
    summary_out = args.summary_out or outputs.get("pseudobulk_genomewide_summary_tsv", "results/tables/pseudobulk_genomewide_summary.tsv")

    result = export_genomewide_pseudobulk_h5ad(
        h5ad,
        config,
        counts_out=counts_out,
        metadata_out=metadata_out,
        genes_out=genes_out,
        groupby=args.groupby,
        chunk_size=args.chunk_size,
        gene_chunk_size=args.gene_chunk_size,
        apply_qc=args.apply_qc,
        min_cells_per_group=args.min_cells_per_group,
        min_donors_per_cell_state=args.min_donors_per_cell_state,
    )
    result.summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)

    row = result.summary.iloc[0].to_dict()
    print(f"Wrote genome-wide pseudobulk counts: {counts_out} ({int(row['n_genes'])} genes x {int(row['n_groups_exported'])} groups)")
    print(f"Wrote genome-wide pseudobulk metadata: {metadata_out}")
    print(f"Wrote genome-wide gene table: {genes_out}")
    print(f"Wrote genome-wide export summary: {summary_out}")


if __name__ == "__main__":
    main()
