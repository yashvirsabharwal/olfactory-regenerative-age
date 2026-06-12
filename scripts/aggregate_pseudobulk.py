#!/usr/bin/env python3
"""Aggregate targeted donor x cell-state pseudobulk counts and run lightweight DE."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.pseudobulk import (
    DEFAULT_CONTRASTS,
    DEFAULT_GROUPBY,
    aggregate_targeted_pseudobulk_h5ad,
    genes_from_gene_sets,
    parse_contrasts,
)
from ora.utils import ensure_parent

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--gene-sets", default="configs/gene_sets.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--genes", nargs="*", default=None)
    parser.add_argument("--groupby", nargs="+", default=list(DEFAULT_GROUPBY))
    parser.add_argument("--contrasts", nargs="+", default=list(DEFAULT_CONTRASTS))
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--min-donors", type=int, default=3)
    parser.add_argument("--apply-qc", action="store_true")
    parser.add_argument("--counts-out", default=None)
    parser.add_argument("--metadata-out", default=None)
    parser.add_argument("--coverage-out", default=None)
    parser.add_argument("--de-out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    gene_set_config = load_config(args.gene_sets)
    outputs = config.get("outputs", {})
    h5ad = args.h5ad or config.get("source", {}).get("h5ad_path", "data/raw/gateway.h5ad")
    genes = args.genes or genes_from_gene_sets(gene_set_config)
    result = aggregate_targeted_pseudobulk_h5ad(
        h5ad,
        config,
        genes,
        groupby=args.groupby,
        chunk_size=args.chunk_size,
        apply_qc=args.apply_qc,
        contrasts=parse_contrasts(args.contrasts),
        min_donors=args.min_donors,
    )

    counts_out = args.counts_out or outputs.get("pseudobulk_counts_tsv", "data/processed/pseudobulk_counts.tsv.gz")
    metadata_out = args.metadata_out or outputs.get("pseudobulk_metadata_tsv", "data/processed/pseudobulk_metadata.tsv")
    coverage_out = args.coverage_out or outputs.get("pseudobulk_gene_coverage_tsv", "results/tables/pseudobulk_gene_coverage.tsv")
    de_out = args.de_out or outputs.get("pseudobulk_de_tsv", "results/tables/pseudobulk_de.tsv")

    result.counts.to_csv(ensure_parent(counts_out), sep="\t", index=False)
    result.metadata.to_csv(ensure_parent(metadata_out), sep="\t", index=False)
    result.coverage.to_csv(ensure_parent(coverage_out), sep="\t", index=False)
    result.de.to_csv(ensure_parent(de_out), sep="\t", index=False)
    ok_tests = int(result.de["status"].eq("ok").sum()) if "status" in result.de else 0
    print(f"Wrote pseudobulk counts: {counts_out} ({result.counts.shape[0]} nonzero rows)")
    print(f"Wrote pseudobulk metadata: {metadata_out} ({result.metadata.shape[0]} groups)")
    print(f"Wrote pseudobulk gene coverage: {coverage_out}")
    print(f"Wrote pseudobulk DE: {de_out} ({ok_tests} ok tests)")


if __name__ == "__main__":
    main()
