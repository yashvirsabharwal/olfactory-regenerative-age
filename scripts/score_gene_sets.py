#!/usr/bin/env python3
"""Score curated gene sets/modules in backed H5AD mode."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.modules import DEFAULT_GROUPBY, score_gene_sets_h5ad
from ora.utils import ensure_parent

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--gene-sets", default="configs/gene_sets.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--coverage-out", default=None)
    parser.add_argument("--donor-features-out", default=None)
    parser.add_argument("--chunk-size", type=int, default=None)
    parser.add_argument("--layer", default=None)
    parser.add_argument("--groupby", nargs="+", default=list(DEFAULT_GROUPBY))
    parser.add_argument("--apply-qc", action="store_true")
    parser.add_argument("--no-log1p", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    gene_set_config = load_config(args.gene_sets)
    outputs = config.get("outputs", {})
    h5ad = args.h5ad or config.get("source", {}).get("h5ad_path", "data/raw/gateway.h5ad")
    summary_out = args.summary_out or outputs.get("module_score_summary_tsv", "results/tables/module_score_summary.tsv")
    coverage_out = args.coverage_out or outputs.get("module_gene_coverage_tsv", "results/tables/module_gene_coverage.tsv")
    donor_features_out = args.donor_features_out or outputs.get(
        "donor_module_features_tsv", "data/processed/donor_module_features.tsv"
    )

    result = score_gene_sets_h5ad(
        h5ad,
        config,
        gene_set_config,
        groupby=args.groupby,
        chunk_size=args.chunk_size,
        layer=args.layer,
        log1p=not args.no_log1p,
        apply_qc=args.apply_qc,
    )
    result.summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.coverage.to_csv(ensure_parent(coverage_out), sep="\t", index=False)
    result.donor_features.to_csv(ensure_parent(donor_features_out), sep="\t", index=False)
    present = int((result.coverage["n_present"] > 0).sum()) if "n_present" in result.coverage else 0
    print(f"Wrote module score summary: {summary_out} ({result.summary.shape[0]} rows)")
    print(f"Wrote module gene coverage: {coverage_out} ({present}/{result.coverage.shape[0]} modules with genes)")
    print(f"Wrote donor module features: {donor_features_out} ({result.donor_features.shape[1] - 1} features)")


if __name__ == "__main__":
    main()
