#!/usr/bin/env python3
"""Score curated gene programs in Milo-style neighborhood memberships."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.neighborhood_programs import score_neighborhood_programs_h5ad
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--memberships", required=True)
    parser.add_argument("--da-table", default=None)
    parser.add_argument("--gene-sets", default="configs/gene_sets.yaml")
    parser.add_argument("--run-name", default="neighborhood_run")
    parser.add_argument("--chunk-neighborhoods", type=int, default=500)
    parser.add_argument("--scores-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--coverage-out", required=True)
    args = parser.parse_args()

    with open(args.gene_sets) as handle:
        gene_set_config = yaml.safe_load(handle) or {}
    memberships = pd.read_csv(args.memberships, sep="\t")
    da_table = pd.read_csv(args.da_table, sep="\t") if args.da_table else None
    scores, summary, coverage = score_neighborhood_programs_h5ad(
        args.h5ad,
        gene_set_config,
        memberships,
        da_table=da_table,
        run_name=args.run_name,
        chunk_neighborhoods=args.chunk_neighborhoods,
    )
    scores.to_csv(ensure_parent(args.scores_out), sep="\t", index=False)
    summary.to_csv(ensure_parent(args.summary_out), sep="\t", index=False)
    coverage.to_csv(ensure_parent(args.coverage_out), sep="\t", index=False)
    print(f"Wrote neighborhood program scores: {args.scores_out} ({scores.shape[0]} rows)")
    print(f"Wrote neighborhood program summary: {args.summary_out} ({summary.shape[0]} rows)")
    print(f"Wrote neighborhood program coverage: {args.coverage_out} ({coverage.shape[0]} rows)")


if __name__ == "__main__":
    main()
