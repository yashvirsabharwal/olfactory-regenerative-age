#!/usr/bin/env python3
"""Annotate full-scale Milo-style neighborhood results with ORA biology themes."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.neighborhood_annotation import build_neighborhood_annotation
from ora.utils import ensure_parent


DEFAULT_TABLES = (
    ("all_full", "results/tables/milo_full_4m_neighborhood_da.tsv"),
    ("lineage_full", "results/tables/milo_full_4m_lineage_neighborhood_da.tsv"),
    ("secretory_full", "results/tables/milo_full_4m_secretory_neighborhood_da.tsv"),
    ("all_matched", "results/tables/milo_full_4m_matched_neighborhood_da.tsv"),
    ("lineage_matched", "results/tables/milo_full_4m_lineage_matched_neighborhood_da.tsv"),
    ("secretory_matched", "results/tables/milo_full_4m_secretory_matched_neighborhood_da.tsv"),
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--table",
        action="append",
        nargs=2,
        metavar=("RUN_NAME", "TSV"),
        help="Run name and Milo-style DA TSV. May be repeated. Defaults to all full 4M tables.",
    )
    parser.add_argument("--fdr-threshold", type=float, default=0.10)
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--top-out", default="results/tables/milo_full_4m_top_neighborhood_themes.tsv")
    parser.add_argument("--summary-out", default="results/tables/milo_full_4m_theme_summary.tsv")
    args = parser.parse_args()

    table_specs = tuple(args.table) if args.table else DEFAULT_TABLES
    tables = {}
    for run_name, path in table_specs:
        candidate = Path(path)
        if not candidate.exists():
            raise SystemExit(f"Missing neighborhood table for `{run_name}`: {path}")
        tables[str(run_name)] = pd.read_csv(candidate, sep="\t")

    top, summary = build_neighborhood_annotation(
        tables,
        fdr_threshold=args.fdr_threshold,
        top_n=args.top_n,
    )
    top.to_csv(ensure_parent(args.top_out), sep="\t", index=False)
    summary.to_csv(ensure_parent(args.summary_out), sep="\t", index=False)
    print(f"Wrote top neighborhood annotations: {args.top_out} ({top.shape[0]} rows)")
    print(f"Wrote neighborhood theme summary: {args.summary_out} ({summary.shape[0]} rows)")


if __name__ == "__main__":
    main()
