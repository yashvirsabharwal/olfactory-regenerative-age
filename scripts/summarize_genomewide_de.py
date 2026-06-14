#!/usr/bin/env python3
"""Create compact report-ready summaries from genome-wide pseudobulk DE results."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.genomewide_de import summarize_genomewide_de
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--de", default=None)
    parser.add_argument("--run-summary", default=None)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--top-hits-out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    de_path = args.de or outputs.get("pseudobulk_genomewide_edger_tsv", "results/tables/pseudobulk_genomewide_edger.tsv.gz")
    run_summary_path = args.run_summary or outputs.get(
        "pseudobulk_genomewide_edger_summary_tsv",
        "results/tables/pseudobulk_genomewide_edger_summary.tsv",
    )
    summary_out = args.summary_out or outputs.get(
        "pseudobulk_genomewide_de_summary_tsv",
        "results/tables/pseudobulk_genomewide_de_summary.tsv",
    )
    top_hits_out = args.top_hits_out or outputs.get(
        "pseudobulk_genomewide_de_top_hits_tsv",
        "results/tables/pseudobulk_genomewide_de_top_hits.tsv",
    )

    summary, top_hits = summarize_genomewide_de(
        de_path,
        run_summary_path,
        fdr_threshold=args.fdr_threshold,
        top_n=args.top_n,
    )
    summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    top_hits.to_csv(ensure_parent(top_hits_out), sep="\t", index=False)
    print(f"Wrote genome-wide DE summary: {summary_out}")
    print(f"Wrote genome-wide DE top hits: {top_hits_out}")


if __name__ == "__main__":
    main()
