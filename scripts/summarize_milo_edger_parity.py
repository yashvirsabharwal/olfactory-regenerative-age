#!/usr/bin/env python3
"""Summarize concordance between Python Milo-style DA and edgeR count DA."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.neighborhood_parity import summarize_edger_parity
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python-da", required=True)
    parser.add_argument("--edger-da", required=True)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--comparison-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--top-n", type=int, default=100)
    args = parser.parse_args()

    python_da = pd.read_csv(args.python_da, sep="\t")
    edger_da = pd.read_csv(args.edger_da, sep="\t")
    comparison, summary = summarize_edger_parity(python_da, edger_da, run_name=args.run_name, top_n=args.top_n)
    comparison.to_csv(ensure_parent(args.comparison_out), sep="\t", index=False)
    summary.to_csv(ensure_parent(args.summary_out), sep="\t", index=False)
    print(f"Wrote edgeR parity comparison: {args.comparison_out} ({comparison.shape[0]} rows)")
    print(f"Wrote edgeR parity summary: {args.summary_out} ({summary.shape[0]} rows)")


if __name__ == "__main__":
    main()
