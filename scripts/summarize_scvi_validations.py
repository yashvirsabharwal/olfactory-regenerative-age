#!/usr/bin/env python3
"""Summarize scVI validation outputs across latent models."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.latent_recompute import summarize_scvi_validation_tables
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--validation",
        action="append",
        nargs=2,
        metavar=("MODEL", "TSV"),
        required=True,
        help="Model label and validation TSV path. Can be repeated.",
    )
    parser.add_argument("--out", default="results/tables/scvi_latent_validation_comparison.tsv")
    args = parser.parse_args()

    validation_paths = {model: path for model, path in args.validation}
    summary = summarize_scvi_validation_tables(validation_paths)
    summary.to_csv(ensure_parent(args.out), sep="\t", index=False)
    print(f"Wrote scVI validation comparison: {args.out} ({summary.shape[0]} rows)")


if __name__ == "__main__":
    main()
