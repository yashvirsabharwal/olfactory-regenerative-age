#!/usr/bin/env python3
"""Compare composition-only and module-augmented ORA repeated-CV summaries."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.model_compare import compare_feature_set_deltas, rank_feature_set_summaries
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-summary", default="results/tables/ora_repeated_cv_summary.tsv")
    parser.add_argument("--augmented-summary", default="results/tables/ora_augmented_repeated_cv_summary.tsv")
    parser.add_argument("--base-label", default="composition")
    parser.add_argument("--augmented-label", default="composition_plus_modules")
    parser.add_argument("--combined-out", default="results/tables/ora_feature_set_model_comparison.tsv")
    parser.add_argument("--delta-out", default="results/tables/ora_feature_set_model_deltas.tsv")
    args = parser.parse_args()

    base = pd.read_csv(args.base_summary, sep="\t")
    augmented = pd.read_csv(args.augmented_summary, sep="\t")
    combined = rank_feature_set_summaries(
        {
            args.base_label: base,
            args.augmented_label: augmented,
        }
    )
    deltas = compare_feature_set_deltas(
        base,
        augmented,
        base_label=args.base_label,
        augmented_label=args.augmented_label,
    )
    combined.to_csv(ensure_parent(args.combined_out), sep="\t", index=False)
    deltas.to_csv(ensure_parent(args.delta_out), sep="\t", index=False)
    best = combined.iloc[0]
    print(
        "Best repeated-CV feature/model pair: "
        f"{best['feature_set']} / {best['model']} "
        f"(MAE {best['mae_mean']:.3f})"
    )
    print(f"Wrote combined model comparison: {args.combined_out} ({combined.shape[0]} rows)")
    print(f"Wrote feature-set deltas: {args.delta_out} ({deltas.shape[0]} rows)")


if __name__ == "__main__":
    main()
