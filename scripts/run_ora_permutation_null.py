#!/usr/bin/env python3
"""Run donor-level shuffled-age null tests for ORA models."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.permutation import run_permutation_null
from ora.utils import ensure_parent


DEFAULT_MODELS = ["random_forest", "xgboost", "catboost", "boosted_ensemble"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", default="data/processed/ora_augmented_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--observed-summary", default=None)
    parser.add_argument("--models", nargs="+", default=DEFAULT_MODELS)
    parser.add_argument("--n-permutations", type=int, default=100)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--random-seed", type=int, default=20260615)
    parser.add_argument("--repeat-performance-out", default=None)
    parser.add_argument("--permutation-summary-out", default=None)
    parser.add_argument("--empirical-summary-out", default=None)
    args = parser.parse_args()

    model_config = load_config(args.model_config)
    model_config["model_names"] = args.models
    gateway_config = load_config(args.gateway_config)
    outputs = gateway_config.get("outputs", {})
    features = pd.read_csv(args.features, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    observed_summary = pd.read_csv(args.observed_summary, sep="\t") if args.observed_summary else None
    result = run_permutation_null(
        features,
        manifest,
        model_config,
        n_permutations=args.n_permutations,
        repeats=args.repeats,
        random_seed=args.random_seed,
        observed_summary=observed_summary,
    )

    repeat_performance_out = args.repeat_performance_out or outputs.get(
        "ora_permutation_repeat_performance_tsv",
        "results/tables/ora_permutation_repeat_performance.tsv",
    )
    permutation_summary_out = args.permutation_summary_out or outputs.get(
        "ora_permutation_summary_tsv",
        "results/tables/ora_permutation_summary.tsv",
    )
    empirical_summary_out = args.empirical_summary_out or outputs.get(
        "ora_permutation_empirical_tsv",
        "results/tables/ora_permutation_empirical.tsv",
    )
    result.repeat_performance.to_csv(ensure_parent(repeat_performance_out), sep="\t", index=False)
    result.permutation_summary.to_csv(ensure_parent(permutation_summary_out), sep="\t", index=False)
    result.empirical_summary.to_csv(ensure_parent(empirical_summary_out), sep="\t", index=False)
    print(
        f"Wrote permutation repeat performance: {repeat_performance_out} "
        f"({result.repeat_performance.shape[0]} rows)"
    )
    print(
        f"Wrote permutation summaries: {permutation_summary_out} "
        f"({result.permutation_summary.shape[0]} rows)"
    )
    print(f"Wrote empirical null summary: {empirical_summary_out} ({result.empirical_summary.shape[0]} rows)")


if __name__ == "__main__":
    main()
