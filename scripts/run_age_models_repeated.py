#!/usr/bin/env python3
"""Run repeated donor-level CV for ORA age prediction models."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.age_model import train_ora_models_repeated
from ora.config import load_config
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", default="data/processed/ora_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--repeats", type=int, default=None)
    parser.add_argument("--repeat-performance-out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--scores-out", default=None)
    parser.add_argument("--feature-stability-out", default=None)
    args = parser.parse_args()

    model_config = load_config(args.model_config)
    gateway_config = load_config(args.gateway_config)
    outputs = gateway_config.get("outputs", {})
    features = pd.read_csv(args.features, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    result = train_ora_models_repeated(features, manifest, model_config, repeats=args.repeats)

    repeat_performance_out = args.repeat_performance_out or outputs.get(
        "ora_repeated_cv_performance_tsv",
        "results/tables/ora_repeated_cv_performance.tsv",
    )
    summary_out = args.summary_out or outputs.get(
        "ora_repeated_cv_summary_tsv",
        "results/tables/ora_repeated_cv_summary.tsv",
    )
    scores_out = args.scores_out or outputs.get(
        "ora_repeated_cv_scores_tsv",
        "results/tables/ora_repeated_cv_scores.tsv",
    )
    feature_stability_out = args.feature_stability_out or outputs.get(
        "ora_repeated_cv_feature_stability_tsv",
        "results/tables/ora_repeated_cv_feature_stability.tsv",
    )
    result.repeat_performance.to_csv(ensure_parent(repeat_performance_out), sep="\t", index=False)
    result.performance_summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.predictions.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    result.feature_stability.to_csv(ensure_parent(feature_stability_out), sep="\t", index=False)
    repeats = int(result.repeat_performance["repeat"].nunique()) if "repeat" in result.repeat_performance else 1
    print(f"Wrote repeated-CV performance: {repeat_performance_out} ({repeats} repeats)")
    print(f"Wrote repeated-CV summary: {summary_out}")
    print(f"Wrote repeated-CV feature stability: {feature_stability_out}")


if __name__ == "__main__":
    main()
