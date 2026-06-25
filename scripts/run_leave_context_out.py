#!/usr/bin/env python3
"""Run leave-one-context-out ORA robustness checks."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config  # noqa: E402
from ora.context_robustness import run_leave_context_out  # noqa: E402
from ora.utils import ensure_parent  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", default="data/processed/ora_scvi_hybrid_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--models", nargs="*", default=["hist_gradient_boosting", "xgboost", "catboost", "boosted_ensemble"])
    parser.add_argument("--contexts", nargs="*", default=["site", "chemistry", "collection_method", "sex", "race_ethnicity", "yield_bin"])
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--min-train-donors", type=int, default=40)
    parser.add_argument("--min-test-donors", type=int, default=5)
    parser.add_argument("--performance-out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--scores-out", default=None)
    parser.add_argument("--feature-stability-out", default=None)
    parser.add_argument("--feasibility-out", default=None)
    args = parser.parse_args()

    model_config = load_config(args.model_config)
    if args.models:
        model_config["model_names"] = args.models
    outputs = load_config(args.gateway_config).get("outputs", {})
    result = run_leave_context_out(
        pd.read_csv(args.features, sep="\t"),
        pd.read_csv(args.manifest, sep="\t"),
        model_config,
        contexts=args.contexts,
        repeats=args.repeats,
        min_train_donors=args.min_train_donors,
        min_test_donors=args.min_test_donors,
    )
    performance_out = args.performance_out or outputs.get("ora_leave_context_out_performance_tsv", "results/tables/ora_leave_context_out_performance.tsv")
    summary_out = args.summary_out or outputs.get("ora_leave_context_out_summary_tsv", "results/tables/ora_leave_context_out_summary.tsv")
    scores_out = args.scores_out or outputs.get("ora_leave_context_out_scores_tsv", "results/tables/ora_leave_context_out_scores.tsv")
    feature_stability_out = args.feature_stability_out or outputs.get("ora_leave_context_out_feature_stability_tsv", "results/tables/ora_leave_context_out_feature_stability.tsv")
    feasibility_out = args.feasibility_out or outputs.get("ora_leave_context_out_feasibility_tsv", "results/tables/ora_leave_context_out_feasibility.tsv")
    result.performance.to_csv(ensure_parent(performance_out), sep="\t", index=False)
    result.summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.scores.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    result.feature_stability.to_csv(ensure_parent(feature_stability_out), sep="\t", index=False)
    result.feasibility.to_csv(ensure_parent(feasibility_out), sep="\t", index=False)
    print(f"Wrote leave-context-out performance: {performance_out} ({result.performance.shape[0]} rows)")
    print(f"Wrote leave-context-out summary: {summary_out} ({result.summary.shape[0]} rows)")
    print(f"Wrote leave-context-out feasibility: {feasibility_out} ({result.feasibility.shape[0]} rows)")


if __name__ == "__main__":
    main()
