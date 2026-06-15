#!/usr/bin/env python3
"""Run leakage-safe donor-level stacked ORA age models."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.stacking import DEFAULT_BASE_MODELS, run_stacked_ora
from ora.utils import ensure_parent


def _require_native_backends(model_names: list[str]) -> None:
    required = {"xgboost": "xgboost", "lightgbm": "lightgbm", "catboost": "catboost"}
    missing = []
    for model_name in model_names:
        package = required.get(model_name)
        if package is None:
            continue
        try:
            __import__(package)
        except Exception as exc:  # pragma: no cover - depends on local native libraries
            missing.append(f"{model_name}: {exc}")
    if missing:
        joined = "\n".join(missing)
        raise SystemExit(
            "Native booster backend check failed. Fix optional runtime dependencies or pass "
            f"--allow-fallback to use sklearn fallback models.\n{joined}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", default="data/processed/ora_augmented_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--base-models", nargs="+", default=DEFAULT_BASE_MODELS)
    parser.add_argument("--allow-fallback", action="store_true")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--inner-folds", type=int, default=3)
    parser.add_argument("--performance-out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--scores-out", default=None)
    parser.add_argument("--weights-out", default=None)
    args = parser.parse_args()

    if not args.allow_fallback:
        _require_native_backends(args.base_models)

    model_config = load_config(args.model_config)
    gateway_config = load_config(args.gateway_config)
    outputs = gateway_config.get("outputs", {})
    features = pd.read_csv(args.features, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    result = run_stacked_ora(
        features,
        manifest,
        model_config,
        base_models=args.base_models,
        repeats=args.repeats,
        inner_folds=args.inner_folds,
    )

    performance_out = args.performance_out or outputs.get(
        "ora_stacking_performance_tsv",
        "results/tables/ora_stacking_performance.tsv",
    )
    summary_out = args.summary_out or outputs.get(
        "ora_stacking_summary_tsv",
        "results/tables/ora_stacking_summary.tsv",
    )
    scores_out = args.scores_out or outputs.get("ora_stacking_scores_tsv", "results/tables/ora_stacking_scores.tsv")
    weights_out = args.weights_out or outputs.get("ora_stacking_weights_tsv", "results/tables/ora_stacking_weights.tsv")
    result.performance.to_csv(ensure_parent(performance_out), sep="\t", index=False)
    result.performance_summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.predictions.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    result.meta_weights.to_csv(ensure_parent(weights_out), sep="\t", index=False)
    print(f"Wrote stacking performance: {performance_out} ({result.performance.shape[0]} rows)")
    print(f"Wrote stacking summary: {summary_out} ({result.performance_summary.shape[0]} rows)")
    print(f"Wrote stacking scores: {scores_out} ({result.predictions.shape[0]} rows)")
    print(f"Wrote stacking weights: {weights_out} ({result.meta_weights.shape[0]} rows)")


if __name__ == "__main__":
    main()
