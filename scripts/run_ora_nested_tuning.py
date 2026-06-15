#!/usr/bin/env python3
"""Run leakage-safe nested hyperparameter tuning for selected ORA models."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.tuning import run_nested_tuning
from ora.utils import ensure_parent


def _require_native_backends(model_names: list[str]) -> None:
    required = {
        "xgboost": "xgboost",
        "lightgbm": "lightgbm",
        "catboost": "catboost",
    }
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
            "Native booster backend check failed. Fix the optional runtime dependencies or pass "
            f"--allow-fallback to use sklearn fallback models.\n{joined}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", default="data/processed/ora_augmented_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--models", nargs="+", default=["xgboost", "catboost"])
    parser.add_argument("--allow-fallback", action="store_true")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--inner-folds", type=int, default=3)
    parser.add_argument("--max-candidates", type=int, default=None)
    parser.add_argument("--performance-out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--scores-out", default=None)
    parser.add_argument("--trace-out", default=None)
    parser.add_argument("--selected-params-out", default=None)
    args = parser.parse_args()

    if not args.allow_fallback:
        _require_native_backends(args.models)

    model_config = load_config(args.model_config)
    model_config["model_names"] = args.models
    tuning_config = load_config(args.tuning_config) if args.tuning_config else {}
    gateway_config = load_config(args.gateway_config)
    outputs = gateway_config.get("outputs", {})
    features = pd.read_csv(args.features, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    result = run_nested_tuning(
        features,
        manifest,
        model_config,
        tuning_config,
        repeats=args.repeats,
        inner_folds=args.inner_folds,
        max_candidates=args.max_candidates,
    )

    performance_out = args.performance_out or outputs.get(
        "ora_nested_tuning_performance_tsv",
        "results/tables/ora_nested_tuning_performance.tsv",
    )
    summary_out = args.summary_out or outputs.get(
        "ora_nested_tuning_summary_tsv",
        "results/tables/ora_nested_tuning_summary.tsv",
    )
    scores_out = args.scores_out or outputs.get(
        "ora_nested_tuning_scores_tsv",
        "results/tables/ora_nested_tuning_scores.tsv",
    )
    trace_out = args.trace_out or outputs.get(
        "ora_nested_tuning_trace_tsv",
        "results/tables/ora_nested_tuning_trace.tsv",
    )
    selected_params_out = args.selected_params_out or outputs.get(
        "ora_nested_tuning_selected_params_tsv",
        "results/tables/ora_nested_tuning_selected_params.tsv",
    )
    result.performance.to_csv(ensure_parent(performance_out), sep="\t", index=False)
    result.performance_summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.predictions.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    result.tuning_trace.to_csv(ensure_parent(trace_out), sep="\t", index=False)
    result.selected_params.to_csv(ensure_parent(selected_params_out), sep="\t", index=False)
    print(f"Wrote nested tuning performance: {performance_out} ({result.performance.shape[0]} rows)")
    print(f"Wrote nested tuning summary: {summary_out} ({result.performance_summary.shape[0]} rows)")
    print(f"Wrote nested tuning scores: {scores_out} ({result.predictions.shape[0]} rows)")
    print(f"Wrote nested tuning trace: {trace_out} ({result.tuning_trace.shape[0]} rows)")
    print(f"Wrote selected params: {selected_params_out} ({result.selected_params.shape[0]} rows)")


if __name__ == "__main__":
    main()
