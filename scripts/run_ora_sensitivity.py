#!/usr/bin/env python3
"""Run ORA model sensitivity scenarios by chemistry, collection method, and cell yield."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.sensitivity import run_ora_sensitivity
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", default="data/processed/ora_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--min-cell-thresholds", nargs="*", type=int, default=[500, 1000, 5000, 10000])
    parser.add_argument("--min-train-donors", type=int, default=20)
    parser.add_argument("--scenarios-out", default=None)
    parser.add_argument("--performance-out", default=None)
    parser.add_argument("--scores-out", default=None)
    args = parser.parse_args()

    model_config = load_config(args.model_config)
    gateway_config = load_config(args.gateway_config)
    outputs = gateway_config.get("outputs", {})
    features = pd.read_csv(args.features, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    result = run_ora_sensitivity(
        features,
        manifest,
        model_config,
        min_cell_thresholds=args.min_cell_thresholds,
        min_train_donors=args.min_train_donors,
    )
    scenarios_out = args.scenarios_out or outputs.get("ora_sensitivity_scenarios_tsv", "results/tables/ora_sensitivity_scenarios.tsv")
    performance_out = args.performance_out or outputs.get("ora_sensitivity_performance_tsv", "results/tables/ora_sensitivity_performance.tsv")
    scores_out = args.scores_out or outputs.get("ora_sensitivity_scores_tsv", "results/tables/ora_sensitivity_scores.tsv")
    result.scenarios.to_csv(ensure_parent(scenarios_out), sep="\t", index=False)
    result.performance.to_csv(ensure_parent(performance_out), sep="\t", index=False)
    result.scores.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    ok = int(result.scenarios["status"].eq("ok").sum()) if "status" in result.scenarios else 0
    print(f"Wrote ORA sensitivity scenarios: {scenarios_out} ({ok} runnable)")
    print(f"Wrote ORA sensitivity performance: {performance_out}")
    print(f"Wrote ORA sensitivity scores: {scores_out}")


if __name__ == "__main__":
    main()
