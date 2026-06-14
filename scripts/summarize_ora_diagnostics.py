#!/usr/bin/env python3
"""Summarize ORA calibration and residual diagnostics."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.diagnostics import summarize_ora_diagnostics
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", default=None)
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--calibration-out", default=None)
    parser.add_argument("--age-bin-out", default=None)
    parser.add_argument("--residuals-out", default=None)
    parser.add_argument("--calibrated-scores-out", default=None)
    args = parser.parse_args()

    model_config = load_config(args.model_config)
    gateway_config = load_config(args.gateway_config)
    outputs = gateway_config.get("outputs", {})
    scores_path = args.scores or outputs.get("donor_ora_scores_tsv", "results/tables/donor_ora_scores.tsv")
    manifest_path = args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv")
    calibration_out = args.calibration_out or outputs.get("ora_calibration_tsv", "results/tables/ora_calibration.tsv")
    age_bin_out = args.age_bin_out or outputs.get("ora_age_bin_errors_tsv", "results/tables/ora_age_bin_errors.tsv")
    residuals_out = args.residuals_out or outputs.get("ora_residual_diagnostics_tsv", "results/tables/ora_residual_diagnostics.tsv")
    calibrated_scores_out = args.calibrated_scores_out or outputs.get(
        "ora_calibrated_scores_tsv",
        "results/tables/ora_calibrated_scores.tsv",
    )

    scores = pd.read_csv(scores_path, sep="\t")
    manifest = pd.read_csv(manifest_path, sep="\t")
    result = summarize_ora_diagnostics(scores, model_config=model_config, manifest=manifest)
    result.calibration.to_csv(ensure_parent(calibration_out), sep="\t", index=False)
    result.age_bin_errors.to_csv(ensure_parent(age_bin_out), sep="\t", index=False)
    result.residual_diagnostics.to_csv(ensure_parent(residuals_out), sep="\t", index=False)
    result.calibrated_scores.to_csv(ensure_parent(calibrated_scores_out), sep="\t", index=False)
    print(f"Wrote ORA calibration: {calibration_out} ({result.calibration.shape[0]} models)")
    print(f"Wrote ORA age-bin errors: {age_bin_out} ({result.age_bin_errors.shape[0]} rows)")
    print(f"Wrote ORA residual diagnostics: {residuals_out} ({result.residual_diagnostics.shape[0]} rows)")
    print(f"Wrote ORA calibrated scores: {calibrated_scores_out} ({result.calibrated_scores.shape[0]} rows)")


if __name__ == "__main__":
    main()
