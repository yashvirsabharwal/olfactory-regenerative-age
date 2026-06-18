#!/usr/bin/env python3
"""Validate an external donor-feature matrix against the ORA feature contract."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.external import validate_external_feature_matrix
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--feature-matrix", default=None)
    parser.add_argument("--gateway-features", default=None)
    parser.add_argument("--dataset-id", default="external_feature_matrix")
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--harmonization-out", default=None)
    args = parser.parse_args()

    external_config = load_config(args.external_config)
    gateway_config = load_config(args.gateway_config)
    outputs = gateway_config.get("outputs", {})
    gateway_features_path = args.gateway_features or outputs.get(
        "ora_augmented_feature_matrix_tsv",
        "data/processed/ora_augmented_feature_matrix.tsv",
    )
    gateway_features = pd.read_csv(gateway_features_path, sep="\t") if Path(gateway_features_path).exists() else None
    summary, harmonization = validate_external_feature_matrix(
        args.feature_matrix,
        external_config,
        gateway_features=gateway_features,
        dataset_id=args.dataset_id,
    )
    summary_out = args.summary_out or outputs.get(
        "external_feature_validation_tsv",
        "results/tables/external_feature_validation.tsv",
    )
    harmonization_out = args.harmonization_out or outputs.get(
        "external_feature_harmonization_tsv",
        "results/tables/external_feature_harmonization.tsv",
    )
    summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    harmonization.to_csv(ensure_parent(harmonization_out), sep="\t", index=False)
    status = summary.loc[0, "status"] if not summary.empty else "empty"
    print(f"Wrote external feature validation summary: {summary_out} ({status})")
    print(f"Wrote external feature harmonization: {harmonization_out} ({harmonization.shape[0]} rows)")


if __name__ == "__main__":
    main()
