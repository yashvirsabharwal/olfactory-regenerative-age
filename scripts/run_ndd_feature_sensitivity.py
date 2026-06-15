#!/usr/bin/env python3
"""Run NDD ORA projection sensitivity across composition and augmented feature sets."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.age_model import project_ora_models
from ora.config import load_config
from ora.ndd import compare_ndd_feature_sets, donor_projection_appendix
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--composition-features", default=None)
    parser.add_argument("--augmented-features", default=None)
    parser.add_argument("--composition-projection-out", default=None)
    parser.add_argument("--composition-summary-out", default=None)
    parser.add_argument("--augmented-projection-out", default=None)
    parser.add_argument("--augmented-summary-out", default=None)
    parser.add_argument("--comparison-out", default=None)
    parser.add_argument("--donor-appendix-out", default=None)
    args = parser.parse_args()

    gateway_config = load_config(args.gateway_config)
    model_config = load_config(args.model_config)
    outputs = gateway_config.get("outputs", {})
    manifest_path = args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv")
    composition_features = args.composition_features or outputs.get(
        "ora_feature_matrix_tsv",
        "data/processed/ora_feature_matrix.tsv",
    )
    augmented_features = args.augmented_features or outputs.get(
        "ora_augmented_feature_matrix_tsv",
        "data/processed/ora_augmented_feature_matrix.tsv",
    )
    composition_projection_out = args.composition_projection_out or outputs.get(
        "ndd_ora_projection_composition_tsv",
        "results/tables/ndd_ora_projection_composition.tsv",
    )
    composition_summary_out = args.composition_summary_out or outputs.get(
        "ndd_ora_projection_composition_summary_tsv",
        "results/tables/ndd_ora_projection_composition_summary.tsv",
    )
    augmented_projection_out = args.augmented_projection_out or outputs.get(
        "ndd_ora_projection_tsv",
        "results/tables/ndd_ora_projection.tsv",
    )
    augmented_summary_out = args.augmented_summary_out or outputs.get(
        "ndd_ora_projection_summary_tsv",
        "results/tables/ndd_ora_projection_summary.tsv",
    )
    comparison_out = args.comparison_out or outputs.get(
        "ndd_ora_projection_feature_comparison_tsv",
        "results/tables/ndd_ora_projection_feature_comparison.tsv",
    )
    donor_appendix_out = args.donor_appendix_out or outputs.get(
        "ndd_ora_projection_donor_appendix_tsv",
        "results/tables/ndd_ora_projection_donor_appendix.tsv",
    )

    manifest = pd.read_csv(manifest_path, sep="\t")
    composition = _project_feature_set(composition_features, manifest, model_config, "composition")
    augmented = _project_feature_set(augmented_features, manifest, model_config, "augmented")
    comparison = compare_ndd_feature_sets(
        {
            "composition": composition.predictions,
            "augmented": augmented.predictions,
        }
    )
    appendix = pd.concat(
        [
            donor_projection_appendix(composition.predictions),
            donor_projection_appendix(augmented.predictions),
        ],
        ignore_index=True,
    )

    composition.predictions.to_csv(ensure_parent(composition_projection_out), sep="\t", index=False)
    composition.summary.to_csv(ensure_parent(composition_summary_out), sep="\t", index=False)
    augmented.predictions.to_csv(ensure_parent(augmented_projection_out), sep="\t", index=False)
    augmented.summary.to_csv(ensure_parent(augmented_summary_out), sep="\t", index=False)
    comparison.to_csv(ensure_parent(comparison_out), sep="\t", index=False)
    appendix.to_csv(ensure_parent(donor_appendix_out), sep="\t", index=False)
    print(f"Wrote composition NDD projection: {composition_projection_out} ({composition.predictions.shape[0]} rows)")
    print(f"Wrote composition NDD summary: {composition_summary_out} ({composition.summary.shape[0]} rows)")
    print(f"Wrote augmented NDD projection: {augmented_projection_out} ({augmented.predictions.shape[0]} rows)")
    print(f"Wrote augmented NDD summary: {augmented_summary_out} ({augmented.summary.shape[0]} rows)")
    print(f"Wrote NDD feature-set comparison: {comparison_out} ({comparison.shape[0]} rows)")
    print(f"Wrote NDD donor appendix: {donor_appendix_out} ({appendix.shape[0]} rows)")


def _project_feature_set(feature_path: str, manifest: pd.DataFrame, model_config: dict, feature_set: str):
    features = pd.read_csv(feature_path, sep="\t")
    result = project_ora_models(features, manifest, model_config)
    result.predictions.insert(0, "feature_set", feature_set)
    result.summary.insert(0, "feature_set", feature_set)
    return result


if __name__ == "__main__":
    main()
