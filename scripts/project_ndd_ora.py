#!/usr/bin/env python3
"""Project healthy-trained ORA models onto all donors, including AD/PD cohorts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.age_model import project_ora_models
from ora.config import load_config
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--features", default=None)
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--scores-out", default=None)
    parser.add_argument("--summary-out", default=None)
    args = parser.parse_args()

    gateway_config = load_config(args.gateway_config)
    model_config = load_config(args.model_config)
    outputs = gateway_config.get("outputs", {})
    features_path = args.features or outputs.get(
        "ora_augmented_feature_matrix_tsv",
        "data/processed/ora_augmented_feature_matrix.tsv",
    )
    manifest_path = args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv")
    scores_out = args.scores_out or outputs.get("ndd_ora_projection_tsv", "results/tables/ndd_ora_projection.tsv")
    summary_out = args.summary_out or outputs.get(
        "ndd_ora_projection_summary_tsv",
        "results/tables/ndd_ora_projection_summary.tsv",
    )

    features = pd.read_csv(features_path, sep="\t")
    manifest = pd.read_csv(manifest_path, sep="\t")
    result = project_ora_models(features, manifest, model_config)
    result.predictions.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    result.summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    ndd_mask = result.predictions["disease_group"].astype(str).isin(["ad", "pd"])
    n_ndd = int(result.predictions[ndd_mask]["donor_id"].nunique())
    print(f"Wrote frozen ORA projections: {scores_out} ({result.predictions['donor_id'].nunique()} donors)")
    print(f"Wrote projection summary: {summary_out} ({result.summary.shape[0]} rows)")
    print(f"Projected NDD donors: {n_ndd}")


if __name__ == "__main__":
    main()
