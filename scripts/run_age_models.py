#!/usr/bin/env python3
"""Train composition-MVP ORA age prediction models."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.age_model import train_ora_models
from ora.config import load_config
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", default="data/processed/ora_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--config", default="configs/models.yaml")
    parser.add_argument("--performance-out", default="results/tables/ora_model_performance.tsv")
    parser.add_argument("--scores-out", default="results/tables/donor_ora_scores.tsv")
    parser.add_argument("--importance-out", default="results/tables/ora_feature_importance.tsv")
    args = parser.parse_args()

    model_config = load_config(args.config)
    features = pd.read_csv(args.features, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    result = train_ora_models(features, manifest, model_config)
    result.performance.to_csv(ensure_parent(args.performance_out), sep="\t", index=False)
    result.predictions.to_csv(ensure_parent(args.scores_out), sep="\t", index=False)
    result.feature_importance.to_csv(ensure_parent(args.importance_out), sep="\t", index=False)
    print(f"Wrote ORA performance: {args.performance_out}")
    print(f"Wrote donor ORA scores: {args.scores_out}")
    print(f"Wrote feature importance: {args.importance_out}")


if __name__ == "__main__":
    main()

