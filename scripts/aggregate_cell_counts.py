#!/usr/bin/env python3
"""Aggregate Gateway cell-state counts and donor-level composition features."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.aggregate import aggregate_cell_counts, build_cell_state_features
from ora.config import load_config, project_path
from ora.io import load_obs
from ora.metadata import resolve_columns
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--counts-out", default=None)
    parser.add_argument("--features-out", default=None)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    args = parser.parse_args()

    config = load_config(args.config)
    h5ad_path = project_path(args.h5ad or config["source"]["h5ad_path"])
    obs = load_obs(h5ad_path)
    columns = resolve_columns(list(obs.columns), config)
    counts = aggregate_cell_counts(obs, config, columns)
    features = build_cell_state_features(counts, config, pseudocount=args.pseudocount)

    outputs = config.get("outputs", {})
    counts_path = args.counts_out or outputs.get("cell_counts_tsv", "data/processed/donor_cell_state_counts.tsv")
    features_path = args.features_out or outputs.get("cell_features_tsv", "data/processed/donor_cell_state_features.tsv")
    counts.to_csv(ensure_parent(counts_path), sep="\t", index=False)
    features.to_csv(ensure_parent(features_path), sep="\t", index=False)
    print(f"Wrote cell counts: {counts_path} ({counts.shape[0]} rows)")
    print(f"Wrote cell-state features: {features_path} ({features.shape[1] - 1} features)")


if __name__ == "__main__":
    main()

