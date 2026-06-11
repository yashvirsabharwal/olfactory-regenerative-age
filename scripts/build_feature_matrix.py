#!/usr/bin/env python3
"""Build the composition-only ORA feature matrix."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--features", default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    feature_path = args.features or outputs.get("cell_features_tsv", "data/processed/donor_cell_state_features.tsv")
    out_path = args.out or outputs.get("ora_feature_matrix_tsv", "data/processed/ora_feature_matrix.tsv")
    features = pd.read_csv(feature_path, sep="\t")
    if "donor_id" not in features.columns:
        raise KeyError("Input feature table must include donor_id.")

    # MVP feature matrix keeps composition-derived biological features. Technical/yield
    # columns remain out of the ORA feature matrix and are read from the manifest.
    keep = ["donor_id"] + [
        col
        for col in features.columns
        if col.startswith("prop__") or col.startswith("clr__") or col.startswith("ratio__")
    ]
    matrix = features[keep].copy()
    matrix.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(f"Wrote ORA feature matrix: {out_path} ({matrix.shape[1] - 1} biological features)")


if __name__ == "__main__":
    main()

