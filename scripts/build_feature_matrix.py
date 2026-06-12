#!/usr/bin/env python3
"""Build ORA feature matrices from composition and optional module features."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.features import build_ora_feature_matrix, feature_kind_counts
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--features", default=None)
    parser.add_argument("--module-features", default=None)
    parser.add_argument(
        "--include-modules",
        action="store_true",
        help="Merge donor-level module_score__ features into the ORA matrix.",
    )
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    feature_path = args.features or outputs.get("cell_features_tsv", "data/processed/donor_cell_state_features.tsv")
    include_modules = args.include_modules or bool(args.module_features)
    if include_modules:
        out_path = args.out or outputs.get(
            "ora_augmented_feature_matrix_tsv",
            "data/processed/ora_augmented_feature_matrix.tsv",
        )
        module_path = args.module_features or outputs.get(
            "donor_module_features_tsv",
            "data/processed/donor_module_features.tsv",
        )
    else:
        out_path = args.out or outputs.get("ora_feature_matrix_tsv", "data/processed/ora_feature_matrix.tsv")
        module_path = None
    features = pd.read_csv(feature_path, sep="\t")
    modules = pd.read_csv(module_path, sep="\t") if module_path else None
    matrix = build_ora_feature_matrix(features, modules)
    matrix.to_csv(ensure_parent(out_path), sep="\t", index=False)
    counts = feature_kind_counts(matrix)
    print(
        f"Wrote ORA feature matrix: {out_path} "
        f"({matrix.shape[1] - 1} biological features; "
        f"{counts['composition']} composition, {counts['module']} module)"
    )


if __name__ == "__main__":
    main()
