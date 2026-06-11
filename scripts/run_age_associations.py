#!/usr/bin/env python3
"""Run covariate-adjusted age association models for donor-level features."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.stats import run_age_associations
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", default="data/processed/donor_cell_state_features.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--out", default="results/tables/age_cell_state_associations.tsv")
    parser.add_argument(
        "--all-donors",
        action="store_true",
        help="Use all donors with age instead of the healthy ORA-training cohort.",
    )
    args = parser.parse_args()

    load_config(args.config)
    features = pd.read_csv(args.features, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    if not args.all_donors and "usable_for_ora_training" in manifest.columns:
        manifest = manifest[manifest["usable_for_ora_training"].astype(bool)].copy()
    covariates = [
        col
        for col in ["sex", "race_ethnicity", "chemistry", "collection_method", "site", "total_cells"]
        if col in manifest.columns
    ]
    feature_columns = [
        col
        for col in features.columns
        if col.startswith("prop__") or col.startswith("clr__") or col.startswith("ratio__")
    ]
    results = run_age_associations(features, manifest, feature_columns=feature_columns, covariates=covariates)
    results.to_csv(ensure_parent(args.out), sep="\t", index=False)
    print(f"Wrote age associations: {args.out} ({results.shape[0]} tests)")


if __name__ == "__main__":
    main()
