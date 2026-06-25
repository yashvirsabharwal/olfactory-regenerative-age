#!/usr/bin/env python3
"""Build regeneration-module age-association and ORA-correlation tables."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.regeneration_modules import (
    build_regeneration_module_age_associations,
    build_regeneration_module_ora_correlations,
    parse_regeneration_module_metadata,
)
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-sets", default="configs/regeneration_gene_sets.yaml")
    parser.add_argument(
        "--donor-module-features",
        default="data/processed/regeneration_donor_module_features.tsv",
    )
    parser.add_argument("--coverage", default="results/tables/regeneration_module_gene_coverage.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--ora-features", default="data/processed/ora_augmented_feature_matrix.tsv")
    parser.add_argument("--feature-map", default="results/tables/regeneration_axis_feature_map.tsv")
    parser.add_argument(
        "--age-out",
        default="results/tables/regeneration_module_age_associations.tsv",
    )
    parser.add_argument(
        "--correlation-out",
        default="results/tables/regeneration_module_ora_correlations.tsv",
    )
    args = parser.parse_args()

    gene_set_config = load_config(args.gene_sets)
    module_metadata = parse_regeneration_module_metadata(gene_set_config)
    donor_module_features = pd.read_csv(args.donor_module_features, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    coverage = _read_optional_table(args.coverage)
    ora_features = pd.read_csv(args.ora_features, sep="\t")
    feature_map = _read_optional_table(args.feature_map)

    age = build_regeneration_module_age_associations(
        donor_module_features=donor_module_features,
        manifest=manifest,
        module_metadata=module_metadata,
        coverage=coverage,
    )
    correlations = build_regeneration_module_ora_correlations(
        donor_module_features=donor_module_features,
        ora_feature_matrix=ora_features,
        module_metadata=module_metadata,
        manifest=manifest,
        feature_map=feature_map,
    )
    age.to_csv(ensure_parent(args.age_out), sep="\t", index=False)
    correlations.to_csv(ensure_parent(args.correlation_out), sep="\t", index=False)
    print(
        "Wrote regeneration module analysis: "
        f"{args.age_out} ({age.shape[0]} rows), "
        f"{args.correlation_out} ({correlations.shape[0]} rows)"
    )


def _read_optional_table(path: str) -> pd.DataFrame | None:
    if not path or not Path(path).exists():
        return None
    return pd.read_csv(path, sep="\t")


if __name__ == "__main__":
    main()
