#!/usr/bin/env python3
"""Compare external mapped GSE184117 features with Gateway age directions."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.external import external_mapped_feature_concordance
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--mapped-features", default=None)
    parser.add_argument("--age-associations", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument(
        "--direct-feature-map",
        action="store_true",
        help="Map external prop__/clr__ features directly to same-named Gateway features.",
    )
    args = parser.parse_args()

    external_config = load_config(args.external_config)
    gateway_config = load_config(args.gateway_config)
    external_outputs = external_config.get("outputs", {})
    gateway_outputs = gateway_config.get("outputs", {})
    mapped_path = args.mapped_features or external_outputs.get(
        "external_10x_mapped_donor_features_tsv",
        "data/processed/gse184117_mapped_donor_features.tsv",
    )
    age_path = args.age_associations or gateway_outputs.get(
        "age_associations_tsv",
        "results/tables/age_cell_state_associations.tsv",
    )
    out_path = args.out or external_outputs.get(
        "external_mapped_feature_concordance_tsv",
        "results/tables/external_mapped_feature_concordance.tsv",
    )
    mapped = pd.read_csv(mapped_path, sep="\t")
    age = pd.read_csv(age_path, sep="\t")
    feature_map = _direct_feature_map(mapped, age) if args.direct_feature_map else None
    concordance = external_mapped_feature_concordance(mapped, age, marker_to_gateway=feature_map)
    concordance.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(
        f"Wrote external mapped-feature concordance: {out_path} "
        f"({concordance.shape[0]} mapped feature rows)"
    )


def _direct_feature_map(mapped_features: pd.DataFrame, age_associations: pd.DataFrame) -> dict[str, tuple[str, ...]]:
    if "feature" not in age_associations:
        return {}
    gateway_features = set(age_associations["feature"].astype(str))
    mapping: dict[str, list[str]] = {}
    for feature in mapped_features.columns:
        if not str(feature).startswith(("prop__", "clr__")) or str(feature) not in gateway_features:
            continue
        panel = str(feature).split("__", 1)[1]
        mapping.setdefault(panel, []).append(str(feature))
    return {panel: tuple(features) for panel, features in mapping.items()}


if __name__ == "__main__":
    main()
