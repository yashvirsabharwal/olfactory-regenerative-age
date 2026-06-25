#!/usr/bin/env python3
"""Build ORA regeneration-axis feature-map artifacts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.regeneration_axis import (
    build_regeneration_axis_feature_map,
    build_regeneration_axis_theme_summary,
    build_regeneration_feature_resource_map,
    write_regeneration_axis_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-matrix", default="data/processed/ora_augmented_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument(
        "--feature-stability",
        default="results/tables/ora_augmented_candidate_repeated_cv_feature_stability.tsv",
    )
    parser.add_argument(
        "--cross-tissue-classification",
        default="results/tables/ora_cross_tissue_feature_classification.tsv",
    )
    parser.add_argument(
        "--resource-out",
        default="resources/feature_maps/regeneration_axis_feature_map.tsv",
    )
    parser.add_argument("--feature-map-out", default="results/tables/regeneration_axis_feature_map.tsv")
    parser.add_argument(
        "--summary-out",
        default="results/tables/regeneration_axis_theme_summary.tsv",
    )
    parser.add_argument(
        "--figure-pdf",
        default="results/figures/manuscript_figure_regeneration_axis.pdf",
    )
    parser.add_argument(
        "--figure-png",
        default="results/figures/manuscript_figure_regeneration_axis.png",
    )
    args = parser.parse_args()

    feature_matrix = pd.read_csv(args.feature_matrix, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    feature_stability = _read_optional_table(args.feature_stability)
    cross_tissue = _read_optional_table(args.cross_tissue_classification)

    resource_map = build_regeneration_feature_resource_map(feature_matrix)
    feature_map = build_regeneration_axis_feature_map(
        feature_matrix=feature_matrix,
        manifest=manifest,
        feature_stability=feature_stability,
        cross_tissue_classification=cross_tissue,
    )
    summary = build_regeneration_axis_theme_summary(feature_map)
    write_regeneration_axis_outputs(
        resource_map=resource_map,
        feature_map=feature_map,
        summary=summary,
        resource_out=args.resource_out,
        feature_map_out=args.feature_map_out,
        summary_out=args.summary_out,
        figure_pdf=args.figure_pdf,
        figure_png=args.figure_png,
    )
    print(
        "Wrote regeneration-axis feature map: "
        f"{args.feature_map_out} ({feature_map.shape[0]} features, {summary.shape[0]} themes)"
    )


def _read_optional_table(path: str) -> pd.DataFrame | None:
    if not path or not Path(path).exists():
        return None
    return pd.read_csv(path, sep="\t")


if __name__ == "__main__":
    main()
