#!/usr/bin/env python3
"""Estimate marker-only composition from an external raw 10x GEO archive."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.external import parse_geo_series_matrix_metadata, score_external_10x_marker_composition
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", default="data/external/GSE184117_RAW.tar")
    parser.add_argument("--metadata", default="data/external/GSE184117_series_matrix.txt.gz")
    parser.add_argument("--dataset-id", default="oliva_2022")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--sample-metadata-out", default=None)
    parser.add_argument("--marker-coverage-out", default=None)
    parser.add_argument("--marker-composition-out", default=None)
    parser.add_argument("--marker-contrasts-out", default=None)
    args = parser.parse_args()

    external_config = load_config(args.external_config)
    outputs = external_config.get("outputs", {})
    sample_metadata_out = args.sample_metadata_out or outputs.get(
        "external_sample_metadata_tsv",
        "results/tables/external_sample_metadata.tsv",
    )
    marker_coverage_out = args.marker_coverage_out or outputs.get(
        "external_10x_marker_coverage_tsv",
        "results/tables/external_10x_marker_coverage.tsv",
    )
    marker_composition_out = args.marker_composition_out or outputs.get(
        "external_10x_marker_composition_tsv",
        "results/tables/external_10x_marker_composition.tsv",
    )
    marker_contrasts_out = args.marker_contrasts_out or outputs.get(
        "external_10x_marker_contrasts_tsv",
        "results/tables/external_10x_marker_contrasts.tsv",
    )

    metadata = parse_geo_series_matrix_metadata(args.metadata, dataset_id=args.dataset_id)
    coverage, composition, contrasts = score_external_10x_marker_composition(
        args.archive,
        metadata,
        external_config.get("marker_panels"),
        dataset_id=args.dataset_id,
    )
    metadata.to_csv(ensure_parent(sample_metadata_out), sep="\t", index=False)
    coverage.to_csv(ensure_parent(marker_coverage_out), sep="\t", index=False)
    composition.to_csv(ensure_parent(marker_composition_out), sep="\t", index=False)
    contrasts.to_csv(ensure_parent(marker_contrasts_out), sep="\t", index=False)
    print(f"Wrote external sample metadata: {sample_metadata_out} ({metadata.shape[0]} samples)")
    print(f"Wrote external marker coverage: {marker_coverage_out} ({coverage.shape[0]} rows)")
    print(f"Wrote external marker composition: {marker_composition_out} ({composition.shape[0]} rows)")
    print(f"Wrote external marker contrasts: {marker_contrasts_out} ({contrasts.shape[0]} marker panels)")


if __name__ == "__main__":
    main()
