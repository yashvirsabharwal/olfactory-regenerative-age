#!/usr/bin/env python3
"""Build a claim-gated external validation evidence ledger."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.external import external_dataset_summary, external_validation_evidence_summary
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--validation-summary", default=None)
    parser.add_argument("--sample-metadata", default=None)
    parser.add_argument("--module-contrasts", default=None)
    parser.add_argument("--marker-contrasts", default=None)
    parser.add_argument("--mapped-features", default=None)
    parser.add_argument("--scanvi-features", default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    config = load_config(args.external_config)
    outputs = config.get("outputs", {})
    summary_path = args.validation_summary or outputs.get(
        "validation_summary_tsv",
        "results/tables/external_validation_summary.tsv",
    )
    evidence_out = args.out or outputs.get(
        "external_validation_evidence_tsv",
        "results/tables/external_validation_evidence.tsv",
    )
    dataset_summary = _read_optional_tsv(summary_path)
    if dataset_summary is None:
        dataset_summary = external_dataset_summary(config)
    evidence = external_validation_evidence_summary(
        config,
        dataset_summary,
        sample_metadata=_read_optional_tsv(
            args.sample_metadata
            or outputs.get("external_sample_metadata_tsv", "results/tables/external_sample_metadata.tsv")
        ),
        module_contrasts=_read_optional_tsv(
            args.module_contrasts
            or outputs.get("external_10x_module_contrasts_tsv", "results/tables/external_10x_module_contrasts.tsv")
        ),
        marker_contrasts=_read_optional_tsv(
            args.marker_contrasts
            or outputs.get("external_10x_marker_contrasts_tsv", "results/tables/external_10x_marker_contrasts.tsv")
        ),
        mapped_features=_read_optional_tsv(
            args.mapped_features
            or outputs.get("external_10x_mapped_donor_features_tsv", "data/processed/gse184117_mapped_donor_features.tsv")
        ),
        scanvi_features=_read_optional_tsv(
            args.scanvi_features
            or outputs.get("external_scanvi_donor_features_tsv", "data/processed/gse184117_scanvi_donor_features.tsv")
        ),
    )
    evidence.to_csv(ensure_parent(evidence_out), sep="\t", index=False)
    print(f"Wrote external validation evidence ledger: {evidence_out} ({evidence.shape[0]} rows)")


def _read_optional_tsv(path: str | None) -> pd.DataFrame | None:
    if not path:
        return None
    candidate = Path(path)
    if not candidate.exists():
        return None
    return pd.read_csv(candidate, sep="\t")


if __name__ == "__main__":
    main()
