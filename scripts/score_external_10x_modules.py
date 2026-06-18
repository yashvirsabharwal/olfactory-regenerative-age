#!/usr/bin/env python3
"""Score sample-level gene modules from an external raw 10x GEO archive."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.external import parse_geo_series_matrix_metadata, score_external_10x_modules
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", default="data/external/GSE184117_RAW.tar")
    parser.add_argument("--metadata", default="data/external/GSE184117_series_matrix.txt.gz")
    parser.add_argument("--dataset-id", default="oliva_2022")
    parser.add_argument("--gene-set-config", default="configs/gene_sets.yaml")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--sample-metadata-out", default=None)
    parser.add_argument("--sample-qc-out", default=None)
    parser.add_argument("--module-scores-out", default=None)
    parser.add_argument("--module-contrasts-out", default=None)
    args = parser.parse_args()

    external_config = load_config(args.external_config)
    gene_set_config = load_config(args.gene_set_config)
    gene_set_config["published_gene_lists"] = external_config.get("published_gene_lists", {})
    outputs = external_config.get("outputs", {})

    sample_metadata_out = args.sample_metadata_out or outputs.get(
        "external_sample_metadata_tsv",
        "results/tables/external_sample_metadata.tsv",
    )
    sample_qc_out = args.sample_qc_out or outputs.get(
        "external_10x_sample_qc_tsv",
        "results/tables/external_10x_sample_qc.tsv",
    )
    module_scores_out = args.module_scores_out or outputs.get(
        "external_10x_module_scores_tsv",
        "results/tables/external_10x_module_scores.tsv",
    )
    module_contrasts_out = args.module_contrasts_out or outputs.get(
        "external_10x_module_contrasts_tsv",
        "results/tables/external_10x_module_contrasts.tsv",
    )

    metadata = parse_geo_series_matrix_metadata(args.metadata, dataset_id=args.dataset_id)
    qc, scores, contrasts = score_external_10x_modules(
        args.archive,
        metadata,
        gene_set_config,
        dataset_id=args.dataset_id,
    )

    metadata.to_csv(ensure_parent(sample_metadata_out), sep="\t", index=False)
    qc.to_csv(ensure_parent(sample_qc_out), sep="\t", index=False)
    scores.to_csv(ensure_parent(module_scores_out), sep="\t", index=False)
    contrasts.to_csv(ensure_parent(module_contrasts_out), sep="\t", index=False)
    usable = (
        int(metadata["usable_for_external_validation"].sum())
        if not metadata.empty and "usable_for_external_validation" in metadata
        else 0
    )
    print(f"Wrote external sample metadata: {sample_metadata_out} ({usable} usable biopsy samples)")
    print(f"Wrote external 10x sample QC: {sample_qc_out} ({qc.shape[0]} samples)")
    print(f"Wrote external 10x module scores: {module_scores_out} ({scores.shape[0]} rows)")
    print(f"Wrote external 10x module contrasts: {module_contrasts_out} ({contrasts.shape[0]} modules)")


if __name__ == "__main__":
    main()
