#!/usr/bin/env python3
"""Build a harmonized marker-mapped AnnData from an external raw 10x GEO archive."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.external import build_external_10x_marker_mapped_anndata, parse_geo_series_matrix_metadata
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", default="data/external/GSE184117_RAW.tar")
    parser.add_argument("--metadata", default="data/external/GSE184117_series_matrix.txt.gz")
    parser.add_argument("--dataset-id", default="oliva_2022")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--out", default=None)
    parser.add_argument("--mapping-qc-out", default=None)
    parser.add_argument("--donor-features-out", default=None)
    args = parser.parse_args()

    external_config = load_config(args.external_config)
    outputs = external_config.get("outputs", {})
    out = args.out or outputs.get("external_10x_mapped_h5ad", "data/processed/gse184117_marker_mapped.h5ad")
    mapping_qc_out = args.mapping_qc_out or outputs.get(
        "external_10x_mapping_qc_tsv",
        "results/tables/external_10x_mapping_qc.tsv",
    )
    donor_features_out = args.donor_features_out or outputs.get(
        "external_10x_mapped_donor_features_tsv",
        "data/processed/gse184117_mapped_donor_features.tsv",
    )

    metadata = parse_geo_series_matrix_metadata(args.metadata, dataset_id=args.dataset_id)
    adata, mapping_qc, donor_features = build_external_10x_marker_mapped_anndata(
        args.archive,
        metadata,
        external_config.get("marker_panels"),
        dataset_id=args.dataset_id,
    )
    ensure_parent(out)
    adata.write_h5ad(out)
    mapping_qc.to_csv(ensure_parent(mapping_qc_out), sep="\t", index=False)
    donor_features.to_csv(ensure_parent(donor_features_out), sep="\t", index=False)
    print(f"Wrote external marker-mapped AnnData: {out} ({adata.n_obs} cells x {adata.n_vars} genes)")
    print(f"Wrote external mapping QC: {mapping_qc_out} ({mapping_qc.shape[0]} samples)")
    print(f"Wrote external mapped donor features: {donor_features_out} ({donor_features.shape[0]} rows)")


if __name__ == "__main__":
    main()
