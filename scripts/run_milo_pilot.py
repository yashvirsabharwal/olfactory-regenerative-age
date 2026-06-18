#!/usr/bin/env python3
"""Run a lightweight Milo-style latent-neighborhood DA pilot."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import anndata as ad
import pandas as pd

from ora.neighborhood import NeighborhoodConfig, run_neighborhood_da
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", default="data/processed/gateway_scvi_stratified_250k.h5ad")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--embedding-key", default="X_scvi")
    parser.add_argument("--out", default="results/tables/milo_pilot_neighborhood_da.tsv")
    parser.add_argument("--summary-out", default="results/tables/milo_pilot_summary.tsv")
    parser.add_argument("--n-neighborhoods", type=int, default=1000)
    parser.add_argument("--n-neighbors", type=int, default=50)
    parser.add_argument("--min-donors", type=int, default=20)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--include-disease", action="store_true", help="Do not restrict donor metadata to healthy donors.")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.h5ad, backed="r")
    try:
        if args.embedding_key not in adata.obsm:
            raise SystemExit(f"Embedding key `{args.embedding_key}` is missing from {args.h5ad}")
        embedding = adata.obsm[args.embedding_key][:]
        cell_metadata = adata.obs.reset_index(drop=True)
    finally:
        close = getattr(adata, "file", None)
        if close is not None:
            close.close()

    donor_metadata = pd.read_csv(args.manifest, sep="\t")
    if not args.include_disease and "is_healthy" in donor_metadata:
        donor_metadata = donor_metadata[donor_metadata["is_healthy"].astype(str).str.lower().eq("true")].copy()
    donor_metadata = donor_metadata[donor_metadata["age"].notna()].copy()

    config = NeighborhoodConfig(
        n_neighborhoods=args.n_neighborhoods,
        n_neighbors=args.n_neighbors,
        min_donors=args.min_donors,
        seed=args.seed,
    )
    neighborhoods, summary = run_neighborhood_da(
        embedding,
        cell_metadata,
        donor_metadata,
        config=config,
    )
    neighborhoods.to_csv(ensure_parent(args.out), sep="\t", index=False)
    summary.to_csv(ensure_parent(args.summary_out), sep="\t", index=False)
    print(f"Wrote Milo-style neighborhood DA: {args.out} ({neighborhoods.shape[0]} rows)")
    print(f"Wrote Milo-style summary: {args.summary_out} ({summary.shape[0]} rows)")


if __name__ == "__main__":
    main()
