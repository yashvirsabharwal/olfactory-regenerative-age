#!/usr/bin/env python3
"""Run a lightweight Milo-style latent-neighborhood DA pilot."""

from __future__ import annotations

import argparse
import re
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
    parser.add_argument("--include-fine-regex", default=None, help="Regex of fine cell types to keep before neighborhood testing.")
    parser.add_argument("--include-coarse-regex", default=None, help="Regex of coarse cell types to keep before neighborhood testing.")
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
    keep_cells = _cell_filter(cell_metadata, fine_regex=args.include_fine_regex, coarse_regex=args.include_coarse_regex)
    if keep_cells is not None:
        embedding = embedding[keep_cells, :]
        cell_metadata = cell_metadata.loc[keep_cells].reset_index(drop=True)

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


def _cell_filter(
    cell_metadata: pd.DataFrame,
    *,
    fine_regex: str | None = None,
    coarse_regex: str | None = None,
) -> pd.Series | None:
    filters = []
    if fine_regex:
        if "fine_celltype" not in cell_metadata:
            raise SystemExit("--include-fine-regex requires `fine_celltype` in AnnData.obs")
        pattern = re.compile(fine_regex, flags=re.IGNORECASE)
        filters.append(cell_metadata["fine_celltype"].astype(str).map(lambda value: bool(pattern.search(value))))
    if coarse_regex:
        if "coarse_celltype" not in cell_metadata:
            raise SystemExit("--include-coarse-regex requires `coarse_celltype` in AnnData.obs")
        pattern = re.compile(coarse_regex, flags=re.IGNORECASE)
        filters.append(cell_metadata["coarse_celltype"].astype(str).map(lambda value: bool(pattern.search(value))))
    if not filters:
        return None
    keep = filters[0].copy()
    for item in filters[1:]:
        keep &= item
    if int(keep.sum()) == 0:
        raise SystemExit("Cell filter kept 0 cells; revise the regex.")
    return keep


if __name__ == "__main__":
    main()
