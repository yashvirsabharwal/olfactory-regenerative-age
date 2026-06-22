#!/usr/bin/env python3
"""Export a stratified latent subset for official MiloR parity testing."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", default="data/processed/gateway_scvi_full_4m_reduced.h5ad")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--embedding-key", default="X_scvi")
    parser.add_argument("--include-coarse-regex", default="Resp_HBC|Olf_INPs|Olf_iOSNs|Olf_mOSNs|Olf_Sus")
    parser.add_argument("--donor-query", default=None)
    parser.add_argument("--include-disease", action="store_true")
    parser.add_argument("--max-cells", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=29)
    parser.add_argument("--metadata-out", required=True)
    parser.add_argument("--embedding-out", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    if not args.include_disease and "is_healthy" in manifest:
        manifest = manifest[manifest["is_healthy"].astype(str).str.lower().eq("true")].copy()
    manifest = manifest[manifest["age"].notna()].copy()
    if args.donor_query:
        try:
            manifest = manifest.query(args.donor_query).copy()
        except Exception as exc:
            raise SystemExit(f"Invalid --donor-query `{args.donor_query}`: {exc}") from exc
    manifest = manifest.drop_duplicates("donor_id").copy()
    if manifest.empty:
        raise SystemExit("No donors remain after filters.")
    keep_donors = set(manifest["donor_id"].astype(str))

    adata = ad.read_h5ad(args.h5ad, backed="r")
    try:
        if args.embedding_key not in adata.obsm:
            raise SystemExit(f"Embedding key `{args.embedding_key}` is missing from {args.h5ad}")
        obs = adata.obs.copy()
        obs["_obs_name"] = obs.index.astype(str)
        obs["_cell_index"] = np.arange(obs.shape[0])
        obs["donor_id"] = obs["donor_id"].astype(str)
        keep = obs["donor_id"].isin(keep_donors)
        if args.include_coarse_regex:
            pattern = re.compile(args.include_coarse_regex, flags=re.IGNORECASE)
            keep &= obs["coarse_celltype"].astype(str).map(lambda value: bool(pattern.search(value)))
        candidate = obs.loc[keep].copy()
        if candidate.empty:
            raise SystemExit("No cells remain after donor and cell-state filters.")
        chosen = _stratified_sample(candidate, max_cells=args.max_cells, seed=args.seed)
        embedding = np.asarray(adata.obsm[args.embedding_key][chosen["_cell_index"].to_numpy(dtype=int), :], dtype=float)
    finally:
        close = getattr(adata, "file", None)
        if close is not None:
            close.close()

    metadata = chosen.merge(manifest, on="donor_id", how="left", suffixes=("", "_manifest"))
    metadata = metadata[
        [
            "_obs_name",
            "_cell_index",
            "donor_id",
            "sample_id",
            "age",
            "sex",
            "chemistry",
            "collection_method",
            "total_cells",
            "fine_celltype",
            "coarse_celltype",
        ]
    ].copy()
    metadata = metadata.rename(columns={"_obs_name": "cell_id", "_cell_index": "cell_index"})
    embedding_out = pd.DataFrame(embedding, columns=[f"X_scvi_{idx + 1}" for idx in range(embedding.shape[1])])
    embedding_out.insert(0, "cell_id", metadata["cell_id"].to_numpy())
    summary = pd.DataFrame(
        [
            {"metric": "candidate_cells", "value": int(candidate.shape[0]), "detail": "cells after donor/cell-state filters"},
            {"metric": "exported_cells", "value": int(metadata.shape[0]), "detail": "cells in MiloR subset"},
            {"metric": "donors", "value": int(metadata["donor_id"].nunique()), "detail": "donors represented"},
            {"metric": "embedding_dimensions", "value": int(embedding.shape[1]), "detail": args.embedding_key},
        ]
    )
    metadata.to_csv(ensure_parent(args.metadata_out), sep="\t", index=False)
    embedding_out.to_csv(ensure_parent(args.embedding_out), sep="\t", index=False)
    summary.to_csv(ensure_parent(args.summary_out), sep="\t", index=False)
    print(f"Wrote MiloR metadata: {args.metadata_out} ({metadata.shape[0]} cells)")
    print(f"Wrote MiloR embedding: {args.embedding_out} ({embedding_out.shape[0]} cells x {embedding.shape[1]} dims)")
    print(f"Wrote MiloR export summary: {args.summary_out}")


def _stratified_sample(frame: pd.DataFrame, *, max_cells: int, seed: int) -> pd.DataFrame:
    if frame.shape[0] <= max_cells:
        return frame.sort_values("_cell_index").reset_index(drop=True)
    rng = np.random.default_rng(seed)
    strata = frame[["donor_id", "coarse_celltype", "fine_celltype"]].astype(str).agg("|".join, axis=1)
    groups = list(frame.groupby(strata, observed=True).indices.values())
    if len(groups) > max_cells:
        group_indices = rng.choice(np.arange(len(groups)), size=max_cells, replace=False)
        chosen = [int(rng.choice(np.asarray(groups[idx], dtype=int), size=1)[0]) for idx in group_indices]
        return frame.iloc[np.sort(np.asarray(chosen, dtype=int))].reset_index(drop=True)
    base = max(1, max_cells // max(len(groups), 1))
    chosen = []
    remainders = []
    for idx in groups:
        idx_array = np.asarray(idx, dtype=int)
        take = min(base, idx_array.size)
        selected = rng.choice(idx_array, size=take, replace=False)
        chosen.extend(selected.tolist())
        if idx_array.size > take:
            remaining = np.setdiff1d(idx_array, selected, assume_unique=False)
            remainders.extend(remaining.tolist())
    if len(chosen) < max_cells and remainders:
        add = min(max_cells - len(chosen), len(remainders))
        chosen.extend(rng.choice(np.asarray(remainders, dtype=int), size=add, replace=False).tolist())
    return frame.iloc[np.sort(np.asarray(chosen[:max_cells], dtype=int))].reset_index(drop=True)


if __name__ == "__main__":
    main()
