#!/usr/bin/env python3
"""Run a guarded scVI latent recomputation workflow."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import re

import numpy as np
import pandas as pd

DEFAULT_PRESERVED_MARKERS = (
    "TP63",
    "KRT5",
    "KRT14",
    "ASCL1",
    "NEUROG1",
    "NEUROD1",
    "GAP43",
    "DCX",
    "TUBB3",
    "OMP",
    "ADCY3",
    "GNAL",
    "CYP2A13",
    "CYP2J2",
    "MUC1",
    "PTPRC",
    "LST1",
    "TYROBP",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--model-dir", default="")
    parser.add_argument("--max-cells", type=int, default=None)
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--gene-list-h5ad", default="")
    parser.add_argument("--gene-list-file", default="")
    parser.add_argument("--batch-key", default="sample_id")
    parser.add_argument("--categorical-covariates", default="chemistry,collection_method,sex")
    parser.add_argument("--sampling-strategy", choices=("random", "stratified"), default="random")
    parser.add_argument("--stratify-keys", default="condition,fine_celltype,sex,flex_version,device_guided")
    parser.add_argument("--include-fine-celltype-regex", default="")
    parser.add_argument("--include-coarse-celltype-regex", default="")
    parser.add_argument("--embedding-key", default="X_scvi")
    parser.add_argument("--max-epochs", type=int, default=200)
    parser.add_argument("--hvg-flavor", default="cell_ranger")
    parser.add_argument("--hvg-batch-key", default="")
    parser.add_argument("--preserve-marker-genes", default=",".join(DEFAULT_PRESERVED_MARKERS))
    parser.add_argument("--accelerator", default="auto")
    parser.add_argument("--devices", default="auto")
    parser.add_argument("--seed", type=int, default=13)
    args = parser.parse_args()

    try:
        import scanpy as sc  # type: ignore
        import scvi  # type: ignore
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "scanpy and scvi-tools are required for latent recomputation. "
            "Install the latent extra first, for example `python -m pip install -e .[latent]`."
        ) from exc

    scvi.settings.seed = args.seed
    selected_var_names = _selected_var_names(sc, h5ad_path=args.gene_list_h5ad, text_path=args.gene_list_file)
    _log(
        "Starting scVI latent run: "
        f"h5ad={args.h5ad}, max_cells={args.max_cells or 'all'}, "
        f"selected_genes={len(selected_var_names) or 'HVG'}, seed={args.seed}"
    )
    adata = _read_pilot_h5ad(
        sc,
        args.h5ad,
        max_cells=args.max_cells,
        seed=args.seed,
        sampling_strategy=args.sampling_strategy,
        stratify_keys=_split_arg(args.stratify_keys),
        include_fine_celltype_regex=args.include_fine_celltype_regex,
        include_coarse_celltype_regex=args.include_coarse_celltype_regex,
        selected_var_names=selected_var_names,
    )
    _log(f"Loaded AnnData into memory: {adata.n_obs} cells x {adata.n_vars} genes")
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        _log("Using counts layer as X")
    batch_key = args.batch_key if args.batch_key in adata.obs else None
    categorical_covariates = [
        key.strip()
        for key in args.categorical_covariates.split(",")
        if key.strip() and key.strip() in adata.obs and key.strip() != batch_key
    ]
    hvg_batch_key = args.hvg_batch_key if args.hvg_batch_key in adata.obs else None
    preserve_marker_genes = _split_arg(args.preserve_marker_genes)
    if selected_var_names:
        adata.var["highly_variable"] = True
        _log("Using fixed gene list; skipping HVG recomputation")
    else:
        _log("Selecting HVGs")
        try:
            _select_hvgs_with_preserved_markers(
                sc,
                adata,
                n_top_genes=args.n_top_genes,
                flavor=args.hvg_flavor,
                batch_key=hvg_batch_key,
                preserve_marker_genes=preserve_marker_genes,
            )
        except ValueError as exc:
            if hvg_batch_key is None or "Bin edges must be unique" not in str(exc):
                raise
            print(f"HVG selection failed for batch key {hvg_batch_key!r}; retrying without an HVG batch key.")
            _select_hvgs_with_preserved_markers(
                sc,
                adata,
                n_top_genes=args.n_top_genes,
                flavor=args.hvg_flavor,
                batch_key=None,
                preserve_marker_genes=preserve_marker_genes,
            )
        _log(f"Selected HVGs/markers: {adata.n_vars} genes")
    _log(
        "Setting up scVI: "
        f"batch_key={batch_key or 'none'}, covariates={','.join(categorical_covariates) or 'none'}"
    )
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key=batch_key,
        categorical_covariate_keys=categorical_covariates or None,
    )
    model = scvi.model.SCVI(adata)
    _log(f"Training scVI for up to {args.max_epochs} epochs")
    model.train(
        max_epochs=args.max_epochs,
        early_stopping=True,
        accelerator=args.accelerator,
        devices=args.devices,
    )
    _log("Writing latent representation")
    adata.obsm[args.embedding_key] = model.get_latent_representation()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    if args.model_dir:
        model_path = Path(args.model_dir)
        model_path.mkdir(parents=True, exist_ok=True)
        model.save(str(model_path), overwrite=True, save_anndata=False)
        _log(f"Wrote scVI model: {args.model_dir}")
    _log(f"Wrote scVI latent H5AD: {args.out} ({adata.n_obs} cells x {adata.n_vars} genes)")


def _log(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def _read_pilot_h5ad(
    sc,
    h5ad_path: str,
    *,
    max_cells: int | None,
    seed: int,
    sampling_strategy: str = "random",
    stratify_keys: tuple[str, ...] = (),
    include_fine_celltype_regex: str = "",
    include_coarse_celltype_regex: str = "",
    selected_var_names: tuple[str, ...] = (),
):
    """Read a bounded pilot subset without materializing the full H5AD first."""

    _log("Opening source H5AD")
    if (
        max_cells is None
        and not include_fine_celltype_regex
        and not include_coarse_celltype_regex
        and not selected_var_names
    ):
        adata = sc.read_h5ad(h5ad_path)
        _log(f"Read full H5AD directly: {adata.n_obs} cells x {adata.n_vars} genes")
        return adata
    backed = sc.read_h5ad(h5ad_path, backed="r")
    try:
        _log(f"Opened backed H5AD: {backed.n_obs} cells x {backed.n_vars} genes")
        var_indices = _selected_var_indices(backed.var_names, selected_var_names)
        selected_gene_count = len(selected_var_names) if selected_var_names else backed.n_vars
        _log(f"Selected {selected_gene_count} genes from backed H5AD")
        candidate_indices = _candidate_indices(
            backed.obs,
            include_fine_celltype_regex=include_fine_celltype_regex,
            include_coarse_celltype_regex=include_coarse_celltype_regex,
        )
        _log(f"Candidate cells after filters: {candidate_indices.shape[0]}")
        if max_cells is None or candidate_indices.shape[0] <= max_cells:
            _log("Materializing selected backed slice into memory")
            adata = backed[np.sort(candidate_indices), var_indices].to_memory()
            _log("Materialized selected backed slice")
            return adata
        rng = np.random.default_rng(seed)
        if sampling_strategy == "stratified":
            indices = _stratified_sample_indices(
                backed.obs.iloc[candidate_indices],
                candidate_indices,
                max_cells=max_cells,
                seed=seed,
                stratify_keys=stratify_keys,
            )
        else:
            indices = np.sort(rng.choice(candidate_indices, size=max_cells, replace=False))
        _log(f"Materializing sampled backed slice into memory: {indices.shape[0]} cells")
        adata = backed[indices, var_indices].to_memory()
        _log("Materialized sampled backed slice")
        return adata
    finally:
        close = getattr(backed, "file", None)
        if close is not None:
            close.close()


def _candidate_indices(
    obs: pd.DataFrame,
    *,
    include_fine_celltype_regex: str = "",
    include_coarse_celltype_regex: str = "",
) -> np.ndarray:
    mask = pd.Series(True, index=obs.index)
    if include_fine_celltype_regex and "fine_celltype" in obs:
        pattern = re.compile(include_fine_celltype_regex, re.IGNORECASE)
        mask &= obs["fine_celltype"].astype(str).map(lambda value: bool(pattern.search(value)))
    if include_coarse_celltype_regex and "coarse_celltype" in obs:
        pattern = re.compile(include_coarse_celltype_regex, re.IGNORECASE)
        mask &= obs["coarse_celltype"].astype(str).map(lambda value: bool(pattern.search(value)))
    return np.flatnonzero(mask.to_numpy())


def _selected_var_names(sc, *, h5ad_path: str, text_path: str) -> tuple[str, ...]:
    if h5ad_path and text_path:
        raise ValueError("Use only one of --gene-list-h5ad or --gene-list-file.")
    if text_path:
        with Path(text_path).open(encoding="utf-8") as handle:
            return tuple(line.strip() for line in handle if line.strip() and not line.startswith("#"))
    if not h5ad_path:
        return ()
    backed = sc.read_h5ad(h5ad_path, backed="r")
    try:
        return tuple(str(name) for name in backed.var_names)
    finally:
        close = getattr(backed, "file", None)
        if close is not None:
            close.close()


def _selected_var_indices(var_names: pd.Index, selected_var_names: tuple[str, ...]) -> slice | np.ndarray:
    if not selected_var_names:
        return slice(None)
    index = pd.Index(var_names.astype(str))
    positions = index.get_indexer(selected_var_names)
    missing = [name for name, position in zip(selected_var_names, positions, strict=False) if position < 0]
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(f"Selected gene list has {len(missing)} genes missing from source H5AD: {preview}")
    return positions


def _stratified_sample_indices(
    obs: pd.DataFrame,
    source_indices: np.ndarray,
    *,
    max_cells: int,
    seed: int,
    stratify_keys: tuple[str, ...],
) -> np.ndarray:
    """Sample approximately equal numbers from observed strata."""

    keys = tuple(key for key in stratify_keys if key in obs)
    if not keys:
        rng = np.random.default_rng(seed)
        return np.sort(rng.choice(source_indices, size=max_cells, replace=False))
    rng = np.random.default_rng(seed)
    strata = obs.loc[:, keys].astype(str).fillna("missing").agg("|".join, axis=1)
    groups = [source_indices[np.flatnonzero(strata.to_numpy() == value)] for value in strata.unique()]
    selected: list[np.ndarray] = []
    remaining: list[np.ndarray] = []
    quota = max(1, max_cells // len(groups))
    for group in groups:
        take = min(quota, group.shape[0])
        chosen = rng.choice(group, size=take, replace=False)
        selected.append(chosen)
        leftovers = np.setdiff1d(group, chosen, assume_unique=False)
        if leftovers.size:
            remaining.append(leftovers)
    selected_indices = np.concatenate(selected) if selected else np.array([], dtype=int)
    slots = max_cells - selected_indices.shape[0]
    if slots > 0 and remaining:
        pool = np.concatenate(remaining)
        fill = rng.choice(pool, size=min(slots, pool.shape[0]), replace=False)
        selected_indices = np.concatenate([selected_indices, fill])
    if selected_indices.shape[0] > max_cells:
        selected_indices = rng.choice(selected_indices, size=max_cells, replace=False)
    return np.sort(selected_indices)


def _select_hvgs_with_preserved_markers(
    sc,
    adata,
    *,
    n_top_genes: int,
    flavor: str,
    batch_key: str | None,
    preserve_marker_genes: tuple[str, ...],
) -> None:
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
        batch_key=batch_key,
        subset=False,
    )
    keep = adata.var["highly_variable"].to_numpy(dtype=bool)
    if preserve_marker_genes:
        marker_lookup = {gene.upper() for gene in preserve_marker_genes}
        marker_keep = np.array([str(gene).upper() in marker_lookup for gene in adata.var_names], dtype=bool)
        if "feature_name" in adata.var:
            feature_keep = adata.var["feature_name"].astype(str).str.upper().isin(marker_lookup).to_numpy()
            marker_keep |= feature_keep
        keep |= marker_keep
    adata._inplace_subset_var(keep)


def _split_arg(value: str) -> tuple[str, ...]:
    return tuple(part.strip() for part in value.split(",") if part.strip())


if __name__ == "__main__":
    main()
