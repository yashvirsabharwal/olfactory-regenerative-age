#!/usr/bin/env python3
"""Build a reduced-gene AnnData object from a large H5AD in row chunks."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from scipy import sparse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, help="Source H5AD path.")
    parser.add_argument("--gene-list-file", required=True, help="One var_name per line to keep.")
    parser.add_argument("--out", required=True, help="Reduced output .h5ad or .zarr path.")
    parser.add_argument("--chunk-dir", required=True, help="Directory for temporary reduced chunk H5ADs.")
    parser.add_argument("--chunk-size", type=int, default=25000)
    parser.add_argument("--max-cells", type=int, default=None, help="Optional smoke-test cell limit.")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--keep-chunks", action="store_true")
    args = parser.parse_args()

    import anndata as ad  # type: ignore

    out_path = Path(args.out)
    chunk_dir = Path(args.chunk_dir)
    if out_path.exists():
        if not args.overwrite:
            raise SystemExit(f"Output already exists: {out_path}. Pass --overwrite to replace it.")
        _remove_path(out_path)
    if chunk_dir.exists():
        if not args.overwrite:
            raise SystemExit(f"Chunk directory already exists: {chunk_dir}. Pass --overwrite to replace it.")
        shutil.rmtree(chunk_dir)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    chunk_dir.mkdir(parents=True, exist_ok=True)

    genes = _read_gene_list(args.gene_list_file)
    _log(f"Opening source H5AD: {args.h5ad}")
    backed = ad.read_h5ad(args.h5ad, backed="r")
    chunk_paths: list[Path] = []
    try:
        n_obs = backed.n_obs if args.max_cells is None else min(args.max_cells, backed.n_obs)
        var_indices = _selected_var_indices(backed.var_names, genes)
        _log(f"Source shape: {backed.n_obs} cells x {backed.n_vars} genes")
        _log(f"Keeping {len(genes)} genes and {n_obs} cells")
        for chunk_id, start in enumerate(range(0, n_obs, args.chunk_size)):
            stop = min(start + args.chunk_size, n_obs)
            chunk_path = chunk_dir / f"chunk_{chunk_id:06d}.h5ad"
            _log(f"Writing chunk {chunk_id}: rows {start}:{stop}")
            chunk = backed[start:stop, :].to_memory()
            chunk = chunk[:, var_indices].copy()
            if not sparse.issparse(chunk.X):
                chunk.X = sparse.csr_matrix(chunk.X)
            else:
                chunk.X = chunk.X.tocsr()
            chunk.write_h5ad(chunk_path)
            chunk_paths.append(chunk_path)
            _log(f"Wrote {chunk_path} ({chunk.n_obs} cells x {chunk.n_vars} genes)")
    finally:
        close = getattr(backed, "file", None)
        if close is not None:
            close.close()

    _log(f"Concatenating {len(chunk_paths)} chunks to {out_path}")
    ad.experimental.concat_on_disk(
        [str(path) for path in chunk_paths],
        str(out_path),
        axis=0,
        join="inner",
        merge="same",
        uns_merge=None,
        max_loaded_elems=25_000_000,
    )
    _log(f"Wrote reduced H5AD: {out_path}")
    if not args.keep_chunks:
        shutil.rmtree(chunk_dir)
        _log(f"Removed chunk directory: {chunk_dir}")


def _read_gene_list(path: str) -> tuple[str, ...]:
    with Path(path).open(encoding="utf-8") as handle:
        genes = tuple(line.strip() for line in handle if line.strip() and not line.startswith("#"))
    if not genes:
        raise ValueError(f"No genes found in {path}")
    return genes


def _selected_var_indices(var_names: pd.Index, selected_var_names: tuple[str, ...]) -> np.ndarray:
    index = pd.Index(var_names.astype(str))
    positions = index.get_indexer(selected_var_names)
    missing = [name for name, position in zip(selected_var_names, positions, strict=False) if position < 0]
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(f"Selected gene list has {len(missing)} genes missing from source H5AD: {preview}")
    return positions


def _remove_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def _log(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


if __name__ == "__main__":
    main()
