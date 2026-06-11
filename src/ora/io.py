"""Memory-safe AnnData IO helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


def read_h5ad_backed(path: str | Path):
    """Open an H5AD in backed read-only mode."""

    try:
        import anndata as ad  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "anndata is required to read H5AD files. Install the project with `pip install -e .`."
        ) from exc
    return ad.read_h5ad(path, backed="r")


def inspect_h5ad(path: str | Path) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """Return schema metadata without loading the full expression matrix."""

    adata = read_h5ad_backed(path)
    try:
        obs_columns = list(adata.obs.columns)
        var_columns = list(adata.var.columns)
        schema = {
            "path": str(path),
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "obs_columns": obs_columns,
            "var_columns": var_columns,
            "obsm_keys": list(adata.obsm.keys()),
            "layers": list(adata.layers.keys()),
            "uns_keys": list(adata.uns.keys()),
        }
        obs_table = _columns_table(adata.obs)
        var_table = _columns_table(adata.var)
        return schema, obs_table, var_table
    finally:
        close = getattr(adata, "file", None)
        if close is not None:
            close.close()


def load_obs(path: str | Path) -> pd.DataFrame:
    """Load only `.obs` from an H5AD file."""

    adata = read_h5ad_backed(path)
    try:
        return adata.obs.copy()
    finally:
        close = getattr(adata, "file", None)
        if close is not None:
            close.close()


def _columns_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col in frame.columns:
        series = frame[col]
        rows.append(
            {
                "column": col,
                "dtype": str(series.dtype),
                "non_null": int(series.notna().sum()),
                "unique": int(series.nunique(dropna=True)),
                "example": _example_value(series),
            }
        )
    return pd.DataFrame(rows)


def _example_value(series: pd.Series) -> str:
    non_null = series.dropna()
    if non_null.empty:
        return ""
    return str(non_null.iloc[0])

