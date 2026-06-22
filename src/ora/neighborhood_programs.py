"""Program scoring for Milo-style neighborhood memberships."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from ora.io import read_h5ad_backed
from ora.modules import DEFAULT_SYMBOL_COLUMNS, parse_gene_sets, resolve_gene_sets


def score_neighborhood_programs_h5ad(
    h5ad_path: str | Path,
    gene_set_config: Mapping[str, Any],
    memberships: pd.DataFrame,
    *,
    da_table: pd.DataFrame | None = None,
    run_name: str = "neighborhood_run",
    chunk_neighborhoods: int = 500,
    log1p: bool | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score configured gene programs in each neighborhood from an H5AD matrix."""

    if memberships.empty:
        empty_scores = _score_columns()
        return pd.DataFrame(columns=empty_scores), _empty_summary(), pd.DataFrame()
    if "neighborhood_id" not in memberships or "cell_index" not in memberships:
        raise ValueError("memberships must contain `neighborhood_id` and `cell_index` columns")

    score_config = dict(gene_set_config.get("score", {}))
    log1p = bool(score_config.get("log1p", True) if log1p is None else log1p)
    symbol_columns = tuple(score_config.get("var_symbol_columns", DEFAULT_SYMBOL_COLUMNS))
    gene_sets = parse_gene_sets(dict(gene_set_config))
    if not gene_sets:
        raise ValueError("No gene sets were configured.")

    adata = read_h5ad_backed(h5ad_path)
    try:
        resolved, coverage = resolve_gene_sets(adata.var, adata.var_names, gene_sets, symbol_columns)
        module_indices = {name: idxs for name, idxs in resolved.items() if idxs}
        if not module_indices:
            return pd.DataFrame(columns=_score_columns()), _empty_summary(), coverage
        scores = _score_memberships(
            adata,
            memberships,
            module_indices=module_indices,
            run_name=run_name,
            chunk_neighborhoods=chunk_neighborhoods,
            log1p=log1p,
        )
    finally:
        close = getattr(adata, "file", None)
        if close is not None:
            close.close()

    scores = _add_program_zscores(scores)
    if da_table is not None and not da_table.empty:
        scores = _merge_da(scores, da_table)
    summary = summarize_neighborhood_programs(scores)
    return scores, summary, coverage


def summarize_neighborhood_programs(scores: pd.DataFrame) -> pd.DataFrame:
    """Summarize program scores across all and age-associated neighborhoods."""

    columns = [
        "run",
        "module",
        "n_neighborhoods",
        "n_significant",
        "median_score",
        "median_z",
        "significant_median_z",
        "negative_significant_median_z",
        "positive_significant_median_z",
        "top_enriched_neighborhoods",
        "top_depleted_neighborhoods",
    ]
    if scores.empty:
        return pd.DataFrame(columns=columns)

    frame = scores.copy()
    if "is_age_associated_fdr_0_10" not in frame:
        frame["is_age_associated_fdr_0_10"] = False
    if "age_direction" not in frame:
        frame["age_direction"] = "not_available"

    rows = []
    for (run, module), group in frame.groupby(["run", "module"], observed=True):
        sig = group[group["is_age_associated_fdr_0_10"]]
        neg = sig[sig["age_direction"].eq("negative")]
        pos = sig[sig["age_direction"].eq("positive")]
        rows.append(
            {
                "run": run,
                "module": module,
                "n_neighborhoods": int(group.shape[0]),
                "n_significant": int(sig.shape[0]),
                "median_score": float(group["program_score"].median()),
                "median_z": float(group["program_z"].median()),
                "significant_median_z": _median_or_nan(sig["program_z"]),
                "negative_significant_median_z": _median_or_nan(neg["program_z"]),
                "positive_significant_median_z": _median_or_nan(pos["program_z"]),
                "top_enriched_neighborhoods": _top_neighborhoods(group, ascending=False),
                "top_depleted_neighborhoods": _top_neighborhoods(group, ascending=True),
            }
        )
    result = pd.DataFrame(rows, columns=columns)
    return result.sort_values(["run", "n_significant", "significant_median_z"], ascending=[True, False, False]).reset_index(drop=True)


def _score_memberships(
    adata: Any,
    memberships: pd.DataFrame,
    *,
    module_indices: dict[str, list[int]],
    run_name: str,
    chunk_neighborhoods: int,
    log1p: bool,
) -> pd.DataFrame:
    union_gene_indices = sorted({idx for indices in module_indices.values() for idx in indices})
    union_lookup = {idx: pos for pos, idx in enumerate(union_gene_indices)}
    module_positions = {
        module: [union_lookup[idx] for idx in indices if idx in union_lookup]
        for module, indices in module_indices.items()
    }
    neighborhood_ids = pd.Index(pd.unique(memberships["neighborhood_id"])).sort_values()
    chunks = [neighborhood_ids[start : start + chunk_neighborhoods] for start in range(0, len(neighborhood_ids), chunk_neighborhoods)]
    rows = []
    for chunk_ids in chunks:
        chunk_members = memberships[memberships["neighborhood_id"].isin(chunk_ids)].copy()
        cell_indices = pd.to_numeric(chunk_members["cell_index"], errors="coerce").astype("Int64")
        chunk_members = chunk_members[cell_indices.notna()].copy()
        chunk_members["cell_index"] = cell_indices[cell_indices.notna()].astype(int).to_numpy()
        if chunk_members.empty:
            continue
        unique_cells, inverse = np.unique(chunk_members["cell_index"].to_numpy(dtype=int), return_inverse=True)
        matrix = _read_dense_block(adata, unique_cells, union_gene_indices)
        if log1p:
            matrix = np.log1p(matrix)
        cell_module_scores = {
            module: matrix[:, positions].mean(axis=1) if positions else np.full(matrix.shape[0], np.nan)
            for module, positions in module_positions.items()
        }
        for module, values in cell_module_scores.items():
            chunk_members["program_score"] = values[inverse]
            grouped = chunk_members.groupby("neighborhood_id", observed=True)["program_score"].agg(["mean", "count"]).reset_index()
            grouped["run"] = run_name
            grouped["module"] = module
            grouped = grouped.rename(columns={"mean": "program_score", "count": "n_membership_cells"})
            rows.append(grouped[["run", "neighborhood_id", "module", "n_membership_cells", "program_score"]])
    if not rows:
        return pd.DataFrame(columns=_score_columns())
    return pd.concat(rows, ignore_index=True)


def _read_dense_block(adata: Any, cell_indices: np.ndarray, gene_indices: list[int]) -> np.ndarray:
    block = adata[cell_indices, gene_indices].X
    if sparse.issparse(block):
        return block.toarray().astype(float, copy=False)
    return np.asarray(block, dtype=float)


def _add_program_zscores(scores: pd.DataFrame) -> pd.DataFrame:
    frame = scores.copy()
    frame["program_z"] = np.nan
    for (_, module), idx in frame.groupby(["run", "module"], observed=True).groups.items():
        values = frame.loc[idx, "program_score"].to_numpy(dtype=float)
        sd = float(np.nanstd(values))
        if sd > 0 and np.isfinite(sd):
            frame.loc[idx, "program_z"] = (values - float(np.nanmean(values))) / sd
        else:
            frame.loc[idx, "program_z"] = 0.0
    return frame


def _merge_da(scores: pd.DataFrame, da_table: pd.DataFrame) -> pd.DataFrame:
    da = da_table.copy()
    keep_cols = [
        "neighborhood_id",
        "top_fine_celltype",
        "top_coarse_celltype",
        "age_coef",
        "age_fdr",
        "status",
    ]
    da = da[[col for col in keep_cols if col in da.columns]].copy()
    merged = scores.merge(da, on="neighborhood_id", how="left")
    merged["age_coef"] = pd.to_numeric(merged.get("age_coef"), errors="coerce")
    merged["age_fdr"] = pd.to_numeric(merged.get("age_fdr"), errors="coerce")
    merged["is_age_associated_fdr_0_10"] = merged["age_fdr"].lt(0.10) & merged.get("status", "").astype(str).eq("tested")
    merged["age_direction"] = np.where(merged["age_coef"].gt(0), "positive", np.where(merged["age_coef"].lt(0), "negative", "zero"))
    return merged


def _top_neighborhoods(group: pd.DataFrame, *, ascending: bool, limit: int = 5) -> str:
    ordered = group.sort_values("program_z", ascending=ascending).head(limit)
    return ";".join(f"{int(row.neighborhood_id)}:{row.program_z:.3g}" for row in ordered.itertuples())


def _median_or_nan(values: pd.Series) -> float:
    if values.empty:
        return float("nan")
    return float(values.median())


def _score_columns() -> list[str]:
    return ["run", "neighborhood_id", "module", "n_membership_cells", "program_score", "program_z"]


def _empty_summary() -> pd.DataFrame:
    return pd.DataFrame(columns=summarize_neighborhood_programs(pd.DataFrame(columns=_score_columns())).columns)
