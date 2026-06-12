"""Targeted pseudobulk aggregation and lightweight DE tests."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .io import read_h5ad_backed
from .metadata import resolve_columns
from .modules import (
    DEFAULT_SYMBOL_COLUMNS,
    GeneSet,
    build_qc_mask,
    build_score_metadata,
    parse_gene_sets,
    resolve_gene_sets,
)
from .stats import bh_fdr
from .utils import normalize_token


DEFAULT_GROUPBY = ("donor_id", "sample_id", "disease_group", "coarse_cell_type", "fine_cell_type")
DEFAULT_CONTRASTS = ("ad:healthy", "pd:healthy")


@dataclass
class PseudobulkResult:
    counts: pd.DataFrame
    metadata: pd.DataFrame
    coverage: pd.DataFrame
    de: pd.DataFrame


def genes_from_gene_sets(gene_set_config: dict[str, Any]) -> list[str]:
    """Return unique genes from configured gene sets, preserving config order."""

    genes: list[str] = []
    seen: set[str] = set()
    for gene_set in parse_gene_sets(gene_set_config):
        for gene in gene_set.genes:
            token = normalize_token(gene)
            if token and token not in seen:
                genes.append(gene)
                seen.add(token)
    return genes


def parse_contrasts(values: Iterable[str]) -> list[tuple[str, str]]:
    """Parse contrast strings in case:control form."""

    contrasts = []
    for value in values:
        if ":" not in value:
            raise ValueError(f"Contrast must use case:control syntax: {value}")
        case, control = value.split(":", 1)
        contrasts.append((normalize_token(case).replace(" ", "_"), normalize_token(control).replace(" ", "_")))
    return contrasts


def aggregate_targeted_pseudobulk_h5ad(
    h5ad_path: str | Path,
    gateway_config: dict[str, Any],
    genes: list[str],
    *,
    groupby: list[str] | tuple[str, ...] = DEFAULT_GROUPBY,
    chunk_size: int = 1_000,
    layer: str | None = None,
    apply_qc: bool = False,
    symbol_columns: list[str] | tuple[str, ...] = DEFAULT_SYMBOL_COLUMNS,
    contrasts: list[tuple[str, str]] | None = None,
    min_donors: int = 3,
) -> PseudobulkResult:
    """Aggregate selected genes to pseudobulk groups and run simple contrasts."""

    if not genes:
        raise ValueError("At least one gene is required for targeted pseudobulk aggregation.")
    adata = read_h5ad_backed(h5ad_path)
    try:
        columns = resolve_columns(list(adata.obs.columns), gateway_config)
        gene_set = GeneSet(name="targeted_pseudobulk", genes=tuple(genes))
        resolved, coverage = resolve_gene_sets(adata.var, adata.var_names, [gene_set], symbol_columns)
        selected_indices = resolved["targeted_pseudobulk"]
        present_genes = [gene for gene in genes if normalize_token(gene) in _present_gene_tokens(coverage)]
        if not selected_indices:
            empty_counts = pd.DataFrame(columns=[*groupby, "gene", "count"])
            empty_meta = pd.DataFrame(columns=[*groupby, "n_cells", "sum_n_counts"])
            empty_de = pd.DataFrame()
            return PseudobulkResult(empty_counts, empty_meta, coverage, empty_de)

        metadata = build_score_metadata(adata.obs, columns, gateway_config)
        metadata["n_counts"] = (
            pd.to_numeric(adata.obs[columns.n_counts], errors="coerce").fillna(0).to_numpy(dtype=float)
            if columns.n_counts
            else 0.0
        )
        missing = [col for col in groupby if col not in metadata.columns]
        if missing:
            raise KeyError(f"Unknown pseudobulk grouping columns: {missing}")
        qc_mask = build_qc_mask(adata.obs, columns, gateway_config) if apply_qc else np.ones(adata.n_obs, dtype=bool)
        group_frame = _group_frame(metadata, list(groupby))
        group_ids, group_meta = _group_ids_and_metadata(group_frame, metadata["n_counts"], list(groupby))
        counts = _aggregate_counts(
            h5ad_path=Path(h5ad_path),
            adata=adata,
            selected_indices=selected_indices,
            group_ids=group_ids,
            qc_mask=qc_mask,
            chunk_size=chunk_size,
            layer=layer,
        )
        gene_symbols = _present_genes_in_index_order(genes, coverage)
        counts_long = counts_to_long(counts, group_meta, gene_symbols, list(groupby))
        de = run_pseudobulk_de(
            counts,
            group_meta,
            gene_symbols,
            contrasts=contrasts or parse_contrasts(DEFAULT_CONTRASTS),
            min_donors=min_donors,
        )
        return PseudobulkResult(counts=counts_long, metadata=group_meta, coverage=coverage, de=de)
    finally:
        close = getattr(adata, "file", None)
        if close is not None:
            close.close()


def counts_to_long(
    counts: np.ndarray,
    group_meta: pd.DataFrame,
    genes: list[str],
    groupby: list[str],
) -> pd.DataFrame:
    """Convert group x gene counts to a compact long table."""

    rows = []
    for gene_idx, gene in enumerate(genes):
        values = counts[:, gene_idx]
        nonzero = values > 0
        if not nonzero.any():
            continue
        frame = group_meta.loc[nonzero, groupby].copy()
        frame["gene"] = gene
        frame["count"] = values[nonzero]
        rows.append(frame)
    if not rows:
        return pd.DataFrame(columns=[*groupby, "gene", "count"])
    return pd.concat(rows, ignore_index=True)


def run_pseudobulk_de(
    counts: np.ndarray,
    group_meta: pd.DataFrame,
    genes: list[str],
    *,
    contrasts: list[tuple[str, str]] | None = None,
    min_donors: int = 3,
    pseudocount: float = 0.5,
) -> pd.DataFrame:
    """Run lightweight Welch tests on donor-level logCPM pseudobulk values."""

    contrasts = contrasts or parse_contrasts(DEFAULT_CONTRASTS)
    if counts.size == 0:
        return pd.DataFrame()
    gene_cols = [f"__gene_{idx}" for idx in range(len(genes))]
    wide = pd.DataFrame(counts, columns=gene_cols)
    meta_cols = ["donor_id", "disease_group", "fine_cell_type", "n_cells", "sum_n_counts"]
    work = pd.concat([group_meta[meta_cols].reset_index(drop=True), wide], axis=1)
    grouped = (
        work.groupby(["donor_id", "disease_group", "fine_cell_type"], observed=True, dropna=False)
        .agg({**{col: "sum" for col in gene_cols}, "n_cells": "sum", "sum_n_counts": "sum"})
        .reset_index()
    )
    library = grouped["sum_n_counts"].to_numpy(dtype=float)
    selected_library = grouped[gene_cols].sum(axis=1).to_numpy(dtype=float)
    library = np.where(library > 0, library, selected_library)
    denominator = library[:, None] + pseudocount * len(genes)
    logcpm = np.log2(((grouped[gene_cols].to_numpy(dtype=float) + pseudocount) / denominator) * 1_000_000 + 1)

    rows = []
    for case, control in contrasts:
        for cell_state, state_idx in grouped.groupby("fine_cell_type", observed=True).groups.items():
            idx = np.asarray(list(state_idx), dtype=int)
            disease = grouped.iloc[idx]["disease_group"].astype(str).map(lambda x: normalize_token(x).replace(" ", "_"))
            case_mask = disease.eq(case).to_numpy()
            control_mask = disease.eq(control).to_numpy()
            n_case = int(grouped.iloc[idx[case_mask]]["donor_id"].nunique())
            n_control = int(grouped.iloc[idx[control_mask]]["donor_id"].nunique())
            for gene_idx, gene in enumerate(genes):
                if n_case < min_donors or n_control < min_donors:
                    rows.append(_empty_de_row(case, control, cell_state, gene, n_case, n_control, "too_few_donors"))
                    continue
                case_values = logcpm[idx[case_mask], gene_idx]
                control_values = logcpm[idx[control_mask], gene_idx]
                t_stat, p_value = _welch_ttest(case_values, control_values)
                rows.append(
                    {
                        "contrast": f"{case}_vs_{control}",
                        "case_group": case,
                        "control_group": control,
                        "fine_cell_type": cell_state,
                        "gene": gene,
                        "n_case": n_case,
                        "n_control": n_control,
                        "mean_logcpm_case": float(np.mean(case_values)),
                        "mean_logcpm_control": float(np.mean(control_values)),
                        "log2fc": float(np.mean(case_values) - np.mean(control_values)),
                        "t_stat": t_stat,
                        "p_value": p_value,
                        "status": "ok",
                    }
                )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr"] = np.nan
    ok = result["status"].eq("ok")
    for _, idx in result[ok].groupby(["contrast", "fine_cell_type"], observed=True).groups.items():
        idx_array = np.asarray(list(idx), dtype=int)
        result.loc[idx_array, "fdr"] = bh_fdr(result.loc[idx_array, "p_value"].to_numpy(dtype=float))
    return result.sort_values(["status", "contrast", "fine_cell_type", "fdr", "p_value", "gene"], na_position="last")


def _aggregate_counts(
    *,
    h5ad_path: Path,
    adata: Any,
    selected_indices: list[int],
    group_ids: np.ndarray,
    qc_mask: np.ndarray,
    chunk_size: int,
    layer: str | None,
) -> np.ndarray:
    n_groups = int(group_ids.max()) + 1 if group_ids.size else 0
    counts = np.zeros((n_groups, len(selected_indices)), dtype=np.float64)
    if layer is None and _h5ad_x_is_csr(h5ad_path):
        _aggregate_csr_counts(h5ad_path, adata.n_obs, adata.n_vars, selected_indices, group_ids, qc_mask, chunk_size, counts)
    else:
        matrix = adata.X if layer is None else adata.layers[layer]
        _aggregate_matrix_counts(matrix, adata.n_obs, selected_indices, group_ids, qc_mask, chunk_size, counts)
    return counts


def _aggregate_csr_counts(
    h5ad_path: Path,
    n_obs: int,
    n_vars: int,
    selected_indices: list[int],
    group_ids: np.ndarray,
    qc_mask: np.ndarray,
    chunk_size: int,
    counts: np.ndarray,
) -> None:
    import h5py

    gene_pos = np.full(n_vars, -1, dtype=np.int64)
    gene_pos[np.asarray(selected_indices, dtype=np.int64)] = np.arange(len(selected_indices), dtype=np.int64)
    with h5py.File(h5ad_path, "r") as handle:
        x_group = handle["X"]
        data_ds = x_group["data"]
        indices_ds = x_group["indices"]
        indptr_ds = x_group["indptr"]
        for start in range(0, n_obs, chunk_size):
            stop = min(start + chunk_size, n_obs)
            keep = qc_mask[start:stop]
            if not keep.any():
                continue
            indptr = np.asarray(indptr_ds[start : stop + 1], dtype=np.int64)
            data_start = int(indptr[0])
            data_stop = int(indptr[-1])
            if data_stop == data_start:
                continue
            local_indptr = indptr - data_start
            values = np.asarray(data_ds[data_start:data_stop], dtype=np.float64)
            gene_indices = np.asarray(indices_ds[data_start:data_stop], dtype=np.int64)
            positions = gene_pos[gene_indices]
            selected = positions >= 0
            if not selected.any():
                continue
            offsets = np.flatnonzero(selected)
            row_ids = np.searchsorted(local_indptr[1:], offsets, side="right")
            kept_rows = keep[row_ids]
            if not kept_rows.any():
                continue
            global_groups = group_ids[start + row_ids[kept_rows]]
            np.add.at(counts, (global_groups, positions[selected][kept_rows]), values[selected][kept_rows])


def _aggregate_matrix_counts(
    matrix: Any,
    n_obs: int,
    selected_indices: list[int],
    group_ids: np.ndarray,
    qc_mask: np.ndarray,
    chunk_size: int,
    counts: np.ndarray,
) -> None:
    for start in range(0, n_obs, chunk_size):
        stop = min(start + chunk_size, n_obs)
        keep = qc_mask[start:stop]
        if not keep.any():
            continue
        x = _as_dense(matrix[start:stop, selected_indices])
        groups = group_ids[start:stop][keep]
        x = x[keep, :]
        for gene_idx in range(x.shape[1]):
            np.add.at(counts[:, gene_idx], groups, x[:, gene_idx])


def _group_frame(metadata: pd.DataFrame, groupby: list[str]) -> pd.DataFrame:
    frame = metadata[groupby].copy()
    for col in groupby:
        frame[col] = frame[col].astype("string").fillna("unknown").replace("", "unknown")
    return frame


def _group_ids_and_metadata(group_frame: pd.DataFrame, n_counts: pd.Series, groupby: list[str]) -> tuple[np.ndarray, pd.DataFrame]:
    grouped = group_frame.groupby(groupby, observed=True, sort=True, dropna=False)
    group_ids = grouped.ngroup().to_numpy(dtype=np.int64)
    meta = grouped.size().rename("n_cells").reset_index()
    sum_n_counts = n_counts.groupby(group_ids).sum().rename("sum_n_counts").reset_index(drop=True)
    meta["sum_n_counts"] = sum_n_counts
    return group_ids, meta


def _present_gene_tokens(coverage: pd.DataFrame) -> set[str]:
    if coverage.empty or "present_genes" not in coverage:
        return set()
    genes = str(coverage.loc[0, "present_genes"]).split(",")
    return {normalize_token(gene) for gene in genes if gene and gene != "nan"}


def _present_genes_in_index_order(genes: list[str], coverage: pd.DataFrame) -> list[str]:
    present = _present_gene_tokens(coverage)
    return [gene for gene in genes if normalize_token(gene) in present]


def _welch_ttest(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    try:
        from scipy import stats

        result = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
        return float(result.statistic), float(result.pvalue)
    except ModuleNotFoundError:
        return np.nan, np.nan


def _empty_de_row(case: str, control: str, cell_state: str, gene: str, n_case: int, n_control: int, status: str) -> dict[str, object]:
    return {
        "contrast": f"{case}_vs_{control}",
        "case_group": case,
        "control_group": control,
        "fine_cell_type": cell_state,
        "gene": gene,
        "n_case": n_case,
        "n_control": n_control,
        "mean_logcpm_case": np.nan,
        "mean_logcpm_control": np.nan,
        "log2fc": np.nan,
        "t_stat": np.nan,
        "p_value": np.nan,
        "status": status,
    }


def _h5ad_x_is_csr(path: str | Path) -> bool:
    try:
        import h5py
    except ModuleNotFoundError:
        return False
    with h5py.File(path, "r") as handle:
        if "X" not in handle:
            return False
        x = handle["X"]
        encoding = x.attrs.get("encoding-type")
        if isinstance(encoding, bytes):
            encoding = encoding.decode("utf-8")
        return str(encoding) == "csr_matrix" and all(key in x for key in ["data", "indices", "indptr"])


def _as_dense(matrix: Any) -> np.ndarray:
    try:
        from scipy import sparse

        if sparse.issparse(matrix):
            matrix = matrix.toarray()
    except ModuleNotFoundError:
        pass
    array = np.asarray(matrix, dtype=float)
    if array.ndim == 1:
        array = array.reshape(-1, 1)
    return array
