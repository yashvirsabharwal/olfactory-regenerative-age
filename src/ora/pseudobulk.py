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
DEFAULT_COVARIATES = ("age", "sex", "chemistry", "collection_method")


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


def run_covariate_pseudobulk_de(
    counts_long: pd.DataFrame,
    metadata: pd.DataFrame,
    manifest: pd.DataFrame,
    *,
    genes: list[str] | None = None,
    contrasts: list[tuple[str, str]] | None = None,
    covariates: list[str] | tuple[str, ...] = DEFAULT_COVARIATES,
    min_donors: int = 3,
    pseudocount: float = 0.5,
) -> pd.DataFrame:
    """Run donor-level covariate-adjusted linear models on targeted pseudobulk logCPM."""

    contrasts = contrasts or parse_contrasts(DEFAULT_CONTRASTS)
    if metadata.empty:
        return pd.DataFrame()
    genes = genes or sorted(counts_long["gene"].dropna().astype(str).unique().tolist())
    if not genes:
        return pd.DataFrame()
    donor_state = _donor_state_table(counts_long, metadata, manifest, genes)
    if donor_state.empty:
        return pd.DataFrame()
    gene_cols = [f"__gene_{gene}" for gene in genes]
    library = donor_state["sum_n_counts"].to_numpy(dtype=float)
    selected_library = donor_state[gene_cols].sum(axis=1).to_numpy(dtype=float)
    library = np.where(library > 0, library, selected_library)
    logcpm = np.log2(((donor_state[gene_cols].to_numpy(dtype=float) + pseudocount) / (library[:, None] + pseudocount * len(genes))) * 1_000_000 + 1)

    rows = []
    for case, control in contrasts:
        for cell_state, state_idx in donor_state.groupby("fine_cell_type", observed=True).groups.items():
            idx = np.asarray(list(state_idx), dtype=int)
            state = donor_state.iloc[idx].reset_index(drop=True)
            disease = state["disease_group"].astype(str).map(lambda x: normalize_token(x).replace(" ", "_"))
            contrast_mask = disease.isin([case, control]).to_numpy()
            if not contrast_mask.any():
                continue
            state = state.loc[contrast_mask].reset_index(drop=True)
            y_state = logcpm[idx[contrast_mask], :]
            disease = disease[contrast_mask].reset_index(drop=True)
            case_mask_all = disease.eq(case).to_numpy()
            control_mask_all = disease.eq(control).to_numpy()
            for gene_idx, gene in enumerate(genes):
                y = y_state[:, gene_idx]
                valid_y = np.isfinite(y)
                frame = state.loc[valid_y].reset_index(drop=True)
                y_valid = y[valid_y]
                disease_valid = disease[valid_y].reset_index(drop=True)
                case_mask = disease_valid.eq(case).to_numpy()
                control_mask = disease_valid.eq(control).to_numpy()
                n_case = int(frame.loc[case_mask, "donor_id"].nunique())
                n_control = int(frame.loc[control_mask, "donor_id"].nunique())
                if n_case < min_donors or n_control < min_donors:
                    rows.append(_empty_adjusted_de_row(case, control, cell_state, gene, n_case, n_control, "too_few_donors"))
                    continue
                design, used_covariates, design_valid = _covariate_design(frame, disease_valid, case, covariates)
                if design.shape[1] < 2 or not design_valid.any():
                    rows.append(_empty_adjusted_de_row(case, control, cell_state, gene, n_case, n_control, "invalid_design"))
                    continue
                y_model = y_valid[design_valid]
                frame_model = frame.loc[design_valid].reset_index(drop=True)
                disease_model = disease_valid[design_valid].reset_index(drop=True)
                n_case_model = int(frame_model.loc[disease_model.eq(case), "donor_id"].nunique())
                n_control_model = int(frame_model.loc[disease_model.eq(control), "donor_id"].nunique())
                if n_case_model < min_donors or n_control_model < min_donors:
                    rows.append(
                        _empty_adjusted_de_row(
                            case,
                            control,
                            cell_state,
                            gene,
                            n_case_model,
                            n_control_model,
                            "too_few_covariate_complete_donors",
                        )
                    )
                    continue
                fit = _ols_case_effect(y_model, design[design_valid, :])
                if fit is None:
                    rows.append(_empty_adjusted_de_row(case, control, cell_state, gene, n_case_model, n_control_model, "insufficient_df"))
                    continue
                case_values = y_model[disease_model.eq(case).to_numpy()]
                control_values = y_model[disease_model.eq(control).to_numpy()]
                rows.append(
                    {
                        "contrast": f"{case}_vs_{control}",
                        "case_group": case,
                        "control_group": control,
                        "fine_cell_type": cell_state,
                        "gene": gene,
                        "n_case": n_case_model,
                        "n_control": n_control_model,
                        "n_total": int(y_model.size),
                        "mean_logcpm_case": float(np.mean(case_values)),
                        "mean_logcpm_control": float(np.mean(control_values)),
                        "log2fc_unadjusted": float(np.mean(case_values) - np.mean(control_values)),
                        "log2fc_adjusted": fit["beta"],
                        "t_stat": fit["t_stat"],
                        "p_value": fit["p_value"],
                        "df_resid": fit["df_resid"],
                        "covariates": ",".join(used_covariates),
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


def _empty_adjusted_de_row(case: str, control: str, cell_state: str, gene: str, n_case: int, n_control: int, status: str) -> dict[str, object]:
    return {
        "contrast": f"{case}_vs_{control}",
        "case_group": case,
        "control_group": control,
        "fine_cell_type": cell_state,
        "gene": gene,
        "n_case": n_case,
        "n_control": n_control,
        "n_total": n_case + n_control,
        "mean_logcpm_case": np.nan,
        "mean_logcpm_control": np.nan,
        "log2fc_unadjusted": np.nan,
        "log2fc_adjusted": np.nan,
        "t_stat": np.nan,
        "p_value": np.nan,
        "df_resid": np.nan,
        "covariates": "",
        "status": status,
    }


def _donor_state_table(
    counts_long: pd.DataFrame,
    metadata: pd.DataFrame,
    manifest: pd.DataFrame,
    genes: list[str],
) -> pd.DataFrame:
    meta_required = {"donor_id", "disease_group", "fine_cell_type", "n_cells", "sum_n_counts"}
    if not meta_required.issubset(metadata.columns):
        missing = sorted(meta_required.difference(metadata.columns))
        raise KeyError(f"Pseudobulk metadata is missing required columns: {missing}")
    donor_state = (
        metadata.groupby(["donor_id", "disease_group", "fine_cell_type"], observed=True, dropna=False)
        .agg(n_cells=("n_cells", "sum"), sum_n_counts=("sum_n_counts", "sum"))
        .reset_index()
    )
    gene_cols = [f"__gene_{gene}" for gene in genes]
    for col in gene_cols:
        donor_state[col] = 0.0
    if not counts_long.empty:
        count_required = {"donor_id", "disease_group", "fine_cell_type", "gene", "count"}
        if not count_required.issubset(counts_long.columns):
            missing = sorted(count_required.difference(counts_long.columns))
            raise KeyError(f"Pseudobulk counts are missing required columns: {missing}")
        counts = (
            counts_long[counts_long["gene"].isin(genes)]
            .groupby(["donor_id", "disease_group", "fine_cell_type", "gene"], observed=True, dropna=False)["count"]
            .sum()
            .reset_index()
        )
        if not counts.empty:
            wide = counts.pivot_table(
                index=["donor_id", "disease_group", "fine_cell_type"],
                columns="gene",
                values="count",
                fill_value=0,
                aggfunc="sum",
            )
            wide.columns = [f"__gene_{gene}" for gene in wide.columns]
            wide = wide.reset_index()
            donor_state = donor_state.drop(columns=gene_cols).merge(
                wide,
                on=["donor_id", "disease_group", "fine_cell_type"],
                how="left",
            )
            for col in gene_cols:
                if col not in donor_state:
                    donor_state[col] = 0.0
                donor_state[col] = pd.to_numeric(donor_state[col], errors="coerce").fillna(0.0)
    manifest_cols = ["donor_id", "age", "sex", "race_ethnicity", "chemistry", "collection_method", "site"]
    available = [col for col in manifest_cols if col in manifest.columns]
    donor_meta = manifest[available].drop_duplicates("donor_id") if available else pd.DataFrame({"donor_id": []})
    if not donor_meta.empty:
        donor_state = donor_state.merge(donor_meta, on="donor_id", how="left", suffixes=("", "_manifest"))
    donor_state["disease_group"] = donor_state["disease_group"].astype(str).map(lambda x: normalize_token(x).replace(" ", "_"))
    return donor_state


def _covariate_design(
    frame: pd.DataFrame,
    disease: pd.Series,
    case: str,
    covariates: list[str] | tuple[str, ...],
) -> tuple[np.ndarray, list[str], np.ndarray]:
    pieces = [
        pd.Series(1.0, index=frame.index, name="intercept"),
        disease.eq(case).astype(float).rename("case"),
    ]
    used: list[str] = []
    for covariate in covariates:
        if covariate not in frame:
            continue
        values = frame[covariate]
        if pd.api.types.is_numeric_dtype(values):
            numeric = pd.to_numeric(values, errors="coerce")
            if numeric.notna().sum() < 2 or numeric.nunique(dropna=True) < 2:
                continue
            pieces.append(numeric.rename(covariate))
            used.append(covariate)
            continue
        categorical = values.fillna("unknown").astype(str)
        if categorical.nunique(dropna=False) < 2:
            continue
        dummies = pd.get_dummies(categorical, prefix=covariate, drop_first=True, dtype=float)
        if dummies.empty:
            continue
        pieces.append(dummies)
        used.append(covariate)
    design_frame = pd.concat(pieces, axis=1)
    design_frame = design_frame.loc[:, design_frame.nunique(dropna=False) > 1]
    if "case" not in design_frame:
        return np.empty((frame.shape[0], 0)), used, np.zeros(frame.shape[0], dtype=bool)
    if "intercept" not in design_frame:
        design_frame.insert(0, "intercept", 1.0)
    valid = np.isfinite(design_frame.to_numpy(dtype=float)).all(axis=1)
    return design_frame.to_numpy(dtype=float), used, valid


def _ols_case_effect(y: np.ndarray, design: np.ndarray) -> dict[str, float] | None:
    if design.shape[0] <= design.shape[1]:
        return None
    rank = int(np.linalg.matrix_rank(design))
    if rank < design.shape[1]:
        return None
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ coef
    df_resid = int(design.shape[0] - rank)
    if df_resid <= 0:
        return None
    rss = float(np.sum(resid**2))
    sigma2 = rss / df_resid
    cov = sigma2 * np.linalg.pinv(design.T @ design)
    se = float(np.sqrt(max(cov[1, 1], 0.0)))
    if se == 0 or not np.isfinite(se):
        return None
    beta = float(coef[1])
    t_stat = beta / se
    p_value = _t_two_sided_p_value(t_stat, df_resid)
    return {"beta": beta, "t_stat": float(t_stat), "p_value": float(p_value), "df_resid": float(df_resid)}


def _t_two_sided_p_value(t_stat: float, df_resid: int) -> float:
    try:
        from scipy import stats

        return float(2 * stats.t.sf(abs(t_stat), df=df_resid))
    except ModuleNotFoundError:
        return float(np.nan)


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
