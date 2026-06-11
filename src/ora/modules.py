"""Chunked gene-set/module scoring for AnnData matrices."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .io import read_h5ad_backed
from .metadata import ColumnMap, collection_method_group, disease_group, parse_age_series, resolve_columns
from .utils import normalize_token


DEFAULT_SYMBOL_COLUMNS = ("feature_name", "gene_symbol", "gene_name", "symbol")
DEFAULT_GROUPBY = ("donor_id", "sample_id", "coarse_cell_type", "fine_cell_type")


@dataclass(frozen=True)
class GeneSet:
    name: str
    genes: tuple[str, ...]
    description: str = ""


@dataclass
class ModuleScoreResult:
    summary: pd.DataFrame
    donor_features: pd.DataFrame
    coverage: pd.DataFrame


def parse_gene_sets(config: dict[str, Any]) -> list[GeneSet]:
    """Parse flexible gene-set config entries."""

    parsed: list[GeneSet] = []
    for name, spec in config.get("gene_sets", {}).items():
        description = ""
        if isinstance(spec, dict):
            genes = spec.get("genes", [])
            description = str(spec.get("description", ""))
        else:
            genes = spec
        genes_tuple = tuple(str(gene).strip() for gene in genes if str(gene).strip())
        if genes_tuple:
            parsed.append(GeneSet(name=str(name), genes=genes_tuple, description=description))
    return parsed


def resolve_gene_sets(
    var: pd.DataFrame,
    var_names: pd.Index,
    gene_sets: list[GeneSet],
    symbol_columns: list[str] | tuple[str, ...] = DEFAULT_SYMBOL_COLUMNS,
) -> tuple[dict[str, list[int]], pd.DataFrame]:
    """Resolve requested gene symbols to matrix column indices."""

    lookup: dict[str, int] = {}
    for idx, value in enumerate(var_names):
        token = normalize_token(value)
        if token and token not in lookup:
            lookup[token] = idx
    for col in symbol_columns:
        if col not in var.columns:
            continue
        for idx, value in enumerate(var[col].astype(str)):
            token = normalize_token(value)
            if token and token not in lookup:
                lookup[token] = idx

    resolved: dict[str, list[int]] = {}
    rows = []
    for gene_set in gene_sets:
        present_genes = []
        missing_genes = []
        indices = []
        seen_indices: set[int] = set()
        for gene in gene_set.genes:
            token = normalize_token(gene)
            if token in lookup:
                idx = lookup[token]
                present_genes.append(gene)
                if idx not in seen_indices:
                    indices.append(idx)
                    seen_indices.add(idx)
            else:
                missing_genes.append(gene)
        resolved[gene_set.name] = indices
        n_requested = len(gene_set.genes)
        n_present = len(present_genes)
        rows.append(
            {
                "module": gene_set.name,
                "description": gene_set.description,
                "n_requested": n_requested,
                "n_present": n_present,
                "coverage_fraction": n_present / n_requested if n_requested else np.nan,
                "present_genes": ",".join(present_genes),
                "missing_genes": ",".join(missing_genes),
            }
        )
    return resolved, pd.DataFrame(rows)


def score_gene_sets_h5ad(
    h5ad_path: str | Path,
    gateway_config: dict[str, Any],
    gene_set_config: dict[str, Any],
    *,
    groupby: list[str] | tuple[str, ...] = DEFAULT_GROUPBY,
    chunk_size: int | None = None,
    layer: str | None = None,
    log1p: bool | None = None,
    apply_qc: bool = False,
) -> ModuleScoreResult:
    """Score configured gene sets and aggregate scores without loading full X."""

    score_config = gene_set_config.get("score", {})
    chunk_size = int(chunk_size or score_config.get("chunk_size", 1_000))
    log1p = bool(score_config.get("log1p", True) if log1p is None else log1p)
    symbol_columns = list(score_config.get("var_symbol_columns", DEFAULT_SYMBOL_COLUMNS))
    gene_sets = parse_gene_sets(gene_set_config)
    if not gene_sets:
        raise ValueError("No gene sets were configured.")

    adata = read_h5ad_backed(h5ad_path)
    try:
        columns = resolve_columns(list(adata.obs.columns), gateway_config)
        resolved, coverage = resolve_gene_sets(adata.var, adata.var_names, gene_sets, symbol_columns)
        module_indices = {name: idxs for name, idxs in resolved.items() if idxs}
        if not module_indices:
            return ModuleScoreResult(
                summary=pd.DataFrame(columns=[*groupby, "module", "n_cells", "mean_score", "sd_score"]),
                donor_features=pd.DataFrame(columns=["donor_id"]),
                coverage=coverage,
            )

        metadata = build_score_metadata(adata.obs, columns, gateway_config)
        missing_group_cols = [col for col in groupby if col not in metadata.columns]
        if missing_group_cols:
            raise KeyError(f"Unknown module-score grouping columns: {missing_group_cols}")
        qc_mask = build_qc_mask(adata.obs, columns, gateway_config) if apply_qc else np.ones(adata.n_obs, dtype=bool)

        if layer is None and _h5ad_x_is_csr(h5ad_path):
            summary, donor_summary = _score_csr_h5ad(
                h5ad_path=Path(h5ad_path),
                n_obs=adata.n_obs,
                n_vars=adata.n_vars,
                metadata=metadata,
                qc_mask=qc_mask,
                module_indices=module_indices,
                groupby=list(groupby),
                chunk_size=chunk_size,
                log1p=log1p,
            )
        else:
            summary, donor_summary = _score_matrix_source(
                matrix=_matrix_source(adata, layer),
                n_obs=adata.n_obs,
                metadata=metadata,
                qc_mask=qc_mask,
                module_indices=module_indices,
                groupby=list(groupby),
                chunk_size=chunk_size,
                log1p=log1p,
            )
        donor_features = donor_module_features(donor_summary)
        return ModuleScoreResult(summary=summary, donor_features=donor_features, coverage=coverage)
    finally:
        close = getattr(adata, "file", None)
        if close is not None:
            close.close()


def build_score_metadata(obs: pd.DataFrame, columns: ColumnMap, config: dict[str, Any]) -> pd.DataFrame:
    """Build canonical per-cell metadata columns for score grouping."""

    frame = pd.DataFrame(index=obs.index)
    frame["donor_id"] = obs[columns.donor_id].astype(str)
    frame["sample_id"] = obs[columns.sample_id].astype(str)
    frame["age"] = parse_age_series(obs[columns.age])
    frame["disease"] = obs[columns.disease].astype(str)
    frame["disease_group"] = obs[columns.disease].map(lambda value: disease_group(value, config))
    frame["coarse_cell_type"] = obs[columns.coarse_cell_type].astype(str)
    frame["fine_cell_type"] = obs[columns.fine_cell_type].astype(str)
    frame["sex"] = obs[columns.sex].astype(str) if columns.sex else "unknown"
    frame["race_ethnicity"] = obs[columns.race_ethnicity].astype(str) if columns.race_ethnicity else "unknown"
    frame["chemistry"] = obs[columns.chemistry].astype(str) if columns.chemistry else "unknown"
    if columns.collection_method:
        frame["collection_method"] = obs[columns.collection_method].map(lambda value: collection_method_group(value, config))
    else:
        frame["collection_method"] = "unknown"
    frame["site"] = obs[columns.site].astype(str) if columns.site else "unknown"
    return frame


def build_qc_mask(obs: pd.DataFrame, columns: ColumnMap, config: dict[str, Any]) -> np.ndarray:
    """Return cells passing configured QC thresholds."""

    thresholds = config.get("qc_thresholds", {})
    mask = pd.Series(True, index=obs.index)
    if columns.n_counts and "n_counts_min" in thresholds:
        mask &= pd.to_numeric(obs[columns.n_counts], errors="coerce") >= float(thresholds["n_counts_min"])
    if columns.n_genes and "n_genes_min" in thresholds:
        mask &= pd.to_numeric(obs[columns.n_genes], errors="coerce") >= float(thresholds["n_genes_min"])
    if columns.percent_mito and "percent_mito_max" in thresholds:
        mask &= pd.to_numeric(obs[columns.percent_mito], errors="coerce") < float(thresholds["percent_mito_max"])
    if columns.coarse_label_confidence and "coarse_label_confidence_min" in thresholds:
        mask &= pd.to_numeric(obs[columns.coarse_label_confidence], errors="coerce") >= float(
            thresholds["coarse_label_confidence_min"]
        )
    return mask.fillna(False).to_numpy(dtype=bool)


def donor_module_features(donor_summary: pd.DataFrame) -> pd.DataFrame:
    """Create a wide donor-level feature table from donor module summaries."""

    if donor_summary.empty:
        return pd.DataFrame(columns=["donor_id"])
    wide = donor_summary.pivot(index="donor_id", columns="module", values="mean_score").reset_index()
    wide.columns = [
        "donor_id" if col == "donor_id" else f"module_score__{_slugify(col)}"
        for col in wide.columns
    ]
    return wide.sort_values("donor_id").reset_index(drop=True)


def _matrix_source(adata: Any, layer: str | None) -> Any:
    if not layer:
        return adata.X
    if layer not in adata.layers:
        raise KeyError(f"Layer not found in H5AD: {layer}")
    return adata.layers[layer]


def _score_matrix_source(
    *,
    matrix: Any,
    n_obs: int,
    metadata: pd.DataFrame,
    qc_mask: np.ndarray,
    module_indices: dict[str, list[int]],
    groupby: list[str],
    chunk_size: int,
    log1p: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    union_indices = sorted({idx for idxs in module_indices.values() for idx in idxs})
    position = {idx: pos for pos, idx in enumerate(union_indices)}
    module_positions = {
        name: [position[idx] for idx in idxs]
        for name, idxs in module_indices.items()
    }
    group_partials = []
    donor_partials = []
    for start in range(0, n_obs, chunk_size):
        stop = min(start + chunk_size, n_obs)
        keep = qc_mask[start:stop]
        if not keep.any():
            continue
        x = _as_dense(matrix[start:stop, union_indices])
        keep_idx = np.flatnonzero(keep)
        x = x[keep_idx, :]
        if log1p:
            x = np.log1p(np.maximum(x, 0))
        scores = _score_chunk(x, module_positions)
        meta_chunk = metadata.iloc[start + keep_idx].reset_index(drop=True)
        group_partials.append(_partial_summary(meta_chunk, scores, groupby))
        donor_partials.append(_partial_summary(meta_chunk, scores, ["donor_id"]))
    return _combine_partials(group_partials, groupby), _combine_partials(donor_partials, ["donor_id"])


def _score_csr_h5ad(
    *,
    h5ad_path: Path,
    n_obs: int,
    n_vars: int,
    metadata: pd.DataFrame,
    qc_mask: np.ndarray,
    module_indices: dict[str, list[int]],
    groupby: list[str],
    chunk_size: int,
    log1p: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    import h5py

    module_items = list(module_indices.items())
    selected_lookup = np.zeros(n_vars, dtype=bool)
    module_lookups = []
    for module, idxs in module_items:
        selected_lookup[np.asarray(idxs, dtype=np.int64)] = True
        lookup = np.zeros(n_vars, dtype=bool)
        lookup[np.asarray(idxs, dtype=np.int64)] = True
        module_lookups.append((module, lookup, len(idxs)))
    group_partials = []
    donor_partials = []
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
            scores = _score_csr_chunk(
                data_ds=data_ds,
                indices_ds=indices_ds,
                indptr_ds=indptr_ds,
                start=start,
                stop=stop,
                module_lookups=module_lookups,
                selected_lookup=selected_lookup,
                log1p=log1p,
            )
            keep_idx = np.flatnonzero(keep)
            scores = scores.iloc[keep_idx].reset_index(drop=True)
            meta_chunk = metadata.iloc[start + keep_idx].reset_index(drop=True)
            group_partials.append(_partial_summary(meta_chunk, scores, groupby))
            donor_partials.append(_partial_summary(meta_chunk, scores, ["donor_id"]))
    return _combine_partials(group_partials, groupby), _combine_partials(donor_partials, ["donor_id"])


def _score_csr_chunk(
    *,
    data_ds: Any,
    indices_ds: Any,
    indptr_ds: Any,
    start: int,
    stop: int,
    module_lookups: list[tuple[str, np.ndarray, int]],
    selected_lookup: np.ndarray,
    log1p: bool,
) -> pd.DataFrame:
    indptr = np.asarray(indptr_ds[start : stop + 1], dtype=np.int64)
    data_start = int(indptr[0])
    data_stop = int(indptr[-1])
    local_indptr = indptr - data_start
    n_rows = stop - start
    module_names = [name for name, _, _ in module_lookups]
    scores = np.zeros((n_rows, len(module_lookups)), dtype=np.float64)
    if data_stop == data_start:
        return pd.DataFrame(scores, columns=module_names)

    values = np.asarray(data_ds[data_start:data_stop], dtype=np.float32)
    gene_indices = np.asarray(indices_ds[data_start:data_stop], dtype=np.int64)
    selected_mask = selected_lookup[gene_indices]
    if not selected_mask.any():
        return pd.DataFrame(scores, columns=module_names)

    selected_offsets = np.flatnonzero(selected_mask)
    selected_values = values[selected_mask]
    if log1p:
        selected_values = np.log1p(np.maximum(selected_values, 0))
    selected_gene_indices = gene_indices[selected_mask]
    row_ids = np.searchsorted(local_indptr[1:], selected_offsets, side="right")

    for module_pos, (_, module_lookup, n_genes) in enumerate(module_lookups):
        module_mask = module_lookup[selected_gene_indices]
        if module_mask.any():
            np.add.at(scores[:, module_pos], row_ids[module_mask], selected_values[module_mask])
        scores[:, module_pos] /= max(n_genes, 1)
    return pd.DataFrame(scores, columns=module_names)


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


def _score_chunk(x: np.ndarray, module_positions: dict[str, list[int]]) -> pd.DataFrame:
    data = {}
    for module, positions in module_positions.items():
        data[module] = x[:, positions].mean(axis=1)
    return pd.DataFrame(data)


def _partial_summary(meta: pd.DataFrame, scores: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    pieces = []
    groups = meta[group_cols].copy()
    for col in group_cols:
        groups[col] = groups[col].astype("string").fillna("unknown").replace("", "unknown")
    for module in scores.columns:
        values = pd.to_numeric(scores[module], errors="coerce")
        valid = np.isfinite(values.to_numpy(dtype=float))
        if not valid.any():
            continue
        work = groups.loc[valid].copy()
        work["_score"] = values.loc[valid].to_numpy(dtype=float)
        work["_score_sq"] = work["_score"] ** 2
        grouped = (
            work.groupby(group_cols, observed=True, dropna=False)
            .agg(n_cells=("_score", "count"), sum_score=("_score", "sum"), sumsq_score=("_score_sq", "sum"))
            .reset_index()
        )
        grouped["module"] = module
        pieces.append(grouped)
    if not pieces:
        return pd.DataFrame(columns=[*group_cols, "module", "n_cells", "sum_score", "sumsq_score"])
    return pd.concat(pieces, ignore_index=True)


def _combine_partials(partials: list[pd.DataFrame], group_cols: list[str]) -> pd.DataFrame:
    columns = [*group_cols, "module", "n_cells", "mean_score", "sd_score"]
    if not partials:
        return pd.DataFrame(columns=columns)
    frame = pd.concat(partials, ignore_index=True)
    if frame.empty:
        return pd.DataFrame(columns=columns)
    combined = (
        frame.groupby([*group_cols, "module"], observed=True, dropna=False)[["n_cells", "sum_score", "sumsq_score"]]
        .sum()
        .reset_index()
    )
    combined["mean_score"] = combined["sum_score"] / combined["n_cells"].replace(0, np.nan)
    variance = (combined["sumsq_score"] - (combined["sum_score"] ** 2) / combined["n_cells"].replace(0, np.nan)) / (
        combined["n_cells"] - 1
    ).replace(0, np.nan)
    combined["sd_score"] = np.sqrt(np.clip(variance, 0, np.inf)).fillna(0)
    combined = combined.drop(columns=["sum_score", "sumsq_score"])
    return combined.sort_values([*group_cols, "module"]).reset_index(drop=True)


def _slugify(value: object) -> str:
    return normalize_token(value).replace("/", " ").replace("+", " plus ").replace(" ", "_")
