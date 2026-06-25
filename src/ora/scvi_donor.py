"""Donor-level scVI embedding baselines."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Iterable

import numpy as np
import pandas as pd

from ora.model_compare import rank_feature_set_summaries


SCVI_FEATURE_COLUMNS = [
    "model",
    "feature_block",
    "cell_state",
    "statistic",
    "latent_dimension",
    "features",
    "selected_features",
    "mean_abs_importance",
    "total_abs_importance",
    "max_selection_fraction",
    "state_rank_within_model",
]


def build_scvi_donor_embedding_features(
    h5ad_path: str | Path,
    *,
    embedding_key: str = "X_scvi",
    donor_col: str = "donor_id",
    cell_state_col: str = "fine_celltype",
    top_cell_states: int = 12,
    min_cells_per_donor: int = 20,
    min_cells_per_state: int = 20,
    chunk_size: int = 250_000,
    include_state_features: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate cell-level scVI coordinates into donor-level model features.

    The feature matrix intentionally excludes donor and state cell counts, because
    those counts are already covered by ORA composition models. Counts are emitted
    in the QC table so reviewers can audit feature coverage and missingness.
    """

    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover - dependency is declared.
        raise RuntimeError("anndata is required to build scVI donor embedding features.") from exc

    h5ad_path = Path(h5ad_path)
    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        if embedding_key not in adata.obsm:
            available = ", ".join(map(str, adata.obsm.keys()))
            raise KeyError(f"Embedding `{embedding_key}` not found in {h5ad_path}. Available: {available}")
        if donor_col not in adata.obs:
            raise KeyError(f"Donor column `{donor_col}` not found in {h5ad_path}.")
        if include_state_features and cell_state_col not in adata.obs:
            raise KeyError(f"Cell-state column `{cell_state_col}` not found in {h5ad_path}.")

        embedding = adata.obsm[embedding_key]
        if len(embedding.shape) != 2:
            raise ValueError(f"Embedding `{embedding_key}` must be two-dimensional.")
        n_cells, n_dims = int(embedding.shape[0]), int(embedding.shape[1])
        if n_cells != int(adata.n_obs):
            raise ValueError(
                f"Embedding `{embedding_key}` has {n_cells} rows but H5AD has {adata.n_obs} cells."
            )

        obs_cols = [donor_col]
        if include_state_features:
            obs_cols.append(cell_state_col)
        obs = adata.obs[obs_cols].copy()

        donor_values = _clean_labels(obs[donor_col])
        donors = sorted(donor_values.dropna().unique().tolist())
        if not donors:
            raise ValueError(f"No usable donor labels found in column `{donor_col}`.")
        donor_codes = pd.Categorical(donor_values, categories=donors).codes.astype(np.int64)

        state_values = pd.Series(pd.NA, index=obs.index, dtype="string")
        states: list[str] = []
        state_codes = np.full(n_cells, -1, dtype=np.int64)
        state_slug_by_name: dict[str, str] = {}
        if include_state_features and top_cell_states > 0:
            state_values = _clean_labels(obs[cell_state_col])
            state_counts = state_values[donor_codes >= 0].value_counts()
            states = state_counts.head(int(top_cell_states)).index.astype(str).tolist()
            state_codes = pd.Categorical(state_values, categories=states).codes.astype(np.int64)
            state_slug_by_name = _unique_slugs(states)

        donor_sums = np.zeros((len(donors), n_dims), dtype=np.float64)
        donor_sumsq = np.zeros((len(donors), n_dims), dtype=np.float64)
        donor_counts = np.zeros(len(donors), dtype=np.int64)

        flat_state_sums = np.zeros((max(len(states), 1) * len(donors), n_dims), dtype=np.float64)
        flat_state_counts = np.zeros(max(len(states), 1) * len(donors), dtype=np.int64)

        for start in range(0, n_cells, max(1, int(chunk_size))):
            stop = min(start + max(1, int(chunk_size)), n_cells)
            chunk = np.asarray(embedding[start:stop], dtype=np.float64)
            chunk_donor_codes = donor_codes[start:stop]
            valid = (chunk_donor_codes >= 0) & np.isfinite(chunk).all(axis=1)
            if not valid.any():
                continue

            valid_donor_codes = chunk_donor_codes[valid]
            valid_chunk = chunk[valid]
            np.add.at(donor_sums, valid_donor_codes, valid_chunk)
            np.add.at(donor_sumsq, valid_donor_codes, valid_chunk * valid_chunk)
            np.add.at(donor_counts, valid_donor_codes, 1)

            if states:
                chunk_state_codes = state_codes[start:stop][valid]
                valid_state = chunk_state_codes >= 0
                if valid_state.any():
                    flat_codes = valid_donor_codes[valid_state] * len(states) + chunk_state_codes[valid_state]
                    np.add.at(flat_state_sums, flat_codes, valid_chunk[valid_state])
                    np.add.at(flat_state_counts, flat_codes, 1)

        features = _feature_frame(
            donors=donors,
            donor_sums=donor_sums,
            donor_sumsq=donor_sumsq,
            donor_counts=donor_counts,
            states=states,
            state_slug_by_name=state_slug_by_name,
            flat_state_sums=flat_state_sums,
            flat_state_counts=flat_state_counts,
            n_dims=n_dims,
            min_cells_per_state=int(min_cells_per_state),
        )
        qc = _qc_frame(
            h5ad_path=h5ad_path,
            embedding_key=embedding_key,
            donor_col=donor_col,
            cell_state_col=cell_state_col,
            n_dims=n_dims,
            n_cells=n_cells,
            donors=donors,
            donor_counts=donor_counts,
            states=states,
            flat_state_counts=flat_state_counts,
            min_cells_per_donor=int(min_cells_per_donor),
            min_cells_per_state=int(min_cells_per_state),
        )
        return features, qc
    finally:
        close = getattr(getattr(adata, "file", None), "close", None)
        if callable(close):
            close()


def summarize_scvi_state_importance(feature_stability: pd.DataFrame) -> pd.DataFrame:
    """Collapse repeated-CV feature stability into global/state latent importance."""

    columns = list(SCVI_FEATURE_COLUMNS)
    if feature_stability.empty:
        return pd.DataFrame(columns=columns)
    required = {"model", "feature", "abs_mean_importance", "selection_fraction"}
    missing = sorted(required.difference(feature_stability.columns))
    if missing:
        raise ValueError(f"Feature stability table is missing required columns: {', '.join(missing)}")

    frame = feature_stability.copy()
    parsed = frame["feature"].map(_parse_scvi_feature)
    parsed_frame = pd.DataFrame(parsed.tolist(), index=frame.index)
    frame = pd.concat([frame, parsed_frame], axis=1)
    frame = frame[frame["feature_block"].notna()].copy()
    if frame.empty:
        return pd.DataFrame(columns=columns)

    frame["abs_mean_importance"] = pd.to_numeric(frame["abs_mean_importance"], errors="coerce").fillna(0.0)
    frame["selection_fraction"] = pd.to_numeric(frame["selection_fraction"], errors="coerce").fillna(0.0)
    summary = (
        frame.groupby(["model", "feature_block", "cell_state", "statistic"], observed=True)
        .agg(
            features=("feature", "nunique"),
            selected_features=("selection_fraction", lambda s: int(np.sum(pd.to_numeric(s, errors="coerce") > 0))),
            mean_abs_importance=("abs_mean_importance", "mean"),
            total_abs_importance=("abs_mean_importance", "sum"),
            max_selection_fraction=("selection_fraction", "max"),
            latent_dimension=("latent_dimension", _join_unique),
        )
        .reset_index()
    )
    summary["state_rank_within_model"] = (
        summary.groupby("model", observed=True)["total_abs_importance"]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    return summary.sort_values(
        ["model", "state_rank_within_model", "total_abs_importance", "cell_state"],
        ascending=[True, True, False, True],
    )[columns].reset_index(drop=True)


def compare_scvi_donor_baseline(summaries: Iterable[tuple[str, str | Path]]) -> pd.DataFrame:
    """Rank ORA and scVI donor-embedding repeated-CV summaries in one table."""

    loaded: dict[str, pd.DataFrame] = {}
    for label, path in summaries:
        if label in loaded:
            raise ValueError(f"Duplicate feature-set label: {label}")
        loaded[str(label)] = pd.read_csv(path, sep="\t")
    if len(loaded) < 2:
        raise ValueError("At least two summary tables are required for a comparison.")
    return rank_feature_set_summaries(loaded)


def _feature_frame(
    *,
    donors: list[str],
    donor_sums: np.ndarray,
    donor_sumsq: np.ndarray,
    donor_counts: np.ndarray,
    states: list[str],
    state_slug_by_name: dict[str, str],
    flat_state_sums: np.ndarray,
    flat_state_counts: np.ndarray,
    n_dims: int,
    min_cells_per_state: int,
) -> pd.DataFrame:
    feature_data: dict[str, object] = {"donor_id": donors}
    valid_donors = donor_counts > 0
    means = np.full_like(donor_sums, np.nan, dtype=np.float64)
    sds = np.full_like(donor_sums, np.nan, dtype=np.float64)
    means[valid_donors] = donor_sums[valid_donors] / donor_counts[valid_donors, None]
    variances = np.maximum(
        donor_sumsq[valid_donors] / donor_counts[valid_donors, None] - means[valid_donors] ** 2,
        0.0,
    )
    sds[valid_donors] = np.sqrt(variances)
    for dim in range(n_dims):
        dim_label = _dim_label(dim)
        feature_data[f"scvi_global_mean__{dim_label}"] = means[:, dim]
        feature_data[f"scvi_global_sd__{dim_label}"] = sds[:, dim]

    if states:
        state_sums = flat_state_sums.reshape(len(donors), len(states), n_dims)
        state_counts = flat_state_counts.reshape(len(donors), len(states))
        state_means = np.full_like(state_sums, np.nan, dtype=np.float64)
        valid_states = state_counts >= max(1, min_cells_per_state)
        row_idx, state_idx = np.where(valid_states)
        if row_idx.size:
            state_means[row_idx, state_idx] = state_sums[row_idx, state_idx] / state_counts[row_idx, state_idx, None]
        for state_pos, state in enumerate(states):
            slug = state_slug_by_name[state]
            for dim in range(n_dims):
                feature_data[f"scvi_state_mean__{slug}__{_dim_label(dim)}"] = state_means[:, state_pos, dim]

    return pd.DataFrame(feature_data)


def _qc_frame(
    *,
    h5ad_path: Path,
    embedding_key: str,
    donor_col: str,
    cell_state_col: str,
    n_dims: int,
    n_cells: int,
    donors: list[str],
    donor_counts: np.ndarray,
    states: list[str],
    flat_state_counts: np.ndarray,
    min_cells_per_donor: int,
    min_cells_per_state: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = [
        _coverage_row(
            h5ad_path=h5ad_path,
            embedding_key=embedding_key,
            donor_col=donor_col,
            cell_state_col=cell_state_col,
            n_dims=n_dims,
            feature_block="global",
            cell_state="all_cells",
            counts=donor_counts,
            total_cells=int(n_cells),
            min_cells=min_cells_per_donor,
        )
    ]
    if states:
        state_counts = flat_state_counts.reshape(len(donors), len(states))
        for idx, state in enumerate(states):
            counts = state_counts[:, idx]
            rows.append(
                _coverage_row(
                    h5ad_path=h5ad_path,
                    embedding_key=embedding_key,
                    donor_col=donor_col,
                    cell_state_col=cell_state_col,
                    n_dims=n_dims,
                    feature_block="cell_state",
                    cell_state=state,
                    counts=counts,
                    total_cells=int(counts.sum()),
                    min_cells=min_cells_per_state,
                )
            )
    return pd.DataFrame(rows)


def _coverage_row(
    *,
    h5ad_path: Path,
    embedding_key: str,
    donor_col: str,
    cell_state_col: str,
    n_dims: int,
    feature_block: str,
    cell_state: str,
    counts: np.ndarray,
    total_cells: int,
    min_cells: int,
) -> dict[str, object]:
    counts = counts.astype(float)
    donors_with_any = int(np.sum(counts > 0))
    donors_with_min = int(np.sum(counts >= max(1, min_cells)))
    return {
        "h5ad_path": str(h5ad_path),
        "embedding_key": embedding_key,
        "donor_col": donor_col,
        "cell_state_col": cell_state_col,
        "feature_block": feature_block,
        "cell_state": cell_state,
        "latent_dimensions": int(n_dims),
        "total_cells": int(total_cells),
        "donors_with_any_cells": donors_with_any,
        "donors_with_min_cells": donors_with_min,
        "min_cells_required": int(min_cells),
        "median_cells_per_donor": float(np.median(counts)) if counts.size else 0.0,
        "min_cells_per_donor": int(np.min(counts)) if counts.size else 0,
        "max_cells_per_donor": int(np.max(counts)) if counts.size else 0,
        "missing_fraction_at_min": float(1.0 - donors_with_min / counts.size) if counts.size else 1.0,
    }


def _parse_scvi_feature(feature: object) -> dict[str, object]:
    text = str(feature)
    if text.startswith("scvi_global_mean__dim"):
        return {
            "feature_block": "global",
            "cell_state": "all_cells",
            "statistic": "mean",
            "latent_dimension": text.rsplit("__", maxsplit=1)[-1],
        }
    if text.startswith("scvi_global_sd__dim"):
        return {
            "feature_block": "global",
            "cell_state": "all_cells",
            "statistic": "sd",
            "latent_dimension": text.rsplit("__", maxsplit=1)[-1],
        }
    if text.startswith("scvi_state_mean__") and "__dim" in text:
        body = text.removeprefix("scvi_state_mean__")
        state, dim = body.rsplit("__", maxsplit=1)
        return {
            "feature_block": "cell_state",
            "cell_state": state,
            "statistic": "mean",
            "latent_dimension": dim,
        }
    return {"feature_block": pd.NA, "cell_state": pd.NA, "statistic": pd.NA, "latent_dimension": pd.NA}


def _clean_labels(series: pd.Series) -> pd.Series:
    cleaned = series.astype("string").str.strip()
    bad = cleaned.str.lower().isin(["", "nan", "none", "na", "<na>"])
    return cleaned.mask(bad)


def _dim_label(dim: int) -> str:
    return f"dim{dim + 1:02d}"


def _unique_slugs(values: list[str]) -> dict[str, str]:
    used: dict[str, int] = {}
    output: dict[str, str] = {}
    for value in values:
        base = _slugify(value)
        count = used.get(base, 0) + 1
        used[base] = count
        output[value] = base if count == 1 else f"{base}_{count}"
    return output


def _slugify(value: object) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "unknown"


def _join_unique(values: pd.Series) -> str:
    unique = []
    for value in values:
        text = str(value)
        if text and text not in unique:
            unique.append(text)
    return ";".join(unique)
