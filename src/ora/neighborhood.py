"""Milo-style latent-neighborhood differential abundance pilots."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests


@dataclass(frozen=True)
class NeighborhoodConfig:
    n_neighborhoods: int = 1000
    n_neighbors: int = 50
    min_donors: int = 20
    seed: int = 13
    age_column: str = "age"
    donor_column: str = "donor_id"
    fine_column: str = "fine_celltype"
    coarse_column: str = "coarse_celltype"
    covariates: tuple[str, ...] = ("sex", "chemistry", "collection_method")
    numeric_covariates: tuple[str, ...] = ("total_cells",)
    seed_stratify_columns: tuple[str, ...] = ()


def run_neighborhood_da(
    embedding: np.ndarray,
    cell_metadata: pd.DataFrame,
    donor_metadata: pd.DataFrame,
    *,
    config: NeighborhoodConfig | None = None,
    return_memberships: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame] | tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run a lightweight Milo-style donor-level DA screen over latent neighborhoods."""

    cfg = config or NeighborhoodConfig()
    if embedding.ndim != 2 or embedding.shape[0] != cell_metadata.shape[0]:
        raise ValueError("embedding must be cells x dimensions and align to cell_metadata")
    cells = cell_metadata.reset_index(drop=True).copy()
    donors = donor_metadata.copy()
    if cfg.donor_column not in cells or cfg.donor_column not in donors:
        raise ValueError(f"`{cfg.donor_column}` must exist in cell and donor metadata")
    if cfg.age_column not in donors:
        raise ValueError(f"`{cfg.age_column}` must exist in donor metadata")

    donor_index = donors.set_index(cfg.donor_column, drop=False)
    keep_cells = cells[cfg.donor_column].isin(donor_index.index).to_numpy()
    if keep_cells.sum() < cfg.n_neighbors:
        raise ValueError("not enough cells with matched donor metadata")
    embedding = np.asarray(embedding[keep_cells, :], dtype=float)
    cells = cells.loc[keep_cells].reset_index(drop=True)

    rng = np.random.default_rng(cfg.seed)
    seed_indices = _seed_indices(cells, n_cells=embedding.shape[0], cfg=cfg, rng=rng)
    n_neighbors = min(cfg.n_neighbors, embedding.shape[0])
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
    nn.fit(embedding)
    neighbor_indices = nn.kneighbors(embedding[seed_indices], return_distance=False)
    total_cells_by_donor = cells[cfg.donor_column].value_counts().rename("donor_cells_in_latent")

    rows = []
    membership_rows = []
    for neighborhood_id, indices in enumerate(neighbor_indices):
        neighborhood_cells = cells.iloc[indices]
        if return_memberships:
            membership_rows.extend(_membership_rows(neighborhood_id, neighborhood_cells, cfg))
        counts = neighborhood_cells[cfg.donor_column].value_counts()
        observed_donors = counts.index.intersection(donor_index.index)
        if len(observed_donors) < cfg.min_donors:
            rows.append(_skipped_row(neighborhood_id, neighborhood_cells, len(observed_donors), "too_few_donors", cfg))
            continue
        model_df = donor_index.loc[observed_donors].copy()
        model_df["neighborhood_cells"] = counts.reindex(observed_donors).fillna(0).astype(float)
        model_df["donor_cells_in_latent"] = total_cells_by_donor.reindex(observed_donors).fillna(0).astype(float)
        result = _fit_neighborhood_model(model_df, cfg)
        rows.append(_neighborhood_row(neighborhood_id, neighborhood_cells, model_df, result, cfg))

    neighborhood_table = pd.DataFrame(rows)
    if "age_pvalue" in neighborhood_table and neighborhood_table["age_pvalue"].notna().any():
        tested = neighborhood_table["age_pvalue"].notna()
        neighborhood_table.loc[tested, "age_fdr"] = multipletests(
            neighborhood_table.loc[tested, "age_pvalue"].to_numpy(),
            method="fdr_bh",
        )[1]
    else:
        neighborhood_table["age_fdr"] = np.nan
    neighborhood_table["age_fdr"] = neighborhood_table["age_fdr"].astype(float)
    summary = summarize_neighborhood_da(neighborhood_table)
    if return_memberships:
        memberships = pd.DataFrame(
            membership_rows,
            columns=["neighborhood_id", "cell_index", "obs_name", cfg.donor_column, cfg.fine_column, cfg.coarse_column],
        )
        return neighborhood_table, summary, memberships
    return neighborhood_table, summary


def summarize_neighborhood_da(neighborhoods: pd.DataFrame) -> pd.DataFrame:
    """Summarize a neighborhood DA result table for reporting."""

    if neighborhoods.empty:
        return pd.DataFrame([{"metric": "neighborhoods", "value": 0, "detail": "no neighborhoods tested"}])
    tested = neighborhoods["status"].eq("tested") if "status" in neighborhoods else pd.Series(False, index=neighborhoods.index)
    sig = neighborhoods["age_fdr"].lt(0.10) if "age_fdr" in neighborhoods else pd.Series(False, index=neighborhoods.index)
    rows = [
        {"metric": "neighborhoods_total", "value": int(neighborhoods.shape[0]), "detail": "seed neighborhoods generated"},
        {"metric": "neighborhoods_tested", "value": int(tested.sum()), "detail": "neighborhoods with enough donors for regression"},
        {"metric": "age_fdr_lt_0_10", "value": int((tested & sig).sum()), "detail": "age-associated neighborhoods at BH FDR < 0.10"},
    ]
    if tested.any():
        top = neighborhoods.loc[tested].sort_values("age_pvalue", na_position="last").head(1).iloc[0]
        rows.append(
            {
                "metric": "top_age_neighborhood",
                "value": int(top["neighborhood_id"]),
                "detail": (
                    f"coef={top['age_coef']:.4g};p={top['age_pvalue']:.3g};"
                    f"fdr={top['age_fdr']:.3g};fine={top['top_fine_celltype']};coarse={top['top_coarse_celltype']}"
                ),
            }
        )
    return pd.DataFrame(rows, columns=["metric", "value", "detail"])


def _fit_neighborhood_model(model_df: pd.DataFrame, cfg: NeighborhoodConfig) -> dict[str, float | str]:
    usable = model_df.copy()
    usable = usable[np.isfinite(pd.to_numeric(usable[cfg.age_column], errors="coerce"))].copy()
    usable[cfg.age_column] = pd.to_numeric(usable[cfg.age_column], errors="coerce")
    usable = usable[usable["donor_cells_in_latent"].gt(0)]
    if usable.shape[0] < cfg.min_donors:
        return {"status": "too_few_donors_after_qc", "age_coef": np.nan, "age_pvalue": np.nan, "n_model_donors": usable.shape[0]}
    proportion = (usable["neighborhood_cells"].to_numpy(dtype=float) + 0.5) / (
        usable["donor_cells_in_latent"].to_numpy(dtype=float) + 1.0
    )
    usable["logit_fraction"] = np.log(proportion / (1.0 - proportion))
    design = pd.DataFrame({"age_scaled": _zscore(usable[cfg.age_column].to_numpy(dtype=float))}, index=usable.index)
    for covariate in cfg.covariates:
        if covariate in usable:
            values = usable[covariate].astype(str).replace({"nan": "unknown", "None": "unknown"})
            if values.nunique() > 1:
                dummies = pd.get_dummies(values, prefix=covariate, drop_first=True, dtype=float)
                design = pd.concat([design, dummies], axis=1)
    for covariate in cfg.numeric_covariates:
        if covariate in usable:
            values = pd.to_numeric(usable[covariate], errors="coerce").to_numpy(dtype=float)
            if np.isfinite(values).sum() >= cfg.min_donors and float(np.nanstd(values)) > 0:
                design[covariate] = _zscore(np.log1p(np.nan_to_num(values, nan=np.nanmedian(values))))
    design = sm.add_constant(design, has_constant="add")
    try:
        fit = sm.OLS(usable["logit_fraction"].astype(float), design.astype(float)).fit()
    except (np.linalg.LinAlgError, ValueError) as exc:
        return {"status": "model_failed", "age_coef": np.nan, "age_pvalue": np.nan, "n_model_donors": usable.shape[0], "error": str(exc)}
    return {
        "status": "tested",
        "age_coef": float(fit.params.get("age_scaled", np.nan)),
        "age_pvalue": float(fit.pvalues.get("age_scaled", np.nan)),
        "n_model_donors": int(usable.shape[0]),
        "model_r2": float(fit.rsquared),
    }


def _neighborhood_row(
    neighborhood_id: int,
    neighborhood_cells: pd.DataFrame,
    model_df: pd.DataFrame,
    result: dict[str, float | str],
    cfg: NeighborhoodConfig,
) -> dict[str, float | int | str]:
    row = _base_neighborhood_annotation(neighborhood_id, neighborhood_cells, cfg)
    row.update(
        {
            "n_donors": int(model_df.shape[0]),
            "mean_age": float(pd.to_numeric(model_df[cfg.age_column], errors="coerce").mean()),
            "age_coef": result.get("age_coef", np.nan),
            "age_pvalue": result.get("age_pvalue", np.nan),
            "model_r2": result.get("model_r2", np.nan),
            "status": str(result.get("status", "tested")),
        }
    )
    return row


def _skipped_row(
    neighborhood_id: int,
    neighborhood_cells: pd.DataFrame,
    n_donors: int,
    status: str,
    cfg: NeighborhoodConfig,
) -> dict[str, float | int | str]:
    row = _base_neighborhood_annotation(neighborhood_id, neighborhood_cells, cfg)
    row.update({"n_donors": int(n_donors), "mean_age": np.nan, "age_coef": np.nan, "age_pvalue": np.nan, "model_r2": np.nan, "status": status})
    return row


def _base_neighborhood_annotation(
    neighborhood_id: int,
    neighborhood_cells: pd.DataFrame,
    cfg: NeighborhoodConfig,
) -> dict[str, float | int | str]:
    fine = _top_label(neighborhood_cells, cfg.fine_column)
    coarse = _top_label(neighborhood_cells, cfg.coarse_column)
    return {
        "neighborhood_id": int(neighborhood_id),
        "n_cells": int(neighborhood_cells.shape[0]),
        "top_fine_celltype": fine[0],
        "top_fine_fraction": fine[1],
        "top_coarse_celltype": coarse[0],
        "top_coarse_fraction": coarse[1],
    }


def _top_label(frame: pd.DataFrame, column: str) -> tuple[str, float]:
    if column not in frame or frame.empty:
        return "unknown", np.nan
    counts = frame[column].astype(str).replace({"nan": "unknown", "None": "unknown"}).value_counts()
    if counts.empty:
        return "unknown", np.nan
    return str(counts.index[0]), float(counts.iloc[0] / counts.sum())


def _membership_rows(
    neighborhood_id: int,
    neighborhood_cells: pd.DataFrame,
    cfg: NeighborhoodConfig,
) -> list[dict[str, int | str]]:
    rows = []
    for local_position, row in neighborhood_cells.iterrows():
        rows.append(
            {
                "neighborhood_id": int(neighborhood_id),
                "cell_index": int(row.get("_cell_index", local_position)),
                "obs_name": str(row.get("_obs_name", local_position)),
                cfg.donor_column: str(row.get(cfg.donor_column, "unknown")),
                cfg.fine_column: str(row.get(cfg.fine_column, "unknown")),
                cfg.coarse_column: str(row.get(cfg.coarse_column, "unknown")),
            }
        )
    return rows


def _seed_indices(
    cells: pd.DataFrame,
    *,
    n_cells: int,
    cfg: NeighborhoodConfig,
    rng: np.random.Generator,
) -> np.ndarray:
    n_seeds = min(cfg.n_neighborhoods, n_cells)
    columns = [column for column in cfg.seed_stratify_columns if column in cells]
    if not columns:
        return np.sort(rng.choice(n_cells, size=n_seeds, replace=False))

    strata = cells[columns].astype(str).agg("|".join, axis=1)
    groups = [index.to_numpy(dtype=int) for _, index in pd.Series(np.arange(n_cells)).groupby(strata, observed=True)]
    if not groups:
        return np.sort(rng.choice(n_cells, size=n_seeds, replace=False))

    per_group = max(1, int(np.ceil(n_seeds / len(groups))))
    selected: list[np.ndarray] = []
    for group in groups:
        take = min(per_group, group.shape[0])
        selected.append(rng.choice(group, size=take, replace=False))
    seeds = np.unique(np.concatenate(selected)) if selected else np.array([], dtype=int)
    if seeds.shape[0] < n_seeds:
        remaining = np.setdiff1d(np.arange(n_cells), seeds, assume_unique=False)
        take = min(n_seeds - seeds.shape[0], remaining.shape[0])
        if take:
            seeds = np.concatenate([seeds, rng.choice(remaining, size=take, replace=False)])
    elif seeds.shape[0] > n_seeds:
        seeds = rng.choice(seeds, size=n_seeds, replace=False)
    return np.sort(seeds.astype(int))


def _zscore(values: np.ndarray) -> np.ndarray:
    sd = float(np.nanstd(values))
    if sd == 0 or not np.isfinite(sd):
        return np.zeros_like(values, dtype=float)
    return (values - float(np.nanmean(values))) / sd
