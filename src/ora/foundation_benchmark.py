"""Foundation-model benchmark planning and input subset helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .age_model import donor_cv_folds
from .config import load_config, project_path
from .metadata import resolve_columns
from .utils import ensure_parent, normalize_token


@dataclass(frozen=True)
class BenchmarkSubsetSpec:
    """Cell-subset configuration for single-cell foundation-model benchmarks."""

    name: str
    max_cells: int
    description: str


DEFAULT_SUBSET_SPECS = (
    BenchmarkSubsetSpec(
        name="lineage",
        max_cells=120_000,
        description="Olfactory basal-to-neuronal lineage cells enriched for HBC, INP, iOSN, and mOSN states.",
    ),
    BenchmarkSubsetSpec(
        name="epithelial",
        max_cells=180_000,
        description="Broad olfactory and respiratory epithelial subset, including sustentacular, secretory, glandular, and neuronal states.",
    ),
    BenchmarkSubsetSpec(
        name="allcell",
        max_cells=250_000,
        description="All-cell donor/fine-cell-type stratified subset for general-purpose foundation embeddings.",
    ),
)


DEFAULT_OUTPUTS = {
    "lineage": "data/processed/foundation_benchmark_lineage_subset.h5ad",
    "epithelial": "data/processed/foundation_benchmark_epithelial_subset.h5ad",
    "allcell": "data/processed/foundation_benchmark_allcell_subset.h5ad",
}


def build_foundation_benchmark_subsets(
    *,
    h5ad_path: str | Path,
    config: dict[str, Any],
    model_config: dict[str, Any],
    output_paths: dict[str, str | Path] | None = None,
    subset_specs: tuple[BenchmarkSubsetSpec, ...] = DEFAULT_SUBSET_SPECS,
    manifest_out: str | Path = "results/tables/foundation_benchmark_subset_manifest.tsv",
    donor_splits_out: str | Path = "results/tables/foundation_benchmark_donor_splits.tsv",
    gene_manifest_out: str | Path = "results/tables/foundation_benchmark_gene_manifest.tsv",
    gene_symbols_out: str | Path = "resources/foundation_benchmark/gateway_gene_symbols.txt",
    gene_ids_out: str | Path = "resources/foundation_benchmark/gateway_gene_ids.txt",
    seed: int = 20260625,
    overwrite: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Write benchmark H5AD subsets plus donor split and gene manifests."""

    try:
        import anndata as ad  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError("anndata is required to build benchmark subsets.") from exc

    output_paths = output_paths or DEFAULT_OUTPUTS
    adata = ad.read_h5ad(project_path(h5ad_path), backed="r")
    try:
        obs = adata.obs.copy()
        obs["_ora_row_position"] = np.arange(adata.n_obs, dtype=np.int64)
        gene_manifest = build_gene_manifest(adata.var.copy(), pd.Index(adata.var_names))
        _write_gene_exports(gene_manifest, gene_manifest_out, gene_symbols_out, gene_ids_out)

        split_table = build_donor_split_table(obs, config, model_config)
        split_table.to_csv(ensure_parent(donor_splits_out), sep="\t", index=False)

        rows: list[dict[str, Any]] = []
        for spec in subset_specs:
            if spec.name not in output_paths:
                raise KeyError(f"Missing output path for benchmark subset `{spec.name}`.")
            out_path = ensure_parent(output_paths[spec.name])
            if out_path.exists() and not overwrite:
                raise FileExistsError(f"Output already exists: {out_path}. Pass overwrite=True to replace it.")

            indices = select_foundation_subset_indices(obs, config, spec.name, spec.max_cells, seed=seed)
            if indices.size == 0:
                raise ValueError(f"No cells selected for benchmark subset `{spec.name}`.")
            subset = adata[indices, :].to_memory()
            subset.uns["ora_foundation_benchmark"] = {
                "subset": spec.name,
                "description": spec.description,
                "max_cells": int(spec.max_cells),
                "sampling_seed": int(seed),
                "source_h5ad": str(h5ad_path),
            }
            subset.write_h5ad(out_path, compression="gzip")
            rows.append(
                summarize_subset(
                    subset.obs,
                    subset_name=spec.name,
                    output_path=out_path,
                    description=spec.description,
                    max_cells=spec.max_cells,
                    seed=seed,
                    n_genes=subset.n_vars,
                )
            )
            del subset

        manifest = pd.DataFrame(rows).sort_values("subset").reset_index(drop=True)
        manifest.to_csv(ensure_parent(manifest_out), sep="\t", index=False)
        return manifest, split_table, gene_manifest
    finally:
        close = getattr(adata, "file", None)
        if close is not None:
            close.close()


def select_foundation_subset_indices(
    obs: pd.DataFrame,
    config: dict[str, Any],
    subset_name: str,
    max_cells: int,
    *,
    seed: int = 20260625,
) -> np.ndarray:
    """Return deterministic row-position indices for a benchmark subset."""

    if max_cells <= 0:
        raise ValueError("max_cells must be positive.")
    work = _subset_columns(obs, config)
    if "_ora_row_position" not in work:
        work["_ora_row_position"] = np.arange(work.shape[0], dtype=np.int64)

    mask = _subset_mask(work, subset_name)
    candidates = work.loc[mask].copy()
    if candidates.empty:
        return np.array([], dtype=np.int64)
    if candidates.shape[0] <= max_cells:
        return np.sort(candidates["_ora_row_position"].to_numpy(dtype=np.int64))

    sampled_positions = _stratified_positions(
        candidates,
        max_cells=max_cells,
        seed=seed + _subset_seed_offset(subset_name),
    )
    return np.sort(sampled_positions.astype(np.int64))


def build_donor_split_table(
    obs_or_manifest: pd.DataFrame,
    config: dict[str, Any],
    model_config: dict[str, Any],
) -> pd.DataFrame:
    """Export the deterministic donor folds used by ORA age models."""

    donor_meta = _donor_manifest_from_obs(obs_or_manifest, config)
    if "usable_for_ora_training" in donor_meta:
        eligible = _boolean_series(donor_meta["usable_for_ora_training"]) & donor_meta["age"].notna()
    else:
        disease = donor_meta.get("disease_group", donor_meta.get("condition", ""))
        eligible = pd.Series([_is_healthy(value, config) for value in disease], index=donor_meta.index) & donor_meta[
            "age"
        ].notna()
    train = donor_meta.loc[eligible].sort_values("donor_id").reset_index(drop=True)
    folds = donor_cv_folds(train, model_config)
    rows: list[dict[str, Any]] = []
    for fold_idx, (train_idx, test_idx) in enumerate(folds, start=1):
        for split, indices in (("train", train_idx), ("test", test_idx)):
            for idx in indices:
                donor = train.iloc[int(idx)]
                rows.append(
                    {
                        "fold": fold_idx,
                        "split": split,
                        "donor_id": donor["donor_id"],
                        "age": donor.get("age", np.nan),
                        "sex": donor.get("sex", ""),
                        "chemistry": donor.get("chemistry", donor.get("flex_version", "")),
                        "collection_method": donor.get("collection_method", donor.get("device_guided", "")),
                        "disease_group": donor.get("disease_group", donor.get("condition", "")),
                    }
                )
    return pd.DataFrame(rows)


def build_gene_manifest(var: pd.DataFrame, var_names: pd.Index) -> pd.DataFrame:
    """Return gene metadata needed by foundation-model tokenizers."""

    output = pd.DataFrame(
        {
            "gene_id": var_names.astype(str),
            "gene_symbol": var["feature_name"].astype(str).to_numpy()
            if "feature_name" in var
            else var_names.astype(str),
        }
    )
    for col in ["feature_biotype", "feature_type", "feature_reference"]:
        output[col] = var[col].astype(str).to_numpy() if col in var else ""
    output["is_protein_coding"] = output["feature_type"].str.lower().eq("protein_coding") | output[
        "feature_biotype"
    ].str.lower().eq("protein_coding")
    return output


def summarize_subset(
    obs: pd.DataFrame,
    *,
    subset_name: str,
    output_path: str | Path,
    description: str,
    max_cells: int,
    seed: int,
    n_genes: int,
) -> dict[str, Any]:
    """Summarize a written benchmark subset."""

    donor_col = _first_existing(obs, ("donor_id",))
    sample_col = _first_existing(obs, ("sample_id",))
    fine_col = _first_existing(obs, ("fine_celltype", "fine_cell_type", "cell_type"))
    coarse_col = _first_existing(obs, ("coarse_celltype", "coarse_cell_type", "cell_type"))
    donor_counts = obs[donor_col].astype(str).value_counts() if donor_col else pd.Series(dtype=int)
    return {
        "subset": subset_name,
        "path": str(output_path),
        "description": description,
        "sampling_seed": int(seed),
        "max_cells": int(max_cells),
        "n_cells": int(obs.shape[0]),
        "n_genes": int(n_genes),
        "n_donors": int(obs[donor_col].nunique()) if donor_col else 0,
        "n_samples": int(obs[sample_col].nunique()) if sample_col else 0,
        "n_fine_cell_types": int(obs[fine_col].nunique()) if fine_col else 0,
        "n_coarse_cell_types": int(obs[coarse_col].nunique()) if coarse_col else 0,
        "min_cells_per_donor": int(donor_counts.min()) if not donor_counts.empty else 0,
        "median_cells_per_donor": float(donor_counts.median()) if not donor_counts.empty else 0.0,
        "max_cells_per_donor": int(donor_counts.max()) if not donor_counts.empty else 0,
    }


def load_configs(config_path: str | Path, model_config_path: str | Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """Load Gateway and model config files for the benchmark CLI."""

    return load_config(config_path), load_config(model_config_path)


def _write_gene_exports(
    gene_manifest: pd.DataFrame,
    gene_manifest_out: str | Path,
    gene_symbols_out: str | Path,
    gene_ids_out: str | Path,
) -> None:
    gene_manifest.to_csv(ensure_parent(gene_manifest_out), sep="\t", index=False)
    ensure_parent(gene_symbols_out).write_text(
        "\n".join(gene_manifest["gene_symbol"].dropna().astype(str).tolist()) + "\n",
        encoding="utf-8",
    )
    ensure_parent(gene_ids_out).write_text(
        "\n".join(gene_manifest["gene_id"].dropna().astype(str).tolist()) + "\n",
        encoding="utf-8",
    )


def _subset_columns(obs: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    columns = {
        "donor_id": _resolve_alias(obs, config, "donor_id"),
        "sample_id": _resolve_alias(obs, config, "sample_id"),
        "coarse_cell_type": _resolve_alias(obs, config, "coarse_cell_type"),
        "fine_cell_type": _resolve_alias(obs, config, "fine_cell_type"),
    }
    missing = [key for key, value in columns.items() if value is None]
    if missing:
        raise KeyError(f"Missing required columns for foundation subset selection: {', '.join(missing)}")
    required = {
        "donor_id": columns["donor_id"],
        "sample_id": columns["sample_id"],
        "coarse_cell_type": columns["coarse_cell_type"],
        "fine_cell_type": columns["fine_cell_type"],
    }
    output = pd.DataFrame({key: obs[col].astype(str).to_numpy() for key, col in required.items()})
    if "_ora_row_position" in obs:
        output["_ora_row_position"] = obs["_ora_row_position"].to_numpy(dtype=np.int64)
    return output


def _subset_mask(work: pd.DataFrame, subset_name: str) -> pd.Series:
    subset = normalize_token(subset_name)
    fine = work["fine_cell_type"].map(normalize_token)
    coarse = work["coarse_cell_type"].map(normalize_token)
    if subset in {"allcell", "all cell", "all cells"}:
        return pd.Series(True, index=work.index)
    if subset == "lineage":
        lineage_terms = (
            "hbc",
            "inp",
            "iosn",
            "mosn",
            "olfactory sensory neuron",
        )
        return _contains_any(fine, lineage_terms) | _contains_any(coarse, ("olf mosn", "olf hbc"))
    if subset == "epithelial":
        epithelial_terms = (
            "olf",
            "resp",
            "hbc",
            "inp",
            "iosn",
            "mosn",
            "sustentacular",
            "sus",
            "microvillar",
            "gland",
            "club",
            "goblet",
            "ciliated",
            "deuterosomal",
            "ionocyte",
        )
        return _contains_any(fine, epithelial_terms) | _contains_any(coarse, epithelial_terms)
    raise ValueError(f"Unknown foundation benchmark subset `{subset_name}`.")


def _stratified_positions(candidates: pd.DataFrame, *, max_cells: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    work = candidates[["donor_id", "fine_cell_type", "_ora_row_position"]].copy()
    work["_random"] = rng.random(work.shape[0])
    group_cols = ["donor_id", "fine_cell_type"]
    n_groups = max(1, work.groupby(group_cols, observed=True).ngroups)
    base_quota = max(1, max_cells // n_groups)
    sampled = (
        work.sort_values("_random")
        .groupby(group_cols, observed=True, group_keys=False)
        .head(base_quota)
        .copy()
    )
    if sampled.shape[0] < max_cells:
        remaining = work.loc[~work["_ora_row_position"].isin(sampled["_ora_row_position"])]
        top_up = remaining.sort_values("_random").head(max_cells - sampled.shape[0])
        sampled = pd.concat([sampled, top_up], ignore_index=True)
    elif sampled.shape[0] > max_cells:
        sampled = sampled.sort_values("_random").head(max_cells)
    return sampled["_ora_row_position"].to_numpy(dtype=np.int64)


def _donor_manifest_from_obs(obs_or_manifest: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    frame = obs_or_manifest.copy()
    if {"donor_id", "age"}.issubset(frame.columns):
        donor = frame.sort_values(["donor_id"]).drop_duplicates("donor_id").copy()
        donor["age"] = pd.to_numeric(donor["age"], errors="coerce")
        return donor

    columns = resolve_columns(list(frame.columns), config)
    donor = pd.DataFrame(
        {
            "donor_id": frame[columns.donor_id].astype(str).to_numpy(),
            "age": _parse_age(frame[columns.age]) if columns.age else np.nan,
            "sex": frame[columns.sex].astype(str).to_numpy() if columns.sex else "",
            "chemistry": frame[columns.chemistry].astype(str).to_numpy() if columns.chemistry else "",
            "collection_method": frame[columns.collection_method].astype(str).to_numpy()
            if columns.collection_method
            else "",
            "condition": frame[columns.disease].astype(str).to_numpy() if columns.disease else "",
        }
    )
    return donor.sort_values("donor_id").drop_duplicates("donor_id").reset_index(drop=True)


def _parse_age(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    missing = numeric.isna()
    if missing.any():
        extracted = values.astype(str).str.extract(r"([0-9]+(?:\.[0-9]+)?)", expand=False)
        numeric = numeric.where(~missing, pd.to_numeric(extracted, errors="coerce"))
    return numeric


def _boolean_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.lower().isin({"true", "1", "yes", "y"})


def _is_healthy(value: object, config: dict[str, Any]) -> bool:
    healthy = {normalize_token(item) for item in config.get("healthy_values", [])}
    return normalize_token(value) in healthy


def _contains_any(series: pd.Series, needles: tuple[str, ...]) -> pd.Series:
    output = pd.Series(False, index=series.index)
    for needle in needles:
        output = output | series.str.contains(needle, regex=False, na=False)
    return output


def _subset_seed_offset(subset_name: str) -> int:
    return sum(ord(ch) for ch in subset_name)


def _first_existing(frame: pd.DataFrame, columns: tuple[str, ...]) -> str | None:
    for col in columns:
        if col in frame:
            return col
    return None


def _resolve_alias(obs: pd.DataFrame, config: dict[str, Any], key: str) -> str | None:
    aliases = config.get("columns", {}).get(key, [])
    if isinstance(aliases, str):
        aliases = [aliases]
    normalized = {normalize_token(col): col for col in obs.columns}
    lowered = {str(col).lower(): col for col in obs.columns}
    for alias in aliases:
        if alias in obs:
            return str(alias)
        if str(alias).lower() in lowered:
            return str(lowered[str(alias).lower()])
        token = normalize_token(alias)
        if token in normalized:
            return str(normalized[token])
    return None
