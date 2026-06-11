"""Metadata harmonization and cohort construction."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any

import numpy as np
import pandas as pd

from .utils import normalize_token


@dataclass(frozen=True)
class ColumnMap:
    donor_id: str
    sample_id: str
    age: str
    sex: str | None
    race_ethnicity: str | None
    disease: str
    chemistry: str | None
    collection_method: str | None
    site: str | None
    coarse_cell_type: str
    fine_cell_type: str
    n_counts: str | None
    n_genes: str | None
    percent_mito: str | None
    coarse_label_confidence: str | None


REQUIRED_CONFIG_COLUMNS = [
    "donor_id",
    "sample_id",
    "age",
    "disease",
    "coarse_cell_type",
    "fine_cell_type",
]


def resolve_columns(obs_columns: list[str], config: dict[str, Any]) -> ColumnMap:
    """Resolve configured metadata aliases to actual AnnData obs column names."""

    column_config = config.get("columns", {})
    resolved: dict[str, str | None] = {}
    lowered = {col.lower(): col for col in obs_columns}
    normalized = {normalize_token(col): col for col in obs_columns}

    for key, aliases in column_config.items():
        if isinstance(aliases, str):
            alias_list = [aliases]
        else:
            alias_list = list(aliases or [])
        resolved[key] = None
        for alias in alias_list:
            if alias in obs_columns:
                resolved[key] = alias
                break
            if alias.lower() in lowered:
                resolved[key] = lowered[alias.lower()]
                break
            token = normalize_token(alias)
            if token in normalized:
                resolved[key] = normalized[token]
                break

    missing = [key for key in REQUIRED_CONFIG_COLUMNS if not resolved.get(key)]
    if missing:
        available = ", ".join(obs_columns)
        raise KeyError(f"Could not resolve required metadata columns {missing}. Available: {available}")

    return ColumnMap(
        donor_id=str(resolved["donor_id"]),
        sample_id=str(resolved["sample_id"]),
        age=str(resolved["age"]),
        sex=resolved.get("sex"),
        race_ethnicity=resolved.get("race_ethnicity"),
        disease=str(resolved["disease"]),
        chemistry=resolved.get("chemistry"),
        collection_method=resolved.get("collection_method"),
        site=resolved.get("site"),
        coarse_cell_type=str(resolved["coarse_cell_type"]),
        fine_cell_type=str(resolved["fine_cell_type"]),
        n_counts=resolved.get("n_counts"),
        n_genes=resolved.get("n_genes"),
        percent_mito=resolved.get("percent_mito"),
        coarse_label_confidence=resolved.get("coarse_label_confidence"),
    )


def build_manifest(obs: pd.DataFrame, config: dict[str, Any], columns: ColumnMap | None = None) -> pd.DataFrame:
    """Build one donor/sample-level manifest from cell metadata."""

    columns = columns or resolve_columns(list(obs.columns), config)
    work = pd.DataFrame(index=obs.index)
    work["_donor_id"] = obs[columns.donor_id].astype(str)
    work["_sample_id"] = obs[columns.sample_id].astype(str)
    work["_age"] = parse_age_series(obs[columns.age])
    work["_disease"] = obs[columns.disease]
    work["_disease_group"] = work["_disease"].map(lambda x: disease_group(x, config))
    work["_fine_cell_type"] = obs[columns.fine_cell_type].astype(str)
    work["_coarse_cell_type"] = obs[columns.coarse_cell_type].astype(str)

    optional_pairs = {
        "sex": columns.sex,
        "race_ethnicity": columns.race_ethnicity,
        "chemistry": columns.chemistry,
        "collection_method": columns.collection_method,
        "site": columns.site,
    }
    for out_col, source_col in optional_pairs.items():
        work[f"_{out_col}"] = obs[source_col] if source_col else np.nan
    if columns.collection_method:
        work["_collection_method"] = work["_collection_method"].map(
            lambda value: collection_method_group(value, config)
        )

    grouped = work.groupby(["_donor_id", "_sample_id"], observed=True, dropna=False)
    manifest = grouped.agg(
        age=("_age", _first_non_null),
        disease=("_disease", _first_non_null),
        disease_group=("_disease_group", _first_non_null),
        sex=("_sex", _first_non_null),
        race_ethnicity=("_race_ethnicity", _first_non_null),
        chemistry=("_chemistry", _first_non_null),
        collection_method=("_collection_method", _first_non_null),
        site=("_site", _first_non_null),
        total_cells=("_fine_cell_type", "size"),
    ).reset_index()
    manifest = manifest.rename(columns={"_donor_id": "donor_id", "_sample_id": "sample_id"})

    lineage_aliases = flatten_lineage_aliases(config)
    mature_aliases = lineage_aliases.get("fully_mature_mosn", set()) | lineage_aliases.get("stressed_mosn", set())
    lineage_tokens = set().union(*lineage_aliases.values()) if lineage_aliases else set()

    per_sample = work.assign(
        _is_lineage=work["_fine_cell_type"].map(lambda x: normalize_token(x) in lineage_tokens),
        _is_mature=work["_fine_cell_type"].map(lambda x: normalize_token(x) in mature_aliases),
    )
    lineage_counts = per_sample.groupby(["_donor_id", "_sample_id"], observed=True).agg(
        lineage_cells=("_is_lineage", "sum"),
        mature_neurons=("_is_mature", "sum"),
    )
    manifest = manifest.merge(
        lineage_counts.reset_index().rename(columns={"_donor_id": "donor_id", "_sample_id": "sample_id"}),
        on=["donor_id", "sample_id"],
        how="left",
    )
    manifest["lineage_cells"] = manifest["lineage_cells"].fillna(0).astype(int)
    manifest["mature_neurons"] = manifest["mature_neurons"].fillna(0).astype(int)
    manifest["has_age"] = manifest["age"].notna()
    manifest["is_healthy"] = manifest["disease_group"].eq("healthy")
    manifest["is_ndd"] = manifest["disease_group"].isin(["ad", "pd", "ndd"])
    manifest["usable_for_ora_training"] = manifest["is_healthy"] & manifest["has_age"]
    return manifest.sort_values(["donor_id", "sample_id"]).reset_index(drop=True)


def parse_age_series(values: pd.Series) -> pd.Series:
    """Parse numeric age from direct numeric or CELLxGENE development-stage strings."""

    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().any():
        return numeric
    extracted = values.astype(str).str.extract(r"(\d+(?:\.\d+)?)\s*-\s*year", flags=re.IGNORECASE)[0]
    return pd.to_numeric(extracted, errors="coerce")


def summarize_cohort(manifest: pd.DataFrame) -> pd.DataFrame:
    """Summarize manifest by disease group."""

    rows = []
    for cohort, frame in [("all", manifest), *list(manifest.groupby("disease_group", dropna=False, observed=True))]:
        label = cohort if isinstance(cohort, str) else str(cohort)
        rows.append(
            {
                "cohort": label,
                "donors": frame["donor_id"].nunique(),
                "samples": frame["sample_id"].nunique(),
                "cells": int(frame["total_cells"].sum()),
                "median_age": frame["age"].median(skipna=True),
                "age_iqr_low": frame["age"].quantile(0.25),
                "age_iqr_high": frame["age"].quantile(0.75),
                "missing_age_samples": int(frame["age"].isna().sum()),
                "lineage_cells": int(frame["lineage_cells"].sum()),
                "mature_neurons": int(frame["mature_neurons"].sum()),
            }
        )
    return pd.DataFrame(rows)


def disease_group(value: object, config: dict[str, Any]) -> str:
    token = normalize_token(value)
    if not token:
        return "unknown"
    healthy = {normalize_token(item) for item in config.get("healthy_values", [])}
    if token in healthy:
        return "healthy"
    ndd_values = config.get("ndd_values", {})
    for group, values in ndd_values.items():
        if token in {normalize_token(item) for item in values}:
            return str(group)
    if "alzheimer" in token:
        return "ad"
    if "parkinson" in token:
        return "pd"
    if token in {"ndd", "neurodegenerative", "neurodegenerative disease"}:
        return "ndd"
    return token


def collection_method_group(value: object, config: dict[str, Any]) -> str:
    token = normalize_token(value)
    if not token:
        return "unknown"
    method_values = config.get("collection_method_values", {})
    for group, values in method_values.items():
        if token in {normalize_token(item) for item in values}:
            return str(group)
    return token


def flatten_lineage_aliases(config: dict[str, Any]) -> dict[str, set[str]]:
    aliases: dict[str, set[str]] = {}
    for key, values in config.get("lineage_cell_types", {}).items():
        aliases[key] = {normalize_token(v) for v in values}
    return aliases


def _first_non_null(values: pd.Series) -> object:
    non_null = values.dropna()
    if non_null.empty:
        return np.nan
    return non_null.iloc[0]
