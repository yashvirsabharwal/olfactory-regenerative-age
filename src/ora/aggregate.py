"""Donor-level cell-state aggregation and feature engineering."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .metadata import ColumnMap, resolve_columns
from .utils import normalize_token


DEFAULT_PSEUDOCOUNT = 0.5


def aggregate_cell_counts(obs: pd.DataFrame, config: dict[str, Any], columns: ColumnMap | None = None) -> pd.DataFrame:
    """Count cells by donor, sample, and annotated cell state."""

    columns = columns or resolve_columns(list(obs.columns), config)
    work = pd.DataFrame(
        {
            "donor_id": obs[columns.donor_id].astype(str).to_numpy(),
            "sample_id": obs[columns.sample_id].astype(str).to_numpy(),
            "coarse_cell_type": obs[columns.coarse_cell_type].astype(str).to_numpy(),
            "fine_cell_type": obs[columns.fine_cell_type].astype(str).to_numpy(),
        }
    )
    counts = (
        work.groupby(["donor_id", "sample_id", "coarse_cell_type", "fine_cell_type"], observed=True)
        .size()
        .rename("cell_count")
        .reset_index()
        .sort_values(["donor_id", "sample_id", "coarse_cell_type", "fine_cell_type"])
    )
    return counts.reset_index(drop=True)


def build_cell_state_features(
    counts: pd.DataFrame,
    config: dict[str, Any],
    pseudocount: float = DEFAULT_PSEUDOCOUNT,
) -> pd.DataFrame:
    """Create donor-level composition, CLR, and lineage-ratio features."""

    required = {"donor_id", "fine_cell_type", "cell_count"}
    missing = sorted(required.difference(counts.columns))
    if missing:
        raise KeyError(f"Cell count table missing columns: {missing}")

    donor_counts = (
        counts.groupby(["donor_id", "fine_cell_type"], observed=True)["cell_count"]
        .sum()
        .unstack(fill_value=0)
        .sort_index(axis=1)
    )
    total_cells = donor_counts.sum(axis=1)
    proportions = donor_counts.div(total_cells.replace(0, np.nan), axis=0).fillna(0)
    safe = proportions + pseudocount / np.maximum(total_cells.to_numpy()[:, None], 1)
    clr = np.log(safe).sub(np.log(safe).mean(axis=1), axis=0)

    feature_frames = [
        pd.DataFrame({"donor_id": donor_counts.index, "total_cells": total_cells.to_numpy()}),
        _prefixed(proportions, "prop__"),
        _prefixed(clr, "clr__"),
    ]
    ratios = compute_lineage_ratios(donor_counts, config, pseudocount=pseudocount)
    if not ratios.empty:
        feature_frames.append(ratios)

    features = feature_frames[0].reset_index(drop=True)
    for frame in feature_frames[1:]:
        features = features.merge(frame.reset_index(drop=True), on="donor_id", how="left")
    return features


def compute_lineage_ratios(
    donor_counts: pd.DataFrame,
    config: dict[str, Any],
    pseudocount: float = DEFAULT_PSEUDOCOUNT,
) -> pd.DataFrame:
    """Compute biologically motivated ORA lineage ratios."""

    alias_map = _alias_to_columns(donor_counts.columns, config)
    if not alias_map:
        return pd.DataFrame({"donor_id": donor_counts.index})

    groups = {
        key: _sum_alias_group(donor_counts, aliases)
        for key, aliases in config.get("lineage_cell_types", {}).items()
    }
    hbc = groups.get("quiescent_hbc", 0) + groups.get("activated_hbc", 0) + groups.get("cycling_hbc", 0)
    inp = groups.get("early_inp", 0) + groups.get("late_inp", 0)
    iosn = groups.get("early_iosn", 0) + groups.get("late_iosn", 0)
    mature = groups.get("early_mature_mosn", 0) + groups.get("fully_mature_mosn", 0)
    stressed = groups.get("stressed_mosn", 0)
    mosn = mature + stressed
    lineage = hbc + inp + iosn + mosn
    total = donor_counts.sum(axis=1)

    output = pd.DataFrame(index=donor_counts.index)
    output.index.name = "donor_id"
    output["ratio__neuronal_fraction"] = (inp + iosn + mosn) / (total + pseudocount)
    output["ratio__mature_neuron_fraction"] = mosn / (total + pseudocount)
    output["ratio__immature_to_mature"] = iosn / (mosn + pseudocount)
    output["ratio__progenitor_to_neuron"] = inp / (iosn + mosn + pseudocount)
    output["ratio__activated_to_quiescent_hbc"] = groups.get("activated_hbc", 0) / (
        groups.get("quiescent_hbc", 0) + pseudocount
    )
    output["ratio__hbc_to_inp"] = inp / (groups.get("activated_hbc", 0) + pseudocount)
    output["ratio__inp_to_iosn"] = iosn / (inp + pseudocount)
    output["ratio__iosn_to_mosn"] = mature / (iosn + pseudocount)
    output["ratio__stressed_to_mature_mosn"] = stressed / (mature + pseudocount)
    output["ratio__lineage_fraction"] = lineage / (total + pseudocount)
    return output.reset_index()


def _prefixed(frame: pd.DataFrame, prefix: str) -> pd.DataFrame:
    clean = frame.copy()
    clean.columns = [prefix + _slugify(col) for col in clean.columns]
    clean.insert(0, "donor_id", clean.index)
    return clean


def _sum_alias_group(donor_counts: pd.DataFrame, aliases: list[str]) -> pd.Series:
    alias_tokens = {normalize_token(alias) for alias in aliases}
    matching = [col for col in donor_counts.columns if normalize_token(col) in alias_tokens]
    if not matching:
        return pd.Series(0.0, index=donor_counts.index)
    return donor_counts[matching].sum(axis=1)


def _alias_to_columns(columns: pd.Index, config: dict[str, Any]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for key, aliases in config.get("lineage_cell_types", {}).items():
        tokens = {normalize_token(alias) for alias in aliases}
        result[key] = [col for col in columns if normalize_token(col) in tokens]
    return result


def _slugify(value: object) -> str:
    text = normalize_token(value)
    return text.replace("/", " ").replace("+", " plus ").replace(" ", "_")
