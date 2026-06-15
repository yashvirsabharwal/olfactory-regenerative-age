"""Donor-level permutation/null tests for ORA age models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .age_model import model_names_from_config, summarize_repeated_performance, train_ora_models_repeated


@dataclass
class PermutationNullResult:
    repeat_performance: pd.DataFrame
    permutation_summary: pd.DataFrame
    empirical_summary: pd.DataFrame


def run_permutation_null(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    model_config: dict[str, Any] | None = None,
    *,
    n_permutations: int = 100,
    repeats: int = 2,
    random_seed: int = 20260615,
    observed_summary: pd.DataFrame | None = None,
) -> PermutationNullResult:
    """Shuffle healthy training-donor ages and compare observed repeated-CV metrics to null."""

    model_config = dict(model_config or {})
    n_permutations = max(1, int(n_permutations))
    repeats = max(1, int(repeats))
    rng = np.random.default_rng(int(random_seed))
    model_names = model_names_from_config(model_config)
    if observed_summary is None or observed_summary.empty:
        observed = train_ora_models_repeated(features, manifest, model_config, repeats=repeats).performance_summary
    else:
        observed = observed_summary.copy()
    observed = observed[observed["model"].astype(str).isin(model_names)].copy()

    performance_rows = []
    summary_rows = []
    for permutation in range(n_permutations):
        permuted = permute_training_ages(manifest, rng)
        perm_config = dict(model_config)
        perm_config["random_seed"] = int(random_seed) + 1000 + permutation
        result = train_ora_models_repeated(features, permuted, perm_config, repeats=repeats)
        repeat_perf = result.repeat_performance.copy()
        repeat_perf.insert(0, "permutation", permutation)
        performance_rows.append(repeat_perf)
        summary = summarize_repeated_performance(repeat_perf)
        summary.insert(0, "permutation", permutation)
        summary_rows.append(summary)

    repeat_performance = pd.concat(performance_rows, ignore_index=True)
    permutation_summary = pd.concat(summary_rows, ignore_index=True)
    empirical = summarize_permutation_null(observed, permutation_summary)
    return PermutationNullResult(
        repeat_performance=repeat_performance,
        permutation_summary=permutation_summary,
        empirical_summary=empirical,
    )


def permute_training_ages(manifest: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Permute ages only among donors eligible for ORA training."""

    output = manifest.copy()
    donor_meta = output.sort_values(["donor_id", "sample_id"] if "sample_id" in output else ["donor_id"])
    donor_meta = donor_meta.drop_duplicates("donor_id")
    eligible = _boolean_series(donor_meta["usable_for_ora_training"]) & donor_meta["age"].notna()
    donor_ids = donor_meta.loc[eligible, "donor_id"].astype(str).to_numpy()
    ages = donor_meta.loc[eligible, "age"].to_numpy()
    if donor_ids.size < 2:
        raise ValueError("At least two eligible training donors are required for age permutation.")
    permuted = ages.copy()
    rng.shuffle(permuted)
    age_map = dict(zip(donor_ids, permuted, strict=True))
    mask = output["donor_id"].astype(str).isin(age_map)
    output.loc[mask, "age"] = output.loc[mask, "donor_id"].astype(str).map(age_map)
    return output


def summarize_permutation_null(observed_summary: pd.DataFrame, permutation_summary: pd.DataFrame) -> pd.DataFrame:
    """Compute empirical p-values for observed metrics against shuffled-label nulls."""

    rows = []
    observed = observed_summary.copy()
    permutations = permutation_summary.copy()
    for metric in ["mae_mean", "rmse_mean", "r2_mean", "spearman_r_mean"]:
        if metric in observed:
            observed[metric] = pd.to_numeric(observed[metric], errors="coerce")
        if metric in permutations:
            permutations[metric] = pd.to_numeric(permutations[metric], errors="coerce")
    for model, obs in observed.groupby("model", observed=True, sort=False):
        null = permutations[permutations["model"].astype(str).eq(str(model))]
        if null.empty:
            continue
        row = {"model": model, "n_permutations": int(null["permutation"].nunique()) if "permutation" in null else int(null.shape[0])}
        obs_row = obs.iloc[0]
        for metric, label, direction in [
            ("mae_mean", "mae", "lower"),
            ("rmse_mean", "rmse", "lower"),
            ("r2_mean", "r2", "higher"),
            ("spearman_r_mean", "spearman_r", "higher"),
        ]:
            observed_value = float(obs_row.get(metric, np.nan))
            null_values = pd.to_numeric(null.get(metric), errors="coerce").dropna().to_numpy(dtype=float)
            row[f"observed_{label}"] = observed_value
            row[f"null_{label}_mean"] = float(np.mean(null_values)) if null_values.size else np.nan
            row[f"null_{label}_sd"] = float(np.std(null_values, ddof=1)) if null_values.size > 1 else 0.0
            row[f"null_{label}_ci_low"] = _quantile(null_values, 0.025)
            row[f"null_{label}_ci_high"] = _quantile(null_values, 0.975)
            row[f"empirical_p_{label}"] = _empirical_p(observed_value, null_values, direction)
            row[f"z_{label}"] = _z_score(observed_value, null_values)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("observed_mae").reset_index(drop=True)


def _empirical_p(observed_value: float, null_values: np.ndarray, direction: str) -> float:
    if not np.isfinite(observed_value) or null_values.size == 0:
        return np.nan
    if direction == "lower":
        hits = np.sum(null_values <= observed_value)
    else:
        hits = np.sum(null_values >= observed_value)
    return float((1 + hits) / (1 + null_values.size))


def _z_score(observed_value: float, null_values: np.ndarray) -> float:
    if not np.isfinite(observed_value) or null_values.size < 2:
        return np.nan
    sd = float(np.std(null_values, ddof=1))
    if sd == 0.0 or not np.isfinite(sd):
        return np.nan
    return float((observed_value - float(np.mean(null_values))) / sd)


def _quantile(values: np.ndarray, q: float) -> float:
    if values.size == 0:
        return np.nan
    return float(np.quantile(values, q))


def _boolean_series(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"true", "1", "yes", "y"})
