"""NDD ORA projection uncertainty and context summaries."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd


@dataclass
class NDDUncertaintyResult:
    uncertainty: pd.DataFrame
    context: pd.DataFrame


def summarize_ndd_projection_uncertainty(
    projection: pd.DataFrame,
    *,
    models: Iterable[str] = (
        "ridge",
        "lasso",
        "elastic_net",
        "random_forest",
        "extra_trees",
        "gradient_boosting",
        "tree_ensemble",
    ),
    diseases: Iterable[str] = ("ad", "pd"),
    n_bootstrap: int = 5000,
    random_seed: int = 42,
) -> NDDUncertaintyResult:
    """Bootstrap NDD ORAA means against all-healthy and matched healthy references."""

    frame = projection.copy()
    frame["disease_group"] = frame["disease_group"].astype(str)
    frame["model"] = frame["model"].astype(str)
    frame["oraa"] = pd.to_numeric(frame["oraa"], errors="coerce")
    rng = np.random.default_rng(int(random_seed))
    rows = []
    for model in models:
        model_frame = frame[frame["model"].eq(str(model))].copy()
        healthy_all = model_frame[model_frame["disease_group"].eq("healthy") & model_frame["oraa"].notna()]
        for disease in diseases:
            disease_frame = model_frame[model_frame["disease_group"].eq(str(disease)) & model_frame["oraa"].notna()]
            if disease_frame.empty:
                continue
            matched = _matched_healthy_reference(healthy_all, disease_frame)
            for label, reference in [("all_healthy", healthy_all), ("matched_healthy", matched)]:
                rows.append(
                    _uncertainty_row(
                        model=str(model),
                        disease=str(disease),
                        reference_label=label,
                        disease_values=disease_frame["oraa"].to_numpy(dtype=float),
                        reference_values=reference["oraa"].to_numpy(dtype=float),
                        n_bootstrap=n_bootstrap,
                        rng=rng,
                    )
                )
    return NDDUncertaintyResult(
        uncertainty=pd.DataFrame(rows),
        context=ndd_projection_context(frame),
    )


def ndd_projection_context(projection: pd.DataFrame) -> pd.DataFrame:
    """Summarize disease, chemistry, and collection-method confounding context."""

    frame = projection.drop_duplicates(["donor_id", "model"]).copy()
    # Context is model-invariant, so keep one row per donor.
    frame = frame.sort_values(["donor_id", "model"]).drop_duplicates("donor_id")
    group_cols = ["disease_group", "chemistry", "collection_method"]
    available = [col for col in group_cols if col in frame]
    if not available:
        return pd.DataFrame()
    return (
        frame.groupby(available, observed=True, dropna=False)
        .agg(
            donors=("donor_id", "nunique"),
            mean_age=("chronological_age", "mean"),
            median_total_cells=("total_cells", "median"),
        )
        .reset_index()
        .sort_values(["disease_group", "donors"], ascending=[True, False])
    )


def _matched_healthy_reference(healthy: pd.DataFrame, disease: pd.DataFrame) -> pd.DataFrame:
    reference = healthy.copy()
    for col in ["chemistry", "collection_method"]:
        if col in reference and col in disease:
            values = set(disease[col].dropna().astype(str))
            if values:
                reference = reference[reference[col].astype(str).isin(values)]
    return reference


def _uncertainty_row(
    *,
    model: str,
    disease: str,
    reference_label: str,
    disease_values: np.ndarray,
    reference_values: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> dict[str, object]:
    disease_mean = float(np.mean(disease_values)) if disease_values.size else np.nan
    reference_mean = float(np.mean(reference_values)) if reference_values.size else np.nan
    disease_boot = _bootstrap_means(disease_values, n_bootstrap, rng)
    reference_boot = _bootstrap_means(reference_values, n_bootstrap, rng)
    diff_boot = disease_boot - reference_boot if disease_boot.size and reference_boot.size else np.array([], dtype=float)
    return {
        "model": model,
        "disease_group": disease,
        "reference": reference_label,
        "n_disease": int(disease_values.size),
        "n_reference": int(reference_values.size),
        "mean_oraa": disease_mean,
        "mean_oraa_ci_low": _quantile(disease_boot, 0.025),
        "mean_oraa_ci_high": _quantile(disease_boot, 0.975),
        "reference_mean_oraa": reference_mean,
        "difference_vs_reference": disease_mean - reference_mean if np.isfinite(reference_mean) else np.nan,
        "difference_ci_low": _quantile(diff_boot, 0.025),
        "difference_ci_high": _quantile(diff_boot, 0.975),
        "p_directional_negative": float(np.mean(diff_boot >= 0)) if diff_boot.size else np.nan,
    }


def _bootstrap_means(values: np.ndarray, n_bootstrap: int, rng: np.random.Generator) -> np.ndarray:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.array([], dtype=float)
    draws = rng.integers(0, values.size, size=(int(n_bootstrap), values.size))
    return values[draws].mean(axis=1)


def _quantile(values: np.ndarray, q: float) -> float:
    if values.size == 0:
        return float("nan")
    return float(np.quantile(values, q))
