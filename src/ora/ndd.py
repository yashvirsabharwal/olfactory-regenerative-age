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
        "xgboost",
        "lightgbm",
        "catboost",
        "boosted_ensemble",
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


def compare_ndd_feature_sets(
    projections: dict[str, pd.DataFrame],
    *,
    baseline: str = "composition",
    augmented: str = "augmented",
    diseases: Iterable[str] = ("ad", "pd"),
    negative_epsilon: float = 1e-6,
) -> pd.DataFrame:
    """Compare disease ORAA means across frozen-projection feature sets."""

    summaries = []
    disease_set = {str(disease) for disease in diseases}
    for feature_set, projection in projections.items():
        if projection is None or projection.empty:
            continue
        frame = projection.copy()
        frame["feature_set"] = str(feature_set)
        frame["disease_group"] = frame["disease_group"].astype(str)
        frame = frame[frame["disease_group"].isin(disease_set)].copy()
        if frame.empty:
            continue
        for col in ["ora", "oraa", "n_features"]:
            if col in frame:
                frame[col] = pd.to_numeric(frame[col], errors="coerce")
        summary = (
            frame.groupby(["feature_set", "model", "disease_group"], observed=True, dropna=False)
            .agg(
                donors=("donor_id", "nunique"),
                mean_ora=("ora", "mean"),
                mean_oraa=("oraa", "mean"),
                sd_oraa=("oraa", lambda s: float(pd.to_numeric(s, errors="coerce").std(ddof=0))),
                n_features=("n_features", "max"),
            )
            .reset_index()
        )
        summaries.append(summary)
    if not summaries:
        return pd.DataFrame()

    summary = pd.concat(summaries, ignore_index=True)
    rows = []
    for (model, disease), group in summary.groupby(["model", "disease_group"], observed=True, dropna=False):
        by_feature = {str(row.feature_set): row for row in group.itertuples(index=False)}
        base = by_feature.get(baseline)
        aug = by_feature.get(augmented)
        base_oraa = _row_value(base, "mean_oraa")
        aug_oraa = _row_value(aug, "mean_oraa")
        rows.append(
            {
                "model": model,
                "disease_group": disease,
                f"{baseline}_donors": _row_value(base, "donors"),
                f"{augmented}_donors": _row_value(aug, "donors"),
                f"{baseline}_n_features": _row_value(base, "n_features"),
                f"{augmented}_n_features": _row_value(aug, "n_features"),
                f"{baseline}_mean_ora": _row_value(base, "mean_ora"),
                f"{augmented}_mean_ora": _row_value(aug, "mean_ora"),
                f"{baseline}_mean_oraa": base_oraa,
                f"{augmented}_mean_oraa": aug_oraa,
                f"{augmented}_minus_{baseline}_oraa": aug_oraa - base_oraa
                if np.isfinite(base_oraa) and np.isfinite(aug_oraa)
                else np.nan,
                "sign_stable_negative": bool(
                    np.isfinite(base_oraa)
                    and np.isfinite(aug_oraa)
                    and base_oraa < -abs(float(negative_epsilon))
                    and aug_oraa < -abs(float(negative_epsilon))
                ),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["disease_group", "model"]).reset_index(drop=True)


def donor_projection_appendix(
    projection: pd.DataFrame,
    *,
    feature_set: str | None = None,
    models: Iterable[str] | None = None,
    diseases: Iterable[str] = ("ad", "pd"),
) -> pd.DataFrame:
    """Return donor-level AD/PD projection rows for appendix export."""

    if projection is None or projection.empty:
        return pd.DataFrame()
    frame = projection.copy()
    if feature_set is not None and "feature_set" not in frame:
        frame.insert(0, "feature_set", str(feature_set))
    elif "feature_set" not in frame:
        frame.insert(0, "feature_set", "unknown")
    disease_set = {str(disease) for disease in diseases}
    frame["disease_group"] = frame["disease_group"].astype(str)
    frame = frame[frame["disease_group"].isin(disease_set)].copy()
    if models is not None:
        model_set = {str(model) for model in models}
        frame = frame[frame["model"].astype(str).isin(model_set)].copy()
    if frame.empty:
        return pd.DataFrame()
    columns = [
        "feature_set",
        "donor_id",
        "disease_group",
        "sex",
        "race_ethnicity",
        "chronological_age",
        "chemistry",
        "collection_method",
        "site",
        "total_cells",
        "model",
        "ora",
        "oraa",
        "training_n",
        "n_features",
    ]
    for col in columns:
        if col not in frame:
            frame[col] = np.nan
    for col in ["chronological_age", "total_cells", "ora", "oraa", "training_n", "n_features"]:
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame[columns].sort_values(["feature_set", "disease_group", "donor_id", "model"]).reset_index(drop=True)


def ndd_projection_diagnostics(
    projection: pd.DataFrame,
    *,
    models: Iterable[str] | None = None,
    diseases: Iterable[str] = ("ad", "pd"),
    min_donors_ok: int = 2,
) -> pd.DataFrame:
    """Summarize AD/PD projected ORAA by key covariate and yield strata."""

    if projection is None or projection.empty:
        return pd.DataFrame()
    frame = projection.copy()
    frame["model"] = frame["model"].astype(str)
    frame["disease_group"] = frame["disease_group"].astype(str)
    disease_set = {str(disease) for disease in diseases}
    frame = frame[frame["disease_group"].isin(disease_set)].copy()
    if models is not None:
        model_set = {str(model) for model in models}
        frame = frame[frame["model"].isin(model_set)].copy()
    if frame.empty:
        return pd.DataFrame()

    for col in ["chronological_age", "total_cells", "ora", "oraa"]:
        if col in frame:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame["age_bin"] = _age_bin(frame.get("chronological_age", pd.Series(np.nan, index=frame.index)))
    frame["cell_yield_quartile"] = _yield_quartile(frame)
    diagnostics = ["sex", "age_bin", "chemistry", "collection_method", "site", "cell_yield_quartile"]
    rows = []
    for diagnostic in diagnostics:
        if diagnostic not in frame:
            continue
        work = frame.copy()
        work[diagnostic] = _clean_level(work[diagnostic])
        grouped = work.groupby(["model", "disease_group", diagnostic], observed=True, dropna=False)
        for (model, disease, level), group in grouped:
            n_donors = int(group["donor_id"].nunique())
            rows.append(
                {
                    "model": model,
                    "disease_group": disease,
                    "diagnostic": diagnostic,
                    "level": level,
                    "n_donors": n_donors,
                    "mean_age": float(group["chronological_age"].mean()) if "chronological_age" in group else np.nan,
                    "min_age": float(group["chronological_age"].min()) if "chronological_age" in group else np.nan,
                    "max_age": float(group["chronological_age"].max()) if "chronological_age" in group else np.nan,
                    "median_total_cells": float(group["total_cells"].median()) if "total_cells" in group else np.nan,
                    "mean_ora": float(group["ora"].mean()) if "ora" in group else np.nan,
                    "mean_oraa": float(group["oraa"].mean()) if "oraa" in group else np.nan,
                    "sd_oraa": float(group["oraa"].std(ddof=0)) if "oraa" in group else np.nan,
                    "donor_ids": ",".join(sorted(group["donor_id"].astype(str).unique())),
                    "status": "ok" if n_donors >= int(min_donors_ok) else "single_donor_stratum",
                }
            )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["disease_group", "diagnostic", "level", "model"]).reset_index(drop=True)


def ndd_label_permutation(
    projection: pd.DataFrame,
    *,
    models: Iterable[str] | None = None,
    diseases: Iterable[str] = ("ad", "pd"),
    strata: Iterable[str] = ("chemistry", "collection_method"),
    n_permutations: int = 5000,
    random_seed: int = 42,
) -> pd.DataFrame:
    """Permutation test on frozen ORAA scores within compatible disease/reference strata."""

    if projection is None or projection.empty:
        return pd.DataFrame()
    frame = projection.copy()
    frame["model"] = frame["model"].astype(str)
    frame["disease_group"] = frame["disease_group"].astype(str)
    frame["oraa"] = pd.to_numeric(frame["oraa"], errors="coerce")
    disease_set = {str(disease) for disease in diseases}
    if models is not None:
        model_set = {str(model) for model in models}
        frame = frame[frame["model"].isin(model_set)].copy()
    rows = []
    rng = np.random.default_rng(int(random_seed))
    strata_cols = [col for col in strata if col in frame]
    for model, model_frame in frame.groupby("model", observed=True):
        healthy = model_frame[model_frame["disease_group"].eq("healthy") & model_frame["oraa"].notna()].copy()
        for disease in sorted(disease_set):
            disease_frame = model_frame[model_frame["disease_group"].eq(disease) & model_frame["oraa"].notna()].copy()
            if disease_frame.empty or healthy.empty:
                continue
            reference = _matched_healthy_reference(healthy, disease_frame)
            eligible = pd.concat([disease_frame, reference], ignore_index=True)
            eligible = eligible.dropna(subset=["oraa"])
            n_disease = int(disease_frame["donor_id"].nunique())
            n_reference = int(reference["donor_id"].nunique())
            if eligible["donor_id"].nunique() < n_disease + 1 or n_reference == 0:
                rows.append(_empty_permutation_row(model, disease, n_disease, n_reference, "too_few_matched_donors"))
                continue
            observed = _disease_reference_difference(disease_frame["oraa"], reference["oraa"])
            null = _permuted_differences(eligible, n_disease, n_permutations, rng)
            rows.append(
                {
                    "model": model,
                    "disease_group": disease,
                    "strata": ",".join(strata_cols) if strata_cols else "none",
                    "n_disease": n_disease,
                    "n_reference": n_reference,
                    "observed_difference_vs_reference": observed,
                    "null_mean": float(np.mean(null)) if null.size else np.nan,
                    "null_ci_low": _quantile(null, 0.025),
                    "null_ci_high": _quantile(null, 0.975),
                    "empirical_p_negative": float((np.sum(null <= observed) + 1) / (null.size + 1)) if null.size else np.nan,
                    "n_permutations": int(null.size),
                    "status": "ok" if null.size else "no_null",
                }
            )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["disease_group", "model"]).reset_index(drop=True)


def _matched_healthy_reference(healthy: pd.DataFrame, disease: pd.DataFrame) -> pd.DataFrame:
    reference = healthy.copy()
    for col in ["chemistry", "collection_method"]:
        if col in reference and col in disease:
            values = set(disease[col].dropna().astype(str))
            if values:
                reference = reference[reference[col].astype(str).isin(values)]
    return reference


def _disease_reference_difference(disease_values: pd.Series, reference_values: pd.Series) -> float:
    disease_arr = pd.to_numeric(disease_values, errors="coerce").dropna().to_numpy(dtype=float)
    reference_arr = pd.to_numeric(reference_values, errors="coerce").dropna().to_numpy(dtype=float)
    if disease_arr.size == 0 or reference_arr.size == 0:
        return float("nan")
    return float(disease_arr.mean() - reference_arr.mean())


def _permuted_differences(
    eligible: pd.DataFrame,
    n_disease: int,
    n_permutations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    donor_frame = eligible.sort_values("donor_id").drop_duplicates("donor_id").copy()
    values = pd.to_numeric(donor_frame["oraa"], errors="coerce").to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    if values.size <= n_disease or n_disease <= 0:
        return np.array([], dtype=float)
    out = np.empty(int(n_permutations), dtype=float)
    for i in range(int(n_permutations)):
        disease_idx = rng.choice(values.size, size=n_disease, replace=False)
        mask = np.zeros(values.size, dtype=bool)
        mask[disease_idx] = True
        out[i] = values[mask].mean() - values[~mask].mean()
    return out


def _empty_permutation_row(model: str, disease: str, n_disease: int, n_reference: int, status: str) -> dict[str, object]:
    return {
        "model": model,
        "disease_group": disease,
        "strata": "chemistry,collection_method",
        "n_disease": n_disease,
        "n_reference": n_reference,
        "observed_difference_vs_reference": np.nan,
        "null_mean": np.nan,
        "null_ci_low": np.nan,
        "null_ci_high": np.nan,
        "empirical_p_negative": np.nan,
        "n_permutations": 0,
        "status": status,
    }


def _row_value(row: object | None, field: str) -> float:
    if row is None:
        return float("nan")
    value = getattr(row, field)
    if pd.isna(value):
        return float("nan")
    return float(value)


def _age_bin(age: pd.Series) -> pd.Series:
    values = pd.to_numeric(age, errors="coerce")
    bins = [-np.inf, 59, 69, 79, np.inf]
    labels = ["lt60", "60_69", "70_79", "80_plus"]
    return pd.cut(values, bins=bins, labels=labels).astype("object").fillna("missing")


def _yield_quartile(frame: pd.DataFrame) -> pd.Series:
    if "total_cells" not in frame:
        return pd.Series("missing", index=frame.index, dtype=object)
    donor_yield = frame[["donor_id", "total_cells"]].drop_duplicates("donor_id").copy()
    values = pd.to_numeric(donor_yield["total_cells"], errors="coerce")
    labels = ["q1_low", "q2", "q3", "q4_high"]
    try:
        donor_yield["cell_yield_quartile"] = pd.qcut(values, q=4, labels=labels, duplicates="drop").astype("object")
    except ValueError:
        donor_yield["cell_yield_quartile"] = "unbinned"
    donor_yield["cell_yield_quartile"] = donor_yield["cell_yield_quartile"].fillna("missing")
    return frame["donor_id"].map(donor_yield.set_index("donor_id")["cell_yield_quartile"]).fillna("missing")


def _clean_level(series: pd.Series) -> pd.Series:
    cleaned = series.fillna("missing").astype(str).str.strip()
    return cleaned.mask(cleaned.eq(""), "missing")


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
