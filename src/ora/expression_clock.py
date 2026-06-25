"""Leakage-safe donor-level expression aging-clock baselines."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .age_model import (
    BackendInfo,
    _boolean_series,
    _combine_backend_info,
    _performance_row,
    donor_cv_folds,
    fit_model_predictions,
    summarize_feature_stability,
    summarize_repeated_performance,
)
from .diagnostics import summarize_ora_diagnostics
from .utils import ensure_parent


@dataclass
class ExpressionClockResult:
    """Container for expression-clock benchmark outputs."""

    repeat_performance: pd.DataFrame
    performance_summary: pd.DataFrame
    predictions: pd.DataFrame
    feature_stability: pd.DataFrame
    feasibility: pd.DataFrame
    comparison: pd.DataFrame


def run_expression_clock_baseline(
    *,
    counts_path: str | Path,
    metadata_path: str | Path,
    manifest: pd.DataFrame,
    model_config: dict[str, Any],
    models: Iterable[str],
    repeats: int = 10,
    n_pcs: int = 20,
    top_variable_genes: int = 5000,
    min_detection_donors: int = 20,
    chunksize: int = 750,
    comparison_summaries: Iterable[tuple[str, str | Path]] = (),
) -> ExpressionClockResult:
    """Run a fold-internal pseudobulk-expression aging-clock-style benchmark."""

    expression, gene_qc = load_donor_logcpm_expression(
        counts_path,
        metadata_path,
        min_detection_donors=min_detection_donors,
        chunksize=chunksize,
    )
    train = _training_donor_frame(manifest, expression.index)
    if train.shape[0] < 10:
        raise ValueError("At least 10 eligible healthy donors are required for expression-clock CV.")
    x = expression.loc[train["donor_id"].astype(str)].to_numpy(dtype=np.float64)
    y = train["age"].to_numpy(dtype=float)
    model_names = [str(model) for model in models]

    repeat_rows: list[pd.DataFrame] = []
    prediction_rows: list[pd.DataFrame] = []
    importance_rows: list[pd.DataFrame] = []
    base_seed = int(model_config.get("random_seed", 42))
    repeats = max(1, int(repeats))
    for repeat in range(repeats):
        repeat_config = dict(model_config)
        repeat_config["random_seed"] = base_seed + repeat
        folds = donor_cv_folds(train, repeat_config)
        fold_cache = [
            _fit_fold_pcs(x, train_idx, test_idx, n_pcs=n_pcs, top_variable_genes=top_variable_genes)
            for train_idx, test_idx in folds
        ]
        for model_name in model_names:
            pred = np.full(train.shape[0], np.nan, dtype=float)
            backends: list[BackendInfo] = []
            model_importance: list[pd.DataFrame] = []
            for fold_idx, ((train_idx, test_idx), fold_data) in enumerate(zip(folds, fold_cache, strict=True)):
                fold_pred, importance, backend = fit_model_predictions(
                    model_name,
                    fold_data["x_train"],
                    y[train_idx],
                    fold_data["x_test"],
                    repeat_config,
                )
                pred[test_idx] = fold_pred
                backends.append(backend)
                if importance is not None:
                    model_importance.append(
                        pd.DataFrame(
                            {
                                "repeat": repeat,
                                "fold": fold_idx,
                                "model": model_name,
                                "feature": [f"fold_internal_pseudobulk_pc__PC{i + 1:02d}" for i in range(len(importance))],
                                "importance": np.asarray(importance, dtype=float),
                                "pca_variance_explained": fold_data["variance_explained"],
                                "n_selected_genes": fold_data["n_selected_genes"],
                            }
                        )
                    )
            row = _performance_row(model_name, y, pred, _combine_backend_info(model_name, backends))
            row.update(
                {
                    "repeat": repeat,
                    "feature_set": "fold_internal_pseudobulk_expression_pcs",
                    "n_pcs": int(n_pcs),
                    "top_variable_genes": int(top_variable_genes),
                    "min_detection_donors": int(min_detection_donors),
                    "median_pca_variance_explained": float(
                        np.nanmedian([fold["variance_explained"] for fold in fold_cache])
                    ),
                    "median_selected_genes": int(np.nanmedian([fold["n_selected_genes"] for fold in fold_cache])),
                }
            )
            repeat_rows.append(pd.DataFrame([row]))
            prediction_rows.append(
                pd.DataFrame(
                    {
                        "repeat": repeat,
                        "feature_set": "fold_internal_pseudobulk_expression_pcs",
                        "donor_id": train["donor_id"].to_numpy(),
                        "model": model_name,
                        "chronological_age": y,
                        "ora": pred,
                    }
                )
            )
            importance_rows.extend(model_importance)

    repeat_performance = pd.concat(repeat_rows, ignore_index=True)
    predictions = pd.concat(prediction_rows, ignore_index=True)
    stability_input = pd.concat(importance_rows, ignore_index=True) if importance_rows else pd.DataFrame()
    feature_stability = summarize_feature_stability(stability_input)
    summary = _summarize_with_calibration(repeat_performance, predictions, manifest, model_config)
    feasibility = expression_clock_feasibility(
        expression=expression,
        gene_qc=gene_qc,
        train=train,
        n_pcs=n_pcs,
        top_variable_genes=top_variable_genes,
        min_detection_donors=min_detection_donors,
    )
    comparison = compare_expression_clock_to_benchmarks(
        ("fold_internal_expression_clock", summary),
        comparison_summaries=comparison_summaries,
    )
    return ExpressionClockResult(
        repeat_performance=repeat_performance,
        performance_summary=summary,
        predictions=predictions,
        feature_stability=feature_stability,
        feasibility=feasibility,
        comparison=comparison,
    )


def load_donor_logcpm_expression(
    counts_path: str | Path,
    metadata_path: str | Path,
    *,
    min_detection_donors: int = 20,
    chunksize: int = 750,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate pseudobulk groups to donors and return donor x gene logCPM."""

    metadata = pd.read_csv(metadata_path, sep="\t")
    required = {"pseudobulk_id", "donor_id"}
    missing = required.difference(metadata.columns)
    if missing:
        raise KeyError(f"Pseudobulk metadata missing columns: {sorted(missing)}")
    metadata = metadata.dropna(subset=["pseudobulk_id", "donor_id"]).copy()
    metadata["pseudobulk_id"] = metadata["pseudobulk_id"].astype(str)
    metadata["donor_id"] = metadata["donor_id"].astype(str)

    header = pd.read_csv(counts_path, sep="\t", nrows=0).columns.tolist()
    pb_cols = [col for col in header if col in set(metadata["pseudobulk_id"])]
    if not pb_cols:
        raise ValueError("No pseudobulk count columns overlap metadata pseudobulk IDs.")
    pb_to_donor = metadata.drop_duplicates("pseudobulk_id").set_index("pseudobulk_id")["donor_id"].to_dict()
    donors = sorted({pb_to_donor[col] for col in pb_cols})
    donor_groups = [(donor, [col for col in pb_cols if pb_to_donor[col] == donor]) for donor in donors]

    gene_ids: list[str] = []
    gene_symbols: list[str] = []
    parts: list[np.ndarray] = []
    for chunk in pd.read_csv(counts_path, sep="\t", chunksize=max(1, int(chunksize))):
        gene_ids.extend(chunk["gene_id"].astype(str).tolist())
        if "gene_symbol" in chunk:
            gene_symbols.extend(chunk["gene_symbol"].astype(str).tolist())
        else:
            gene_symbols.extend(chunk["gene_id"].astype(str).tolist())
        donor_counts = np.zeros((chunk.shape[0], len(donors)), dtype=np.float64)
        for donor_idx, (_, cols) in enumerate(donor_groups):
            donor_counts[:, donor_idx] = chunk[cols].sum(axis=1).to_numpy(dtype=np.float64)
        parts.append(donor_counts)

    counts = np.vstack(parts).T
    detected = (counts > 0).sum(axis=0)
    keep = detected >= int(min_detection_donors)
    if not np.any(keep):
        raise ValueError("No genes pass the minimum donor-detection threshold.")
    gene_ids_array = np.asarray(gene_ids, dtype=object)[keep]
    gene_symbols_array = np.asarray(gene_symbols, dtype=object)[keep]
    counts = counts[:, keep]
    library_sizes = counts.sum(axis=1)
    safe_library_sizes = np.where(library_sizes > 0, library_sizes, 1.0)
    logcpm = np.log1p(counts / safe_library_sizes[:, None] * 1_000_000.0)
    expression = pd.DataFrame(logcpm, index=pd.Index(donors, name="donor_id"), columns=gene_ids_array)
    gene_qc = pd.DataFrame(
        {
            "gene_id": gene_ids_array,
            "gene_symbol": gene_symbols_array,
            "detected_donors": detected[keep],
            "mean_logcpm": logcpm.mean(axis=0),
            "variance_logcpm": logcpm.var(axis=0),
        }
    ).sort_values("variance_logcpm", ascending=False)
    return expression, gene_qc.reset_index(drop=True)


def expression_clock_feasibility(
    *,
    expression: pd.DataFrame,
    gene_qc: pd.DataFrame,
    train: pd.DataFrame,
    n_pcs: int,
    top_variable_genes: int,
    min_detection_donors: int,
) -> pd.DataFrame:
    """Summarize direct clock feasibility and the implemented fallback."""

    return pd.DataFrame(
        [
            {
                "baseline": "public_scageclock_direct",
                "status": "not_run_no_locked_gateway_compatible_public_implementation",
                "n_training_donors": int(train.shape[0]),
                "n_expression_genes": int(expression.shape[1]),
                "n_pcs": "",
                "top_variable_genes": "",
                "min_detection_donors": "",
                "notes": "No clean public scAgeClock-style implementation/checkpoint was identified for direct Gateway inference; use tissue-specific expression baseline as the prespecified fallback.",
            },
            {
                "baseline": "fold_internal_pseudobulk_expression_pcs",
                "status": "ok",
                "n_training_donors": int(train.shape[0]),
                "n_expression_genes": int(expression.shape[1]),
                "n_pcs": int(n_pcs),
                "top_variable_genes": int(min(top_variable_genes, expression.shape[1])),
                "min_detection_donors": int(min_detection_donors),
                "notes": f"Fold-internal variance selection, scaling, and PCA; {gene_qc.shape[0]} genes pass detection filter.",
            },
        ]
    )


def write_expression_clock_outputs(
    result: ExpressionClockResult,
    *,
    summary_out: str | Path,
    repeat_performance_out: str | Path,
    scores_out: str | Path,
    feature_stability_out: str | Path,
    feasibility_out: str | Path,
    comparison_out: str | Path,
) -> None:
    """Write expression-clock benchmark result tables."""

    result.performance_summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.repeat_performance.to_csv(ensure_parent(repeat_performance_out), sep="\t", index=False)
    result.predictions.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    result.feature_stability.to_csv(ensure_parent(feature_stability_out), sep="\t", index=False)
    result.feasibility.to_csv(ensure_parent(feasibility_out), sep="\t", index=False)
    result.comparison.to_csv(ensure_parent(comparison_out), sep="\t", index=False)


def _fit_fold_pcs(
    expression: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    *,
    n_pcs: int,
    top_variable_genes: int,
) -> dict[str, Any]:
    from sklearn.decomposition import PCA  # type: ignore

    x_train_raw = expression[train_idx]
    x_test_raw = expression[test_idx]
    variances = np.nanvar(x_train_raw, axis=0)
    finite = np.isfinite(variances)
    if not np.any(finite):
        raise ValueError("No finite expression variances in training fold.")
    candidate = np.flatnonzero(finite)
    order = candidate[np.argsort(variances[candidate])[::-1]]
    selected = order[: max(1, min(int(top_variable_genes), order.size))]
    train_selected = x_train_raw[:, selected]
    test_selected = x_test_raw[:, selected]
    means = np.nanmean(train_selected, axis=0)
    stds = np.nanstd(train_selected, axis=0)
    stds = np.where(stds > 0, stds, 1.0)
    train_scaled = (train_selected - means) / stds
    test_scaled = (test_selected - means) / stds
    components = max(1, min(int(n_pcs), train_scaled.shape[0] - 1, train_scaled.shape[1]))
    pca = PCA(n_components=components, random_state=0)
    x_train = pca.fit_transform(train_scaled)
    x_test = pca.transform(test_scaled)
    return {
        "x_train": x_train,
        "x_test": x_test,
        "n_selected_genes": int(selected.size),
        "variance_explained": float(np.sum(pca.explained_variance_ratio_)),
    }


def _training_donor_frame(manifest: pd.DataFrame, available_donors: pd.Index) -> pd.DataFrame:
    donor_meta = (
        manifest.sort_values(["donor_id", "sample_id"] if "sample_id" in manifest else ["donor_id"])
        .drop_duplicates("donor_id")
        .copy()
    )
    donor_meta["donor_id"] = donor_meta["donor_id"].astype(str)
    donor_meta["age"] = pd.to_numeric(donor_meta["age"], errors="coerce")
    mask = donor_meta["donor_id"].isin(set(available_donors.astype(str))) & donor_meta["age"].notna()
    if "usable_for_ora_training" in donor_meta:
        mask = mask & _boolean_series(donor_meta["usable_for_ora_training"])
    elif "disease_group" in donor_meta:
        mask = mask & donor_meta["disease_group"].astype(str).str.lower().eq("healthy")
    return donor_meta.loc[mask].sort_values("donor_id").reset_index(drop=True)


def _summarize_with_calibration(
    repeat_performance: pd.DataFrame,
    predictions: pd.DataFrame,
    manifest: pd.DataFrame,
    model_config: dict[str, Any],
) -> pd.DataFrame:
    summary = summarize_repeated_performance(repeat_performance)
    calibration_rows = []
    for repeat, scores in predictions.groupby("repeat", observed=True):
        cal = summarize_ora_diagnostics(scores, model_config=model_config, manifest=manifest).calibration
        cal.insert(0, "repeat", repeat)
        calibration_rows.append(cal)
    if calibration_rows:
        cal_frame = pd.concat(calibration_rows, ignore_index=True)
        cal_summary = (
            cal_frame.groupby("model", observed=True)
            .agg(
                calibration_slope_mean=("calibration_slope_ora_on_age", "mean"),
                calibration_slope_sd=("calibration_slope_ora_on_age", "std"),
                calibration_slope_ci_low=(
                    "calibration_slope_ora_on_age",
                    lambda s: float(np.nanquantile(s, 0.025)),
                ),
                calibration_slope_ci_high=(
                    "calibration_slope_ora_on_age",
                    lambda s: float(np.nanquantile(s, 0.975)),
                ),
            )
            .reset_index()
        )
        summary = summary.merge(cal_summary, on="model", how="left")
    if "feature_set" in repeat_performance:
        summary.insert(0, "feature_set", repeat_performance["feature_set"].dropna().astype(str).iloc[0])
    for col in ["n_pcs", "top_variable_genes", "min_detection_donors", "median_pca_variance_explained"]:
        if col in repeat_performance:
            summary[col] = repeat_performance.groupby("model", observed=True)[col].median().reindex(
                summary["model"]
            ).to_numpy()
    return summary


def compare_expression_clock_to_benchmarks(
    expression_summary: tuple[str, pd.DataFrame],
    *,
    comparison_summaries: Iterable[tuple[str, str | Path]] = (),
) -> pd.DataFrame:
    """Rank expression-clock summary against existing benchmark summaries."""

    rows = [_best_row(expression_summary[0], expression_summary[1])]
    for label, path in comparison_summaries:
        table = pd.read_csv(path, sep="\t")
        rows.append(_best_row(label, table))
    return pd.DataFrame(rows).sort_values("mae_mean").reset_index(drop=True)


def _best_row(label: str, table: pd.DataFrame) -> dict[str, Any]:
    if table.empty or "mae_mean" not in table:
        return {"benchmark": label, "status": "missing_mae", "model": "", "mae_mean": np.nan}
    best = table.sort_values("mae_mean").iloc[0]
    return {
        "benchmark": label,
        "status": "ok",
        "model": best.get("model", ""),
        "feature_set": best.get("feature_set", label),
        "mae_mean": best.get("mae_mean", np.nan),
        "rmse_mean": best.get("rmse_mean", np.nan),
        "r2_mean": best.get("r2_mean", np.nan),
        "spearman_r_mean": best.get("spearman_r_mean", np.nan),
        "calibration_slope_mean": best.get("calibration_slope_mean", np.nan),
    }
