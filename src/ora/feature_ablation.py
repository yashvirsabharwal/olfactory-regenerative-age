"""Feature-family ablation benchmarks for donor-level ORA models."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .age_model import train_ora_models_repeated
from .diagnostics import summarize_ora_diagnostics
from .features import (
    COMPOSITION_PREFIXES,
    MODULE_PREFIXES,
    SCVI_GLOBAL_PREFIXES,
    SCVI_STATE_PREFIXES,
)
from .permutation import run_permutation_null
from .utils import ensure_parent


FEATURE_SET_ORDER = [
    "null_model",
    "technical_covariates_only",
    "proportions_only",
    "clr_only",
    "lineage_ratios_only",
    "composition_without_ratios",
    "modules_only",
    "composition_plus_modules",
    "scvi_global_only",
    "scvi_cell_state_only",
    "scvi_donor_embedding",
    "pseudobulk_expression_pcs",
    "ora_scvi_hybrid",
]


@dataclass
class FeatureFamilyAblationResult:
    feature_matrices: dict[str, pd.DataFrame]
    feasibility: pd.DataFrame
    summary: pd.DataFrame
    deltas: pd.DataFrame
    repeat_performance: pd.DataFrame
    scores: pd.DataFrame


def run_feature_family_ablation(
    *,
    feature_matrix: pd.DataFrame,
    manifest: pd.DataFrame,
    model_config: dict[str, Any],
    output_dir: str | Path = "data/processed/feature_family_ablation",
    pseudobulk_counts_path: str | Path | None = None,
    pseudobulk_metadata_path: str | Path | None = None,
    models: list[str] | None = None,
    repeats: int | None = 10,
    n_permutations: int = 0,
    permutation_repeats: int = 1,
    random_seed: int = 20260624,
    baseline_feature_set: str = "ora_scvi_hybrid",
) -> FeatureFamilyAblationResult:
    """Build matched feature-family matrices and run repeated donor-level CV."""

    output_dir = Path(output_dir)
    families, feasibility = build_feature_family_matrices(
        feature_matrix,
        manifest,
        pseudobulk_counts_path=pseudobulk_counts_path,
        pseudobulk_metadata_path=pseudobulk_metadata_path,
        random_seed=random_seed,
    )
    model_config = dict(model_config or {})
    if models:
        model_config["model_names"] = models

    summary_rows = []
    repeat_rows = []
    score_rows = []
    for feature_set, matrix in families.items():
        matrix_path = output_dir / f"{feature_set}.tsv"
        matrix.to_csv(ensure_parent(matrix_path), sep="\t", index=False)
        config = dict(model_config)
        if feature_set == "null_model":
            config["model_names"] = ["null_model"]
        result = train_ora_models_repeated(matrix, manifest, config, repeats=repeats)
        repeated = result.repeat_performance.copy()
        repeated.insert(0, "feature_set", feature_set)
        repeat_rows.append(repeated)

        scores = result.predictions.copy()
        scores.insert(0, "feature_set", feature_set)
        score_rows.append(scores)

        summary = result.performance_summary.copy()
        summary.insert(0, "feature_set", feature_set)
        calibration = summarize_repeated_calibration(scores, manifest, model_config)
        summary = summary.merge(calibration, on=["feature_set", "model"], how="left")
        summary["n_features"] = matrix.shape[1] - 1
        summary_rows.append(summary)

    summary = pd.concat(summary_rows, ignore_index=True)
    summary = _add_permutation_p_values(
        summary,
        families=families,
        manifest=manifest,
        model_config=model_config,
        n_permutations=n_permutations,
        permutation_repeats=permutation_repeats,
        random_seed=random_seed,
    )
    summary = _ordered_summary(summary)
    deltas = feature_family_deltas(summary, baseline_feature_set=baseline_feature_set)
    return FeatureFamilyAblationResult(
        feature_matrices=families,
        feasibility=feasibility,
        summary=summary,
        deltas=deltas,
        repeat_performance=pd.concat(repeat_rows, ignore_index=True),
        scores=pd.concat(score_rows, ignore_index=True),
    )


def build_feature_family_matrices(
    feature_matrix: pd.DataFrame,
    manifest: pd.DataFrame,
    *,
    pseudobulk_counts_path: str | Path | None = None,
    pseudobulk_metadata_path: str | Path | None = None,
    random_seed: int = 20260624,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Return donor-level matrices for each available ablation feature family."""

    if "donor_id" not in feature_matrix:
        raise KeyError("Feature matrix must include donor_id.")
    rows: list[dict[str, object]] = []
    matrices: dict[str, pd.DataFrame] = {}

    def add(label: str, columns: list[str], note: str = "") -> None:
        if not columns:
            rows.append(_feasibility_row(label, "skipped_no_features", 0, note or "No matching columns."))
            return
        matrix = feature_matrix[["donor_id", *columns]].copy()
        matrices[label] = matrix
        rows.append(_feasibility_row(label, "ok", len(columns), note))

    columns = [col for col in feature_matrix.columns if col != "donor_id"]
    prop_cols = [col for col in columns if col.startswith("prop__")]
    clr_cols = [col for col in columns if col.startswith("clr__")]
    ratio_cols = [col for col in columns if col.startswith("ratio__")]
    module_cols = [col for col in columns if col.startswith(MODULE_PREFIXES)]
    composition_cols = [col for col in columns if col.startswith(COMPOSITION_PREFIXES)]
    scvi_global_cols = [col for col in columns if col.startswith(SCVI_GLOBAL_PREFIXES)]
    scvi_state_cols = [col for col in columns if col.startswith(SCVI_STATE_PREFIXES)]
    scvi_cols = [*scvi_global_cols, *scvi_state_cols]

    null_matrix = feature_matrix[["donor_id"]].copy()
    null_matrix["null_feature"] = 0.0
    matrices["null_model"] = null_matrix
    rows.append(_feasibility_row("null_model", "ok", 1, "Mean-age null model with a dummy feature."))

    technical = build_technical_covariate_matrix(manifest)
    matrices["technical_covariates_only"] = technical
    rows.append(_feasibility_row("technical_covariates_only", "ok", technical.shape[1] - 1, "Manifest-derived technical and demographic covariates."))

    add("proportions_only", prop_cols)
    add("clr_only", clr_cols)
    add("lineage_ratios_only", ratio_cols)
    add("composition_without_ratios", [*prop_cols, *clr_cols])
    add("modules_only", module_cols)
    add("composition_plus_modules", [*composition_cols, *module_cols])
    add("scvi_global_only", scvi_global_cols)
    add("scvi_cell_state_only", scvi_state_cols)
    add("scvi_donor_embedding", scvi_cols)

    pseudobulk = maybe_build_pseudobulk_pc_features(
        pseudobulk_counts_path,
        pseudobulk_metadata_path,
        random_seed=random_seed,
    )
    if pseudobulk is None:
        rows.append(
            _feasibility_row(
                "pseudobulk_expression_pcs",
                "skipped_missing_inputs",
                0,
                "Pseudobulk count/metadata inputs were not provided or are unavailable.",
            )
        )
    else:
        matrices["pseudobulk_expression_pcs"] = pseudobulk
        rows.append(_feasibility_row("pseudobulk_expression_pcs", "ok", pseudobulk.shape[1] - 1, "Donor-level PCA on aggregated genomewide pseudobulk logCPM."))

    if composition_cols and module_cols and scvi_cols:
        matrices["ora_scvi_hybrid"] = feature_matrix[["donor_id", *composition_cols, *module_cols, *scvi_cols]].copy()
        rows.append(_feasibility_row("ora_scvi_hybrid", "ok", len(composition_cols) + len(module_cols) + len(scvi_cols), "Composition, modules, and donor scVI embeddings."))
    else:
        rows.append(_feasibility_row("ora_scvi_hybrid", "skipped_no_features", 0, "Requires composition, module, and scVI columns."))

    ordered = {label: matrices[label] for label in FEATURE_SET_ORDER if label in matrices}
    return ordered, pd.DataFrame(rows)


def build_technical_covariate_matrix(manifest: pd.DataFrame) -> pd.DataFrame:
    """Build numeric donor-level technical/demographic covariates."""

    if "donor_id" not in manifest:
        raise KeyError("Manifest must include donor_id.")
    donor_meta = manifest.sort_values(["donor_id", "sample_id"] if "sample_id" in manifest else ["donor_id"])
    donor_meta = donor_meta.drop_duplicates("donor_id").copy()
    output = donor_meta[["donor_id"]].copy()
    numeric_cols = ["total_cells", "lineage_cells", "mature_neurons"]
    for col in numeric_cols:
        if col in donor_meta:
            values = pd.to_numeric(donor_meta[col], errors="coerce").fillna(0.0)
            output[f"technical__log10_{col}"] = np.log10(values + 1.0)
    categorical_cols = ["sex", "race_ethnicity", "chemistry", "collection_method", "site"]
    for col in categorical_cols:
        if col not in donor_meta:
            continue
        values = donor_meta[col].fillna("unknown").astype(str)
        dummies = pd.get_dummies(values, prefix=f"technical__{col}", dtype=float)
        output = pd.concat([output, dummies], axis=1)
    if output.shape[1] == 1:
        output["technical__dummy"] = 0.0
    return output


def maybe_build_pseudobulk_pc_features(
    counts_path: str | Path | None,
    metadata_path: str | Path | None,
    *,
    n_components: int = 20,
    n_variable_genes: int = 2000,
    random_seed: int = 20260624,
) -> pd.DataFrame | None:
    """Build donor-level expression PCs from genomewide pseudobulk counts if available."""

    if counts_path in {None, ""} or metadata_path in {None, ""}:
        return None
    counts_path = Path(counts_path)
    metadata_path = Path(metadata_path)
    if not counts_path.exists() or not metadata_path.exists():
        return None
    counts = pd.read_csv(counts_path, sep="\t")
    metadata = pd.read_csv(metadata_path, sep="\t")
    if "pseudobulk_id" not in metadata or "donor_id" not in metadata:
        return None
    pb_cols = [col for col in counts.columns if col in set(metadata["pseudobulk_id"].astype(str))]
    if len(pb_cols) < 2:
        return None
    metadata = metadata.drop_duplicates("pseudobulk_id").set_index("pseudobulk_id")
    pb_cols = [col for col in pb_cols if col in metadata.index]
    donors = sorted(metadata.loc[pb_cols, "donor_id"].astype(str).unique())
    if len(donors) < 3:
        return None

    matrix = counts[pb_cols].to_numpy(dtype=np.float32, copy=True)
    donor_counts = np.zeros((matrix.shape[0], len(donors)), dtype=np.float32)
    donor_index = {donor: idx for idx, donor in enumerate(donors)}
    pb_donor_codes = metadata.loc[pb_cols, "donor_id"].astype(str).map(donor_index).to_numpy(dtype=int)
    for pb_idx, donor_code in enumerate(pb_donor_codes):
        donor_counts[:, donor_code] += matrix[:, pb_idx]

    gene_totals = donor_counts.sum(axis=1)
    keep = gene_totals > 0
    donor_counts = donor_counts[keep]
    if donor_counts.shape[0] < 2:
        return None
    library_sizes = donor_counts.sum(axis=0)
    library_sizes[library_sizes <= 0] = 1.0
    log_cpm = np.log1p((donor_counts / library_sizes[None, :]) * 1_000_000.0).T
    variances = np.var(log_cpm, axis=0)
    keep_genes = min(int(n_variable_genes), log_cpm.shape[1])
    variable_idx = np.argsort(variances)[-keep_genes:]
    x = log_cpm[:, variable_idx]
    x = (x - x.mean(axis=0)) / np.where(x.std(axis=0) == 0, 1.0, x.std(axis=0))

    from sklearn.decomposition import PCA

    n_components = min(int(n_components), x.shape[0] - 1, x.shape[1])
    if n_components < 1:
        return None
    pcs = PCA(n_components=n_components, random_state=int(random_seed)).fit_transform(x)
    output = pd.DataFrame({"donor_id": donors})
    for idx in range(n_components):
        output[f"pseudobulk_pc__PC{idx + 1:02d}"] = pcs[:, idx]
    return output


def summarize_repeated_calibration(scores: pd.DataFrame, manifest: pd.DataFrame, model_config: dict[str, Any]) -> pd.DataFrame:
    """Summarize calibration slope over repeated out-of-fold predictions."""

    rows = []
    for (feature_set, repeat), repeat_scores in scores.groupby(["feature_set", "repeat"], observed=True):
        calibration = summarize_ora_diagnostics(repeat_scores, model_config=model_config, manifest=manifest).calibration
        calibration.insert(0, "repeat", repeat)
        calibration.insert(0, "feature_set", feature_set)
        rows.append(calibration)
    if not rows:
        return pd.DataFrame(columns=["feature_set", "model", "calibration_slope_mean"])
    frame = pd.concat(rows, ignore_index=True)
    summary = (
        frame.groupby(["feature_set", "model"], observed=True)
        .agg(
            calibration_slope_mean=("calibration_slope_ora_on_age", "mean"),
            calibration_slope_sd=("calibration_slope_ora_on_age", "std"),
            calibration_slope_ci_low=("calibration_slope_ora_on_age", lambda s: float(np.nanquantile(s, 0.025))),
            calibration_slope_ci_high=("calibration_slope_ora_on_age", lambda s: float(np.nanquantile(s, 0.975))),
        )
        .reset_index()
    )
    return summary


def feature_family_deltas(summary: pd.DataFrame, *, baseline_feature_set: str = "ora_scvi_hybrid") -> pd.DataFrame:
    """Compute best-model deltas against the selected baseline and overall best family."""

    best = (
        summary.sort_values(["mae_mean", "rmse_mean", "feature_set", "model"])
        .groupby("feature_set", observed=True, as_index=False)
        .head(1)
        .copy()
    )
    overall = best.sort_values("mae_mean").iloc[0]
    baseline = best[best["feature_set"].eq(baseline_feature_set)]
    baseline_row = baseline.iloc[0] if not baseline.empty else overall
    rows = []
    for _, row in best.iterrows():
        rows.append(
            {
                "feature_set": row["feature_set"],
                "best_model": row["model"],
                "best_mae_mean": row["mae_mean"],
                "best_mae_ci_low": row.get("mae_ci_low"),
                "best_mae_ci_high": row.get("mae_ci_high"),
                "best_rmse_mean": row.get("rmse_mean"),
                "best_r2_mean": row.get("r2_mean"),
                "best_spearman_r_mean": row.get("spearman_r_mean"),
                "baseline_feature_set": baseline_row["feature_set"],
                "delta_mae_vs_baseline": float(row["mae_mean"] - baseline_row["mae_mean"]),
                "overall_best_feature_set": overall["feature_set"],
                "delta_mae_vs_overall_best": float(row["mae_mean"] - overall["mae_mean"]),
                "n_features": row.get("n_features"),
                "permutation_p_mae": row.get("permutation_p_mae"),
            }
        )
    return pd.DataFrame(rows).sort_values(["best_mae_mean", "feature_set"]).reset_index(drop=True)


def write_ablation_figure(summary: pd.DataFrame, output_pdf: str | Path, output_png: str | Path | None = None) -> None:
    """Write a compact performance figure for best model in each feature family."""

    import matplotlib.pyplot as plt

    best = (
        summary.sort_values(["mae_mean", "rmse_mean", "feature_set", "model"])
        .groupby("feature_set", observed=True, as_index=False)
        .head(1)
        .copy()
        .sort_values("mae_mean", ascending=True)
    )
    best["err_low"] = best["mae_mean"] - pd.to_numeric(best["mae_ci_low"], errors="coerce")
    best["err_high"] = pd.to_numeric(best["mae_ci_high"], errors="coerce") - best["mae_mean"]
    fig_height = max(4.5, 0.34 * best.shape[0] + 1.4)
    fig, ax = plt.subplots(figsize=(8, fig_height))
    y = np.arange(best.shape[0])
    ax.barh(y, best["mae_mean"], color="#4C78A8")
    ax.errorbar(best["mae_mean"], y, xerr=[best["err_low"], best["err_high"]], fmt="none", ecolor="#222222", capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels([_label_feature_set(label) for label in best["feature_set"]])
    ax.invert_yaxis()
    ax.set_xlabel("Repeated-CV MAE (years)")
    ax.set_title("Feature-family ablation")
    ax.grid(axis="x", alpha=0.25)
    for idx, (_, row) in enumerate(best.iterrows()):
        ax.text(float(row["mae_mean"]) + 0.08, idx, str(row["model"]).replace("_", " "), va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(ensure_parent(output_pdf))
    if output_png:
        fig.savefig(ensure_parent(output_png), dpi=200)
    plt.close(fig)


def _add_permutation_p_values(
    summary: pd.DataFrame,
    *,
    families: dict[str, pd.DataFrame],
    manifest: pd.DataFrame,
    model_config: dict[str, Any],
    n_permutations: int,
    permutation_repeats: int,
    random_seed: int,
) -> pd.DataFrame:
    output = summary.copy()
    output["permutation_p_mae"] = np.nan
    output["permutation_n"] = 0
    output["permutation_scope"] = "not_run"
    if n_permutations <= 0:
        return output

    best_rows = (
        output.sort_values(["feature_set", "mae_mean", "rmse_mean"])
        .groupby("feature_set", observed=True, as_index=False)
        .head(1)
    )
    for offset, (_, best) in enumerate(best_rows.iterrows()):
        feature_set = str(best["feature_set"])
        if feature_set == "null_model" or feature_set not in families:
            continue
        model_name = str(best["model"])
        perm_config = dict(model_config)
        perm_config["model_names"] = [model_name]
        observed = output[output["feature_set"].eq(feature_set) & output["model"].astype(str).eq(model_name)]
        result = run_permutation_null(
            families[feature_set],
            manifest,
            perm_config,
            n_permutations=n_permutations,
            repeats=permutation_repeats,
            random_seed=random_seed + 100 * (offset + 1),
            observed_summary=observed,
        )
        empirical = result.empirical_summary
        if empirical.empty:
            continue
        p = float(empirical["empirical_p_mae"].iloc[0])
        mask = output["feature_set"].eq(feature_set) & output["model"].astype(str).eq(model_name)
        output.loc[mask, "permutation_p_mae"] = p
        output.loc[mask, "permutation_n"] = int(n_permutations)
        output.loc[mask, "permutation_scope"] = "best_model_only"
    return output


def _ordered_summary(summary: pd.DataFrame) -> pd.DataFrame:
    output = summary.copy()
    order = {label: idx for idx, label in enumerate(FEATURE_SET_ORDER)}
    output["_feature_set_order"] = output["feature_set"].map(order).fillna(len(order))
    output["mae_rank_overall"] = output["mae_mean"].rank(method="min", ascending=True).astype(int)
    output["mae_rank_within_feature_set"] = output.groupby("feature_set", observed=True)["mae_mean"].rank(method="min", ascending=True).astype(int)
    output["is_best_within_feature_set"] = output["mae_rank_within_feature_set"].eq(1)
    output["is_best_overall"] = output["mae_rank_overall"].eq(1)
    return output.sort_values(["mae_mean", "_feature_set_order", "model"]).drop(columns=["_feature_set_order"]).reset_index(drop=True)


def _feasibility_row(feature_set: str, status: str, n_features: int, note: str) -> dict[str, object]:
    return {"feature_set": feature_set, "status": status, "n_features": int(n_features), "note": note}


def _label_feature_set(label: object) -> str:
    return str(label).replace("_", " ")
