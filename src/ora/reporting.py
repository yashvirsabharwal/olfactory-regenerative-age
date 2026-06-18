"""MVP report generation for ORA analysis outputs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ora.utils import ensure_parent


FIGURE_NAMES = {
    "cohort": "mvp_cohort_overview.png",
    "performance": "mvp_model_performance.png",
    "performance_comparison": "mvp_model_comparison.png",
    "associations": "mvp_top_age_associations.png",
    "predictions": "mvp_predicted_vs_age.png",
    "calibrated_predictions": "mvp_calibrated_predicted_vs_age.png",
    "calibration_age_bins": "mvp_calibration_age_bins.png",
    "importance": "mvp_feature_importance.png",
    "ndd_projection": "mvp_ndd_projection.png",
    "ndd_matched_reference": "mvp_ndd_matched_reference.png",
    "module_scores": "mvp_module_scores.png",
    "pseudobulk_de": "mvp_pseudobulk_de.png",
    "pseudobulk_covariate_de": "mvp_pseudobulk_covariate_de.png",
}


def generate_mvp_report(
    *,
    manifest: pd.DataFrame,
    cohort_summary: pd.DataFrame,
    associations: pd.DataFrame,
    performance: pd.DataFrame,
    scores: pd.DataFrame,
    importance: pd.DataFrame,
    out_md: str | Path,
    figure_dir: str | Path,
    augmented_performance: pd.DataFrame | None = None,
    augmented_scores: pd.DataFrame | None = None,
    augmented_importance: pd.DataFrame | None = None,
    ora_calibration: pd.DataFrame | None = None,
    ora_age_bin_errors: pd.DataFrame | None = None,
    ora_residual_diagnostics: pd.DataFrame | None = None,
    ora_calibrated_scores: pd.DataFrame | None = None,
    ndd_projection: pd.DataFrame | None = None,
    ndd_projection_summary: pd.DataFrame | None = None,
    ndd_projection_uncertainty: pd.DataFrame | None = None,
    ndd_projection_context: pd.DataFrame | None = None,
    ndd_projection_feature_comparison: pd.DataFrame | None = None,
    ndd_projection_donor_appendix: pd.DataFrame | None = None,
    ndd_projection_diagnostics: pd.DataFrame | None = None,
    module_summary: pd.DataFrame | None = None,
    module_coverage: pd.DataFrame | None = None,
    donor_module_features: pd.DataFrame | None = None,
    external_validation_summary: pd.DataFrame | None = None,
    external_gene_list_coverage: pd.DataFrame | None = None,
    external_feature_contract: pd.DataFrame | None = None,
    external_sample_metadata: pd.DataFrame | None = None,
    external_10x_sample_qc: pd.DataFrame | None = None,
    external_10x_module_contrasts: pd.DataFrame | None = None,
    external_10x_marker_composition: pd.DataFrame | None = None,
    external_10x_marker_contrasts: pd.DataFrame | None = None,
    external_marker_age_concordance: pd.DataFrame | None = None,
    external_validation_evidence: pd.DataFrame | None = None,
    pseudobulk_de: pd.DataFrame | None = None,
    pseudobulk_coverage: pd.DataFrame | None = None,
    pseudobulk_metadata: pd.DataFrame | None = None,
    pseudobulk_covariate_de: pd.DataFrame | None = None,
    pseudobulk_genomewide_summary: pd.DataFrame | None = None,
    pseudobulk_genomewide_qc_summary: pd.DataFrame | None = None,
    pseudobulk_genomewide_gene_qc: pd.DataFrame | None = None,
    pseudobulk_genomewide_disease_summary: pd.DataFrame | None = None,
    pseudobulk_genomewide_de_summary: pd.DataFrame | None = None,
    pseudobulk_genomewide_de_top_hits: pd.DataFrame | None = None,
    pseudobulk_genomewide_de_audit: pd.DataFrame | None = None,
    pseudobulk_genomewide_donor_balance: pd.DataFrame | None = None,
    pseudobulk_genomewide_matched_feasibility: pd.DataFrame | None = None,
    pseudobulk_genomewide_de_summary_matched: pd.DataFrame | None = None,
    pseudobulk_genomewide_de_top_hits_matched: pd.DataFrame | None = None,
    pseudobulk_genomewide_de_audit_matched: pd.DataFrame | None = None,
    pseudobulk_genomewide_limma_de_summary: pd.DataFrame | None = None,
    pseudobulk_genomewide_limma_de_top_hits: pd.DataFrame | None = None,
    pseudobulk_genomewide_limma_de_audit: pd.DataFrame | None = None,
    pseudobulk_genomewide_limma_de_summary_matched: pd.DataFrame | None = None,
    pseudobulk_genomewide_limma_de_top_hits_matched: pd.DataFrame | None = None,
    pseudobulk_genomewide_limma_de_audit_matched: pd.DataFrame | None = None,
    ora_sensitivity_scenarios: pd.DataFrame | None = None,
    ora_sensitivity_performance: pd.DataFrame | None = None,
    ora_repeated_cv_summary: pd.DataFrame | None = None,
    ora_repeated_cv_feature_stability: pd.DataFrame | None = None,
    ora_augmented_candidate_repeated_cv_summary: pd.DataFrame | None = None,
    ora_feature_interpretation: pd.DataFrame | None = None,
    ora_feature_set_model_comparison: pd.DataFrame | None = None,
    ora_permutation_empirical: pd.DataFrame | None = None,
    ora_nested_tuning_summary: pd.DataFrame | None = None,
    ora_stacking_summary: pd.DataFrame | None = None,
    ora_model_card: pd.DataFrame | None = None,
    ndd_label_permutation: pd.DataFrame | None = None,
    latent_space_readiness: pd.DataFrame | None = None,
    latent_space_local_audit: pd.DataFrame | None = None,
    latent_space_portal_assets: pd.DataFrame | None = None,
    latent_recompute_feasibility: pd.DataFrame | None = None,
    scvi_pilot_validation: pd.DataFrame | None = None,
    output_provenance: pd.DataFrame | None = None,
    source: dict[str, Any] | None = None,
    paper_defaults: dict[str, Any] | None = None,
    schema: dict[str, Any] | None = None,
    top_n: int = 12,
) -> list[Path]:
    """Write a markdown report and supporting PNG figures."""

    figure_path = Path(figure_dir)
    figure_path.mkdir(parents=True, exist_ok=True)
    figure_paths = _write_figures(
        cohort_summary=cohort_summary,
        associations=associations,
        performance=performance,
        scores=scores,
        importance=importance,
        ora_calibrated_scores=ora_calibrated_scores,
        ora_age_bin_errors=ora_age_bin_errors,
        ora_calibration=ora_calibration,
        augmented_performance=augmented_performance,
        ndd_projection=ndd_projection,
        module_summary=module_summary,
        pseudobulk_de=pseudobulk_de,
        pseudobulk_covariate_de=pseudobulk_covariate_de,
        figure_dir=figure_path,
        top_n=top_n,
    )
    markdown = render_mvp_markdown(
        manifest=manifest,
        cohort_summary=cohort_summary,
        associations=associations,
        performance=performance,
        importance=importance,
        ora_calibration=ora_calibration,
        ora_age_bin_errors=ora_age_bin_errors,
        ora_residual_diagnostics=ora_residual_diagnostics,
        ora_calibrated_scores=ora_calibrated_scores,
        augmented_performance=augmented_performance,
        augmented_importance=augmented_importance,
        ndd_projection=ndd_projection,
        ndd_projection_summary=ndd_projection_summary,
        ndd_projection_uncertainty=ndd_projection_uncertainty,
        ndd_projection_context=ndd_projection_context,
        ndd_projection_feature_comparison=ndd_projection_feature_comparison,
        ndd_projection_donor_appendix=ndd_projection_donor_appendix,
        ndd_projection_diagnostics=ndd_projection_diagnostics,
        module_summary=module_summary,
        module_coverage=module_coverage,
        donor_module_features=donor_module_features,
        external_validation_summary=external_validation_summary,
        external_gene_list_coverage=external_gene_list_coverage,
        external_feature_contract=external_feature_contract,
        external_sample_metadata=external_sample_metadata,
        external_10x_sample_qc=external_10x_sample_qc,
        external_10x_module_contrasts=external_10x_module_contrasts,
        external_10x_marker_composition=external_10x_marker_composition,
        external_10x_marker_contrasts=external_10x_marker_contrasts,
        external_marker_age_concordance=external_marker_age_concordance,
        external_validation_evidence=external_validation_evidence,
        pseudobulk_de=pseudobulk_de,
        pseudobulk_coverage=pseudobulk_coverage,
        pseudobulk_metadata=pseudobulk_metadata,
        pseudobulk_covariate_de=pseudobulk_covariate_de,
        pseudobulk_genomewide_summary=pseudobulk_genomewide_summary,
        pseudobulk_genomewide_qc_summary=pseudobulk_genomewide_qc_summary,
        pseudobulk_genomewide_gene_qc=pseudobulk_genomewide_gene_qc,
        pseudobulk_genomewide_disease_summary=pseudobulk_genomewide_disease_summary,
        pseudobulk_genomewide_de_summary=pseudobulk_genomewide_de_summary,
        pseudobulk_genomewide_de_top_hits=pseudobulk_genomewide_de_top_hits,
        pseudobulk_genomewide_de_audit=pseudobulk_genomewide_de_audit,
        pseudobulk_genomewide_donor_balance=pseudobulk_genomewide_donor_balance,
        pseudobulk_genomewide_matched_feasibility=pseudobulk_genomewide_matched_feasibility,
        pseudobulk_genomewide_de_summary_matched=pseudobulk_genomewide_de_summary_matched,
        pseudobulk_genomewide_de_top_hits_matched=pseudobulk_genomewide_de_top_hits_matched,
        pseudobulk_genomewide_de_audit_matched=pseudobulk_genomewide_de_audit_matched,
        pseudobulk_genomewide_limma_de_summary=pseudobulk_genomewide_limma_de_summary,
        pseudobulk_genomewide_limma_de_top_hits=pseudobulk_genomewide_limma_de_top_hits,
        pseudobulk_genomewide_limma_de_audit=pseudobulk_genomewide_limma_de_audit,
        pseudobulk_genomewide_limma_de_summary_matched=pseudobulk_genomewide_limma_de_summary_matched,
        pseudobulk_genomewide_limma_de_top_hits_matched=pseudobulk_genomewide_limma_de_top_hits_matched,
        pseudobulk_genomewide_limma_de_audit_matched=pseudobulk_genomewide_limma_de_audit_matched,
        ora_sensitivity_scenarios=ora_sensitivity_scenarios,
        ora_sensitivity_performance=ora_sensitivity_performance,
        ora_repeated_cv_summary=ora_repeated_cv_summary,
        ora_repeated_cv_feature_stability=ora_repeated_cv_feature_stability,
        ora_augmented_candidate_repeated_cv_summary=ora_augmented_candidate_repeated_cv_summary,
        ora_feature_interpretation=ora_feature_interpretation,
        ora_feature_set_model_comparison=ora_feature_set_model_comparison,
        ora_permutation_empirical=ora_permutation_empirical,
        ora_nested_tuning_summary=ora_nested_tuning_summary,
        ora_stacking_summary=ora_stacking_summary,
        ora_model_card=ora_model_card,
        ndd_label_permutation=ndd_label_permutation,
        latent_space_readiness=latent_space_readiness,
        latent_space_local_audit=latent_space_local_audit,
        latent_space_portal_assets=latent_space_portal_assets,
        latent_recompute_feasibility=latent_recompute_feasibility,
        scvi_pilot_validation=scvi_pilot_validation,
        output_provenance=output_provenance,
        out_md=out_md,
        figure_paths=figure_paths,
        source=source or {},
        paper_defaults=paper_defaults or {},
        schema=schema or {},
        top_n=top_n,
    )
    report_path = ensure_parent(out_md)
    report_path.write_text(markdown, encoding="utf-8")
    return [report_path, *figure_paths.values()]


def load_schema(path: str | Path | None) -> dict[str, Any]:
    if not path:
        return {}
    schema_path = Path(path)
    if not schema_path.exists():
        return {}
    return json.loads(schema_path.read_text(encoding="utf-8"))


def rank_associations(associations: pd.DataFrame, top_n: int = 12) -> pd.DataFrame:
    """Return the top age associations by FDR, then p-value."""

    required = {"feature", "status", "p_value", "fdr", "beta_per_10_years"}
    if associations.empty or not required.issubset(associations.columns):
        return pd.DataFrame(columns=list(required))
    frame = associations[associations["status"].eq("ok")].copy()
    frame["p_value"] = pd.to_numeric(frame["p_value"], errors="coerce")
    frame["fdr"] = pd.to_numeric(frame["fdr"], errors="coerce")
    frame["beta_per_10_years"] = pd.to_numeric(frame["beta_per_10_years"], errors="coerce")
    frame = frame[np.isfinite(frame["p_value"]) & np.isfinite(frame["fdr"])]
    return frame.sort_values(["fdr", "p_value", "feature"]).head(top_n).reset_index(drop=True)


def rank_pseudobulk_de(pseudobulk_de: pd.DataFrame | None, top_n: int = 12) -> pd.DataFrame:
    """Return the top pseudobulk DE rows by FDR, then p-value."""

    required = {
        "contrast",
        "fine_cell_type",
        "gene",
        "n_case",
        "n_control",
        "log2fc",
        "p_value",
        "fdr",
        "status",
    }
    if pseudobulk_de is None or pseudobulk_de.empty or not required.issubset(pseudobulk_de.columns):
        return pd.DataFrame(columns=list(required))
    frame = pseudobulk_de[pseudobulk_de["status"].eq("ok")].copy()
    for col in ["n_case", "n_control", "log2fc", "p_value", "fdr"]:
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame = frame[np.isfinite(frame["p_value"]) & np.isfinite(frame["fdr"])]
    return frame.sort_values(["fdr", "p_value", "contrast", "fine_cell_type", "gene"]).head(top_n).reset_index(drop=True)


def rank_pseudobulk_covariate_de(pseudobulk_de: pd.DataFrame | None, top_n: int = 12) -> pd.DataFrame:
    """Return top covariate-adjusted pseudobulk DE rows by FDR, then p-value."""

    required = {
        "contrast",
        "fine_cell_type",
        "gene",
        "n_case",
        "n_control",
        "log2fc_adjusted",
        "p_value",
        "fdr",
        "status",
    }
    if pseudobulk_de is None or pseudobulk_de.empty or not required.issubset(pseudobulk_de.columns):
        return pd.DataFrame(columns=list(required))
    frame = pseudobulk_de[pseudobulk_de["status"].eq("ok")].copy()
    for col in ["n_case", "n_control", "log2fc_adjusted", "p_value", "fdr"]:
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame = frame[np.isfinite(frame["p_value"]) & np.isfinite(frame["fdr"])]
    return frame.sort_values(["fdr", "p_value", "contrast", "fine_cell_type", "gene"]).head(top_n).reset_index(drop=True)


def best_predictive_model(performance: pd.DataFrame) -> pd.Series:
    """Return the best non-null model row by MAE, falling back to the overall best row."""

    if performance.empty:
        return pd.Series(dtype=object)
    frame = performance.copy()
    frame["mae"] = pd.to_numeric(frame["mae"], errors="coerce")
    candidates = frame[~frame["model"].eq("null_model") & frame["mae"].notna()]
    if candidates.empty:
        candidates = frame[frame["mae"].notna()]
    if candidates.empty:
        return pd.Series(dtype=object)
    return candidates.sort_values(["mae", "model"]).iloc[0]


def combined_performance(
    performance: pd.DataFrame,
    augmented_performance: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Combine baseline and module-augmented performance tables."""

    pieces = []
    if performance is not None and not performance.empty:
        base = performance.copy()
        base.insert(0, "feature_set", "composition")
        pieces.append(base)
    if augmented_performance is not None and not augmented_performance.empty:
        aug = augmented_performance.copy()
        aug.insert(0, "feature_set", "composition_plus_modules")
        pieces.append(aug)
    if not pieces:
        return pd.DataFrame()
    return pd.concat(pieces, ignore_index=True)


def _sort_metric_table(frame: pd.DataFrame | None, metric: str, *, ascending: bool = True) -> pd.DataFrame:
    if frame is None:
        return pd.DataFrame()
    if frame.empty or metric not in frame:
        return frame
    out = frame.copy()
    out[metric] = pd.to_numeric(out[metric], errors="coerce")
    sort_cols = [metric]
    if "feature_set" in out:
        sort_cols.append("feature_set")
    if "model" in out:
        sort_cols.append("model")
    return out.sort_values(sort_cols, ascending=[ascending] + [True] * (len(sort_cols) - 1)).reset_index(drop=True)


def render_mvp_markdown(
    *,
    manifest: pd.DataFrame,
    cohort_summary: pd.DataFrame,
    associations: pd.DataFrame,
    performance: pd.DataFrame,
    importance: pd.DataFrame,
    ora_calibration: pd.DataFrame | None,
    ora_age_bin_errors: pd.DataFrame | None,
    ora_residual_diagnostics: pd.DataFrame | None,
    ora_calibrated_scores: pd.DataFrame | None,
    augmented_performance: pd.DataFrame | None,
    augmented_importance: pd.DataFrame | None,
    ndd_projection: pd.DataFrame | None,
    ndd_projection_summary: pd.DataFrame | None,
    ndd_projection_uncertainty: pd.DataFrame | None,
    ndd_projection_context: pd.DataFrame | None,
    ndd_projection_feature_comparison: pd.DataFrame | None,
    ndd_projection_donor_appendix: pd.DataFrame | None,
    ndd_projection_diagnostics: pd.DataFrame | None,
    module_summary: pd.DataFrame | None,
    module_coverage: pd.DataFrame | None,
    donor_module_features: pd.DataFrame | None,
    external_validation_summary: pd.DataFrame | None,
    external_gene_list_coverage: pd.DataFrame | None,
    external_feature_contract: pd.DataFrame | None,
    external_sample_metadata: pd.DataFrame | None,
    external_10x_sample_qc: pd.DataFrame | None,
    external_10x_module_contrasts: pd.DataFrame | None,
    external_10x_marker_composition: pd.DataFrame | None,
    external_10x_marker_contrasts: pd.DataFrame | None,
    external_marker_age_concordance: pd.DataFrame | None,
    external_validation_evidence: pd.DataFrame | None,
    pseudobulk_de: pd.DataFrame | None,
    pseudobulk_coverage: pd.DataFrame | None,
    pseudobulk_metadata: pd.DataFrame | None,
    pseudobulk_covariate_de: pd.DataFrame | None,
    pseudobulk_genomewide_summary: pd.DataFrame | None,
    pseudobulk_genomewide_qc_summary: pd.DataFrame | None,
    pseudobulk_genomewide_gene_qc: pd.DataFrame | None,
    pseudobulk_genomewide_disease_summary: pd.DataFrame | None,
    pseudobulk_genomewide_de_summary: pd.DataFrame | None,
    pseudobulk_genomewide_de_top_hits: pd.DataFrame | None,
    pseudobulk_genomewide_de_audit: pd.DataFrame | None,
    pseudobulk_genomewide_donor_balance: pd.DataFrame | None,
    pseudobulk_genomewide_matched_feasibility: pd.DataFrame | None,
    pseudobulk_genomewide_de_summary_matched: pd.DataFrame | None,
    pseudobulk_genomewide_de_top_hits_matched: pd.DataFrame | None,
    pseudobulk_genomewide_de_audit_matched: pd.DataFrame | None,
    pseudobulk_genomewide_limma_de_summary: pd.DataFrame | None,
    pseudobulk_genomewide_limma_de_top_hits: pd.DataFrame | None,
    pseudobulk_genomewide_limma_de_audit: pd.DataFrame | None,
    pseudobulk_genomewide_limma_de_summary_matched: pd.DataFrame | None,
    pseudobulk_genomewide_limma_de_top_hits_matched: pd.DataFrame | None,
    pseudobulk_genomewide_limma_de_audit_matched: pd.DataFrame | None,
    ora_sensitivity_scenarios: pd.DataFrame | None,
    ora_sensitivity_performance: pd.DataFrame | None,
    ora_repeated_cv_summary: pd.DataFrame | None,
    ora_repeated_cv_feature_stability: pd.DataFrame | None,
    ora_augmented_candidate_repeated_cv_summary: pd.DataFrame | None,
    ora_feature_interpretation: pd.DataFrame | None,
    ora_feature_set_model_comparison: pd.DataFrame | None,
    ora_permutation_empirical: pd.DataFrame | None,
    ora_nested_tuning_summary: pd.DataFrame | None,
    ora_stacking_summary: pd.DataFrame | None,
    ora_model_card: pd.DataFrame | None,
    ndd_label_permutation: pd.DataFrame | None,
    latent_space_readiness: pd.DataFrame | None,
    latent_space_local_audit: pd.DataFrame | None,
    latent_space_portal_assets: pd.DataFrame | None,
    latent_recompute_feasibility: pd.DataFrame | None,
    scvi_pilot_validation: pd.DataFrame | None,
    output_provenance: pd.DataFrame | None,
    out_md: str | Path,
    figure_paths: dict[str, Path],
    source: dict[str, Any],
    paper_defaults: dict[str, Any],
    schema: dict[str, Any],
    top_n: int = 12,
) -> str:
    """Render the markdown body for the MVP report."""

    report_path = Path(out_md)
    top_assoc = rank_associations(associations, top_n=top_n)
    top_pseudobulk = rank_pseudobulk_de(pseudobulk_de, top_n=top_n)
    top_pseudobulk_adjusted = rank_pseudobulk_covariate_de(pseudobulk_covariate_de, top_n=top_n)
    best_model = best_predictive_model(performance)
    combined_perf = combined_performance(performance, augmented_performance)
    healthy_train = _usable_training_donors(manifest)
    lines = [
        "# Gateway ORA MVP Report",
        "",
        "## Source",
        "",
        f"- Dataset: {source.get('name', 'Gateway human olfactory epithelium atlas')}",
        f"- DOI: {source.get('doi', 'not recorded')}",
        f"- CELLxGENE collection: {source.get('collection_id', 'not recorded')}",
        f"- Matrix inspected: {_format_int(schema.get('n_obs'))} cells x {_format_int(schema.get('n_vars'))} genes",
        f"- Paper default cells: {_format_int(paper_defaults.get('cells_total'))}",
        f"- Embeddings available in export: {_format_list(schema.get('obsm_keys', []))}",
        "",
        "## Cohort",
        "",
        f"- Usable healthy training donors: {_format_int(healthy_train)}",
        f"- Donors with any age value: {_format_int(_donors_with_age(manifest))}",
        "",
        _markdown_table(
            cohort_summary,
            ["cohort", "donors", "samples", "cells", "median_age", "missing_age_samples", "lineage_cells", "mature_neurons"],
            max_rows=20,
        ),
        "",
        _figure_link(report_path, figure_paths.get("cohort"), "Cohort overview"),
        "",
        "## ORA Age Prediction",
        "",
        _model_summary_sentence(best_model),
        "",
        _markdown_table(_sort_metric_table(performance, "mae"), ["model", "n", "mae", "rmse", "r2", "spearman_r"], max_rows=20),
        "",
        _figure_link(report_path, figure_paths.get("performance"), "Model performance"),
        "",
        _figure_link(report_path, figure_paths.get("predictions"), "Predicted age versus chronological age"),
        "",
    ]
    if _has_ora_diagnostics(ora_calibration, ora_age_bin_errors, ora_residual_diagnostics):
        lines.extend(
            [
                "## ORA Calibration Diagnostics",
                "",
                _ora_calibration_summary_sentence(ora_calibration),
                "",
                _markdown_table(
                    ora_calibration if ora_calibration is not None else pd.DataFrame(),
                    [
                        "model",
                        "n",
                        "calibration_slope_ora_on_age",
                        "calibration_intercept_ora_on_age",
                        "mae",
                        "calibrated_mae",
                        "spearman_r",
                    ],
                    max_rows=20,
                ),
                "",
                _markdown_table(
                    _diagnostic_model_subset(ora_age_bin_errors, group_col="group", levels=["young", "middle", "old"]),
                    ["model", "group", "level", "n", "mean_error", "mae", "calibrated_mean_error", "calibrated_mae"],
                    max_rows=30,
                ),
                "",
                _markdown_table(
                    _top_residual_diagnostics(ora_residual_diagnostics, top_n=20),
                    ["model", "group", "level", "n", "mean_error", "mae", "mean_oraa", "calibrated_mean_error"],
                    max_rows=20,
                ),
                "",
                _figure_link(report_path, figure_paths.get("calibrated_predictions"), "Calibrated predicted age versus chronological age"),
                "",
                _figure_link(report_path, figure_paths.get("calibration_age_bins"), "Raw and calibrated age-bin error"),
                "",
            ]
        )
    if augmented_performance is not None and not augmented_performance.empty:
        lines.extend(
            [
                "## Module-Augmented ORA",
                "",
                _augmented_summary_sentence(performance, augmented_performance),
                "",
                _markdown_table(
                    _sort_metric_table(combined_perf, "mae"),
                    ["feature_set", "model", "n", "mae", "rmse", "r2", "spearman_r"],
                    max_rows=24,
                ),
                "",
                _figure_link(report_path, figure_paths.get("performance_comparison"), "ORA model comparison"),
                "",
            ]
        )
        if augmented_importance is not None and not augmented_importance.empty:
            lines.extend(
                [
                    "### Augmented Feature Importance",
                    "",
                    _markdown_table(
                        _top_importance(augmented_importance, top_n=10),
                        ["model", "feature", "importance", "stability"],
                        max_rows=20,
                    ),
                    "",
                ]
            )
    if ora_repeated_cv_summary is not None and not ora_repeated_cv_summary.empty:
        lines.extend(
            [
                "## Repeated-CV ORA Stability",
                "",
                _ora_repeated_cv_summary_sentence(ora_repeated_cv_summary),
                "",
                _markdown_table(
                    _sort_metric_table(ora_repeated_cv_summary, "mae_mean"),
                    [
                        "model",
                        "repeats",
                        "n",
                        "mae_mean",
                        "mae_ci_low",
                        "mae_ci_high",
                        "spearman_r_mean",
                        "spearman_r_ci_low",
                        "spearman_r_ci_high",
                    ],
                    max_rows=12,
                ),
                "",
                _markdown_table(
                    _top_repeated_cv_features(ora_repeated_cv_feature_stability, top_n=10),
                    ["model", "feature", "mean_importance", "selection_fraction"],
                    max_rows=10,
                ),
                "",
            ]
        )
        if _has_feature_interpretation(ora_feature_interpretation):
            lines.extend(
                [
                    "### Top Feature Biological Interpretation",
                    "",
                    _feature_interpretation_summary_sentence(ora_feature_interpretation),
                    "",
                    _markdown_table(
                        ora_feature_interpretation if ora_feature_interpretation is not None else pd.DataFrame(),
                        [
                            "feature",
                            "biology_theme",
                            "supporting_models",
                            "age_direction",
                            "interpretation",
                            "caution",
                        ],
                        max_rows=16,
                    ),
                    "",
                ]
            )
        candidate_repeats = _top_candidate_repeated_cv(ora_augmented_candidate_repeated_cv_summary, top_n=10)
        if not candidate_repeats.empty:
            lines.extend(
                [
                    "### Augmented Candidate Model Repeated CV",
                    "",
                    "This targeted benchmark reruns the strongest boosted candidates on the module-augmented feature set without rerunning every baseline model.",
                    "",
                    _markdown_table(
                        candidate_repeats,
                        [
                            "model",
                            "repeats",
                            "n",
                            "mae_mean",
                            "mae_ci_low",
                            "mae_ci_high",
                            "spearman_r_mean",
                            "spearman_r_ci_low",
                            "spearman_r_ci_high",
                        ],
                        max_rows=10,
                    ),
                    "",
                ]
            )
        comparison = _top_feature_set_model_comparison(ora_feature_set_model_comparison, top_n=10)
        if not comparison.empty:
            lines.extend(
                [
                    "### Feature-Set Model Comparison",
                    "",
                    _markdown_table(
                        comparison,
                        [
                            "feature_set",
                            "model",
                            "mae_mean",
                            "mae_ci_low",
                            "mae_ci_high",
                            "rmse_mean",
                            "r2_mean",
                            "spearman_r_mean",
                        ],
                        max_rows=10,
                    ),
                    "",
                ]
            )
        permutation = _top_permutation_empirical(ora_permutation_empirical, top_n=10)
        if not permutation.empty:
            lines.extend(
                [
                    "### Shuffled-Age Null Test",
                    "",
                    _markdown_table(
                        permutation,
                        [
                            "model",
                            "n_permutations",
                            "observed_mae",
                            "null_mae_mean",
                            "empirical_p_mae",
                            "observed_spearman_r",
                            "null_spearman_r_mean",
                            "empirical_p_spearman_r",
                        ],
                        max_rows=10,
                    ),
                    "",
                ]
            )
        nested = _top_nested_tuning_summary(ora_nested_tuning_summary, top_n=10)
        if not nested.empty:
            lines.extend(
                [
                    "### Nested Booster Tuning",
                    "",
                    _markdown_table(
                        nested,
                        [
                            "model",
                            "repeats",
                            "n",
                            "mae_mean",
                            "mae_ci_low",
                            "mae_ci_high",
                            "rmse_mean",
                            "r2_mean",
                            "spearman_r_mean",
                        ],
                        max_rows=10,
                    ),
                    "",
                ]
            )
        stacking = _top_stacking_summary(ora_stacking_summary, top_n=10)
        if not stacking.empty:
            lines.extend(
                [
                    "### Leakage-Safe OOF Stacking",
                    "",
                    _markdown_table(
                        stacking,
                        [
                            "model",
                            "repeats",
                            "n",
                            "mae_mean",
                            "mae_ci_low",
                            "mae_ci_high",
                            "rmse_mean",
                            "r2_mean",
                            "spearman_r_mean",
                        ],
                        max_rows=10,
                    ),
                    "",
                ]
            )
    if ora_model_card is not None and not ora_model_card.empty:
        lines.extend(
            [
                "## ORA Model Card",
                "",
                _markdown_table(
                    ora_model_card,
                    [
                        "model",
                        "feature_set",
                        "role",
                        "n",
                        "repeats",
                        "mae_mean",
                        "spearman_r_mean",
                        "calibration_slope",
                        "permutation_p_mae",
                        "limitations",
                    ],
                    max_rows=16,
                ),
                "",
            ]
        )
    if _has_ndd_projection(ndd_projection, ndd_projection_summary):
        lines.extend(
            [
                "## NDD ORA Projection",
                "",
                _ndd_projection_summary_sentence(ndd_projection),
                "",
                _markdown_table(
                    ndd_projection_summary if ndd_projection_summary is not None else pd.DataFrame(),
                    ["model", "disease_group", "donors", "training_donors", "ndd_donors", "mean_age", "mean_ora", "mean_oraa", "sd_oraa"],
                    max_rows=30,
                ),
                "",
                _markdown_table(
                    _top_ndd_projection_rows(ndd_projection, top_n=10),
                    ["donor_id", "model", "disease_group", "chronological_age", "ora", "oraa"],
                    max_rows=10,
                ),
                "",
                _markdown_table(
                    ndd_projection_uncertainty if ndd_projection_uncertainty is not None else pd.DataFrame(),
                    [
                        "model",
                        "disease_group",
                        "reference",
                        "n_disease",
                        "n_reference",
                        "mean_oraa",
                        "mean_oraa_ci_low",
                        "mean_oraa_ci_high",
                        "difference_vs_reference",
                        "difference_ci_low",
                        "difference_ci_high",
                    ],
                    max_rows=20,
                ),
                "",
                _markdown_table(
                    ndd_projection_context if ndd_projection_context is not None else pd.DataFrame(),
                    ["disease_group", "chemistry", "collection_method", "donors", "mean_age", "median_total_cells"],
                    max_rows=20,
                ),
                "",
                _figure_link(report_path, figure_paths.get("ndd_projection"), "NDD ORA projection"),
                "",
                _figure_link(report_path, figure_paths.get("ndd_matched_reference"), "Matched NDD ORA reference"),
                "",
            ]
        )
        ndd_feature_comparison = _top_ndd_feature_comparison(ndd_projection_feature_comparison, top_n=24)
        if not ndd_feature_comparison.empty:
            lines.extend(
                [
                    "### Feature-Set Projection Sensitivity",
                    "",
                    _markdown_table(
                        ndd_feature_comparison,
                        [
                            "model",
                            "disease_group",
                            "composition_mean_oraa",
                            "augmented_mean_oraa",
                            "augmented_minus_composition_oraa",
                            "sign_stable_negative",
                        ],
                        max_rows=24,
                    ),
                    "",
                ]
            )
        ndd_donor_appendix = _top_ndd_donor_appendix(ndd_projection_donor_appendix, top_n=12)
        if not ndd_donor_appendix.empty:
            lines.extend(
                [
                    "### Donor Appendix Preview",
                    "",
                    _markdown_table(
                        ndd_donor_appendix,
                        [
                            "feature_set",
                            "donor_id",
                            "disease_group",
                            "chronological_age",
                            "chemistry",
                            "collection_method",
                            "model",
                            "ora",
                            "oraa",
                        ],
                        max_rows=12,
                    ),
                    "",
                ]
            )
        ndd_diagnostics = _top_ndd_projection_diagnostics(ndd_projection_diagnostics, top_n=30)
        if not ndd_diagnostics.empty:
            lines.extend(
                [
                    "### Projection Diagnostics",
                    "",
                    _markdown_table(
                        ndd_diagnostics,
                        [
                            "model",
                            "disease_group",
                            "diagnostic",
                            "level",
                            "n_donors",
                            "mean_age",
                            "median_total_cells",
                            "mean_oraa",
                            "status",
                        ],
                        max_rows=30,
                    ),
                    "",
                ]
            )
        if ndd_label_permutation is not None and not ndd_label_permutation.empty:
            ndd_label_display = ndd_label_permutation
            if "model" in ndd_label_display:
                ndd_label_display = ndd_label_display[~ndd_label_display["model"].eq("null_model")]
            lines.extend(
                [
                    "### NDD Label Permutation",
                    "",
                    "Frozen-score label permutations test whether AD/PD ORAA differences are unusually negative within matched chemistry and collection-method context.",
                    "",
                    _markdown_table(
                        ndd_label_display,
                        [
                            "model",
                            "disease_group",
                            "n_disease",
                            "n_reference",
                            "observed_difference_vs_reference",
                            "null_mean",
                            "empirical_p_negative",
                            "status",
                        ],
                        max_rows=24,
                    ),
                    "",
                ]
            )
    if _has_module_tables(module_summary, module_coverage, donor_module_features):
        lines.extend(
            [
                "## Module Scores",
                "",
                _markdown_table(
                    module_coverage if module_coverage is not None else pd.DataFrame(),
                    ["module", "n_requested", "n_present", "coverage_fraction", "missing_genes"],
                    max_rows=20,
                ),
                "",
                _markdown_table(
                    _module_feature_summary(donor_module_features),
                    ["module", "mean", "sd", "min", "max"],
                    max_rows=20,
                ),
                "",
                _figure_link(report_path, figure_paths.get("module_scores"), "Module scores by cell state"),
                "",
            ]
        )
    if _has_external_validation_tables(external_validation_summary, external_gene_list_coverage, external_feature_contract):
        lines.extend(
            [
                "## External Validation Readiness",
                "",
                _external_validation_summary_sentence(external_validation_summary, external_gene_list_coverage),
                "",
                _markdown_table(
                    external_validation_summary if external_validation_summary is not None else pd.DataFrame(),
                    [
                        "dataset_id",
                        "status",
                        "validation_use",
                        "expected_level",
                        "readiness_class",
                        "ready_for_feature_validation",
                        "ready_for_raw_adapter",
                        "files_missing",
                    ],
                    max_rows=20,
                ),
                "",
                _markdown_table(
                    external_gene_list_coverage if external_gene_list_coverage is not None else pd.DataFrame(),
                    ["gene_list", "n_requested", "n_present", "coverage_fraction", "missing_genes"],
                    max_rows=20,
                ),
                "",
                _markdown_table(
                    external_feature_contract if external_feature_contract is not None else pd.DataFrame(),
                    ["field", "kind"],
                    max_rows=20,
                ),
                "",
            ]
        )
        if _has_external_validation_evidence(external_validation_evidence):
            lines.extend(
                [
                    "### External Evidence Ledger",
                    "",
                    _external_evidence_summary_sentence(external_validation_evidence),
                    "",
                    _markdown_table(
                        external_validation_evidence if external_validation_evidence is not None else pd.DataFrame(),
                        [
                            "dataset_id",
                            "evidence_type",
                            "feature_level",
                            "validation_strength",
                            "n_samples",
                            "n_donors",
                            "n_features",
                            "supports_primary_claim",
                            "supports_ndd_claim",
                            "next_action",
                        ],
                        max_rows=20,
                    ),
                    "",
                ]
            )
        if _has_external_10x_tables(external_sample_metadata, external_10x_sample_qc, external_10x_module_contrasts):
            lines.extend(
                [
                    "### GSE184117 Sample-Level Module Sanity Check",
                    "",
                    _external_10x_summary_sentence(
                        external_sample_metadata,
                        external_10x_sample_qc,
                        external_10x_module_contrasts,
                    ),
                    "",
                    _markdown_table(
                        external_sample_metadata if external_sample_metadata is not None else pd.DataFrame(),
                        [
                            "sample_id",
                            "donor_id",
                            "age",
                            "disease_state",
                            "disease_group",
                            "sample_class",
                            "usable_for_external_validation",
                        ],
                        max_rows=12,
                    ),
                    "",
                    _markdown_table(
                        external_10x_sample_qc if external_10x_sample_qc is not None else pd.DataFrame(),
                        ["sample_id", "disease_group", "sample_class", "n_cells", "n_genes", "detected_genes", "total_counts"],
                        max_rows=12,
                    ),
                    "",
                    _markdown_table(
                        _top_external_module_contrasts(external_10x_module_contrasts),
                        [
                            "module",
                            "n_healthy",
                            "n_presbyosmia",
                            "mean_healthy",
                            "mean_presbyosmia",
                            "presbyosmia_minus_healthy",
                            "p_value",
                            "direction",
                            "status",
                        ],
                        max_rows=16,
                    ),
                    "",
                ]
            )
        if _has_external_marker_tables(external_10x_marker_composition, external_10x_marker_contrasts):
            lines.extend(
                [
                    "### GSE184117 Marker-Only Composition Sanity Check",
                    "",
                    _external_marker_summary_sentence(external_10x_marker_composition, external_10x_marker_contrasts),
                    "",
                    _markdown_table(
                        _top_external_marker_contrasts(external_10x_marker_contrasts),
                        [
                            "marker_panel",
                            "n_healthy",
                            "n_presbyosmia",
                            "mean_fraction_healthy",
                            "mean_fraction_presbyosmia",
                            "presbyosmia_minus_healthy",
                            "p_value",
                            "direction",
                            "status",
                        ],
                        max_rows=16,
                    ),
                    "",
                ]
            )
            if _has_external_marker_age_concordance(external_marker_age_concordance):
                lines.extend(
                    [
                        "### GSE184117 Marker-Age Concordance",
                        "",
                        _external_marker_age_concordance_sentence(external_marker_age_concordance),
                        "",
                        _markdown_table(
                            _top_external_marker_age_concordance(external_marker_age_concordance),
                            [
                                "marker_panel",
                                "gateway_feature",
                                "external_direction",
                                "gateway_age_direction",
                                "concordance",
                                "external_delta",
                                "gateway_beta_per_10_years",
                                "gateway_fdr",
                                "status",
                            ],
                            max_rows=20,
                        ),
                        "",
                    ]
                )
    if _has_latent_space_audit(
        latent_space_readiness,
        latent_space_local_audit,
        latent_space_portal_assets,
    ) or _has_scvi_pilot_validation(scvi_pilot_validation):
        lines.extend(
            [
                "## Latent-Space Readiness",
                "",
                _latent_space_summary_sentence(latent_space_readiness),
                "",
                _markdown_table(
                    latent_space_readiness if latent_space_readiness is not None else pd.DataFrame(),
                    [
                        "status",
                        "local_embeddings",
                        "portal_embeddings",
                        "usable_local_embeddings",
                        "visualization_only_embeddings",
                        "recommendation",
                    ],
                    max_rows=5,
                ),
                "",
                _markdown_table(
                    latent_space_local_audit if latent_space_local_audit is not None else pd.DataFrame(),
                    [
                        "embedding_key",
                        "n_cells",
                        "n_dimensions",
                        "readiness",
                        "recommended_use",
                        "notes",
                    ],
                    max_rows=20,
                ),
                "",
                _markdown_table(
                    latent_space_portal_assets if latent_space_portal_assets is not None else pd.DataFrame(),
                    [
                        "dataset_id",
                        "dataset_version_id",
                        "asset_filetype",
                        "asset_filesize_bytes",
                        "portal_embeddings",
                        "default_embedding",
                        "status",
                    ],
                    max_rows=10,
                ),
                "",
                _markdown_table(
                    latent_recompute_feasibility if latent_recompute_feasibility is not None else pd.DataFrame(),
                    ["check", "status", "detail", "recommendation"],
                    max_rows=20,
                ),
                "",
                _scvi_pilot_summary_sentence(scvi_pilot_validation),
                "",
                _markdown_table(
                    scvi_pilot_validation if scvi_pilot_validation is not None else pd.DataFrame(),
                    ["check", "status", "detail", "recommendation"],
                    max_rows=20,
                ),
                "",
            ]
        )
    if _has_pseudobulk_tables(pseudobulk_de, pseudobulk_coverage, pseudobulk_metadata):
        lines.extend(
            [
                "## Pseudobulk Differential Expression",
                "",
                _pseudobulk_summary_sentence(pseudobulk_de, pseudobulk_metadata),
                "",
                _markdown_table(
                    pseudobulk_coverage if pseudobulk_coverage is not None else pd.DataFrame(),
                    ["module", "n_requested", "n_present", "coverage_fraction", "missing_genes"],
                    max_rows=20,
                ),
                "",
                _markdown_table(
                    _pseudobulk_metadata_summary(pseudobulk_metadata),
                    ["disease_group", "groups", "cells"],
                    max_rows=10,
                ),
                "",
                _markdown_table(
                    top_pseudobulk,
                    ["contrast", "fine_cell_type", "gene", "n_case", "n_control", "log2fc", "p_value", "fdr"],
                    max_rows=top_n,
                ),
                "",
                _figure_link(report_path, figure_paths.get("pseudobulk_de"), "Top pseudobulk DE hits"),
                "",
            ]
        )
    if pseudobulk_covariate_de is not None and not pseudobulk_covariate_de.empty:
        lines.extend(
            [
                "## Covariate-Adjusted Pseudobulk DE",
                "",
                _pseudobulk_adjusted_summary_sentence(pseudobulk_covariate_de),
                "",
                _markdown_table(
                    top_pseudobulk_adjusted,
                    [
                        "contrast",
                        "fine_cell_type",
                        "gene",
                        "n_case",
                        "n_control",
                        "log2fc_adjusted",
                        "p_value",
                        "fdr",
                        "covariates",
                    ],
                    max_rows=top_n,
                ),
                "",
                _figure_link(report_path, figure_paths.get("pseudobulk_covariate_de"), "Top covariate-adjusted pseudobulk DE hits"),
                "",
            ]
        )
    if pseudobulk_genomewide_summary is not None and not pseudobulk_genomewide_summary.empty:
        lines.extend(
            [
                "## Genome-Wide Pseudobulk Export",
                "",
                _pseudobulk_genomewide_summary_sentence(pseudobulk_genomewide_summary),
                "",
                _pseudobulk_genomewide_qc_sentence(pseudobulk_genomewide_qc_summary),
                "",
                _markdown_table(
                    pseudobulk_genomewide_summary,
                    [
                        "n_cells",
                        "n_genes",
                        "n_groups_total",
                        "n_groups_exported",
                        "n_groups_failed_min_cells",
                        "n_groups_failed_min_donors",
                        "min_cells_per_group",
                        "min_donors_per_cell_state",
                    ],
                    max_rows=5,
                ),
                "",
                _markdown_table(
                    pseudobulk_genomewide_disease_summary if pseudobulk_genomewide_disease_summary is not None else pd.DataFrame(),
                    ["disease_group", "groups", "donors", "cells", "matrix_total_count", "median_detected_genes"],
                    max_rows=10,
                ),
                "",
                _markdown_table(
                    _top_genomewide_variable_genes(pseudobulk_genomewide_gene_qc, top_n=10),
                    ["gene_symbol", "total_count", "detected_group_fraction", "variance_log1p"],
                    max_rows=10,
                ),
                "",
            ]
        )
    if pseudobulk_genomewide_de_summary is not None and not pseudobulk_genomewide_de_summary.empty:
        lines.extend(
            [
                "## Genome-Wide edgeR DE",
                "",
                _pseudobulk_genomewide_de_sentence(pseudobulk_genomewide_de_summary),
                "",
                _markdown_table(
                    pseudobulk_genomewide_de_summary,
                    [
                        "contrast",
                        "tested_rows",
                        "tested_genes",
                        "tested_cell_states",
                        "ok_cell_state_models",
                        "significant_rows",
                        "significant_genes",
                        "significant_cell_states",
                        "sex_linked_significant_rows",
                    ],
                    max_rows=10,
                ),
                "",
                _markdown_table(
                    pseudobulk_genomewide_de_top_hits if pseudobulk_genomewide_de_top_hits is not None else pd.DataFrame(),
                    ["contrast", "fine_cell_type", "gene_symbol", "log2fc", "p_value", "fdr", "is_sex_linked_initial"],
                    max_rows=top_n,
                ),
                "",
                "Top non-sex-linked sentinel hits:",
                "",
                _markdown_table(
                    _non_sex_linked_genomewide_de_hits(pseudobulk_genomewide_de_top_hits, top_n=top_n),
                    ["contrast", "fine_cell_type", "gene_symbol", "log2fc", "p_value", "fdr", "is_sex_linked_initial"],
                    max_rows=top_n,
                ),
                "",
            ]
        )
        if pseudobulk_genomewide_de_audit is not None and not pseudobulk_genomewide_de_audit.empty:
            lines.extend(
                [
                    "### Genome-Wide DE Audit",
                    "",
                    _markdown_table(
                        pseudobulk_genomewide_de_audit,
                        [
                            "contrast",
                            "significant_rows",
                            "is_sex_linked_initial_significant_rows",
                            "is_mitochondrial_significant_rows",
                            "is_ribosomal_significant_rows",
                            "is_hemoglobin_significant_rows",
                            "is_immunoglobulin_significant_rows",
                        ],
                        max_rows=10,
                    ),
                    "",
                ]
            )
        donor_balance = _top_donor_balance_issues(pseudobulk_genomewide_donor_balance, top_n=20)
        if not donor_balance.empty:
            lines.extend(
                [
                    "### Genome-Wide Donor Balance",
                    "",
                    _markdown_table(
                        donor_balance,
                        ["contrast", "fine_cell_type", "n_case", "n_control", "status", "balance_status"],
                        max_rows=20,
                    ),
                    "",
                ]
            )
        matched_feasibility = _top_matched_de_feasibility(pseudobulk_genomewide_matched_feasibility, top_n=20)
        if not matched_feasibility.empty:
            lines.extend(
                [
                    "### Matched FLEX v2/Device DE Feasibility",
                    "",
                    _markdown_table(
                        matched_feasibility,
                        [
                            "contrast",
                            "fine_cell_type",
                            "n_case",
                            "n_matched_healthy",
                            "ready_for_matched_de",
                        ],
                        max_rows=20,
                    ),
                    "",
                ]
            )
        if pseudobulk_genomewide_de_summary_matched is not None and not pseudobulk_genomewide_de_summary_matched.empty:
            lines.extend(
                [
                    "### Matched FLEX v2/Device edgeR DE Sensitivity",
                    "",
                    _matched_genomewide_de_sentence(pseudobulk_genomewide_de_summary_matched),
                    "",
                    _markdown_table(
                        pseudobulk_genomewide_de_summary_matched,
                        [
                            "contrast",
                            "tested_rows",
                            "tested_genes",
                            "tested_cell_states",
                            "ok_cell_state_models",
                            "significant_rows",
                            "significant_genes",
                            "significant_cell_states",
                            "sex_linked_significant_rows",
                        ],
                        max_rows=10,
                    ),
                    "",
                    _markdown_table(
                        _non_sex_linked_genomewide_de_hits(pseudobulk_genomewide_de_top_hits_matched, top_n=top_n),
                        ["contrast", "fine_cell_type", "gene_symbol", "log2fc", "p_value", "fdr", "is_sex_linked_initial"],
                        max_rows=top_n,
                    ),
                    "",
                ]
            )
            if pseudobulk_genomewide_de_audit_matched is not None and not pseudobulk_genomewide_de_audit_matched.empty:
                lines.extend(
                    [
                        "Matched sentinel audit:",
                        "",
                        _markdown_table(
                            pseudobulk_genomewide_de_audit_matched,
                            [
                                "contrast",
                                "significant_rows",
                                "is_sex_linked_initial_significant_rows",
                                "is_mitochondrial_significant_rows",
                                "is_ribosomal_significant_rows",
                                "is_hemoglobin_significant_rows",
                                "is_immunoglobulin_significant_rows",
                            ],
                            max_rows=10,
                        ),
                        "",
                    ]
                )
    if pseudobulk_genomewide_limma_de_summary is not None and not pseudobulk_genomewide_limma_de_summary.empty:
        lines.extend(
            [
                "## Genome-Wide limma-voom DE Parity",
                "",
                _genomewide_de_method_sentence(pseudobulk_genomewide_limma_de_summary, "limma-voom"),
                "",
                _markdown_table(
                    pseudobulk_genomewide_limma_de_summary,
                    [
                        "contrast",
                        "tested_rows",
                        "tested_genes",
                        "tested_cell_states",
                        "ok_cell_state_models",
                        "significant_rows",
                        "significant_genes",
                        "significant_cell_states",
                        "sex_linked_significant_rows",
                    ],
                    max_rows=10,
                ),
                "",
                "Top non-sex-linked limma-voom hits:",
                "",
                _markdown_table(
                    _non_sex_linked_genomewide_de_hits(pseudobulk_genomewide_limma_de_top_hits, top_n=top_n),
                    ["contrast", "fine_cell_type", "gene_symbol", "log2fc", "p_value", "fdr", "is_sex_linked_initial"],
                    max_rows=top_n,
                ),
                "",
            ]
        )
        if pseudobulk_genomewide_limma_de_audit is not None and not pseudobulk_genomewide_limma_de_audit.empty:
            lines.extend(
                [
                    "### limma-voom Sentinel Audit",
                    "",
                    _markdown_table(
                        pseudobulk_genomewide_limma_de_audit,
                        [
                            "contrast",
                            "significant_rows",
                            "is_sex_linked_initial_significant_rows",
                            "is_mitochondrial_significant_rows",
                            "is_ribosomal_significant_rows",
                            "is_hemoglobin_significant_rows",
                            "is_immunoglobulin_significant_rows",
                        ],
                        max_rows=10,
                    ),
                    "",
                ]
            )
        if pseudobulk_genomewide_limma_de_summary_matched is not None and not pseudobulk_genomewide_limma_de_summary_matched.empty:
            lines.extend(
                [
                    "### Matched FLEX v2/Device limma-voom DE Sensitivity",
                    "",
                    _genomewide_de_method_sentence(
                        pseudobulk_genomewide_limma_de_summary_matched,
                        "Matched FLEX v2/device limma-voom",
                    ),
                    "",
                    _markdown_table(
                        pseudobulk_genomewide_limma_de_summary_matched,
                        [
                            "contrast",
                            "tested_rows",
                            "tested_genes",
                            "tested_cell_states",
                            "ok_cell_state_models",
                            "significant_rows",
                            "significant_genes",
                            "significant_cell_states",
                            "sex_linked_significant_rows",
                        ],
                        max_rows=10,
                    ),
                    "",
                    _markdown_table(
                        _non_sex_linked_genomewide_de_hits(pseudobulk_genomewide_limma_de_top_hits_matched, top_n=top_n),
                        ["contrast", "fine_cell_type", "gene_symbol", "log2fc", "p_value", "fdr", "is_sex_linked_initial"],
                        max_rows=top_n,
                    ),
                    "",
                ]
            )
            if (
                pseudobulk_genomewide_limma_de_audit_matched is not None
                and not pseudobulk_genomewide_limma_de_audit_matched.empty
            ):
                lines.extend(
                    [
                        "Matched limma-voom sentinel audit:",
                        "",
                        _markdown_table(
                            pseudobulk_genomewide_limma_de_audit_matched,
                            [
                                "contrast",
                                "significant_rows",
                                "is_sex_linked_initial_significant_rows",
                                "is_mitochondrial_significant_rows",
                                "is_ribosomal_significant_rows",
                                "is_hemoglobin_significant_rows",
                                "is_immunoglobulin_significant_rows",
                            ],
                            max_rows=10,
                        ),
                        "",
                    ]
                )
    if ora_sensitivity_performance is not None and not ora_sensitivity_performance.empty:
        lines.extend(
            [
                "## ORA Sensitivity",
                "",
                _ora_sensitivity_summary_sentence(ora_sensitivity_scenarios, ora_sensitivity_performance),
                "",
                _markdown_table(
                    _top_ora_sensitivity_performance(ora_sensitivity_performance, model="random_forest"),
                    ["scenario", "model", "n", "mae", "rmse", "r2", "spearman_r", "healthy_train_donors"],
                    max_rows=20,
                ),
                "",
            ]
        )
    lines.extend(
        [
        "## Age Associations",
        "",
        f"- Association tests with status ok: {_format_int(int(associations['status'].eq('ok').sum())) if 'status' in associations else 'not recorded'}",
        "- Top table is ranked by FDR, then p-value; beta is per 10 years.",
        "",
        _markdown_table(
            top_assoc,
            ["feature", "n", "beta_per_10_years", "p_value", "fdr", "direction"],
            max_rows=top_n,
        ),
        "",
        _figure_link(report_path, figure_paths.get("associations"), "Top age associations"),
        "",
        "## Feature Importance",
        "",
        _markdown_table(_top_importance(importance, top_n=10), ["model", "feature", "importance", "stability"], max_rows=20),
        "",
        _figure_link(report_path, figure_paths.get("importance"), "Top model features"),
        "",
        "## Reproducibility Snapshot",
        "",
        _output_provenance_sentence(output_provenance),
        "",
        "## Interpretation Notes",
        "",
        "- The composition baseline and module-augmented ORA models are trained only on healthy donors with valid age.",
        "- NDD ORA projections use frozen healthy-trained models; projected AD/PD donors are not included in training or cross-validation.",
        "- Module scores are average log1p expression over curated marker sets, summarized at donor and cell-state levels.",
        "- Pseudobulk DE includes both unadjusted donor-level logCPM Welch contrasts and targeted covariate-adjusted linear models.",
        "- ORA predictions are under-dispersed across chronological age; calibration diagnostics support using ORA as a relative tissue-state axis rather than an absolute biological-age estimator.",
        "- Genome-wide pseudobulk counts now have local edgeR quasi-likelihood and limma-voom workflows; DESeq2 remains an adapter hook.",
        "- Genome-wide NDD DE is discovery-oriented only: AD/PD sample sizes are five donors each and sex/chemistry/collection imbalance can dominate top hits.",
        "- Trajectory density, Milo, and cNMF remain deferred commands.",
        "- Chemistry, collection method, site, and yield are treated as covariates or sensitivity variables rather than biological ORA features.",
        "- AD/PD donors are excluded from ORA training and reserved for later frozen-model projection.",
        "",
        ]
    )
    return "\n".join(line for line in lines if line is not None)


def _write_figures(
    *,
    cohort_summary: pd.DataFrame,
    associations: pd.DataFrame,
    performance: pd.DataFrame,
    scores: pd.DataFrame,
    importance: pd.DataFrame,
    ora_calibrated_scores: pd.DataFrame | None,
    ora_age_bin_errors: pd.DataFrame | None,
    ora_calibration: pd.DataFrame | None,
    augmented_performance: pd.DataFrame | None,
    ndd_projection: pd.DataFrame | None,
    module_summary: pd.DataFrame | None,
    pseudobulk_de: pd.DataFrame | None,
    pseudobulk_covariate_de: pd.DataFrame | None,
    figure_dir: Path,
    top_n: int,
) -> dict[str, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = {key: figure_dir / name for key, name in FIGURE_NAMES.items()}
    _plot_cohort(cohort_summary, paths["cohort"], plt)
    _plot_performance(performance, paths["performance"], plt)
    if augmented_performance is not None and not augmented_performance.empty:
        _plot_performance_comparison(
            combined_performance(performance, augmented_performance),
            paths["performance_comparison"],
            plt,
        )
    else:
        paths.pop("performance_comparison", None)
    _plot_associations(rank_associations(associations, top_n=top_n), paths["associations"], plt)
    _plot_predictions(scores, paths["predictions"], plt)
    if ora_calibrated_scores is not None and not ora_calibrated_scores.empty:
        _plot_calibrated_predictions(ora_calibrated_scores, ora_calibration, paths["calibrated_predictions"], plt)
    else:
        paths.pop("calibrated_predictions", None)
    if ora_age_bin_errors is not None and not ora_age_bin_errors.empty:
        _plot_calibration_age_bins(ora_age_bin_errors, ora_calibration, paths["calibration_age_bins"], plt)
    else:
        paths.pop("calibration_age_bins", None)
    _plot_importance(importance, paths["importance"], plt)
    if ndd_projection is not None and not ndd_projection.empty:
        _plot_ndd_projection(ndd_projection, paths["ndd_projection"], plt)
        _plot_ndd_matched_reference(ndd_projection, paths["ndd_matched_reference"], plt)
    else:
        paths.pop("ndd_projection", None)
        paths.pop("ndd_matched_reference", None)
    if module_summary is not None and not module_summary.empty:
        _plot_module_scores(module_summary, paths["module_scores"], plt)
    else:
        paths.pop("module_scores", None)
    if pseudobulk_de is not None and not pseudobulk_de.empty:
        _plot_pseudobulk_de(pseudobulk_de, paths["pseudobulk_de"], plt, top_n=top_n)
    else:
        paths.pop("pseudobulk_de", None)
    if pseudobulk_covariate_de is not None and not pseudobulk_covariate_de.empty:
        _plot_pseudobulk_covariate_de(pseudobulk_covariate_de, paths["pseudobulk_covariate_de"], plt, top_n=top_n)
    else:
        paths.pop("pseudobulk_covariate_de", None)
    return paths


def _plot_cohort(cohort_summary: pd.DataFrame, path: Path, plt: Any) -> None:
    frame = cohort_summary.copy()
    if frame.empty or "cohort" not in frame:
        _blank_figure(path, plt, "No cohort summary available")
        return
    cohorts = frame["cohort"].astype(str).tolist()
    donors = pd.to_numeric(frame.get("donors"), errors="coerce").fillna(0)
    cells = pd.to_numeric(frame.get("cells"), errors="coerce").fillna(0) / 1_000_000
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    axes[0].bar(cohorts, donors, color="#287c8e")
    axes[0].set_title("Donors")
    axes[0].set_ylabel("Count")
    axes[1].bar(cohorts, cells, color="#c2674f")
    axes[1].set_title("Cells")
    axes[1].set_ylabel("Millions")
    for ax in axes:
        ax.tick_params(axis="x", rotation=25)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Gateway Cohort Overview")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_performance(performance: pd.DataFrame, path: Path, plt: Any) -> None:
    if performance.empty or "model" not in performance:
        _blank_figure(path, plt, "No model performance available")
        return
    frame = performance.copy()
    frame["mae"] = pd.to_numeric(frame.get("mae"), errors="coerce")
    frame["spearman_r"] = pd.to_numeric(frame.get("spearman_r"), errors="coerce")
    labels = frame["model"].astype(str).str.replace("_", " ").tolist()
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    axes[0].bar(labels, frame["mae"], color="#4d7fb8")
    axes[0].set_title("MAE")
    axes[0].set_ylabel("Years")
    axes[1].bar(labels, frame["spearman_r"], color="#7a9f45")
    axes[1].set_title("Spearman r")
    axes[1].axhline(0, color="#555555", linewidth=0.8)
    for ax in axes:
        ax.tick_params(axis="x", rotation=25)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("ORA Model Performance")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_performance_comparison(performance: pd.DataFrame, path: Path, plt: Any) -> None:
    if performance.empty or not {"feature_set", "model", "mae", "spearman_r"}.issubset(performance.columns):
        _blank_figure(path, plt, "No augmented model performance available")
        return
    frame = performance[~performance["model"].eq("null_model")].copy()
    frame["mae"] = pd.to_numeric(frame["mae"], errors="coerce")
    frame["spearman_r"] = pd.to_numeric(frame["spearman_r"], errors="coerce")
    frame = frame.sort_values(["model", "feature_set"])
    labels = [
        f"{str(row.model).replace('_', ' ')}\n{str(row.feature_set).replace('_', ' ')}"
        for row in frame.itertuples()
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    axes[0].bar(labels, frame["mae"], color="#4d7fb8")
    axes[0].set_title("MAE")
    axes[0].set_ylabel("Years")
    axes[1].bar(labels, frame["spearman_r"], color="#7a9f45")
    axes[1].set_title("Spearman r")
    axes[1].axhline(0, color="#555555", linewidth=0.8)
    for ax in axes:
        ax.tick_params(axis="x", rotation=25, labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Composition Versus Module-Augmented ORA")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_associations(top_assoc: pd.DataFrame, path: Path, plt: Any) -> None:
    if top_assoc.empty:
        _blank_figure(path, plt, "No age associations available")
        return
    frame = top_assoc.sort_values("beta_per_10_years")
    labels = [str(item).replace("__", ": ") for item in frame["feature"]]
    beta = pd.to_numeric(frame["beta_per_10_years"], errors="coerce")
    colors = np.where(beta >= 0, "#287c8e", "#c2674f")
    fig, ax = plt.subplots(figsize=(8, max(4, 0.34 * len(frame))), constrained_layout=True)
    ax.barh(labels, beta, color=colors)
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_xlabel("Beta per 10 years")
    ax.set_title("Top Age Associations")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_predictions(scores: pd.DataFrame, path: Path, plt: Any) -> None:
    if scores.empty or not {"model", "chronological_age", "ora"}.issubset(scores.columns):
        _blank_figure(path, plt, "No ORA predictions available")
        return
    frame = scores[scores["model"].isin(_display_models())].copy()
    if frame.empty:
        frame = scores.copy()
    models = list(frame["model"].drop_duplicates())
    fig, axes = plt.subplots(1, len(models), figsize=(4.4 * len(models), 4), squeeze=False, constrained_layout=True)
    for ax, model in zip(axes.ravel(), models):
        sub = frame[frame["model"].eq(model)]
        x = pd.to_numeric(sub["chronological_age"], errors="coerce")
        y = pd.to_numeric(sub["ora"], errors="coerce")
        ax.scatter(x, y, s=20, alpha=0.72, color="#4d7fb8", edgecolor="none")
        low = float(np.nanmin([x.min(), y.min()]))
        high = float(np.nanmax([x.max(), y.max()]))
        ax.plot([low, high], [low, high], color="#555555", linewidth=0.9, linestyle="--")
        ax.set_title(str(model).replace("_", " "))
        ax.set_xlabel("Chronological age")
        ax.set_ylabel("Predicted age")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Donor-Level ORA Predictions")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_calibrated_predictions(
    calibrated_scores: pd.DataFrame,
    calibration: pd.DataFrame | None,
    path: Path,
    plt: Any,
) -> None:
    required = {"model", "chronological_age", "ora", "calibrated_ora"}
    if calibrated_scores.empty or not required.issubset(calibrated_scores.columns):
        _blank_figure(path, plt, "No calibrated ORA predictions available")
        return
    models = _calibration_display_models(calibration, calibrated_scores, max_models=4)
    frame = calibrated_scores[calibrated_scores["model"].astype(str).isin(models)].copy()
    if frame.empty:
        _blank_figure(path, plt, "No calibrated display models available")
        return
    fig, axes = plt.subplots(1, len(models), figsize=(4.4 * len(models), 4), squeeze=False, constrained_layout=True)
    for ax, model in zip(axes.ravel(), models):
        sub = frame[frame["model"].astype(str).eq(model)].copy()
        x = pd.to_numeric(sub["chronological_age"], errors="coerce")
        raw = pd.to_numeric(sub["ora"], errors="coerce")
        calibrated = pd.to_numeric(sub["calibrated_ora"], errors="coerce")
        valid = x.notna() & raw.notna() & calibrated.notna()
        if not valid.any():
            continue
        ax.scatter(x[valid], raw[valid], s=16, alpha=0.34, color="#8a8f98", edgecolor="none", label="Raw")
        ax.scatter(x[valid], calibrated[valid], s=20, alpha=0.72, color="#287c8e", edgecolor="none", label="Calibrated")
        low = float(np.nanmin([x[valid].min(), raw[valid].min(), calibrated[valid].min()]))
        high = float(np.nanmax([x[valid].max(), raw[valid].max(), calibrated[valid].max()]))
        ax.plot([low, high], [low, high], color="#555555", linewidth=0.9, linestyle="--")
        ax.set_title(str(model).replace("_", " "))
        ax.set_xlabel("Chronological age")
        ax.set_ylabel("Predicted age")
        ax.spines[["top", "right"]].set_visible(False)
    axes.ravel()[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("Raw Versus Recalibrated ORA Predictions")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_calibration_age_bins(
    age_bin_errors: pd.DataFrame,
    calibration: pd.DataFrame | None,
    path: Path,
    plt: Any,
) -> None:
    required = {"model", "group", "level", "mae", "calibrated_mae"}
    if age_bin_errors.empty or not required.issubset(age_bin_errors.columns):
        _blank_figure(path, plt, "No age-bin calibration errors available")
        return
    models = _calibration_display_models(calibration, age_bin_errors, max_models=4)
    frame = age_bin_errors[
        age_bin_errors["model"].astype(str).isin(models)
        & age_bin_errors["group"].astype(str).eq("age_bin")
    ].copy()
    if frame.empty:
        _blank_figure(path, plt, "No age-bin calibration display models available")
        return
    frame["mae"] = pd.to_numeric(frame["mae"], errors="coerce")
    frame["calibrated_mae"] = pd.to_numeric(frame["calibrated_mae"], errors="coerce")
    age_order = [label for label in ["young", "middle", "old"] if label in set(frame["level"].astype(str))]
    if not age_order:
        age_order = sorted(frame["level"].astype(str).unique())
    fig, axes = plt.subplots(1, len(models), figsize=(4.2 * len(models), 3.8), squeeze=False, constrained_layout=True)
    width = 0.34
    for ax, model in zip(axes.ravel(), models):
        sub = frame[frame["model"].astype(str).eq(model)].set_index("level")
        raw = [float(sub.loc[level, "mae"]) if level in sub.index else np.nan for level in age_order]
        calibrated = [float(sub.loc[level, "calibrated_mae"]) if level in sub.index else np.nan for level in age_order]
        x = np.arange(len(age_order), dtype=float)
        ax.bar(x - width / 2, raw, width=width, color="#c2674f", label="Raw")
        ax.bar(x + width / 2, calibrated, width=width, color="#287c8e", label="Calibrated")
        ax.set_xticks(x, [label.title() for label in age_order])
        ax.set_title(str(model).replace("_", " "))
        ax.set_ylabel("MAE (years)")
        ax.spines[["top", "right"]].set_visible(False)
    axes.ravel()[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("ORA Error By Age Bin")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_importance(importance: pd.DataFrame, path: Path, plt: Any) -> None:
    top = _top_importance(importance, top_n=10)
    if top.empty:
        _blank_figure(path, plt, "No feature importance available")
        return
    models = list(top["model"].drop_duplicates())
    fig, axes = plt.subplots(1, len(models), figsize=(4.8 * len(models), 4.2), squeeze=False, constrained_layout=True)
    for ax, model in zip(axes.ravel(), models):
        sub = top[top["model"].eq(model)].copy()
        sub["abs_importance"] = pd.to_numeric(sub["importance"], errors="coerce").abs()
        sub = sub.sort_values("abs_importance")
        labels = [str(item).replace("__", ": ") for item in sub["feature"]]
        ax.barh(labels, sub["abs_importance"], color="#7a9f45")
        ax.set_title(str(model).replace("_", " "))
        ax.set_xlabel("Absolute importance")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Top ORA Model Features")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_ndd_projection(projection: pd.DataFrame, path: Path, plt: Any) -> None:
    required = {"model", "disease_group", "oraa"}
    if projection.empty or not required.issubset(projection.columns):
        _blank_figure(path, plt, "No NDD projection available")
        return
    frame = projection[projection["model"].isin(_display_models())].copy()
    if frame.empty:
        frame = projection.copy()
    frame["oraa"] = pd.to_numeric(frame["oraa"], errors="coerce")
    frame = frame[frame["oraa"].notna()]
    if frame.empty:
        _blank_figure(path, plt, "No NDD projection ORAA available")
        return
    order = [group for group in ["healthy", "ad", "pd"] if group in set(frame["disease_group"].astype(str))]
    if not order:
        order = sorted(frame["disease_group"].astype(str).unique())
    models = list(frame["model"].drop_duplicates())
    fig, axes = plt.subplots(1, len(models), figsize=(4.5 * len(models), 4), squeeze=False, constrained_layout=True)
    colors = {"healthy": "#4d7fb8", "ad": "#c2674f", "pd": "#7a9f45"}
    for ax, model in zip(axes.ravel(), models):
        sub = frame[frame["model"].eq(model)].copy()
        values = [sub[sub["disease_group"].astype(str).eq(group)]["oraa"].to_numpy(dtype=float) for group in order]
        positions = np.arange(len(order))
        for pos, group_values, group in zip(positions, values, order):
            if group_values.size == 0:
                continue
            jitter = np.linspace(-0.08, 0.08, group_values.size) if group_values.size > 1 else np.array([0.0])
            ax.scatter(
                np.full(group_values.size, pos) + jitter,
                group_values,
                s=24,
                alpha=0.72,
                color=colors.get(group, "#666666"),
                edgecolor="none",
            )
            ax.hlines(np.nanmean(group_values), pos - 0.22, pos + 0.22, color="#222222", linewidth=1.1)
        ax.axhline(0, color="#555555", linewidth=0.8, linestyle="--")
        ax.set_xticks(positions, [group.upper() if group in {"ad", "pd"} else group for group in order])
        ax.set_title(str(model).replace("_", " "))
        ax.set_ylabel("ORA acceleration")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Frozen Healthy-Trained ORA Projection")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_ndd_matched_reference(projection: pd.DataFrame, path: Path, plt: Any) -> None:
    required = {"model", "disease_group", "chemistry", "collection_method", "oraa"}
    if projection.empty or not required.issubset(projection.columns):
        _blank_figure(path, plt, "No matched NDD reference available")
        return
    frame = projection.copy()
    frame["disease_group"] = frame["disease_group"].astype(str)
    frame["model"] = frame["model"].astype(str)
    frame["oraa"] = pd.to_numeric(frame["oraa"], errors="coerce")
    disease = frame[frame["disease_group"].isin(["ad", "pd"]) & frame["oraa"].notna()].copy()
    healthy = frame[frame["disease_group"].eq("healthy") & frame["oraa"].notna()].copy()
    if disease.empty or healthy.empty:
        _blank_figure(path, plt, "No disease or healthy matched reference donors")
        return
    for col in ["chemistry", "collection_method"]:
        values = set(disease[col].dropna().astype(str))
        if values:
            healthy = healthy[healthy[col].astype(str).isin(values)].copy()
    if healthy.empty:
        _blank_figure(path, plt, "No matched healthy reference donors")
        return
    models = [model for model in ["random_forest", "xgboost", "catboost", "boosted_ensemble"] if model in set(frame["model"])]
    if not models:
        models = list(frame["model"].drop_duplicates())[:4]
    groups = ["matched_healthy", "ad", "pd"]
    colors = {"matched_healthy": "#4d7fb8", "ad": "#c2674f", "pd": "#7a9f45"}
    fig, axes = plt.subplots(1, len(models), figsize=(4.4 * len(models), 4), squeeze=False, constrained_layout=True)
    for ax, model in zip(axes.ravel(), models):
        pieces = [
            healthy[healthy["model"].eq(model)].assign(reference_group="matched_healthy"),
            disease[disease["model"].eq(model) & disease["disease_group"].eq("ad")].assign(reference_group="ad"),
            disease[disease["model"].eq(model) & disease["disease_group"].eq("pd")].assign(reference_group="pd"),
        ]
        sub = pd.concat(pieces, ignore_index=True)
        for pos, group in enumerate(groups):
            values = sub[sub["reference_group"].eq(group)]["oraa"].to_numpy(dtype=float)
            if values.size == 0:
                continue
            jitter = np.linspace(-0.09, 0.09, values.size) if values.size > 1 else np.array([0.0])
            ax.scatter(
                np.full(values.size, pos) + jitter,
                values,
                s=26,
                alpha=0.75,
                color=colors[group],
                edgecolor="none",
            )
            ax.hlines(np.nanmean(values), pos - 0.23, pos + 0.23, color="#222222", linewidth=1.1)
        ax.axhline(0, color="#555555", linewidth=0.8, linestyle="--")
        ax.set_xticks(np.arange(len(groups)), ["Matched\nhealthy", "AD", "PD"])
        ax.set_title(str(model).replace("_", " "))
        ax.set_ylabel("ORA acceleration")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("FLEX v2/Device-Matched NDD ORA Reference")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_module_scores(module_summary: pd.DataFrame, path: Path, plt: Any) -> None:
    required = {"donor_id", "sample_id", "coarse_cell_type", "fine_cell_type", "module", "n_cells", "mean_score"}
    if module_summary.empty or not required.issubset(module_summary.columns):
        _blank_figure(path, plt, "No module scores available")
        return
    work = module_summary.copy()
    work["n_cells"] = pd.to_numeric(work["n_cells"], errors="coerce").fillna(0)
    work["mean_score"] = pd.to_numeric(work["mean_score"], errors="coerce").fillna(0)
    totals = (
        work.drop_duplicates(["donor_id", "sample_id", "coarse_cell_type", "fine_cell_type"], keep="first")
        .groupby("fine_cell_type", observed=True)["n_cells"]
        .sum()
        .sort_values(ascending=False)
        .head(18)
    )
    top_states = totals.index.tolist()
    work = work[work["fine_cell_type"].isin(top_states)].copy()
    work["_weighted"] = work["mean_score"] * work["n_cells"]
    weighted = (
        work.groupby(["module", "fine_cell_type"], observed=True)[["_weighted", "n_cells"]]
        .sum()
        .reset_index()
    )
    weighted["mean_score"] = weighted["_weighted"] / weighted["n_cells"].replace(0, np.nan)
    pivot = weighted.pivot(index="module", columns="fine_cell_type", values="mean_score").fillna(0)
    pivot = pivot.reindex(columns=top_states)
    scaled = pivot.sub(pivot.mean(axis=1), axis=0).div(pivot.std(axis=1).replace(0, 1), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    image = ax.imshow(scaled.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(scaled.shape[0]), [str(item).replace("_", " ") for item in scaled.index], fontsize=8)
    ax.set_xticks(range(scaled.shape[1]), [str(item).replace("_", " ") for item in scaled.columns], rotation=45, ha="right", fontsize=8)
    ax.set_title("Module Scores Across Common Cell States")
    ax.set_xlabel("Fine cell state")
    ax.set_ylabel("Module")
    fig.colorbar(image, ax=ax, label="Row-scaled mean score")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_pseudobulk_de(pseudobulk_de: pd.DataFrame, path: Path, plt: Any, top_n: int) -> None:
    top = rank_pseudobulk_de(pseudobulk_de, top_n=top_n)
    if top.empty:
        _blank_figure(path, plt, "No pseudobulk DE hits available")
        return
    frame = top.copy()
    frame["fdr"] = pd.to_numeric(frame["fdr"], errors="coerce").clip(lower=np.finfo(float).tiny)
    frame["log2fc"] = pd.to_numeric(frame["log2fc"], errors="coerce")
    frame["neg_log10_fdr"] = -np.log10(frame["fdr"])
    frame = frame.sort_values("neg_log10_fdr")
    labels = [
        f"{row.gene} | {str(row.fine_cell_type).replace('_', ' ')} | {str(row.contrast).replace('_', ' ')}"
        for row in frame.itertuples()
    ]
    colors = np.where(frame["log2fc"] >= 0, "#287c8e", "#c2674f")
    fig, ax = plt.subplots(figsize=(9.5, max(4, 0.34 * len(frame))), constrained_layout=True)
    ax.barh(labels, frame["neg_log10_fdr"], color=colors)
    ax.set_xlabel("-log10 FDR")
    ax.set_title("Top Targeted Pseudobulk DE Hits")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_pseudobulk_covariate_de(pseudobulk_de: pd.DataFrame, path: Path, plt: Any, top_n: int) -> None:
    top = rank_pseudobulk_covariate_de(pseudobulk_de, top_n=top_n)
    if top.empty:
        _blank_figure(path, plt, "No covariate-adjusted pseudobulk DE hits available")
        return
    frame = top.copy()
    frame["fdr"] = pd.to_numeric(frame["fdr"], errors="coerce").clip(lower=np.finfo(float).tiny)
    frame["log2fc_adjusted"] = pd.to_numeric(frame["log2fc_adjusted"], errors="coerce")
    frame["neg_log10_fdr"] = -np.log10(frame["fdr"])
    frame = frame.sort_values("neg_log10_fdr")
    labels = [
        f"{row.gene} | {str(row.fine_cell_type).replace('_', ' ')} | {str(row.contrast).replace('_', ' ')}"
        for row in frame.itertuples()
    ]
    colors = np.where(frame["log2fc_adjusted"] >= 0, "#287c8e", "#c2674f")
    fig, ax = plt.subplots(figsize=(9.5, max(4, 0.34 * len(frame))), constrained_layout=True)
    ax.barh(labels, frame["neg_log10_fdr"], color=colors)
    ax.set_xlabel("-log10 FDR")
    ax.set_title("Top Covariate-Adjusted Pseudobulk DE Hits")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _top_importance(importance: pd.DataFrame, top_n: int = 10) -> pd.DataFrame:
    if importance.empty or not {"model", "feature", "importance"}.issubset(importance.columns):
        return pd.DataFrame(columns=["model", "feature", "importance", "stability"])
    frame = importance[importance["model"].isin(_display_models())].copy()
    if frame.empty:
        frame = importance.copy()
    frame["abs_importance"] = pd.to_numeric(frame["importance"], errors="coerce").abs()
    frame = frame[np.isfinite(frame["abs_importance"])]
    pieces = []
    for model, sub in frame.groupby("model", sort=False):
        pieces.append(sub.sort_values("abs_importance", ascending=False).head(top_n))
    if not pieces:
        return pd.DataFrame(columns=["model", "feature", "importance", "stability"])
    return pd.concat(pieces, ignore_index=True).drop(columns=["abs_importance"], errors="ignore")


def _module_feature_summary(donor_module_features: pd.DataFrame | None) -> pd.DataFrame:
    if donor_module_features is None or donor_module_features.empty:
        return pd.DataFrame(columns=["module", "mean", "sd", "min", "max"])
    module_cols = [col for col in donor_module_features.columns if col.startswith("module_score__")]
    rows = []
    for col in module_cols:
        values = pd.to_numeric(donor_module_features[col], errors="coerce")
        rows.append(
            {
                "module": col.replace("module_score__", ""),
                "mean": values.mean(),
                "sd": values.std(ddof=0),
                "min": values.min(),
                "max": values.max(),
            }
        )
    return pd.DataFrame(rows).sort_values("sd", ascending=False).reset_index(drop=True)


def _has_module_tables(
    module_summary: pd.DataFrame | None,
    module_coverage: pd.DataFrame | None,
    donor_module_features: pd.DataFrame | None,
) -> bool:
    return any(frame is not None and not frame.empty for frame in [module_summary, module_coverage, donor_module_features])


def _has_ora_diagnostics(
    ora_calibration: pd.DataFrame | None,
    ora_age_bin_errors: pd.DataFrame | None,
    ora_residual_diagnostics: pd.DataFrame | None,
) -> bool:
    return any(frame is not None and not frame.empty for frame in [ora_calibration, ora_age_bin_errors, ora_residual_diagnostics])


def _ora_calibration_summary_sentence(ora_calibration: pd.DataFrame | None) -> str:
    if ora_calibration is None or ora_calibration.empty:
        return "_No ORA calibration diagnostics available._"
    frame = ora_calibration.copy()
    frame["mae"] = pd.to_numeric(frame.get("mae"), errors="coerce")
    frame["calibrated_mae"] = pd.to_numeric(frame.get("calibrated_mae"), errors="coerce")
    non_null = frame[~frame["model"].eq("null_model")] if "model" in frame else frame
    if non_null.empty:
        non_null = frame
    best = non_null.sort_values("mae").iloc[0]
    best_calibrated = non_null.sort_values("calibrated_mae").iloc[0]
    return (
        f"Out-of-fold calibration diagnostics compare raw ORA to a simple linear recalibration of predicted age. "
        f"Best raw MAE is {_format_table_value(best.get('mae'))} for `{best.get('model')}`; "
        f"best recalibrated MAE is {_format_table_value(best_calibrated.get('calibrated_mae'))} "
        f"for `{best_calibrated.get('model')}`."
    )


def _diagnostic_model_subset(
    diagnostics: pd.DataFrame | None,
    *,
    group_col: str,
    levels: list[str],
) -> pd.DataFrame:
    columns = ["model", "group", "level", "n", "mean_error", "mae", "calibrated_mean_error", "calibrated_mae"]
    if diagnostics is None or diagnostics.empty or not set(columns).issubset(diagnostics.columns):
        return pd.DataFrame(columns=columns)
    frame = diagnostics.copy()
    frame = frame[frame[group_col].eq("age_bin") & frame["level"].astype(str).isin(levels)]
    frame["_level_order"] = frame["level"].astype(str).map({level: idx for idx, level in enumerate(levels)}).fillna(len(levels))
    frame["_model_order"] = frame["model"].astype(str).map(_model_order_map()).fillna(99)
    return (
        frame.sort_values(["_model_order", "_level_order"])[columns]
        .reset_index(drop=True)
    )


def _top_residual_diagnostics(diagnostics: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = ["model", "group", "level", "n", "mean_error", "mae", "mean_oraa", "calibrated_mean_error"]
    if diagnostics is None or diagnostics.empty or not set(columns).issubset(diagnostics.columns):
        return pd.DataFrame(columns=columns)
    frame = diagnostics.copy()
    frame = frame[~frame["group"].eq("age_bin")].copy()
    frame["n"] = pd.to_numeric(frame["n"], errors="coerce")
    frame["mean_error"] = pd.to_numeric(frame["mean_error"], errors="coerce")
    frame["abs_mean_error"] = frame["mean_error"].abs()
    return (
        frame[frame["n"].ge(5)]
        .sort_values(["abs_mean_error", "mae"], ascending=[False, False])
        .head(top_n)[columns]
        .reset_index(drop=True)
    )


def _model_order_map() -> dict[str, int]:
    return {
        "null_model": 0,
        "ridge": 1,
        "lasso": 2,
        "elastic_net": 3,
        "random_forest": 4,
        "extra_trees": 5,
        "gradient_boosting": 6,
        "tree_ensemble": 7,
        "xgboost": 8,
        "lightgbm": 9,
        "catboost": 10,
        "boosted_ensemble": 11,
    }


def _display_models() -> list[str]:
    return [
        "ridge",
        "lasso",
        "elastic_net",
        "random_forest",
        "extra_trees",
        "gradient_boosting",
        "hist_gradient_boosting",
        "tree_ensemble",
        "xgboost",
        "lightgbm",
        "catboost",
        "boosted_ensemble",
    ]


def _calibration_display_models(
    calibration: pd.DataFrame | None,
    fallback: pd.DataFrame,
    *,
    max_models: int = 4,
) -> list[str]:
    if calibration is not None and not calibration.empty and {"model", "calibrated_mae"}.issubset(calibration.columns):
        frame = calibration[~calibration["model"].astype(str).eq("null_model")].copy()
        frame["calibrated_mae"] = pd.to_numeric(frame["calibrated_mae"], errors="coerce")
        frame = frame.dropna(subset=["calibrated_mae"]).sort_values(["calibrated_mae", "model"])
        models = frame["model"].astype(str).head(max_models).tolist()
        if models:
            return models
    if "model" not in fallback:
        return []
    models = [
        model
        for model in _display_models()
        if model in set(fallback["model"].astype(str)) and model != "null_model"
    ]
    return models[:max_models]


def _has_external_validation_tables(
    external_validation_summary: pd.DataFrame | None,
    external_gene_list_coverage: pd.DataFrame | None,
    external_feature_contract: pd.DataFrame | None,
) -> bool:
    return any(
        frame is not None and not frame.empty
        for frame in [external_validation_summary, external_gene_list_coverage, external_feature_contract]
    )


def _has_external_10x_tables(
    external_sample_metadata: pd.DataFrame | None,
    external_10x_sample_qc: pd.DataFrame | None,
    external_10x_module_contrasts: pd.DataFrame | None,
) -> bool:
    return any(
        frame is not None and not frame.empty
        for frame in [external_sample_metadata, external_10x_sample_qc, external_10x_module_contrasts]
    )


def _has_external_marker_tables(
    external_10x_marker_composition: pd.DataFrame | None,
    external_10x_marker_contrasts: pd.DataFrame | None,
) -> bool:
    return any(
        frame is not None and not frame.empty
        for frame in [external_10x_marker_composition, external_10x_marker_contrasts]
    )


def _has_external_marker_age_concordance(external_marker_age_concordance: pd.DataFrame | None) -> bool:
    return external_marker_age_concordance is not None and not external_marker_age_concordance.empty


def _has_external_validation_evidence(external_validation_evidence: pd.DataFrame | None) -> bool:
    return external_validation_evidence is not None and not external_validation_evidence.empty


def _external_evidence_summary_sentence(external_validation_evidence: pd.DataFrame | None) -> str:
    if external_validation_evidence is None or external_validation_evidence.empty:
        return "_No external validation evidence ledger available._"
    frame = external_validation_evidence
    n_rows = frame.shape[0]
    strengths = (
        frame["validation_strength"].dropna().astype(str).value_counts().to_dict()
        if "validation_strength" in frame
        else {}
    )
    sanity_rows = int(
        frame["supports_primary_claim"].astype(str).eq("sanity_only").sum()
        if "supports_primary_claim" in frame
        else 0
    )
    strength_text = ", ".join(f"{key}: {_format_int(value)}" for key, value in sorted(strengths.items()))
    return (
        f"The external evidence ledger has {_format_int(n_rows)} rows across configured sources and generated "
        f"sanity checks. Validation-strength classes: {strength_text or 'none'}. "
        f"{_format_int(sanity_rows)} {_plural(sanity_rows, 'row is', 'rows are')} explicitly gated as sanity-only."
    )


def _external_validation_summary_sentence(
    external_validation_summary: pd.DataFrame | None,
    external_gene_list_coverage: pd.DataFrame | None,
) -> str:
    dataset_count = (
        _format_int(external_validation_summary["dataset_id"].nunique())
        if external_validation_summary is not None
        and not external_validation_summary.empty
        and "dataset_id" in external_validation_summary
        else "0"
    )
    feature_ready = (
        int(external_validation_summary["ready_for_feature_validation"].astype(bool).sum())
        if external_validation_summary is not None
        and not external_validation_summary.empty
        and "ready_for_feature_validation" in external_validation_summary
        else 0
    )
    raw_ready = (
        int(external_validation_summary["ready_for_raw_adapter"].astype(bool).sum())
        if external_validation_summary is not None
        and not external_validation_summary.empty
        and "ready_for_raw_adapter" in external_validation_summary
        else 0
    )
    gene_lists = (
        _format_int(external_gene_list_coverage["gene_list"].nunique())
        if external_gene_list_coverage is not None
        and not external_gene_list_coverage.empty
        and "gene_list" in external_gene_list_coverage
        else "0"
    )
    return (
        f"External validation registry currently tracks {dataset_count} candidate datasets; "
        f"{_format_int(feature_ready)} have all files needed for donor-level feature replication and "
        f"{_format_int(raw_ready)} {_plural(raw_ready, 'is', 'are')} ready for raw-adapter work. "
        f"Published validation gene-list coverage is available for {gene_lists} lists and can now be scored "
        "through the module-scoring command by passing the external dataset config."
    )


def _external_10x_summary_sentence(
    sample_metadata: pd.DataFrame | None,
    sample_qc: pd.DataFrame | None,
    module_contrasts: pd.DataFrame | None,
) -> str:
    usable = (
        int(sample_metadata["usable_for_external_validation"].astype(bool).sum())
        if sample_metadata is not None
        and not sample_metadata.empty
        and "usable_for_external_validation" in sample_metadata
        else 0
    )
    samples = (
        _format_int(sample_qc["sample_id"].nunique())
        if sample_qc is not None and not sample_qc.empty and "sample_id" in sample_qc
        else "0"
    )
    modules = (
        _format_int(module_contrasts["module"].nunique())
        if module_contrasts is not None and not module_contrasts.empty and "module" in module_contrasts
        else "0"
    )
    return (
        f"GSE184117 public metadata resolves {_format_int(usable)} usable biopsy samples with age and "
        f"olfaction status. Raw 10x module scoring ran for {samples} samples and {modules} modules. "
        "These sample-level contrasts are descriptive only because donor n is 3 versus 3 and cell labels "
        "are still absent from the public archive."
    )


def _top_external_module_contrasts(module_contrasts: pd.DataFrame | None) -> pd.DataFrame:
    columns = [
        "module",
        "n_healthy",
        "n_presbyosmia",
        "mean_healthy",
        "mean_presbyosmia",
        "presbyosmia_minus_healthy",
        "p_value",
        "direction",
        "status",
    ]
    if module_contrasts is None or module_contrasts.empty:
        return pd.DataFrame(columns=columns)
    frame = module_contrasts.copy()
    frame["abs_delta"] = pd.to_numeric(frame.get("presbyosmia_minus_healthy"), errors="coerce").abs()
    frame["p_value"] = pd.to_numeric(frame.get("p_value"), errors="coerce")
    return frame.sort_values(["p_value", "abs_delta", "module"], ascending=[True, False, True])[columns].reset_index(drop=True)


def _external_marker_summary_sentence(
    marker_composition: pd.DataFrame | None,
    marker_contrasts: pd.DataFrame | None,
) -> str:
    samples = (
        _format_int(marker_composition["sample_id"].nunique())
        if marker_composition is not None
        and not marker_composition.empty
        and "sample_id" in marker_composition
        else "0"
    )
    panels = (
        _format_int(marker_contrasts["marker_panel"].nunique())
        if marker_contrasts is not None
        and not marker_contrasts.empty
        and "marker_panel" in marker_contrasts
        else "0"
    )
    return (
        f"Marker-only composition scored {panels} coarse olfactory/epithelial/immune panels across {samples} "
        "GSE184117 samples. These are conservative marker assignments for external sanity checking, not a replacement "
        "for reference-mapped cell labels."
    )


def _top_external_marker_contrasts(marker_contrasts: pd.DataFrame | None) -> pd.DataFrame:
    columns = [
        "marker_panel",
        "n_healthy",
        "n_presbyosmia",
        "mean_fraction_healthy",
        "mean_fraction_presbyosmia",
        "presbyosmia_minus_healthy",
        "p_value",
        "direction",
        "status",
    ]
    if marker_contrasts is None or marker_contrasts.empty:
        return pd.DataFrame(columns=columns)
    frame = marker_contrasts.copy()
    frame["abs_delta"] = pd.to_numeric(frame.get("presbyosmia_minus_healthy"), errors="coerce").abs()
    frame["p_value"] = pd.to_numeric(frame.get("p_value"), errors="coerce")
    return frame.sort_values(["p_value", "abs_delta", "marker_panel"], ascending=[True, False, True])[columns].reset_index(drop=True)


def _external_marker_age_concordance_sentence(concordance: pd.DataFrame | None) -> str:
    if concordance is None or concordance.empty:
        return "_No marker-age concordance table available._"
    rows = len(concordance)
    concordant = int(concordance["concordance"].astype(str).eq("concordant").sum()) if "concordance" in concordance else 0
    discordant = int(concordance["concordance"].astype(str).eq("discordant").sum()) if "concordance" in concordance else 0
    panels = _format_int(concordance["marker_panel"].nunique()) if "marker_panel" in concordance else "0"
    return (
        f"Marker-only GSE184117 presbyosmia shifts were mapped to Gateway age-associated features for {panels} panels, "
        f"yielding {_format_int(concordant)} concordant and {_format_int(discordant)} discordant marker-feature rows "
        f"out of {_format_int(rows)}. This remains a sanity check because the external cohort is 3 versus 3 and lacks cell labels."
    )


def _top_external_marker_age_concordance(concordance: pd.DataFrame | None) -> pd.DataFrame:
    columns = [
        "marker_panel",
        "gateway_feature",
        "external_direction",
        "gateway_age_direction",
        "concordance",
        "external_delta",
        "gateway_beta_per_10_years",
        "gateway_fdr",
        "status",
    ]
    if concordance is None or concordance.empty:
        return pd.DataFrame(columns=columns)
    frame = concordance.copy()
    frame["gateway_fdr"] = pd.to_numeric(frame.get("gateway_fdr"), errors="coerce")
    frame["external_abs_delta"] = pd.to_numeric(frame.get("external_delta"), errors="coerce").abs()
    order = {"concordant": 0, "discordant": 1, "not_evaluable": 2}
    frame["_order"] = frame["concordance"].astype(str).map(order).fillna(3)
    return (
        frame.sort_values(["_order", "gateway_fdr", "external_abs_delta", "marker_panel"], ascending=[True, True, False, True])
        .head(20)[columns]
        .reset_index(drop=True)
    )


def _has_latent_space_audit(
    latent_space_readiness: pd.DataFrame | None,
    latent_space_local_audit: pd.DataFrame | None,
    latent_space_portal_assets: pd.DataFrame | None,
) -> bool:
    return any(
        frame is not None and not frame.empty
        for frame in [latent_space_readiness, latent_space_local_audit, latent_space_portal_assets]
    )


def _latent_space_summary_sentence(latent_space_readiness: pd.DataFrame | None) -> str:
    if latent_space_readiness is None or latent_space_readiness.empty:
        return "_No latent-space readiness audit available._"
    row = latent_space_readiness.iloc[0]
    status = str(row.get("status", "unknown"))
    local = str(row.get("local_embeddings", "")) or "none"
    portal = str(row.get("portal_embeddings", "")) or "none"
    recommendation = str(row.get("recommendation", "Review latent-space audit before trajectory work."))
    return (
        f"Latent-space audit status is `{status}`. Local embeddings: {local}; CELLxGENE portal embeddings: {portal}. "
        f"{recommendation}"
    )


def _has_scvi_pilot_validation(scvi_pilot_validation: pd.DataFrame | None) -> bool:
    return scvi_pilot_validation is not None and not scvi_pilot_validation.empty


def _scvi_pilot_summary_sentence(scvi_pilot_validation: pd.DataFrame | None) -> str:
    if scvi_pilot_validation is None or scvi_pilot_validation.empty:
        return "_No scVI pilot validation table available._"
    pilot = scvi_pilot_validation[scvi_pilot_validation["check"].eq("pilot_h5ad")]
    embedding = scvi_pilot_validation[scvi_pilot_validation["check"].eq("embedding_dimensions")]
    gate = scvi_pilot_validation[scvi_pilot_validation["check"].eq("claim_gate")]
    pilot_detail = str(pilot["detail"].iloc[0]) if not pilot.empty else "pilot H5AD not summarized"
    embedding_status = str(embedding["status"].iloc[0]) if not embedding.empty else "unknown"
    gate_status = str(gate["status"].iloc[0]) if not gate.empty else "unknown"
    return (
        f"A recomputed scVI latent run is available ({pilot_detail}); embedding validation is `{embedding_status}`. "
        f"Downstream latent claims remain `{gate_status}` until marker-continuity, mixing, and seed-stability diagnostics mature."
    )


def _has_ndd_projection(
    ndd_projection: pd.DataFrame | None,
    ndd_projection_summary: pd.DataFrame | None,
) -> bool:
    return any(frame is not None and not frame.empty for frame in [ndd_projection, ndd_projection_summary])


def _ndd_projection_summary_sentence(ndd_projection: pd.DataFrame | None) -> str:
    if ndd_projection is None or ndd_projection.empty:
        return "_No frozen ORA projection rows available._"
    frame = ndd_projection.copy()
    n_donors = frame["donor_id"].nunique() if "donor_id" in frame else None
    n_train = (
        frame.loc[frame["is_training_donor"].astype(bool), "donor_id"].nunique()
        if {"donor_id", "is_training_donor"}.issubset(frame.columns)
        else None
    )
    n_ndd = (
        frame.loc[frame["disease_group"].astype(str).isin(["ad", "pd"]), "donor_id"].nunique()
        if {"donor_id", "disease_group"}.issubset(frame.columns)
        else None
    )
    return (
        f"Frozen ORA models were trained on {_format_int(n_train)} healthy age-known donors "
        f"and projected onto {_format_int(n_donors)} donors total, including {_format_int(n_ndd)} AD/PD donors."
    )


def _top_ndd_projection_rows(ndd_projection: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = ["donor_id", "model", "disease_group", "chronological_age", "ora", "oraa"]
    if ndd_projection is None or ndd_projection.empty or not set(columns).issubset(ndd_projection.columns):
        return pd.DataFrame(columns=columns)
    frame = ndd_projection[
        ndd_projection["disease_group"].astype(str).isin(["ad", "pd"])
        & ndd_projection["model"].isin(["elastic_net", "random_forest"])
    ].copy()
    if frame.empty:
        return pd.DataFrame(columns=columns)
    frame["abs_oraa"] = pd.to_numeric(frame["oraa"], errors="coerce").abs()
    frame = frame[frame["abs_oraa"].notna()]
    return frame.sort_values(["abs_oraa", "model"], ascending=[False, True]).head(top_n)[columns].reset_index(drop=True)


def _top_ndd_feature_comparison(comparison: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = [
        "model",
        "disease_group",
        "composition_mean_oraa",
        "augmented_mean_oraa",
        "augmented_minus_composition_oraa",
        "sign_stable_negative",
    ]
    if comparison is None or comparison.empty or not set(columns).issubset(comparison.columns):
        return pd.DataFrame(columns=columns)
    frame = comparison[comparison["model"].isin(_display_models())].copy()
    if frame.empty:
        frame = comparison.copy()
    frame["abs_augmented_minus_composition_oraa"] = pd.to_numeric(
        frame["augmented_minus_composition_oraa"],
        errors="coerce",
    ).abs()
    return (
        frame.sort_values(["disease_group", "abs_augmented_minus_composition_oraa", "model"], ascending=[True, False, True])
        .head(top_n)[columns]
        .reset_index(drop=True)
    )


def _top_ndd_donor_appendix(appendix: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = [
        "feature_set",
        "donor_id",
        "disease_group",
        "chronological_age",
        "chemistry",
        "collection_method",
        "model",
        "ora",
        "oraa",
    ]
    if appendix is None or appendix.empty or not set(columns).issubset(appendix.columns):
        return pd.DataFrame(columns=columns)
    frame = appendix[appendix["model"].isin(["random_forest", "xgboost", "catboost", "boosted_ensemble"])].copy()
    if frame.empty:
        frame = appendix.copy()
    frame["abs_oraa"] = pd.to_numeric(frame["oraa"], errors="coerce").abs()
    return frame.sort_values(["abs_oraa", "feature_set", "model"], ascending=[False, True, True]).head(top_n)[columns].reset_index(drop=True)


def _top_ndd_projection_diagnostics(diagnostics: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = [
        "model",
        "disease_group",
        "diagnostic",
        "level",
        "n_donors",
        "mean_age",
        "median_total_cells",
        "mean_oraa",
        "status",
    ]
    if diagnostics is None or diagnostics.empty or not set(columns).issubset(diagnostics.columns):
        return pd.DataFrame(columns=columns)
    frame = diagnostics[diagnostics["model"].isin(["random_forest", "xgboost", "catboost", "boosted_ensemble"])].copy()
    if frame.empty:
        frame = diagnostics.copy()
    frame["abs_mean_oraa"] = pd.to_numeric(frame["mean_oraa"], errors="coerce").abs()
    return (
        frame.sort_values(["status", "disease_group", "diagnostic", "abs_mean_oraa"], ascending=[True, True, True, False])
        .head(top_n)[columns]
        .reset_index(drop=True)
    )


def _has_pseudobulk_tables(
    pseudobulk_de: pd.DataFrame | None,
    pseudobulk_coverage: pd.DataFrame | None,
    pseudobulk_metadata: pd.DataFrame | None,
) -> bool:
    return any(frame is not None and not frame.empty for frame in [pseudobulk_de, pseudobulk_coverage, pseudobulk_metadata])


def _pseudobulk_summary_sentence(
    pseudobulk_de: pd.DataFrame | None,
    pseudobulk_metadata: pd.DataFrame | None,
) -> str:
    if pseudobulk_de is None or pseudobulk_de.empty:
        return "_No pseudobulk DE rows available._"
    ok = pseudobulk_de[pseudobulk_de["status"].eq("ok")].copy() if "status" in pseudobulk_de else pseudobulk_de.copy()
    contrasts = _format_int(ok["contrast"].nunique()) if "contrast" in ok else "not recorded"
    genes = _format_int(ok["gene"].nunique()) if "gene" in ok else "not recorded"
    states = _format_int(ok["fine_cell_type"].nunique()) if "fine_cell_type" in ok else "not recorded"
    groups = _format_int(len(pseudobulk_metadata)) if pseudobulk_metadata is not None and not pseudobulk_metadata.empty else "not recorded"
    cells = (
        _format_int(pd.to_numeric(pseudobulk_metadata.get("n_cells"), errors="coerce").sum())
        if pseudobulk_metadata is not None and not pseudobulk_metadata.empty and "n_cells" in pseudobulk_metadata
        else "not recorded"
    )
    return (
        f"Targeted pseudobulk produced {_format_int(len(ok))} valid tests across {contrasts} contrasts, "
        f"{states} fine cell states, and {genes} genes. Aggregated groups: {groups}; cells represented: {cells}."
    )


def _pseudobulk_adjusted_summary_sentence(pseudobulk_de: pd.DataFrame | None) -> str:
    if pseudobulk_de is None or pseudobulk_de.empty:
        return "_No covariate-adjusted pseudobulk DE rows available._"
    ok = pseudobulk_de[pseudobulk_de["status"].eq("ok")].copy() if "status" in pseudobulk_de else pseudobulk_de.copy()
    contrasts = _format_int(ok["contrast"].nunique()) if "contrast" in ok else "not recorded"
    genes = _format_int(ok["gene"].nunique()) if "gene" in ok else "not recorded"
    states = _format_int(ok["fine_cell_type"].nunique()) if "fine_cell_type" in ok else "not recorded"
    covariates = _format_list(sorted({item for value in ok.get("covariates", pd.Series(dtype=str)).dropna().astype(str) for item in value.split(",") if item}))
    return (
        f"Covariate-adjusted targeted pseudobulk produced {_format_int(len(ok))} valid tests across "
        f"{contrasts} contrasts, {states} fine cell states, and {genes} genes. Covariates used where variable: {covariates}."
    )


def _pseudobulk_genomewide_summary_sentence(summary: pd.DataFrame | None) -> str:
    if summary is None or summary.empty:
        return "_No genome-wide pseudobulk export summary available._"
    row = summary.iloc[0]
    return (
        f"Genome-wide pseudobulk export scanned {_format_int(row.get('n_cells'))} cells and "
        f"{_format_int(row.get('n_genes'))} genes. It retained {_format_int(row.get('n_groups_exported'))} "
        f"of {_format_int(row.get('n_groups_total'))} donor/sample/cell-state groups for downstream DE after "
        f"minimum cell and donor filters."
    )


def _pseudobulk_genomewide_qc_sentence(summary: pd.DataFrame | None) -> str:
    if summary is None or summary.empty:
        return "_No genome-wide pseudobulk QC summary available._"
    row = summary.iloc[0]
    metadata_total = float(row.get("metadata_total_counts", np.nan))
    matrix_total = float(row.get("matrix_total_counts", np.nan))
    ratio = matrix_total / metadata_total if metadata_total else np.nan
    return (
        f"QC confirmed matrix/metadata column alignment: {row.get('matrix_columns_match_metadata')}. "
        f"Median detected genes per pseudobulk group: {_format_int(row.get('median_group_detected_genes'))}; "
        f"median gene detected-group fraction: {_format_table_value(row.get('median_gene_detected_group_fraction'))}. "
        f"Matrix total counts / metadata total counts: {_format_table_value(ratio)}."
    )


def _top_genomewide_variable_genes(gene_qc: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = ["gene_symbol", "total_count", "detected_group_fraction", "variance_log1p"]
    if gene_qc is None or gene_qc.empty or not set(columns).issubset(gene_qc.columns):
        return pd.DataFrame(columns=columns)
    frame = gene_qc.copy()
    frame["variance_log1p"] = pd.to_numeric(frame["variance_log1p"], errors="coerce")
    frame["total_count"] = pd.to_numeric(frame["total_count"], errors="coerce")
    return frame.sort_values(["variance_log1p", "total_count"], ascending=[False, False]).head(top_n)[columns].reset_index(drop=True)


def _pseudobulk_genomewide_de_sentence(summary: pd.DataFrame | None) -> str:
    if summary is None or summary.empty:
        return "_No genome-wide edgeR DE summary available._"
    tested = int(pd.to_numeric(summary.get("tested_rows"), errors="coerce").fillna(0).sum())
    significant = int(pd.to_numeric(summary.get("significant_rows"), errors="coerce").fillna(0).sum())
    sex_linked = int(pd.to_numeric(summary.get("sex_linked_significant_rows"), errors="coerce").fillna(0).sum())
    states = int(pd.to_numeric(summary.get("ok_cell_state_models"), errors="coerce").fillna(0).sum())
    threshold = pd.to_numeric(summary.get("fdr_threshold"), errors="coerce").dropna()
    fdr_threshold = threshold.iloc[0] if not threshold.empty else 0.05
    return (
        f"edgeR quasi-likelihood models tested {_format_int(tested)} gene/cell-state/contrast rows across "
        f"{_format_int(states)} successful cell-state contrast models. At FDR < {_format_table_value(fdr_threshold)}, "
        f"{_format_int(significant)} rows were significant; {_format_int(sex_linked)} significant rows were in the "
        "initial sex-linked sentinel list, so top hits require sex-balance sensitivity checks."
    )


def _matched_genomewide_de_sentence(summary: pd.DataFrame | None) -> str:
    if summary is None or summary.empty:
        return "_No matched FLEX v2/device edgeR DE summary available._"
    tested = int(pd.to_numeric(summary.get("tested_rows"), errors="coerce").fillna(0).sum())
    significant = int(pd.to_numeric(summary.get("significant_rows"), errors="coerce").fillna(0).sum())
    states = int(pd.to_numeric(summary.get("ok_cell_state_models"), errors="coerce").fillna(0).sum())
    sex_linked = int(pd.to_numeric(summary.get("sex_linked_significant_rows"), errors="coerce").fillna(0).sum())
    return (
        f"The matched FLEX v2/device subset reran edgeR on {_format_int(tested)} tested rows across "
        f"{_format_int(states)} successful cell-state contrast models. It found {_format_int(significant)} "
        f"FDR-significant rows, including {_format_int(sex_linked)} sex-linked sentinel rows."
    )


def _genomewide_de_method_sentence(summary: pd.DataFrame | None, method_label: str) -> str:
    if summary is None or summary.empty:
        return f"_No {method_label} genome-wide DE summary available._"
    tested = int(pd.to_numeric(summary.get("tested_rows"), errors="coerce").fillna(0).sum())
    significant = int(pd.to_numeric(summary.get("significant_rows"), errors="coerce").fillna(0).sum())
    states = int(pd.to_numeric(summary.get("ok_cell_state_models"), errors="coerce").fillna(0).sum())
    sex_linked = int(pd.to_numeric(summary.get("sex_linked_significant_rows"), errors="coerce").fillna(0).sum())
    threshold = pd.to_numeric(summary.get("fdr_threshold"), errors="coerce").dropna()
    fdr_threshold = threshold.iloc[0] if not threshold.empty else 0.05
    return (
        f"{method_label} models tested {_format_int(tested)} gene/cell-state/contrast rows across "
        f"{_format_int(states)} successful cell-state contrast models. At FDR < {_format_table_value(fdr_threshold)}, "
        f"{_format_int(significant)} rows were significant, including {_format_int(sex_linked)} sex-linked sentinel "
        f"{_plural(sex_linked, 'row', 'rows')}."
    )


def _non_sex_linked_genomewide_de_hits(de_hits: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = ["contrast", "fine_cell_type", "gene_symbol", "log2fc", "p_value", "fdr", "is_sex_linked_initial"]
    if de_hits is None or de_hits.empty or not set(columns).issubset(de_hits.columns):
        return pd.DataFrame(columns=columns)
    frame = de_hits.copy()
    is_sex = frame["is_sex_linked_initial"].astype(str).str.lower().isin({"true", "1", "yes"})
    frame = frame.loc[~is_sex, columns].copy()
    frame["fdr"] = pd.to_numeric(frame["fdr"], errors="coerce")
    frame["p_value"] = pd.to_numeric(frame["p_value"], errors="coerce")
    return frame.sort_values(["fdr", "p_value", "contrast", "fine_cell_type", "gene_symbol"]).head(top_n).reset_index(drop=True)


def _top_donor_balance_issues(balance: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = ["contrast", "fine_cell_type", "n_case", "n_control", "status", "balance_status"]
    if balance is None or balance.empty or not set(columns).issubset(balance.columns):
        return pd.DataFrame(columns=columns)
    frame = balance.copy()
    frame["n_case"] = pd.to_numeric(frame["n_case"], errors="coerce")
    frame["n_control"] = pd.to_numeric(frame["n_control"], errors="coerce")
    issue_order = {
        "low_case_and_control_donors": 0,
        "low_case_donors": 1,
        "low_control_donors": 2,
        "ok": 3,
    }
    frame["_issue_order"] = frame["balance_status"].map(issue_order).fillna(4)
    return frame.sort_values(["_issue_order", "contrast", "fine_cell_type"]).head(top_n)[columns].reset_index(drop=True)


def _top_matched_de_feasibility(feasibility: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = ["contrast", "fine_cell_type", "n_case", "n_matched_healthy", "ready_for_matched_de"]
    if feasibility is None or feasibility.empty or not set(columns).issubset(feasibility.columns):
        return pd.DataFrame(columns=columns)
    frame = feasibility.copy()
    frame["ready_for_matched_de"] = frame["ready_for_matched_de"].astype(bool)
    frame["n_case"] = pd.to_numeric(frame["n_case"], errors="coerce")
    frame["n_matched_healthy"] = pd.to_numeric(frame["n_matched_healthy"], errors="coerce")
    return (
        frame.sort_values(["ready_for_matched_de", "contrast", "n_case", "n_matched_healthy"], ascending=[False, True, False, False])
        .head(top_n)[columns]
        .reset_index(drop=True)
    )


def _ora_sensitivity_summary_sentence(
    scenarios: pd.DataFrame | None,
    performance: pd.DataFrame | None,
) -> str:
    if performance is None or performance.empty:
        return "_No ORA sensitivity performance rows available._"
    runnable = int(scenarios["status"].eq("ok").sum()) if scenarios is not None and "status" in scenarios else performance["scenario"].nunique()
    rf = performance[performance["model"].eq("random_forest")].copy() if "model" in performance else performance.copy()
    if rf.empty:
        return f"ORA sensitivity generated performance rows for {performance['scenario'].nunique()} scenarios."
    best = rf.sort_values("mae").iloc[0]
    worst = rf.sort_values("mae", ascending=False).iloc[0]
    return (
        f"ORA sensitivity reran age models across {_format_int(runnable)} runnable strata. "
        f"Best random-forest MAE was {_format_table_value(best.get('mae'))} in `{best.get('scenario')}`; "
        f"weakest was {_format_table_value(worst.get('mae'))} in `{worst.get('scenario')}`."
    )


def _top_ora_sensitivity_performance(
    performance: pd.DataFrame | None,
    model: str = "random_forest",
) -> pd.DataFrame:
    columns = ["scenario", "model", "n", "mae", "rmse", "r2", "spearman_r", "healthy_train_donors"]
    if performance is None or performance.empty or not set(columns).issubset(performance.columns):
        return pd.DataFrame(columns=columns)
    frame = performance[performance["model"].eq(model)].copy()
    if frame.empty:
        frame = performance.copy()
    frame["mae"] = pd.to_numeric(frame["mae"], errors="coerce")
    return frame.sort_values(["mae", "scenario"]).head(20)[columns].reset_index(drop=True)


def _ora_repeated_cv_summary_sentence(summary: pd.DataFrame | None) -> str:
    if summary is None or summary.empty:
        return "_No repeated-CV ORA summary available._"
    frame = summary.copy()
    if "mae_mean" in frame:
        frame["mae_mean"] = pd.to_numeric(frame["mae_mean"], errors="coerce")
        frame = frame.sort_values(["mae_mean", "model"])
    row = frame.iloc[0]
    model = row.get("model", "best model")
    return (
        f"Repeated donor-level CV used {_format_int(row.get('repeats'))} repeats. "
        f"Best mean-MAE model was {model} with MAE {_format_table_value(row.get('mae_mean'))} "
        f"({_format_table_value(row.get('mae_ci_low'))}-{_format_table_value(row.get('mae_ci_high'))}); "
        f"Spearman r mean was {_format_table_value(row.get('spearman_r_mean'))} "
        f"({_format_table_value(row.get('spearman_r_ci_low'))}-{_format_table_value(row.get('spearman_r_ci_high'))})."
    )


def _top_feature_set_model_comparison(comparison: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = [
        "feature_set",
        "model",
        "mae_mean",
        "mae_ci_low",
        "mae_ci_high",
        "rmse_mean",
        "r2_mean",
        "spearman_r_mean",
    ]
    if comparison is None or comparison.empty or not set(columns).issubset(comparison.columns):
        return pd.DataFrame(columns=columns)
    frame = comparison.copy()
    frame["mae_mean"] = pd.to_numeric(frame["mae_mean"], errors="coerce")
    return frame.sort_values(["mae_mean", "rmse_mean", "model"]).head(top_n)[columns].reset_index(drop=True)


def _top_candidate_repeated_cv(summary: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = [
        "model",
        "repeats",
        "n",
        "mae_mean",
        "mae_ci_low",
        "mae_ci_high",
        "spearman_r_mean",
        "spearman_r_ci_low",
        "spearman_r_ci_high",
    ]
    if summary is None or summary.empty or not set(columns).issubset(summary.columns):
        return pd.DataFrame(columns=columns)
    return _sort_metric_table(summary, "mae_mean").head(top_n)[columns].reset_index(drop=True)


def _top_permutation_empirical(permutation: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = [
        "model",
        "n_permutations",
        "observed_mae",
        "null_mae_mean",
        "empirical_p_mae",
        "observed_spearman_r",
        "null_spearman_r_mean",
        "empirical_p_spearman_r",
    ]
    if permutation is None or permutation.empty or not set(columns).issubset(permutation.columns):
        return pd.DataFrame(columns=columns)
    frame = permutation.copy()
    frame["observed_mae"] = pd.to_numeric(frame["observed_mae"], errors="coerce")
    return frame.sort_values(["observed_mae", "model"]).head(top_n)[columns].reset_index(drop=True)


def _top_nested_tuning_summary(summary: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = [
        "model",
        "repeats",
        "n",
        "mae_mean",
        "mae_ci_low",
        "mae_ci_high",
        "rmse_mean",
        "r2_mean",
        "spearman_r_mean",
    ]
    if summary is None or summary.empty or not set(columns).issubset(summary.columns):
        return pd.DataFrame(columns=columns)
    frame = summary.copy()
    frame["mae_mean"] = pd.to_numeric(frame["mae_mean"], errors="coerce")
    return frame.sort_values(["mae_mean", "model"]).head(top_n)[columns].reset_index(drop=True)


def _output_provenance_sentence(provenance: pd.DataFrame | None) -> str:
    if provenance is None or provenance.empty:
        return "_No output provenance table available._"
    total = len(provenance)
    exists = int(provenance["exists"].astype(bool).sum()) if "exists" in provenance else 0
    missing = total - exists
    hashed = int(provenance["checksum_status"].astype(str).eq("sha256").sum()) if "checksum_status" in provenance else 0
    return (
        f"Output provenance tracks {_format_int(total)} expected artifacts; {_format_int(exists)} are present, "
        f"{_format_int(missing)} are missing, and {_format_int(hashed)} small artifacts have SHA-256 checksums."
    )


def _top_stacking_summary(summary: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = [
        "model",
        "repeats",
        "n",
        "mae_mean",
        "mae_ci_low",
        "mae_ci_high",
        "rmse_mean",
        "r2_mean",
        "spearman_r_mean",
    ]
    if summary is None or summary.empty or not set(columns).issubset(summary.columns):
        return pd.DataFrame(columns=columns)
    frame = summary.copy()
    frame["mae_mean"] = pd.to_numeric(frame["mae_mean"], errors="coerce")
    return frame.sort_values(["mae_mean", "model"]).head(top_n)[columns].reset_index(drop=True)


def _top_repeated_cv_features(feature_stability: pd.DataFrame | None, top_n: int) -> pd.DataFrame:
    columns = ["model", "feature", "mean_importance", "selection_fraction"]
    if feature_stability is None or feature_stability.empty or not set(columns).issubset(feature_stability.columns):
        return pd.DataFrame(columns=columns)
    frame = feature_stability.copy()
    frame["selection_fraction"] = pd.to_numeric(frame["selection_fraction"], errors="coerce")
    frame["mean_importance"] = pd.to_numeric(frame["mean_importance"], errors="coerce")
    frame["abs_mean_importance"] = frame["mean_importance"].abs()
    return (
        frame.sort_values(["selection_fraction", "abs_mean_importance"], ascending=[False, False])
        .head(top_n)[columns]
        .reset_index(drop=True)
    )


def _has_feature_interpretation(feature_interpretation: pd.DataFrame | None) -> bool:
    return feature_interpretation is not None and not feature_interpretation.empty


def _feature_interpretation_summary_sentence(feature_interpretation: pd.DataFrame | None) -> str:
    if feature_interpretation is None or feature_interpretation.empty:
        return "_No feature interpretation table available._"
    frame = feature_interpretation
    themes = (
        frame["biology_theme"].dropna().astype(str).value_counts().head(4).to_dict()
        if "biology_theme" in frame
        else {}
    )
    theme_text = ", ".join(f"{theme}: {_format_int(count)}" for theme, count in themes.items())
    supported = 0
    if "n_supporting_models" in frame:
        supported = int(pd.to_numeric(frame["n_supporting_models"], errors="coerce").fillna(0).ge(2).sum())
    return (
        f"Top ORA features are mapped to broad biological themes for manuscript triage; "
        f"{_format_int(supported)} appear among top features for at least two model families. "
        f"Leading themes: {theme_text or 'not classified'}."
    )


def _pseudobulk_metadata_summary(pseudobulk_metadata: pd.DataFrame | None) -> pd.DataFrame:
    if pseudobulk_metadata is None or pseudobulk_metadata.empty or "disease_group" not in pseudobulk_metadata:
        return pd.DataFrame(columns=["disease_group", "groups", "cells"])
    frame = pseudobulk_metadata.copy()
    frame["n_cells"] = pd.to_numeric(frame.get("n_cells"), errors="coerce").fillna(0)
    summary = (
        frame.groupby("disease_group", observed=True)["n_cells"]
        .agg(groups="size", cells="sum")
        .reset_index()
        .sort_values("cells", ascending=False)
    )
    return summary


def _blank_figure(path: Path, plt: Any, message: str) -> None:
    fig, ax = plt.subplots(figsize=(5, 2.5), constrained_layout=True)
    ax.text(0.5, 0.5, message, ha="center", va="center")
    ax.set_axis_off()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _markdown_table(frame: pd.DataFrame, columns: list[str], max_rows: int) -> str:
    available = [col for col in columns if col in frame.columns]
    if frame.empty or not available:
        return "_No rows available._"
    table = frame[available].head(max_rows).copy()
    for col in table.columns:
        if pd.api.types.is_bool_dtype(table[col]):
            table[col] = table[col].map(lambda value: "True" if bool(value) else "False")
        elif pd.api.types.is_numeric_dtype(table[col]):
            table[col] = table[col].map(_format_table_value)
        else:
            table[col] = table[col].fillna("").astype(str)
    return _render_markdown_table(table)


def _render_markdown_table(table: pd.DataFrame) -> str:
    headers = [str(col) for col in table.columns]
    rows = [[_escape_table_cell(value) for value in row] for row in table.to_numpy()]
    header = "| " + " | ".join(_escape_table_cell(col) for col in headers) + " |"
    divider = "| " + " | ".join("---" for _ in headers) + " |"
    body = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([header, divider, *body])


def _escape_table_cell(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _figure_link(report_path: Path, figure_path: Path | None, alt: str) -> str:
    if figure_path is None:
        return ""
    try:
        rel = figure_path.relative_to(report_path.parent)
    except ValueError:
        rel = Path("../figures") / figure_path.name
    return f"![{alt}]({rel.as_posix()})"


def _model_summary_sentence(best_model: pd.Series) -> str:
    if best_model.empty:
        return "_No model performance rows available._"
    return (
        f"Best non-null model by MAE: **{best_model.get('model')}** "
        f"(MAE {_format_float(best_model.get('mae'))} years, "
        f"Spearman r {_format_float(best_model.get('spearman_r'))})."
    )


def _augmented_summary_sentence(performance: pd.DataFrame, augmented_performance: pd.DataFrame) -> str:
    base = best_predictive_model(performance)
    augmented = best_predictive_model(augmented_performance)
    if base.empty or augmented.empty:
        return "_Module-augmented model performance is available, but the baseline comparison is incomplete._"
    delta = float(augmented.get("mae")) - float(base.get("mae"))
    direction = "lower" if delta < 0 else "higher" if delta > 0 else "unchanged"
    return (
        f"Best module-augmented model by MAE: **{augmented.get('model')}** "
        f"(MAE {_format_float(augmented.get('mae'))} years, Spearman r {_format_float(augmented.get('spearman_r'))}); "
        f"this is {_format_float(abs(delta))} years {direction} than the best composition-only model."
    )


def _usable_training_donors(manifest: pd.DataFrame) -> int | None:
    if "donor_id" not in manifest or "usable_for_ora_training" not in manifest:
        return None
    frame = manifest[manifest["usable_for_ora_training"].astype(bool)]
    if "age" in frame:
        frame = frame[frame["age"].notna()]
    return int(frame["donor_id"].nunique())


def _donors_with_age(manifest: pd.DataFrame) -> int | None:
    if "donor_id" not in manifest or "age" not in manifest:
        return None
    return int(manifest[manifest["age"].notna()]["donor_id"].nunique())


def _format_table_value(value: object) -> str:
    if pd.isna(value):
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number.is_integer():
        return f"{int(number):,}"
    if abs(number) >= 1000 and number.is_integer():
        return f"{int(number):,}"
    if abs(number) >= 10:
        return f"{number:.2f}"
    if abs(number) >= 0.01:
        return f"{number:.4f}"
    return f"{number:.3g}"


def _format_float(value: object) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "not recorded"
    if not np.isfinite(number):
        return "not recorded"
    return f"{number:.2f}"


def _format_int(value: object) -> str:
    if value is None:
        return "not recorded"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "not recorded"
    if not np.isfinite(number):
        return "not recorded"
    return f"{int(number):,}"


def _plural(count: int, singular: str, plural: str) -> str:
    return singular if int(count) == 1 else plural


def _format_list(values: object) -> str:
    if not values:
        return "none recorded"
    if isinstance(values, (list, tuple)):
        return ", ".join(str(value) for value in values)
    return str(values)
