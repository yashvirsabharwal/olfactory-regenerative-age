#!/usr/bin/env python3
"""Generate the composition-MVP ORA report and figures."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.reporting import generate_mvp_report, load_schema


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--cohort-summary", default=None)
    parser.add_argument("--associations", default=None)
    parser.add_argument("--performance", default=None)
    parser.add_argument("--scores", default=None)
    parser.add_argument("--importance", default=None)
    parser.add_argument("--ora-calibration", default=None)
    parser.add_argument("--ora-age-bin-errors", default=None)
    parser.add_argument("--ora-residual-diagnostics", default=None)
    parser.add_argument("--ora-calibrated-scores", default=None)
    parser.add_argument("--augmented-performance", default=None)
    parser.add_argument("--augmented-scores", default=None)
    parser.add_argument("--augmented-importance", default=None)
    parser.add_argument("--ndd-projection", default=None)
    parser.add_argument("--ndd-projection-summary", default=None)
    parser.add_argument("--ndd-projection-uncertainty", default=None)
    parser.add_argument("--ndd-projection-context", default=None)
    parser.add_argument("--ndd-projection-feature-comparison", default=None)
    parser.add_argument("--ndd-projection-donor-appendix", default=None)
    parser.add_argument("--ndd-projection-diagnostics", default=None)
    parser.add_argument("--module-summary", default=None)
    parser.add_argument("--module-coverage", default=None)
    parser.add_argument("--donor-module-features", default=None)
    parser.add_argument("--external-validation-summary", default=None)
    parser.add_argument("--external-gene-list-coverage", default=None)
    parser.add_argument("--external-feature-contract", default=None)
    parser.add_argument("--external-sample-metadata", default=None)
    parser.add_argument("--external-10x-sample-qc", default=None)
    parser.add_argument("--external-10x-module-contrasts", default=None)
    parser.add_argument("--external-10x-marker-composition", default=None)
    parser.add_argument("--external-10x-marker-contrasts", default=None)
    parser.add_argument("--external-marker-age-concordance", default=None)
    parser.add_argument("--external-validation-evidence", default=None)
    parser.add_argument("--pseudobulk-de", default=None)
    parser.add_argument("--pseudobulk-coverage", default=None)
    parser.add_argument("--pseudobulk-metadata", default=None)
    parser.add_argument("--pseudobulk-covariate-de", default=None)
    parser.add_argument("--pseudobulk-genomewide-summary", default=None)
    parser.add_argument("--pseudobulk-genomewide-qc-summary", default=None)
    parser.add_argument("--pseudobulk-genomewide-gene-qc", default=None)
    parser.add_argument("--pseudobulk-genomewide-disease-summary", default=None)
    parser.add_argument("--pseudobulk-genomewide-de-summary", default=None)
    parser.add_argument("--pseudobulk-genomewide-de-top-hits", default=None)
    parser.add_argument("--pseudobulk-genomewide-de-audit", default=None)
    parser.add_argument("--pseudobulk-genomewide-donor-balance", default=None)
    parser.add_argument("--pseudobulk-genomewide-matched-feasibility", default=None)
    parser.add_argument("--pseudobulk-genomewide-de-summary-matched", default=None)
    parser.add_argument("--pseudobulk-genomewide-de-top-hits-matched", default=None)
    parser.add_argument("--pseudobulk-genomewide-de-audit-matched", default=None)
    parser.add_argument("--pseudobulk-genomewide-limma-de-summary", default=None)
    parser.add_argument("--pseudobulk-genomewide-limma-de-top-hits", default=None)
    parser.add_argument("--pseudobulk-genomewide-limma-de-audit", default=None)
    parser.add_argument("--pseudobulk-genomewide-limma-de-summary-matched", default=None)
    parser.add_argument("--pseudobulk-genomewide-limma-de-top-hits-matched", default=None)
    parser.add_argument("--pseudobulk-genomewide-limma-de-audit-matched", default=None)
    parser.add_argument("--ora-sensitivity-scenarios", default=None)
    parser.add_argument("--ora-sensitivity-performance", default=None)
    parser.add_argument("--ora-repeated-cv-summary", default=None)
    parser.add_argument("--ora-repeated-cv-feature-stability", default=None)
    parser.add_argument("--ora-augmented-candidate-repeated-cv-summary", default=None)
    parser.add_argument("--ora-feature-interpretation", default=None)
    parser.add_argument("--ora-feature-set-model-comparison", default=None)
    parser.add_argument("--ora-permutation-empirical", default=None)
    parser.add_argument("--ora-nested-tuning-summary", default=None)
    parser.add_argument("--ora-stacking-summary", default=None)
    parser.add_argument("--ora-model-card", default=None)
    parser.add_argument("--ndd-label-permutation", default=None)
    parser.add_argument("--latent-space-readiness", default=None)
    parser.add_argument("--latent-space-local-audit", default=None)
    parser.add_argument("--latent-space-portal-assets", default=None)
    parser.add_argument("--latent-recompute-feasibility", default=None)
    parser.add_argument("--scvi-pilot-validation", default=None)
    parser.add_argument("--output-provenance", default=None)
    parser.add_argument("--schema", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--figure-dir", default=None)
    parser.add_argument("--top-n", type=int, default=12)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    paths = {
        "manifest": args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv"),
        "cohort_summary": args.cohort_summary or outputs.get("cohort_summary_tsv", "results/tables/cohort_summary.tsv"),
        "associations": args.associations or outputs.get("age_associations_tsv", "results/tables/age_cell_state_associations.tsv"),
        "performance": args.performance or outputs.get("model_performance_tsv", "results/tables/ora_model_performance.tsv"),
        "scores": args.scores or outputs.get("donor_ora_scores_tsv", "results/tables/donor_ora_scores.tsv"),
        "importance": args.importance or outputs.get("feature_importance_tsv", "results/tables/ora_feature_importance.tsv"),
        "ora_calibration": args.ora_calibration or outputs.get("ora_calibration_tsv", "results/tables/ora_calibration.tsv"),
        "ora_age_bin_errors": args.ora_age_bin_errors
        or outputs.get("ora_age_bin_errors_tsv", "results/tables/ora_age_bin_errors.tsv"),
        "ora_residual_diagnostics": args.ora_residual_diagnostics
        or outputs.get("ora_residual_diagnostics_tsv", "results/tables/ora_residual_diagnostics.tsv"),
        "ora_calibrated_scores": args.ora_calibrated_scores
        or outputs.get("ora_calibrated_scores_tsv", "results/tables/ora_calibrated_scores.tsv"),
        "augmented_performance": args.augmented_performance
        or outputs.get("augmented_model_performance_tsv", "results/tables/ora_augmented_model_performance.tsv"),
        "augmented_scores": args.augmented_scores
        or outputs.get("augmented_donor_ora_scores_tsv", "results/tables/augmented_donor_ora_scores.tsv"),
        "augmented_importance": args.augmented_importance
        or outputs.get("augmented_feature_importance_tsv", "results/tables/ora_augmented_feature_importance.tsv"),
        "ndd_projection": args.ndd_projection or outputs.get("ndd_ora_projection_tsv", "results/tables/ndd_ora_projection.tsv"),
        "ndd_projection_summary": args.ndd_projection_summary
        or outputs.get("ndd_ora_projection_summary_tsv", "results/tables/ndd_ora_projection_summary.tsv"),
        "ndd_projection_uncertainty": args.ndd_projection_uncertainty
        or outputs.get("ndd_ora_projection_uncertainty_tsv", "results/tables/ndd_ora_projection_uncertainty.tsv"),
        "ndd_projection_context": args.ndd_projection_context
        or outputs.get("ndd_ora_projection_context_tsv", "results/tables/ndd_ora_projection_context.tsv"),
        "ndd_projection_feature_comparison": args.ndd_projection_feature_comparison
        or outputs.get("ndd_ora_projection_feature_comparison_tsv", "results/tables/ndd_ora_projection_feature_comparison.tsv"),
        "ndd_projection_donor_appendix": args.ndd_projection_donor_appendix
        or outputs.get("ndd_ora_projection_donor_appendix_tsv", "results/tables/ndd_ora_projection_donor_appendix.tsv"),
        "ndd_projection_diagnostics": args.ndd_projection_diagnostics
        or outputs.get("ndd_ora_projection_diagnostics_tsv", "results/tables/ndd_ora_projection_diagnostics.tsv"),
        "module_summary": args.module_summary
        or outputs.get("module_score_summary_tsv", "results/tables/module_score_summary.tsv"),
        "module_coverage": args.module_coverage
        or outputs.get("module_gene_coverage_tsv", "results/tables/module_gene_coverage.tsv"),
        "donor_module_features": args.donor_module_features
        or outputs.get("donor_module_features_tsv", "data/processed/donor_module_features.tsv"),
        "external_validation_summary": args.external_validation_summary
        or outputs.get("external_validation_summary_tsv", "results/tables/external_validation_summary.tsv"),
        "external_gene_list_coverage": args.external_gene_list_coverage
        or outputs.get("external_gene_list_coverage_tsv", "results/tables/external_gene_list_coverage.tsv"),
        "external_feature_contract": args.external_feature_contract
        or outputs.get("external_feature_contract_tsv", "results/tables/external_feature_contract.tsv"),
        "external_sample_metadata": args.external_sample_metadata
        or outputs.get("external_sample_metadata_tsv", "results/tables/external_sample_metadata.tsv"),
        "external_10x_sample_qc": args.external_10x_sample_qc
        or outputs.get("external_10x_sample_qc_tsv", "results/tables/external_10x_sample_qc.tsv"),
        "external_10x_module_contrasts": args.external_10x_module_contrasts
        or outputs.get("external_10x_module_contrasts_tsv", "results/tables/external_10x_module_contrasts.tsv"),
        "external_10x_marker_composition": args.external_10x_marker_composition
        or outputs.get("external_10x_marker_composition_tsv", "results/tables/external_10x_marker_composition.tsv"),
        "external_10x_marker_contrasts": args.external_10x_marker_contrasts
        or outputs.get("external_10x_marker_contrasts_tsv", "results/tables/external_10x_marker_contrasts.tsv"),
        "external_marker_age_concordance": args.external_marker_age_concordance
        or outputs.get("external_marker_age_concordance_tsv", "results/tables/external_marker_age_concordance.tsv"),
        "external_validation_evidence": args.external_validation_evidence
        or outputs.get("external_validation_evidence_tsv", "results/tables/external_validation_evidence.tsv"),
        "pseudobulk_de": args.pseudobulk_de or outputs.get("pseudobulk_de_tsv", "results/tables/pseudobulk_de.tsv"),
        "pseudobulk_coverage": args.pseudobulk_coverage
        or outputs.get("pseudobulk_gene_coverage_tsv", "results/tables/pseudobulk_gene_coverage.tsv"),
        "pseudobulk_metadata": args.pseudobulk_metadata
        or outputs.get("pseudobulk_metadata_tsv", "data/processed/pseudobulk_metadata.tsv"),
        "pseudobulk_covariate_de": args.pseudobulk_covariate_de
        or outputs.get("pseudobulk_covariate_de_tsv", "results/tables/pseudobulk_covariate_de.tsv"),
        "pseudobulk_genomewide_summary": args.pseudobulk_genomewide_summary
        or outputs.get("pseudobulk_genomewide_summary_tsv", "results/tables/pseudobulk_genomewide_summary.tsv"),
        "pseudobulk_genomewide_qc_summary": args.pseudobulk_genomewide_qc_summary
        or outputs.get("pseudobulk_genomewide_qc_summary_tsv", "results/tables/pseudobulk_genomewide_qc_summary.tsv"),
        "pseudobulk_genomewide_gene_qc": args.pseudobulk_genomewide_gene_qc
        or outputs.get("pseudobulk_genomewide_gene_qc_tsv", "results/tables/pseudobulk_genomewide_gene_qc.tsv"),
        "pseudobulk_genomewide_disease_summary": args.pseudobulk_genomewide_disease_summary
        or outputs.get("pseudobulk_genomewide_disease_summary_tsv", "results/tables/pseudobulk_genomewide_disease_summary.tsv"),
        "pseudobulk_genomewide_de_summary": args.pseudobulk_genomewide_de_summary
        or outputs.get("pseudobulk_genomewide_de_summary_tsv", "results/tables/pseudobulk_genomewide_de_summary.tsv"),
        "pseudobulk_genomewide_de_top_hits": args.pseudobulk_genomewide_de_top_hits
        or outputs.get("pseudobulk_genomewide_de_top_hits_tsv", "results/tables/pseudobulk_genomewide_de_top_hits.tsv"),
        "pseudobulk_genomewide_de_audit": args.pseudobulk_genomewide_de_audit
        or outputs.get("pseudobulk_genomewide_de_audit_tsv", "results/tables/pseudobulk_genomewide_de_audit.tsv"),
        "pseudobulk_genomewide_donor_balance": args.pseudobulk_genomewide_donor_balance
        or outputs.get("pseudobulk_genomewide_donor_balance_tsv", "results/tables/pseudobulk_genomewide_donor_balance.tsv"),
        "pseudobulk_genomewide_matched_feasibility": args.pseudobulk_genomewide_matched_feasibility
        or outputs.get(
            "pseudobulk_genomewide_matched_feasibility_tsv",
            "results/tables/pseudobulk_genomewide_matched_feasibility.tsv",
        ),
        "pseudobulk_genomewide_de_summary_matched": args.pseudobulk_genomewide_de_summary_matched
        or outputs.get(
            "pseudobulk_genomewide_de_summary_matched_tsv",
            "results/tables/pseudobulk_genomewide_de_summary_matched_flex_v2_device.tsv",
        ),
        "pseudobulk_genomewide_de_top_hits_matched": args.pseudobulk_genomewide_de_top_hits_matched
        or outputs.get(
            "pseudobulk_genomewide_de_top_hits_matched_tsv",
            "results/tables/pseudobulk_genomewide_de_top_hits_matched_flex_v2_device.tsv",
        ),
        "pseudobulk_genomewide_de_audit_matched": args.pseudobulk_genomewide_de_audit_matched
        or outputs.get(
            "pseudobulk_genomewide_de_audit_matched_tsv",
            "results/tables/pseudobulk_genomewide_de_audit_matched_flex_v2_device.tsv",
        ),
        "pseudobulk_genomewide_limma_de_summary": args.pseudobulk_genomewide_limma_de_summary
        or outputs.get(
            "pseudobulk_genomewide_limma_de_summary_tsv",
            "results/tables/pseudobulk_genomewide_limma_voom_de_summary.tsv",
        ),
        "pseudobulk_genomewide_limma_de_top_hits": args.pseudobulk_genomewide_limma_de_top_hits
        or outputs.get(
            "pseudobulk_genomewide_limma_de_top_hits_tsv",
            "results/tables/pseudobulk_genomewide_limma_voom_de_top_hits.tsv",
        ),
        "pseudobulk_genomewide_limma_de_audit": args.pseudobulk_genomewide_limma_de_audit
        or outputs.get(
            "pseudobulk_genomewide_limma_de_audit_tsv",
            "results/tables/pseudobulk_genomewide_limma_voom_de_audit.tsv",
        ),
        "pseudobulk_genomewide_limma_de_summary_matched": args.pseudobulk_genomewide_limma_de_summary_matched
        or outputs.get(
            "pseudobulk_genomewide_limma_de_summary_matched_tsv",
            "results/tables/pseudobulk_genomewide_limma_voom_de_summary_matched_flex_v2_device.tsv",
        ),
        "pseudobulk_genomewide_limma_de_top_hits_matched": args.pseudobulk_genomewide_limma_de_top_hits_matched
        or outputs.get(
            "pseudobulk_genomewide_limma_de_top_hits_matched_tsv",
            "results/tables/pseudobulk_genomewide_limma_voom_de_top_hits_matched_flex_v2_device.tsv",
        ),
        "pseudobulk_genomewide_limma_de_audit_matched": args.pseudobulk_genomewide_limma_de_audit_matched
        or outputs.get(
            "pseudobulk_genomewide_limma_de_audit_matched_tsv",
            "results/tables/pseudobulk_genomewide_limma_voom_de_audit_matched_flex_v2_device.tsv",
        ),
        "ora_sensitivity_scenarios": args.ora_sensitivity_scenarios
        or outputs.get("ora_sensitivity_scenarios_tsv", "results/tables/ora_sensitivity_scenarios.tsv"),
        "ora_sensitivity_performance": args.ora_sensitivity_performance
        or outputs.get("ora_sensitivity_performance_tsv", "results/tables/ora_sensitivity_performance.tsv"),
        "ora_repeated_cv_summary": args.ora_repeated_cv_summary
        or outputs.get("ora_repeated_cv_summary_tsv", "results/tables/ora_repeated_cv_summary.tsv"),
        "ora_repeated_cv_feature_stability": args.ora_repeated_cv_feature_stability
        or outputs.get("ora_repeated_cv_feature_stability_tsv", "results/tables/ora_repeated_cv_feature_stability.tsv"),
        "ora_augmented_candidate_repeated_cv_summary": args.ora_augmented_candidate_repeated_cv_summary
        or outputs.get(
            "ora_augmented_candidate_repeated_cv_summary_tsv",
            "results/tables/ora_augmented_candidate_repeated_cv_summary.tsv",
        ),
        "ora_feature_interpretation": args.ora_feature_interpretation
        or outputs.get("ora_feature_interpretation_tsv", "results/tables/ora_feature_interpretation.tsv"),
        "ora_feature_set_model_comparison": args.ora_feature_set_model_comparison
        or outputs.get("ora_feature_set_model_comparison_tsv", "results/tables/ora_feature_set_model_comparison.tsv"),
        "ora_permutation_empirical": args.ora_permutation_empirical
        or outputs.get("ora_permutation_empirical_tsv", "results/tables/ora_permutation_empirical.tsv"),
        "ora_nested_tuning_summary": args.ora_nested_tuning_summary
        or outputs.get("ora_nested_tuning_summary_tsv", "results/tables/ora_nested_tuning_summary.tsv"),
        "ora_stacking_summary": args.ora_stacking_summary
        or outputs.get("ora_stacking_summary_tsv", "results/tables/ora_stacking_summary.tsv"),
        "ora_model_card": args.ora_model_card or outputs.get("ora_model_card_tsv", "results/tables/ora_model_card.tsv"),
        "ndd_label_permutation": args.ndd_label_permutation
        or outputs.get("ndd_label_permutation_tsv", "results/tables/ndd_label_permutation.tsv"),
        "latent_space_readiness": args.latent_space_readiness
        or outputs.get("latent_space_readiness_tsv", "results/tables/latent_space_readiness.tsv"),
        "latent_space_local_audit": args.latent_space_local_audit
        or outputs.get("latent_space_local_audit_tsv", "results/tables/latent_space_local_audit.tsv"),
        "latent_space_portal_assets": args.latent_space_portal_assets
        or outputs.get("latent_space_portal_assets_tsv", "results/tables/latent_space_portal_assets.tsv"),
        "latent_recompute_feasibility": args.latent_recompute_feasibility
        or outputs.get("latent_recompute_feasibility_tsv", "results/tables/latent_recompute_feasibility.tsv"),
        "scvi_pilot_validation": args.scvi_pilot_validation
        or outputs.get("scvi_scaled_validation_tsv")
        or outputs.get("scvi_pilot_validation_tsv", "results/tables/scvi_pilot_validation.tsv"),
        "output_provenance": args.output_provenance
        or outputs.get("output_provenance_tsv", "results/reports/output_provenance.tsv"),
        "schema": args.schema or outputs.get("schema_json", "results/reports/h5ad_schema.json"),
        "out": args.out or outputs.get("mvp_report_md", "results/reports/mvp_report.md"),
        "figure_dir": args.figure_dir or outputs.get("figure_dir", "results/figures"),
    }

    written = generate_mvp_report(
        manifest=pd.read_csv(paths["manifest"], sep="\t"),
        cohort_summary=pd.read_csv(paths["cohort_summary"], sep="\t"),
        associations=pd.read_csv(paths["associations"], sep="\t"),
        performance=pd.read_csv(paths["performance"], sep="\t"),
        scores=pd.read_csv(paths["scores"], sep="\t"),
        importance=pd.read_csv(paths["importance"], sep="\t"),
        ora_calibration=_read_optional_tsv(paths["ora_calibration"]),
        ora_age_bin_errors=_read_optional_tsv(paths["ora_age_bin_errors"]),
        ora_residual_diagnostics=_read_optional_tsv(paths["ora_residual_diagnostics"]),
        ora_calibrated_scores=_read_optional_tsv(paths["ora_calibrated_scores"]),
        augmented_performance=_read_optional_tsv(paths["augmented_performance"]),
        augmented_scores=_read_optional_tsv(paths["augmented_scores"]),
        augmented_importance=_read_optional_tsv(paths["augmented_importance"]),
        ndd_projection=_read_optional_tsv(paths["ndd_projection"]),
        ndd_projection_summary=_read_optional_tsv(paths["ndd_projection_summary"]),
        ndd_projection_uncertainty=_read_optional_tsv(paths["ndd_projection_uncertainty"]),
        ndd_projection_context=_read_optional_tsv(paths["ndd_projection_context"]),
        ndd_projection_feature_comparison=_read_optional_tsv(paths["ndd_projection_feature_comparison"]),
        ndd_projection_donor_appendix=_read_optional_tsv(paths["ndd_projection_donor_appendix"]),
        ndd_projection_diagnostics=_read_optional_tsv(paths["ndd_projection_diagnostics"]),
        module_summary=_read_optional_tsv(paths["module_summary"]),
        module_coverage=_read_optional_tsv(paths["module_coverage"]),
        donor_module_features=_read_optional_tsv(paths["donor_module_features"]),
        external_validation_summary=_read_optional_tsv(paths["external_validation_summary"]),
        external_gene_list_coverage=_read_optional_tsv(paths["external_gene_list_coverage"]),
        external_feature_contract=_read_optional_tsv(paths["external_feature_contract"]),
        external_sample_metadata=_read_optional_tsv(paths["external_sample_metadata"]),
        external_10x_sample_qc=_read_optional_tsv(paths["external_10x_sample_qc"]),
        external_10x_module_contrasts=_read_optional_tsv(paths["external_10x_module_contrasts"]),
        external_10x_marker_composition=_read_optional_tsv(paths["external_10x_marker_composition"]),
        external_10x_marker_contrasts=_read_optional_tsv(paths["external_10x_marker_contrasts"]),
        external_marker_age_concordance=_read_optional_tsv(paths["external_marker_age_concordance"]),
        external_validation_evidence=_read_optional_tsv(paths["external_validation_evidence"]),
        pseudobulk_de=_read_optional_tsv(paths["pseudobulk_de"]),
        pseudobulk_coverage=_read_optional_tsv(paths["pseudobulk_coverage"]),
        pseudobulk_metadata=_read_optional_tsv(paths["pseudobulk_metadata"]),
        pseudobulk_covariate_de=_read_optional_tsv(paths["pseudobulk_covariate_de"]),
        pseudobulk_genomewide_summary=_read_optional_tsv(paths["pseudobulk_genomewide_summary"]),
        pseudobulk_genomewide_qc_summary=_read_optional_tsv(paths["pseudobulk_genomewide_qc_summary"]),
        pseudobulk_genomewide_gene_qc=_read_optional_tsv(paths["pseudobulk_genomewide_gene_qc"]),
        pseudobulk_genomewide_disease_summary=_read_optional_tsv(paths["pseudobulk_genomewide_disease_summary"]),
        pseudobulk_genomewide_de_summary=_read_optional_tsv(paths["pseudobulk_genomewide_de_summary"]),
        pseudobulk_genomewide_de_top_hits=_read_optional_tsv(paths["pseudobulk_genomewide_de_top_hits"]),
        pseudobulk_genomewide_de_audit=_read_optional_tsv(paths["pseudobulk_genomewide_de_audit"]),
        pseudobulk_genomewide_donor_balance=_read_optional_tsv(paths["pseudobulk_genomewide_donor_balance"]),
        pseudobulk_genomewide_matched_feasibility=_read_optional_tsv(paths["pseudobulk_genomewide_matched_feasibility"]),
        pseudobulk_genomewide_de_summary_matched=_read_optional_tsv(paths["pseudobulk_genomewide_de_summary_matched"]),
        pseudobulk_genomewide_de_top_hits_matched=_read_optional_tsv(paths["pseudobulk_genomewide_de_top_hits_matched"]),
        pseudobulk_genomewide_de_audit_matched=_read_optional_tsv(paths["pseudobulk_genomewide_de_audit_matched"]),
        pseudobulk_genomewide_limma_de_summary=_read_optional_tsv(paths["pseudobulk_genomewide_limma_de_summary"]),
        pseudobulk_genomewide_limma_de_top_hits=_read_optional_tsv(paths["pseudobulk_genomewide_limma_de_top_hits"]),
        pseudobulk_genomewide_limma_de_audit=_read_optional_tsv(paths["pseudobulk_genomewide_limma_de_audit"]),
        pseudobulk_genomewide_limma_de_summary_matched=_read_optional_tsv(
            paths["pseudobulk_genomewide_limma_de_summary_matched"]
        ),
        pseudobulk_genomewide_limma_de_top_hits_matched=_read_optional_tsv(
            paths["pseudobulk_genomewide_limma_de_top_hits_matched"]
        ),
        pseudobulk_genomewide_limma_de_audit_matched=_read_optional_tsv(
            paths["pseudobulk_genomewide_limma_de_audit_matched"]
        ),
        ora_sensitivity_scenarios=_read_optional_tsv(paths["ora_sensitivity_scenarios"]),
        ora_sensitivity_performance=_read_optional_tsv(paths["ora_sensitivity_performance"]),
        ora_repeated_cv_summary=_read_optional_tsv(paths["ora_repeated_cv_summary"]),
        ora_repeated_cv_feature_stability=_read_optional_tsv(paths["ora_repeated_cv_feature_stability"]),
        ora_augmented_candidate_repeated_cv_summary=_read_optional_tsv(
            paths["ora_augmented_candidate_repeated_cv_summary"]
        ),
        ora_feature_interpretation=_read_optional_tsv(paths["ora_feature_interpretation"]),
        ora_feature_set_model_comparison=_read_optional_tsv(paths["ora_feature_set_model_comparison"]),
        ora_permutation_empirical=_read_optional_tsv(paths["ora_permutation_empirical"]),
        ora_nested_tuning_summary=_read_optional_tsv(paths["ora_nested_tuning_summary"]),
        ora_stacking_summary=_read_optional_tsv(paths["ora_stacking_summary"]),
        ora_model_card=_read_optional_tsv(paths["ora_model_card"]),
        ndd_label_permutation=_read_optional_tsv(paths["ndd_label_permutation"]),
        latent_space_readiness=_read_optional_tsv(paths["latent_space_readiness"]),
        latent_space_local_audit=_read_optional_tsv(paths["latent_space_local_audit"]),
        latent_space_portal_assets=_read_optional_tsv(paths["latent_space_portal_assets"]),
        latent_recompute_feasibility=_read_optional_tsv(paths["latent_recompute_feasibility"]),
        scvi_pilot_validation=_read_optional_tsv(paths["scvi_pilot_validation"]),
        output_provenance=_read_optional_tsv(paths["output_provenance"]),
        schema=load_schema(paths["schema"]),
        source=config.get("source", {}),
        paper_defaults=config.get("paper_defaults", {}),
        out_md=paths["out"],
        figure_dir=paths["figure_dir"],
        top_n=args.top_n,
    )
    print(f"Wrote MVP report: {paths['out']}")
    print(f"Wrote {len(written) - 1} figures: {paths['figure_dir']}")


def _read_optional_tsv(path: str | None) -> pd.DataFrame | None:
    if not path:
        return None
    candidate = Path(path)
    if not candidate.exists():
        return None
    return pd.read_csv(candidate, sep="\t")


if __name__ == "__main__":
    main()
