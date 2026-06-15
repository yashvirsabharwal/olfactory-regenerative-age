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
    parser.add_argument("--augmented-performance", default=None)
    parser.add_argument("--augmented-scores", default=None)
    parser.add_argument("--augmented-importance", default=None)
    parser.add_argument("--ndd-projection", default=None)
    parser.add_argument("--ndd-projection-summary", default=None)
    parser.add_argument("--ndd-projection-uncertainty", default=None)
    parser.add_argument("--ndd-projection-context", default=None)
    parser.add_argument("--module-summary", default=None)
    parser.add_argument("--module-coverage", default=None)
    parser.add_argument("--donor-module-features", default=None)
    parser.add_argument("--external-validation-summary", default=None)
    parser.add_argument("--external-gene-list-coverage", default=None)
    parser.add_argument("--external-feature-contract", default=None)
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
    parser.add_argument("--ora-sensitivity-scenarios", default=None)
    parser.add_argument("--ora-sensitivity-performance", default=None)
    parser.add_argument("--ora-repeated-cv-summary", default=None)
    parser.add_argument("--ora-repeated-cv-feature-stability", default=None)
    parser.add_argument("--ora-feature-set-model-comparison", default=None)
    parser.add_argument("--ora-permutation-empirical", default=None)
    parser.add_argument("--ora-nested-tuning-summary", default=None)
    parser.add_argument("--ora-stacking-summary", default=None)
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
        "ora_sensitivity_scenarios": args.ora_sensitivity_scenarios
        or outputs.get("ora_sensitivity_scenarios_tsv", "results/tables/ora_sensitivity_scenarios.tsv"),
        "ora_sensitivity_performance": args.ora_sensitivity_performance
        or outputs.get("ora_sensitivity_performance_tsv", "results/tables/ora_sensitivity_performance.tsv"),
        "ora_repeated_cv_summary": args.ora_repeated_cv_summary
        or outputs.get("ora_repeated_cv_summary_tsv", "results/tables/ora_repeated_cv_summary.tsv"),
        "ora_repeated_cv_feature_stability": args.ora_repeated_cv_feature_stability
        or outputs.get("ora_repeated_cv_feature_stability_tsv", "results/tables/ora_repeated_cv_feature_stability.tsv"),
        "ora_feature_set_model_comparison": args.ora_feature_set_model_comparison
        or outputs.get("ora_feature_set_model_comparison_tsv", "results/tables/ora_feature_set_model_comparison.tsv"),
        "ora_permutation_empirical": args.ora_permutation_empirical
        or outputs.get("ora_permutation_empirical_tsv", "results/tables/ora_permutation_empirical.tsv"),
        "ora_nested_tuning_summary": args.ora_nested_tuning_summary
        or outputs.get("ora_nested_tuning_summary_tsv", "results/tables/ora_nested_tuning_summary.tsv"),
        "ora_stacking_summary": args.ora_stacking_summary
        or outputs.get("ora_stacking_summary_tsv", "results/tables/ora_stacking_summary.tsv"),
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
        augmented_performance=_read_optional_tsv(paths["augmented_performance"]),
        augmented_scores=_read_optional_tsv(paths["augmented_scores"]),
        augmented_importance=_read_optional_tsv(paths["augmented_importance"]),
        ndd_projection=_read_optional_tsv(paths["ndd_projection"]),
        ndd_projection_summary=_read_optional_tsv(paths["ndd_projection_summary"]),
        ndd_projection_uncertainty=_read_optional_tsv(paths["ndd_projection_uncertainty"]),
        ndd_projection_context=_read_optional_tsv(paths["ndd_projection_context"]),
        module_summary=_read_optional_tsv(paths["module_summary"]),
        module_coverage=_read_optional_tsv(paths["module_coverage"]),
        donor_module_features=_read_optional_tsv(paths["donor_module_features"]),
        external_validation_summary=_read_optional_tsv(paths["external_validation_summary"]),
        external_gene_list_coverage=_read_optional_tsv(paths["external_gene_list_coverage"]),
        external_feature_contract=_read_optional_tsv(paths["external_feature_contract"]),
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
        ora_sensitivity_scenarios=_read_optional_tsv(paths["ora_sensitivity_scenarios"]),
        ora_sensitivity_performance=_read_optional_tsv(paths["ora_sensitivity_performance"]),
        ora_repeated_cv_summary=_read_optional_tsv(paths["ora_repeated_cv_summary"]),
        ora_repeated_cv_feature_stability=_read_optional_tsv(paths["ora_repeated_cv_feature_stability"]),
        ora_feature_set_model_comparison=_read_optional_tsv(paths["ora_feature_set_model_comparison"]),
        ora_permutation_empirical=_read_optional_tsv(paths["ora_permutation_empirical"]),
        ora_nested_tuning_summary=_read_optional_tsv(paths["ora_nested_tuning_summary"]),
        ora_stacking_summary=_read_optional_tsv(paths["ora_stacking_summary"]),
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
