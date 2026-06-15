from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.reporting import generate_mvp_report, rank_associations, rank_pseudobulk_covariate_de, rank_pseudobulk_de


class ReportingTests(unittest.TestCase):
    def test_rank_associations_sorts_by_fdr_then_p_value(self):
        associations = pd.DataFrame(
            {
                "feature": ["b", "a", "c", "bad"],
                "status": ["ok", "ok", "ok", "too_few_samples"],
                "p_value": [0.03, 0.01, 0.02, None],
                "fdr": [0.05, 0.05, 0.04, None],
                "beta_per_10_years": [1.0, 2.0, -1.0, None],
            }
        )

        ranked = rank_associations(associations, top_n=3)

        self.assertEqual(ranked["feature"].tolist(), ["c", "a", "b"])

    def test_rank_pseudobulk_de_sorts_ok_rows(self):
        de = pd.DataFrame(
            {
                "contrast": ["ad_vs_healthy", "pd_vs_healthy", "ad_vs_healthy", "ad_vs_healthy"],
                "fine_cell_type": ["qHBC", "mOSN", "qHBC", "qHBC"],
                "gene": ["TP63", "SNCA", "OMP", "BAD"],
                "n_case": [5, 5, 5, 1],
                "n_control": [10, 10, 10, 10],
                "log2fc": [1.0, -2.0, 0.5, 0.0],
                "p_value": [0.01, 0.001, 0.02, None],
                "fdr": [0.03, 0.01, 0.02, None],
                "status": ["ok", "ok", "ok", "too_few_donors"],
            }
        )

        ranked = rank_pseudobulk_de(de, top_n=3)

        self.assertEqual(ranked["gene"].tolist(), ["SNCA", "OMP", "TP63"])

    def test_rank_pseudobulk_covariate_de_sorts_ok_rows(self):
        de = pd.DataFrame(
            {
                "contrast": ["ad_vs_healthy", "pd_vs_healthy", "ad_vs_healthy", "ad_vs_healthy"],
                "fine_cell_type": ["qHBC", "mOSN", "qHBC", "qHBC"],
                "gene": ["TP63", "SNCA", "OMP", "BAD"],
                "n_case": [5, 5, 5, 1],
                "n_control": [10, 10, 10, 10],
                "log2fc_adjusted": [1.0, -2.0, 0.5, 0.0],
                "p_value": [0.01, 0.001, 0.02, None],
                "fdr": [0.03, 0.01, 0.02, None],
                "status": ["ok", "ok", "ok", "too_few_donors"],
            }
        )

        ranked = rank_pseudobulk_covariate_de(de, top_n=3)

        self.assertEqual(ranked["gene"].tolist(), ["SNCA", "OMP", "TP63"])

    def test_generate_mvp_report_writes_markdown_and_figures(self):
        manifest = pd.DataFrame(
            {
                "donor_id": ["d1", "d2", "d3", "d4"],
                "sample_id": ["s1", "s2", "s3", "s4"],
                "age": [40, 50, 60, 70],
                "usable_for_ora_training": [True, True, True, True],
            }
        )
        cohort_summary = pd.DataFrame(
            {
                "cohort": ["all", "healthy"],
                "donors": [4, 4],
                "samples": [4, 4],
                "cells": [4000, 4000],
                "median_age": [55, 55],
                "missing_age_samples": [0, 0],
                "lineage_cells": [1000, 1000],
                "mature_neurons": [100, 100],
            }
        )
        associations = pd.DataFrame(
            {
                "feature": ["prop__hbc", "clr__mosn"],
                "n": [4, 4],
                "beta_per_10_years": [0.1, -0.2],
                "p_value": [0.01, 0.02],
                "fdr": [0.02, 0.03],
                "direction": ["positive", "negative"],
                "status": ["ok", "ok"],
            }
        )
        performance = pd.DataFrame(
            {
                "model": ["null_model", "elastic_net", "random_forest"],
                "n": [4, 4, 4],
                "mae": [12.0, 8.0, 9.0],
                "rmse": [14.0, 10.0, 11.0],
                "r2": [0.0, 0.2, 0.1],
                "spearman_r": [0.0, 0.7, 0.6],
            }
        )
        scores = pd.DataFrame(
            {
                "donor_id": ["d1", "d2", "d3", "d4"] * 2,
                "model": ["elastic_net"] * 4 + ["random_forest"] * 4,
                "chronological_age": [40, 50, 60, 70] * 2,
                "ora": [42, 49, 61, 69, 41, 51, 59, 71],
                "oraa": [1, -1, 1, -1, 1, 1, -1, -1],
            }
        )
        importance = pd.DataFrame(
            {
                "model": ["elastic_net", "elastic_net", "random_forest", "random_forest"],
                "feature": ["prop__hbc", "clr__mosn", "prop__hbc", "clr__mosn"],
                "importance": [0.3, -0.5, 0.2, 0.1],
                "stability": [1.0, 0.8, 1.0, 1.0],
            }
        )
        pseudobulk_de = pd.DataFrame(
            {
                "contrast": ["ad_vs_healthy", "pd_vs_healthy"],
                "case_group": ["ad", "pd"],
                "control_group": ["healthy", "healthy"],
                "fine_cell_type": ["qHBC", "mOSN"],
                "gene": ["TP63", "SNCA"],
                "n_case": [5, 5],
                "n_control": [4, 4],
                "mean_logcpm_case": [4.0, 2.0],
                "mean_logcpm_control": [2.0, 4.0],
                "log2fc": [2.0, -2.0],
                "t_stat": [3.0, -3.0],
                "p_value": [0.01, 0.02],
                "status": ["ok", "ok"],
                "fdr": [0.02, 0.03],
            }
        )
        pseudobulk_coverage = pd.DataFrame(
            {
                "module": ["targeted_pseudobulk"],
                "n_requested": [2],
                "n_present": [2],
                "coverage_fraction": [1.0],
                "missing_genes": [""],
            }
        )
        pseudobulk_metadata = pd.DataFrame(
            {
                "donor_id": ["d1", "d2"],
                "sample_id": ["s1", "s2"],
                "disease_group": ["healthy", "ad"],
                "coarse_cell_type": ["HBC", "HBC"],
                "fine_cell_type": ["qHBC", "qHBC"],
                "n_cells": [100, 80],
                "sum_n_counts": [1000, 800],
            }
        )
        pseudobulk_covariate_de = pd.DataFrame(
            {
                "contrast": ["ad_vs_healthy", "pd_vs_healthy"],
                "case_group": ["ad", "pd"],
                "control_group": ["healthy", "healthy"],
                "fine_cell_type": ["qHBC", "mOSN"],
                "gene": ["TP63", "SNCA"],
                "n_case": [5, 5],
                "n_control": [4, 4],
                "n_total": [9, 9],
                "mean_logcpm_case": [4.0, 2.0],
                "mean_logcpm_control": [2.0, 4.0],
                "log2fc_unadjusted": [2.0, -2.0],
                "log2fc_adjusted": [1.5, -1.5],
                "t_stat": [3.0, -3.0],
                "p_value": [0.01, 0.02],
                "df_resid": [5, 5],
                "covariates": ["age,sex", "age,sex"],
                "status": ["ok", "ok"],
                "fdr": [0.02, 0.03],
            }
        )
        pseudobulk_genomewide_summary = pd.DataFrame(
            {
                "n_cells": [4000],
                "n_genes": [6],
                "n_groups_total": [4],
                "n_groups_exported": [2],
                "n_groups_failed_min_cells": [1],
                "n_groups_failed_min_donors": [1],
                "min_cells_per_group": [10],
                "min_donors_per_cell_state": [3],
            }
        )
        pseudobulk_genomewide_qc_summary = pd.DataFrame(
            {
                "n_genes": [6],
                "n_groups": [2],
                "metadata_rows": [2],
                "matrix_columns_match_metadata": [True],
                "matrix_total_counts": [1000],
                "metadata_total_counts": [1000],
                "median_group_detected_genes": [5],
                "median_gene_detected_group_fraction": [0.8],
            }
        )
        pseudobulk_genomewide_gene_qc = pd.DataFrame(
            {
                "gene_symbol": ["BPIFA1", "OMP"],
                "total_count": [1000, 500],
                "detected_group_fraction": [1.0, 0.5],
                "variance_log1p": [8.0, 2.0],
            }
        )
        pseudobulk_genomewide_disease_summary = pd.DataFrame(
            {
                "disease_group": ["healthy", "ad"],
                "groups": [10, 2],
                "donors": [4, 1],
                "cells": [3000, 1000],
                "matrix_total_count": [700, 300],
                "median_detected_genes": [5, 6],
            }
        )
        pseudobulk_genomewide_de_summary = pd.DataFrame(
            {
                "contrast": ["ad_vs_healthy"],
                "tested_rows": [100],
                "tested_genes": [50],
                "tested_cell_states": [2],
                "ok_cell_state_models": [2],
                "fdr_threshold": [0.05],
                "significant_rows": [3],
                "significant_genes": [2],
                "significant_cell_states": [1],
                "sex_linked_significant_rows": [1],
            }
        )
        pseudobulk_genomewide_de_top_hits = pd.DataFrame(
            {
                "contrast": ["ad_vs_healthy", "ad_vs_healthy"],
                "fine_cell_type": ["qHBC", "mOSN"],
                "gene_symbol": ["USP9Y", "MAFF"],
                "log2fc": [2.0, 1.0],
                "p_value": [0.001, 0.002],
                "fdr": [0.01, 0.02],
                "is_sex_linked_initial": [True, False],
            }
        )
        ora_sensitivity_scenarios = pd.DataFrame(
            {
                "scenario": ["baseline", "collection_method__brush"],
                "status": ["ok", "ok"],
            }
        )
        ora_sensitivity_performance = pd.DataFrame(
            {
                "scenario": ["baseline", "collection_method__brush"],
                "model": ["random_forest", "random_forest"],
                "n": [4, 4],
                "mae": [9.0, 12.0],
                "rmse": [11.0, 14.0],
                "r2": [0.1, -0.1],
                "spearman_r": [0.6, 0.2],
                "healthy_train_donors": [4, 4],
            }
        )
        ora_repeated_cv_summary = pd.DataFrame(
            {
                "model": ["random_forest"],
                "repeats": [2],
                "n": [4],
                "mae_mean": [9.5],
                "mae_ci_low": [9.0],
                "mae_ci_high": [10.0],
                "spearman_r_mean": [0.5],
                "spearman_r_ci_low": [0.4],
                "spearman_r_ci_high": [0.6],
            }
        )
        ora_repeated_cv_feature_stability = pd.DataFrame(
            {
                "model": ["elastic_net"],
                "feature": ["clr__cdc1"],
                "mean_importance": [-1.2],
                "selection_fraction": [1.0],
            }
        )
        ora_permutation_empirical = pd.DataFrame(
            {
                "model": ["random_forest"],
                "n_permutations": [10],
                "observed_mae": [9.5],
                "null_mae_mean": [14.0],
                "empirical_p_mae": [0.09],
                "observed_spearman_r": [0.5],
                "null_spearman_r_mean": [0.02],
                "empirical_p_spearman_r": [0.09],
            }
        )
        ora_calibration = pd.DataFrame(
            {
                "model": ["ridge", "random_forest"],
                "n": [4, 4],
                "calibration_slope_ora_on_age": [0.8, 0.9],
                "calibration_intercept_ora_on_age": [10.0, 5.0],
                "mae": [9.0, 8.0],
                "calibrated_mae": [7.0, 6.5],
                "spearman_r": [0.5, 0.6],
            }
        )
        ora_age_bin_errors = pd.DataFrame(
            {
                "model": ["ridge", "ridge", "ridge"],
                "group": ["age_bin", "age_bin", "age_bin"],
                "level": ["young", "middle", "old"],
                "n": [1, 2, 1],
                "mean_error": [2.0, -1.0, 3.0],
                "mae": [2.0, 1.5, 3.0],
                "calibrated_mean_error": [0.5, -0.2, 0.8],
                "calibrated_mae": [0.5, 0.8, 0.8],
            }
        )
        ora_residual_diagnostics = pd.DataFrame(
            {
                "model": ["ridge", "random_forest"],
                "group": ["chemistry", "collection_method"],
                "level": ["flex_v2", "device"],
                "n": [6, 6],
                "mean_error": [4.0, -3.0],
                "mae": [4.5, 3.5],
                "mean_oraa": [0.2, -0.3],
                "calibrated_mean_error": [1.0, -1.2],
            }
        )
        ndd_projection = pd.DataFrame(
            {
                "donor_id": ["d1", "d2", "d3", "d4"] * 2,
                "model": ["elastic_net"] * 4 + ["random_forest"] * 4,
                "disease_group": ["healthy", "healthy", "ad", "pd"] * 2,
                "chronological_age": [40, 50, 70, 75] * 2,
                "ora": [41, 49, 80, 82, 42, 48, 78, 81],
                "oraa": [0.1, -0.2, 5.0, 4.0, 0.2, -0.1, 4.5, 3.5],
                "is_training_donor": [True, True, False, False] * 2,
            }
        )
        ndd_projection_summary = pd.DataFrame(
            {
                "model": ["elastic_net", "elastic_net"],
                "disease_group": ["healthy", "ad"],
                "donors": [2, 1],
                "training_donors": [2, 0],
                "ndd_donors": [0, 1],
                "mean_age": [45, 70],
                "mean_ora": [45, 80],
                "mean_oraa": [0, 5],
                "sd_oraa": [0.2, 0],
            }
        )
        ndd_projection_uncertainty = pd.DataFrame(
            {
                "model": ["random_forest"],
                "disease_group": ["ad"],
                "reference": ["matched_healthy"],
                "n_disease": [1],
                "n_reference": [2],
                "mean_oraa": [4.5],
                "mean_oraa_ci_low": [4.0],
                "mean_oraa_ci_high": [5.0],
                "difference_vs_reference": [4.0],
                "difference_ci_low": [3.0],
                "difference_ci_high": [5.0],
            }
        )
        ndd_projection_context = pd.DataFrame(
            {
                "disease_group": ["ad", "healthy"],
                "chemistry": ["flex_v2", "flex_v2"],
                "collection_method": ["device", "device"],
                "donors": [1, 2],
                "mean_age": [70, 45],
                "median_total_cells": [1000, 800],
            }
        )
        external_validation_summary = pd.DataFrame(
            {
                "dataset_id": ["toy_external"],
                "status": ["metadata_pending"],
                "validation_use": ["aging replication"],
                "expected_level": ["donor_feature_matrix"],
                "ready_for_feature_validation": [False],
                "ready_for_raw_adapter": [False],
                "files_missing": ["feature_matrix,expression,metadata"],
            }
        )
        external_gene_list_coverage = pd.DataFrame(
            {
                "gene_list": ["aging"],
                "n_requested": [2],
                "n_present": [2],
                "coverage_fraction": [1.0],
                "missing_genes": [""],
            }
        )
        external_feature_contract = pd.DataFrame(
            {
                "field": ["donor_id", "age", "module_score__"],
                "kind": ["required_column", "required_column", "accepted_feature_prefix"],
            }
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "reports" / "mvp.md"
            figures = Path(tmpdir) / "figures"
            written = generate_mvp_report(
                manifest=manifest,
                cohort_summary=cohort_summary,
                associations=associations,
                performance=performance,
                scores=scores,
                importance=importance,
                ndd_projection=ndd_projection,
                ndd_projection_summary=ndd_projection_summary,
                ndd_projection_uncertainty=ndd_projection_uncertainty,
                ndd_projection_context=ndd_projection_context,
                ora_calibration=ora_calibration,
                ora_age_bin_errors=ora_age_bin_errors,
                ora_residual_diagnostics=ora_residual_diagnostics,
                external_validation_summary=external_validation_summary,
                external_gene_list_coverage=external_gene_list_coverage,
                external_feature_contract=external_feature_contract,
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
                ora_sensitivity_scenarios=ora_sensitivity_scenarios,
                ora_sensitivity_performance=ora_sensitivity_performance,
                ora_repeated_cv_summary=ora_repeated_cv_summary,
                ora_repeated_cv_feature_stability=ora_repeated_cv_feature_stability,
                ora_permutation_empirical=ora_permutation_empirical,
                out_md=out,
                figure_dir=figures,
                source={"name": "test", "doi": "doi"},
                paper_defaults={"cells_total": 4000},
                schema={"n_obs": 4000, "n_vars": 6, "obsm_keys": ["X_umap"]},
            )

            self.assertTrue(out.exists())
            report_text = out.read_text(encoding="utf-8")
            self.assertIn("Gateway ORA MVP Report", report_text)
            self.assertIn("NDD ORA Projection", report_text)
            self.assertIn("matched_healthy", report_text)
            self.assertIn("flex_v2", report_text)
            self.assertIn("Pseudobulk Differential Expression", report_text)
            self.assertIn("Covariate-Adjusted Pseudobulk DE", report_text)
            self.assertIn("Genome-Wide Pseudobulk Export", report_text)
            self.assertIn("Genome-Wide edgeR DE", report_text)
            self.assertIn("sex-linked sentinel", report_text)
            self.assertIn("Top non-sex-linked sentinel hits", report_text)
            self.assertIn("USP9Y", report_text)
            self.assertIn("MAFF", report_text)
            self.assertIn("BPIFA1", report_text)
            self.assertIn("Matrix total counts", report_text)
            self.assertIn("ORA Sensitivity", report_text)
            self.assertIn("collection_method__brush", report_text)
            self.assertIn("Repeated-CV ORA Stability", report_text)
            self.assertIn("Shuffled-Age Null Test", report_text)
            self.assertIn("clr__cdc1", report_text)
            self.assertIn("ORA Calibration Diagnostics", report_text)
            self.assertIn("calibrated MAE", report_text)
            self.assertIn("External Validation Readiness", report_text)
            self.assertIn("toy_external", report_text)
            self.assertIn("module_score__", report_text)
            self.assertIn("ad_vs_healthy", report_text)
            self.assertGreaterEqual(len(written), 6)
            self.assertTrue((figures / "mvp_model_performance.png").exists())
            self.assertTrue((figures / "mvp_ndd_projection.png").exists())
            self.assertTrue((figures / "mvp_pseudobulk_de.png").exists())
            self.assertTrue((figures / "mvp_pseudobulk_covariate_de.png").exists())


if __name__ == "__main__":
    unittest.main()
