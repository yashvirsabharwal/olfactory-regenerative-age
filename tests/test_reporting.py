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
                pseudobulk_de=pseudobulk_de,
                pseudobulk_coverage=pseudobulk_coverage,
                pseudobulk_metadata=pseudobulk_metadata,
                pseudobulk_covariate_de=pseudobulk_covariate_de,
                pseudobulk_genomewide_summary=pseudobulk_genomewide_summary,
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
            self.assertIn("Pseudobulk Differential Expression", report_text)
            self.assertIn("Covariate-Adjusted Pseudobulk DE", report_text)
            self.assertIn("Genome-Wide Pseudobulk Export", report_text)
            self.assertIn("ad_vs_healthy", report_text)
            self.assertGreaterEqual(len(written), 6)
            self.assertTrue((figures / "mvp_model_performance.png").exists())
            self.assertTrue((figures / "mvp_ndd_projection.png").exists())
            self.assertTrue((figures / "mvp_pseudobulk_de.png").exists())
            self.assertTrue((figures / "mvp_pseudobulk_covariate_de.png").exists())


if __name__ == "__main__":
    unittest.main()
