from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.reporting import generate_mvp_report, rank_associations


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
                out_md=out,
                figure_dir=figures,
                source={"name": "test", "doi": "doi"},
                paper_defaults={"cells_total": 4000},
                schema={"n_obs": 4000, "n_vars": 6, "obsm_keys": ["X_umap"]},
            )

            self.assertTrue(out.exists())
            self.assertIn("Gateway ORA MVP Report", out.read_text(encoding="utf-8"))
            self.assertGreaterEqual(len(written), 6)
            self.assertTrue((figures / "mvp_model_performance.png").exists())


if __name__ == "__main__":
    unittest.main()
