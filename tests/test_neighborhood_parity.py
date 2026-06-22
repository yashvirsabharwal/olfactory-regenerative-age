import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ora.neighborhood_parity import export_neighborhood_count_inputs, summarize_edger_parity


class NeighborhoodParityTests(unittest.TestCase):
    def test_export_neighborhood_count_inputs_filters_and_scales_donors(self):
        memberships = pd.DataFrame(
            {
                "neighborhood_id": [0, 0, 0, 1, 1, 2],
                "donor_id": ["d1", "d1", "d2", "d2", "d3", "blocked"],
            }
        )
        donors = pd.DataFrame(
            {
                "donor_id": ["d1", "d2", "d3"],
                "age": [30, 50, 70],
                "sex": ["female", "male", "female"],
                "total_cells": [100, 200, 300],
            }
        )

        counts, design, summary = export_neighborhood_count_inputs(memberships, donors)

        self.assertEqual(counts.columns.tolist(), ["neighborhood_id", "d1", "d2", "d3"])
        self.assertEqual(counts.loc[counts["neighborhood_id"].eq(0), "d1"].iloc[0], 2)
        self.assertEqual(counts.loc[counts["neighborhood_id"].eq(1), "d3"].iloc[0], 1)
        self.assertIn("age_scaled", design.columns)
        self.assertAlmostEqual(float(design["age_scaled"].mean()), 0.0)
        self.assertEqual(int(summary.loc[summary["metric"].eq("donors"), "value"].iloc[0]), 3)

    def test_summarize_edger_parity_reports_overlap_and_direction(self):
        python_da = pd.DataFrame(
            {
                "neighborhood_id": [0, 1, 2],
                "age_coef": [-1.0, 0.5, -0.2],
                "age_pvalue": [0.001, 0.02, 0.8],
                "age_fdr": [0.01, 0.08, 0.9],
                "top_fine_celltype": ["A", "B", "C"],
                "top_coarse_celltype": ["X", "Y", "Z"],
                "status": ["tested", "tested", "tested"],
            }
        )
        edger_da = pd.DataFrame(
            {
                "neighborhood_id": [0, 1, 2],
                "logFC": [-2.0, -0.4, -0.1],
                "logCPM": [1.0, 1.0, 1.0],
                "F": [10, 5, 1],
                "PValue": [0.001, 0.04, 0.9],
                "FDR": [0.01, 0.12, 0.95],
            }
        )

        comparison, summary = summarize_edger_parity(python_da, edger_da, run_name="toy", top_n=2)

        self.assertEqual(comparison.shape[0], 3)
        self.assertEqual(int(summary.loc[summary["metric"].eq("significant_overlap"), "value"].iloc[0]), 1)
        self.assertEqual(float(summary.loc[summary["metric"].eq("python_sig_direction_agreement"), "value"].iloc[0]), 0.5)


if __name__ == "__main__":
    unittest.main()
