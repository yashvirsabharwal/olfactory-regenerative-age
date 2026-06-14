from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.ndd import summarize_ndd_projection_uncertainty


class NDDUncertaintyTests(unittest.TestCase):
    def test_uncertainty_uses_matched_healthy_reference(self):
        projection = pd.DataFrame(
            {
                "donor_id": ["h1", "h2", "h3", "a1", "a2"] * 2,
                "model": ["random_forest"] * 5 + ["elastic_net"] * 5,
                "disease_group": ["healthy", "healthy", "healthy", "ad", "ad"] * 2,
                "chemistry": ["v1", "v2", "v2", "v2", "v2"] * 2,
                "collection_method": ["brush", "device", "device", "device", "device"] * 2,
                "chronological_age": [50, 60, 70, 75, 80] * 2,
                "total_cells": [100, 200, 300, 400, 500] * 2,
                "oraa": [1.0, 2.0, 3.0, -5.0, -7.0, 1.0, 2.0, 3.0, -4.0, -6.0],
            }
        )

        result = summarize_ndd_projection_uncertainty(
            projection,
            n_bootstrap=100,
            random_seed=1,
        )
        rf_matched = result.uncertainty[
            result.uncertainty["model"].eq("random_forest")
            & result.uncertainty["disease_group"].eq("ad")
            & result.uncertainty["reference"].eq("matched_healthy")
        ].iloc[0]

        self.assertEqual(int(rf_matched["n_disease"]), 2)
        self.assertEqual(int(rf_matched["n_reference"]), 2)
        self.assertLess(rf_matched["difference_vs_reference"], 0)
        self.assertIn("collection_method", result.context.columns)


if __name__ == "__main__":
    unittest.main()
