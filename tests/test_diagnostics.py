from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.diagnostics import summarize_ora_diagnostics


class DiagnosticsTests(unittest.TestCase):
    def test_summarize_ora_diagnostics_reports_calibration_and_residual_strata(self):
        donors = [f"d{i}" for i in range(8)]
        ages = np.array([30, 40, 50, 60, 70, 80, 55, 65], dtype=float)
        scores = pd.DataFrame(
            {
                "donor_id": donors * 2,
                "model": ["ridge"] * 8 + ["random_forest"] * 8,
                "chronological_age": np.concatenate([ages, ages]),
                "ora": np.concatenate([10 + 0.8 * ages, 5 + 0.9 * ages]),
                "oraa": np.concatenate([np.linspace(-1, 1, 8), np.linspace(1, -1, 8)]),
                "sex": ["female", "male"] * 8,
                "chemistry": ["flex_v1"] * 8 + ["flex_v2"] * 8,
                "collection_method": ["brush", "device"] * 8,
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "sample_id": [f"s{i}" for i in range(8)],
                "total_cells": np.linspace(1000, 8000, 8),
                "race_ethnicity": ["reported"] * 8,
            }
        )

        result = summarize_ora_diagnostics(
            scores,
            model_config={"age_bins": {"young": [0, 49], "middle": [50, 69], "old": [70, 120]}},
            manifest=manifest,
        )

        self.assertEqual(set(result.calibration["model"]), {"ridge", "random_forest"})
        ridge = result.calibration[result.calibration["model"].eq("ridge")].iloc[0]
        self.assertAlmostEqual(ridge["calibration_slope_ora_on_age"], 0.8)
        self.assertLess(ridge["calibrated_mae"], ridge["mae"])
        self.assertEqual(set(result.age_bin_errors["level"]), {"young", "middle", "old"})
        self.assertIn("total_cells_bin", set(result.residual_diagnostics["group"]))
        self.assertTrue({"calibrated_ora", "calibrated_error"}.issubset(result.calibrated_scores.columns))


if __name__ == "__main__":
    unittest.main()
