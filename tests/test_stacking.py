from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.stacking import run_stacked_ora


class StackingTests(unittest.TestCase):
    def test_run_stacked_ora_returns_predictions_and_meta_weights(self):
        donors = [f"d{i}" for i in range(12)]
        features = pd.DataFrame(
            {
                "donor_id": donors,
                "prop__age_signal": np.linspace(0, 1, 12),
                "clr__age_signal": np.linspace(1, 0, 12),
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "sample_id": [f"s{i}" for i in range(12)],
                "age": np.linspace(35, 80, 12),
                "sex": ["F", "M"] * 6,
                "race_ethnicity": ["reported"] * 12,
                "disease_group": ["healthy"] * 12,
                "chemistry": ["v2"] * 12,
                "collection_method": ["device"] * 12,
                "site": ["site1"] * 12,
                "total_cells": np.arange(12) + 100,
                "usable_for_ora_training": [True] * 12,
            }
        )

        result = run_stacked_ora(
            features,
            manifest,
            {"outer_cv_folds": 3, "random_seed": 1},
            base_models=["ridge", "null_model"],
            repeats=2,
            inner_folds=2,
        )

        self.assertEqual(set(result.performance["model"]), {"stacked_ensemble"})
        self.assertEqual(result.performance["repeat"].nunique(), 2)
        self.assertEqual(set(result.predictions["donor_id"]), set(donors))
        self.assertTrue(result.predictions["oraa"].notna().all())
        self.assertTrue({"ridge", "null_model", "intercept"}.issubset(set(result.meta_weights["base_model"])))
        self.assertIn("mae_mean", result.performance_summary.columns)

    def test_run_stacked_ora_requires_base_models(self):
        with self.assertRaises(ValueError):
            run_stacked_ora(pd.DataFrame(), pd.DataFrame(), {}, base_models=[])


if __name__ == "__main__":
    unittest.main()
