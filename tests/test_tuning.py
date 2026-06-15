from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.tuning import candidate_params_for_model, run_nested_tuning


class NestedTuningTests(unittest.TestCase):
    def test_candidate_params_include_base_and_deduplicate(self):
        candidates = candidate_params_for_model(
            "xgboost",
            {"models": {"xgboost": {"enabled": True, "max_depth": 2}}},
            {"search_spaces": {"xgboost": [{"max_depth": 2}, {"max_depth": 3}]}},
        )

        self.assertEqual(candidates[0], {"max_depth": 2})
        self.assertEqual(candidates[-1], {"max_depth": 3})
        self.assertEqual(len(candidates), 2)

    def test_run_nested_tuning_returns_held_out_predictions_and_trace(self):
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

        result = run_nested_tuning(
            features,
            manifest,
            {"outer_cv_folds": 3, "random_seed": 1, "model_names": ["ridge"]},
            repeats=1,
            inner_folds=2,
            max_candidates=1,
        )

        self.assertEqual(set(result.performance["model"]), {"ridge"})
        self.assertEqual(set(result.predictions["donor_id"]), set(donors))
        self.assertFalse(result.tuning_trace.empty)
        self.assertFalse(result.selected_params.empty)
        self.assertIn("mae_mean", result.performance_summary.columns)


if __name__ == "__main__":
    unittest.main()
