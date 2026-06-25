from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.negative_controls import run_negative_controls


class NegativeControlTests(unittest.TestCase):
    def test_negative_controls_emit_baselines_and_explainability(self):
        donors = [f"d{i:02d}" for i in range(36)]
        ages = np.linspace(25, 80, len(donors))
        features = pd.DataFrame(
            {
                "donor_id": donors,
                "clr__age_signal": ages / 10.0,
                "prop__weak_noise": np.sin(ages),
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "sample_id": donors,
                "age": ages,
                "sex": ["female" if i % 2 else "male" for i in range(len(donors))],
                "chemistry": ["flex_v1"] * 24 + ["flex_v2"] * 12,
                "collection_method": ["brush", "device"] * 18,
                "site": [""] * len(donors),
                "total_cells": np.linspace(5000, 8000, len(donors)),
                "usable_for_ora_training": [True] * len(donors),
            }
        )
        scores = pd.DataFrame(
            {
                "donor_id": donors,
                "model": ["ridge"] * len(donors),
                "chronological_age": ages,
                "ora": ages + np.sin(ages),
                "oraa": np.sin(ages),
            }
        )
        model_config = {
            "outer_cv_folds": 4,
            "random_seed": 7,
            "age_bins": {"young": [0, 49], "middle": [50, 69], "old": [70, 120]},
        }

        result = run_negative_controls(
            features,
            manifest,
            scores,
            model_config,
            n_shuffles=5,
        )

        self.assertIn("biological_ridge_cv", set(result.performance["control"]))
        self.assertIn("technical_only_ridge_cv", set(result.performance["control"]))
        self.assertIn("age_shuffle_within_technical_strata", set(result.performance["control"]))
        self.assertIn("disease_label_negative_control", set(result.performance["control"]))
        self.assertIn("technical_only_ridge_cv", set(result.baseline_comparison["comparison"]))
        self.assertIn("total_cell_yield", set(result.covariate_explainability["covariate_set"]))


if __name__ == "__main__":
    unittest.main()
