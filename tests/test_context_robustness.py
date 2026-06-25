from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.context_robustness import build_context_splits, run_leave_context_out


class ContextRobustnessTests(unittest.TestCase):
    def test_build_context_splits_records_single_level_and_small_contexts(self):
        manifest = pd.DataFrame(
            {
                "donor_id": [f"d{i}" for i in range(8)],
                "age": np.linspace(30, 70, 8),
                "usable_for_ora_training": [True] * 8,
                "site": ["missing"] * 8,
                "chemistry": ["v1"] * 6 + ["v2"] * 2,
            }
        )

        specs = build_context_splits(manifest, contexts=["site", "chemistry"], min_train_donors=4, min_test_donors=3)
        statuses = {(spec.context, spec.level): spec.status for spec in specs}

        self.assertEqual(statuses[("site", "missing")], "skipped_single_level")
        self.assertEqual(statuses[("chemistry", "v2")], "too_few_test_donors")
        self.assertEqual(statuses[("chemistry", "v1")], "too_few_train_donors")

    def test_run_leave_context_out_returns_held_out_metrics(self):
        donors = [f"d{i}" for i in range(12)]
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "sample_id": donors,
                "age": np.linspace(30, 80, 12),
                "usable_for_ora_training": [True] * 12,
                "disease_group": ["healthy"] * 12,
                "sex": ["F"] * 6 + ["M"] * 6,
                "chemistry": ["v2"] * 12,
                "collection_method": ["device"] * 12,
                "site": ["site"] * 12,
                "total_cells": np.arange(12) + 100,
            }
        )
        features = pd.DataFrame(
            {
                "donor_id": donors,
                "prop__signal": np.linspace(0.0, 1.0, 12),
                "scvi_global_mean__dim01": np.linspace(1.0, 0.0, 12),
            }
        )

        result = run_leave_context_out(
            features,
            manifest,
            {"outer_cv_folds": 2, "random_seed": 1, "model_names": ["ridge"]},
            contexts=["sex"],
            repeats=2,
            min_train_donors=4,
            min_test_donors=4,
        )

        self.assertEqual(set(result.feasibility["status"]), {"ok"})
        self.assertFalse(result.performance.empty)
        self.assertFalse(result.summary.empty)
        self.assertIn("selection_fraction", result.feature_stability.columns)


if __name__ == "__main__":
    unittest.main()
