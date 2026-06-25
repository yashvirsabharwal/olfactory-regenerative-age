from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.sensitivity import filter_manifest_for_scenario, run_ora_sensitivity


class SensitivityTests(unittest.TestCase):
    def test_filter_manifest_for_scenario_handles_cell_threshold(self):
        manifest = pd.DataFrame({"donor_id": ["d1", "d2"], "total_cells": [100, 1000]})
        filtered = filter_manifest_for_scenario(
            manifest,
            {"scenario": "min_total_cells__500", "filter_type": "min_total_cells", "filter_value": 500},
        )
        self.assertEqual(filtered["donor_id"].tolist(), ["d2"])

    def test_filter_manifest_for_scenario_handles_lineage_and_training_rules(self):
        manifest = pd.DataFrame(
            {
                "donor_id": ["d1", "d2", "d3"],
                "lineage_cells": [10, 100, 1000],
                "mature_neurons": [2, 60, 10],
                "passes_strict_ora_training_rule": [False, True, "False"],
            }
        )

        lineage = filter_manifest_for_scenario(
            manifest,
            {"scenario": "min_lineage_cells__100", "filter_type": "min_lineage_cells", "filter_value": 100},
        )
        mature = filter_manifest_for_scenario(
            manifest,
            {"scenario": "min_mature_neurons__50", "filter_type": "min_mature_neurons", "filter_value": 50},
        )
        strict = filter_manifest_for_scenario(
            manifest,
            {
                "scenario": "strict_ora_training_rule",
                "filter_type": "training_rule_column",
                "filter_value": "passes_strict_ora_training_rule",
            },
        )

        self.assertEqual(lineage["donor_id"].tolist(), ["d2", "d3"])
        self.assertEqual(mature["donor_id"].tolist(), ["d2"])
        self.assertEqual(strict["donor_id"].tolist(), ["d2"])

    def test_filter_manifest_handles_compound_and_yield_extremes(self):
        manifest = pd.DataFrame(
            {
                "donor_id": ["d1", "d2", "d3", "d4", "d5"],
                "chemistry": ["flex_v2", "flex_v2", "flex_v1", "flex_v2", "flex_v2"],
                "collection_method": ["device", "brush", "device", "device", "device"],
                "total_cells": [10, 100, 200, 300, 1000],
            }
        )
        compound = filter_manifest_for_scenario(
            manifest,
            {
                "scenario": "matched_flex_v2_device",
                "filter_type": "compound",
                "filter_value": "chemistry=flex_v2;collection_method=device",
            },
        )
        yield_filtered = filter_manifest_for_scenario(
            manifest,
            {
                "scenario": "exclude_yield_extremes__20pct",
                "filter_type": "exclude_yield_extremes",
                "filter_value": 0.20,
            },
        )

        self.assertEqual(compound["donor_id"].tolist(), ["d1", "d4", "d5"])
        self.assertEqual(yield_filtered["donor_id"].tolist(), ["d2", "d3", "d4"])

    def test_run_ora_sensitivity_records_runnable_and_skipped_scenarios(self):
        donors = [f"d{i}" for i in range(8)]
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "sample_id": donors,
                "age": [30, 35, 40, 45, 50, 55, 60, 65],
                "disease_group": ["healthy"] * 8,
                "sex": ["female", "male"] * 4,
                "usable_for_ora_training": [True] * 8,
                "chemistry": ["v1"] * 6 + ["v2"] * 2,
                "collection_method": ["brush", "device"] * 4,
                "site": ["toy"] * 8,
                "total_cells": [1000] * 8,
                "lineage_cells": [100] * 8,
                "mature_neurons": [60] * 8,
                "passes_strict_ora_training_rule": [True] * 8,
            }
        )
        features = pd.DataFrame(
            {
                "donor_id": donors,
                "prop__qHBC": [0.1, 0.2, 0.2, 0.4, 0.5, 0.6, 0.6, 0.8],
                "clr__qHBC": [-2, -1, -0.5, 0, 0.3, 0.6, 0.8, 1.0],
            }
        )
        result = run_ora_sensitivity(
            features,
            manifest,
            {
                "outer_cv_folds": 2,
                "random_seed": 1,
                "missingness_max_fraction": 1.0,
                "model_names": ["ridge", "random_forest"],
            },
            min_cell_thresholds=[500],
            min_train_donors=4,
        )

        self.assertIn("baseline", set(result.scenarios["scenario"]))
        self.assertIn("chemistry__v2", set(result.scenarios["scenario"]))
        self.assertIn("strict_ora_training_rule", set(result.scenarios["scenario"]))
        skipped = result.scenarios.set_index("scenario").loc["chemistry__v2", "status"]
        self.assertEqual(skipped, "too_few_training_donors")
        self.assertFalse(result.performance.empty)
        self.assertIn("scenario", result.performance.columns)
        self.assertTrue({"backend", "backend_package", "fallback_used"}.issubset(result.performance.columns))


if __name__ == "__main__":
    unittest.main()
