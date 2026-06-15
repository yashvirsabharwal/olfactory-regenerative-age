from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.age_model import (
    MODEL_ORDER,
    biological_feature_columns,
    donor_cv_folds,
    model_names_from_config,
    project_ora_models,
    train_ora_models,
    train_ora_models_repeated,
)


class AgeModelTests(unittest.TestCase):
    def test_biological_features_exclude_yield_and_covariates(self):
        features = pd.DataFrame(
            {
                "donor_id": ["d1"],
                "total_cells": [100],
                "lineage_cells": [20],
                "mature_neurons": [5],
                "has_age": [True],
                "is_training_donor": [True],
                "chemistry": ["v2"],
                "prop__hbc": [0.2],
                "clr__mosn": [0.1],
                "ratio__iosn_to_mosn": [0.3],
            }
        )

        cols = biological_feature_columns(features, {"exclude_from_biological_features": ["chemistry"]})

        self.assertEqual(cols, ["prop__hbc", "clr__mosn", "ratio__iosn_to_mosn"])

    def test_cv_folds_are_donor_level_and_cover_all_rows(self):
        data = pd.DataFrame({"donor_id": [f"d{i}" for i in range(10)], "age": np.linspace(20, 80, 10)})
        folds = donor_cv_folds(data, {"outer_cv_folds": 5, "random_seed": 1, "age_bins": {"a": [0, 49], "b": [50, 120]}})
        test_rows = sorted(int(i) for _, test in folds for i in test)

        self.assertEqual(test_rows, list(range(10)))
        for train, test in folds:
            self.assertTrue(set(train).isdisjoint(set(test)))

    def test_train_models_excludes_ndd_and_returns_scores(self):
        donors = [f"d{i}" for i in range(12)]
        ages = np.linspace(35, 80, 12)
        features = pd.DataFrame(
            {
                "donor_id": donors,
                "prop__young_feature": np.linspace(1, 0, 12),
                "clr__old_feature": np.linspace(0, 1, 12),
                "ratio__lineage_fraction": np.linspace(0.2, 0.8, 12),
                "total_cells": np.arange(12) + 100,
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "sample_id": [f"s{i}" for i in range(12)],
                "age": ages,
                "sex": ["F", "M"] * 6,
                "race_ethnicity": ["reported"] * 12,
                "disease_group": ["healthy"] * 10 + ["ad", "pd"],
                "chemistry": ["v2"] * 12,
                "collection_method": ["device"] * 12,
                "site": ["site1"] * 12,
                "total_cells": np.arange(12) + 100,
                "usable_for_ora_training": [True] * 10 + [False, False],
            }
        )

        result = train_ora_models(features, manifest, {"outer_cv_folds": 5, "random_seed": 1})

        self.assertEqual(set(result.performance["model"]), set(MODEL_ORDER))
        self.assertEqual(set(result.predictions["donor_id"]), set(donors[:10]))
        self.assertTrue(result.predictions["oraa"].notna().all())

    def test_model_names_from_config_supports_subset_and_enabled_flags(self):
        self.assertEqual(model_names_from_config({"model_names": ["xgboost", "random_forest"]}), ["random_forest", "xgboost"])
        self.assertNotIn(
            "catboost",
            model_names_from_config({"models": {"catboost": {"enabled": False}}}),
        )
        with self.assertRaises(ValueError):
            model_names_from_config({"model_names": ["not_a_model"]})
        with self.assertRaises(ValueError):
            model_names_from_config({"model_names": []})
        with self.assertRaises(ValueError):
            model_names_from_config({"models": {name: {"enabled": False} for name in MODEL_ORDER}})

    def test_repeated_models_return_intervals_and_feature_stability(self):
        donors = [f"d{i}" for i in range(12)]
        ages = np.linspace(35, 80, 12)
        features = pd.DataFrame(
            {
                "donor_id": donors,
                "prop__young_feature": np.linspace(1, 0, 12),
                "clr__old_feature": np.linspace(0, 1, 12),
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "sample_id": [f"s{i}" for i in range(12)],
                "age": ages,
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

        result = train_ora_models_repeated(
            features,
            manifest,
            {"outer_cv_folds": 3, "outer_cv_repeats": 2, "random_seed": 1},
        )

        self.assertEqual(set(result.performance_summary["model"]), set(MODEL_ORDER))
        self.assertTrue({"mae_mean", "mae_ci_low", "mae_ci_high"}.issubset(result.performance_summary.columns))
        self.assertEqual(result.repeat_performance["repeat"].nunique(), 2)
        self.assertFalse(result.feature_stability.empty)
        self.assertTrue({"selection_fraction", "mean_importance"}.issubset(result.feature_stability.columns))

    def test_project_models_scores_ndd_without_training_on_them(self):
        donors = [f"d{i}" for i in range(14)]
        ages = np.linspace(35, 82, 14)
        features = pd.DataFrame(
            {
                "donor_id": donors,
                "prop__young_feature": np.linspace(1, 0, 14),
                "clr__old_feature": np.linspace(0, 1, 14),
                "module_score__stress": np.linspace(0.2, 1.2, 14),
                "total_cells": np.arange(14) + 100,
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "sample_id": [f"s{i}" for i in range(14)],
                "age": ages,
                "sex": ["F", "M"] * 7,
                "race_ethnicity": ["reported"] * 14,
                "disease_group": ["healthy"] * 11 + ["ad", "pd", "healthy"],
                "disease": ["Healthy"] * 11 + ["AD", "PD", "Healthy"],
                "chemistry": ["v2"] * 14,
                "collection_method": ["device"] * 14,
                "site": ["site1"] * 14,
                "total_cells": np.arange(14) + 100,
                "is_ndd": [False] * 11 + [True, True, False],
                "usable_for_ora_training": [True] * 11 + [False, False, "False"],
            }
        )

        result = project_ora_models(features, manifest, {"outer_cv_folds": 5, "random_seed": 1})
        ndd = result.predictions[result.predictions["disease_group"].isin(["ad", "pd"])]
        missing_age_healthy = result.predictions[result.predictions["donor_id"].eq("d13")]

        self.assertEqual(set(ndd["donor_id"]), {"d11", "d12"})
        self.assertFalse(ndd["is_training_donor"].astype(bool).any())
        self.assertFalse(missing_age_healthy["is_training_donor"].astype(bool).any())
        self.assertTrue(ndd["ora"].notna().all())
        self.assertTrue(ndd["oraa"].notna().all())
        self.assertIn("ad", set(result.summary["disease_group"]))


if __name__ == "__main__":
    unittest.main()
