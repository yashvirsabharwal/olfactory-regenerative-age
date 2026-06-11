from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.age_model import biological_feature_columns, donor_cv_folds, train_ora_models


class AgeModelTests(unittest.TestCase):
    def test_biological_features_exclude_yield_and_covariates(self):
        features = pd.DataFrame(
            {
                "donor_id": ["d1"],
                "total_cells": [100],
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

        self.assertEqual(set(result.performance["model"]), {"null_model", "elastic_net", "random_forest"})
        self.assertEqual(set(result.predictions["donor_id"]), set(donors[:10]))
        self.assertTrue(result.predictions["oraa"].notna().all())


if __name__ == "__main__":
    unittest.main()
