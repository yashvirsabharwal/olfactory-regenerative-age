from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.permutation import permute_training_ages, run_permutation_null


class PermutationTests(unittest.TestCase):
    def test_permute_training_ages_only_changes_eligible_donors(self):
        manifest = _manifest()
        rng = np.random.default_rng(4)

        permuted = permute_training_ages(manifest, rng)

        eligible_original = manifest.loc[:7, "age"].sort_values().to_list()
        eligible_permuted = permuted.loc[:7, "age"].sort_values().to_list()
        self.assertEqual(eligible_permuted, eligible_original)
        self.assertEqual(permuted.loc[8, "age"], manifest.loc[8, "age"])
        self.assertEqual(permuted.loc[9, "age"], manifest.loc[9, "age"])

    def test_run_permutation_null_returns_empirical_summary(self):
        manifest = _manifest(n=10)
        features = pd.DataFrame(
            {
                "donor_id": manifest["donor_id"],
                "prop__age_signal": np.linspace(0.0, 1.0, 10),
                "clr__age_signal": np.linspace(1.0, 0.0, 10),
            }
        )

        result = run_permutation_null(
            features,
            manifest,
            {"outer_cv_folds": 2, "random_seed": 1, "model_names": ["null_model", "ridge"]},
            n_permutations=3,
            repeats=1,
            random_seed=3,
        )

        self.assertEqual(set(result.empirical_summary["model"]), {"null_model", "ridge"})
        self.assertEqual(result.permutation_summary["permutation"].nunique(), 3)
        self.assertTrue({"empirical_p_mae", "observed_mae", "null_mae_mean"}.issubset(result.empirical_summary.columns))


def _manifest(n: int = 10) -> pd.DataFrame:
    donors = [f"d{i}" for i in range(n)]
    return pd.DataFrame(
        {
            "donor_id": donors,
            "sample_id": [f"s{i}" for i in range(n)],
            "age": np.linspace(30, 80, n),
            "sex": ["F", "M"] * (n // 2),
            "race_ethnicity": ["reported"] * n,
            "disease_group": ["healthy"] * (n - 2) + ["ad", "pd"],
            "chemistry": ["v2"] * n,
            "collection_method": ["device"] * n,
            "site": ["site1"] * n,
            "total_cells": np.arange(n) + 100,
            "usable_for_ora_training": [True] * (n - 2) + [False, False],
        }
    )


if __name__ == "__main__":
    unittest.main()
