from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.stats import run_age_associations


class StatsTests(unittest.TestCase):
    def test_age_associations_do_not_drop_all_missing_covariate(self):
        donors = [f"d{i}" for i in range(8)]
        features = pd.DataFrame(
            {
                "donor_id": donors,
                "prop__state": np.linspace(0.1, 0.8, 8),
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "sample_id": [f"s{i}" for i in range(8)],
                "age": np.linspace(30, 80, 8),
                "sex": ["F", "M"] * 4,
                "site": [np.nan] * 8,
            }
        )

        result = run_age_associations(features, manifest, feature_columns=["prop__state"], covariates=["sex", "site"])

        self.assertEqual(result.loc[0, "status"], "ok")
        self.assertEqual(result.loc[0, "n"], 8)
        self.assertTrue(np.isfinite(result.loc[0, "beta_per_10_years"]))


if __name__ == "__main__":
    unittest.main()
