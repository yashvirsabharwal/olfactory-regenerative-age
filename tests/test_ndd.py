from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.ndd import (
    compare_ndd_feature_sets,
    donor_projection_appendix,
    ndd_projection_diagnostics,
    summarize_ndd_projection_uncertainty,
)


class NDDUncertaintyTests(unittest.TestCase):
    def test_uncertainty_uses_matched_healthy_reference(self):
        projection = pd.DataFrame(
            {
                "donor_id": ["h1", "h2", "h3", "a1", "a2"] * 2,
                "model": ["random_forest"] * 5 + ["elastic_net"] * 5,
                "disease_group": ["healthy", "healthy", "healthy", "ad", "ad"] * 2,
                "chemistry": ["v1", "v2", "v2", "v2", "v2"] * 2,
                "collection_method": ["brush", "device", "device", "device", "device"] * 2,
                "chronological_age": [50, 60, 70, 75, 80] * 2,
                "total_cells": [100, 200, 300, 400, 500] * 2,
                "oraa": [1.0, 2.0, 3.0, -5.0, -7.0, 1.0, 2.0, 3.0, -4.0, -6.0],
            }
        )

        result = summarize_ndd_projection_uncertainty(
            projection,
            n_bootstrap=100,
            random_seed=1,
        )
        rf_matched = result.uncertainty[
            result.uncertainty["model"].eq("random_forest")
            & result.uncertainty["disease_group"].eq("ad")
            & result.uncertainty["reference"].eq("matched_healthy")
        ].iloc[0]

        self.assertEqual(int(rf_matched["n_disease"]), 2)
        self.assertEqual(int(rf_matched["n_reference"]), 2)
        self.assertLess(rf_matched["difference_vs_reference"], 0)
        self.assertIn("collection_method", result.context.columns)

    def test_compare_feature_sets_and_appendix(self):
        base = pd.DataFrame(
            {
                "donor_id": ["h1", "a1", "p1"] * 2,
                "model": ["random_forest"] * 3 + ["xgboost"] * 3,
                "disease_group": ["healthy", "ad", "pd"] * 2,
                "sex": ["F", "F", "M"] * 2,
                "race_ethnicity": ["reported"] * 6,
                "chemistry": ["v2"] * 6,
                "collection_method": ["device"] * 6,
                "site": ["s1"] * 6,
                "chronological_age": [55, 70, 72] * 2,
                "total_cells": [100, 200, 300] * 2,
                "ora": [55, 62, 65, 55, 61, 64],
                "oraa": [0, -8, -7, 0, -9, -8],
                "training_n": [10] * 6,
                "n_features": [100] * 6,
            }
        )
        augmented = base.copy()
        augmented["ora"] = augmented["ora"] - 1
        augmented["oraa"] = augmented["oraa"] - 1
        augmented["n_features"] = 110

        comparison = compare_ndd_feature_sets({"composition": base, "augmented": augmented})
        appendix = donor_projection_appendix(augmented, feature_set="augmented", models=["random_forest"])

        rf_ad = comparison[
            comparison["model"].eq("random_forest")
            & comparison["disease_group"].eq("ad")
        ].iloc[0]
        self.assertEqual(rf_ad["composition_n_features"], 100)
        self.assertEqual(rf_ad["augmented_n_features"], 110)
        self.assertEqual(rf_ad["augmented_minus_composition_oraa"], -1)
        self.assertTrue(rf_ad["sign_stable_negative"])
        self.assertEqual(set(appendix["donor_id"]), {"a1", "p1"})
        self.assertEqual(set(appendix["feature_set"]), {"augmented"})

    def test_projection_diagnostics_flags_single_donor_strata(self):
        projection = pd.DataFrame(
            {
                "donor_id": ["a1", "a2", "p1"],
                "model": ["random_forest", "random_forest", "random_forest"],
                "disease_group": ["ad", "ad", "pd"],
                "sex": ["F", "M", "F"],
                "chemistry": ["v2", "v2", "v2"],
                "collection_method": ["device", "device", "device"],
                "site": ["", "", ""],
                "chronological_age": [70, 80, 75],
                "total_cells": [100, 300, 200],
                "ora": [60, 62, 58],
                "oraa": [-8, -10, -12],
            }
        )

        diagnostics = ndd_projection_diagnostics(projection)

        self.assertTrue({"sex", "age_bin", "cell_yield_quartile"}.issubset(set(diagnostics["diagnostic"])))
        sex_rows = diagnostics[diagnostics["diagnostic"].eq("sex")]
        self.assertIn("single_donor_stratum", set(sex_rows["status"]))
        ad_chem = diagnostics[
            diagnostics["disease_group"].eq("ad")
            & diagnostics["diagnostic"].eq("chemistry")
            & diagnostics["level"].eq("v2")
        ].iloc[0]
        self.assertEqual(ad_chem["n_donors"], 2)
        self.assertEqual(ad_chem["status"], "ok")


if __name__ == "__main__":
    unittest.main()
