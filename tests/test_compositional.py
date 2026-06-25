from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.compositional import run_compositional_age_model


class CompositionalModelTests(unittest.TestCase):
    def test_clr_age_model_recovers_direction_and_sensitivities(self):
        donors = [f"d{i:02d}" for i in range(36)]
        ages = np.linspace(25, 80, len(donors))
        count_rows = []
        manifest_rows = []
        for idx, (donor, age) in enumerate(zip(donors, ages, strict=True)):
            old_count = int(20 + age)
            young_count = int(130 - age)
            support_count = 25 + (idx % 3)
            for state, count in [
                ("Old State", old_count),
                ("Young State", young_count),
                ("Support State", support_count),
            ]:
                count_rows.append(
                    {
                        "donor_id": donor,
                        "sample_id": donor,
                        "coarse_cell_type": "epithelium",
                        "fine_cell_type": state,
                        "cell_count": count,
                    }
                )
            manifest_rows.append(
                {
                    "donor_id": donor,
                    "sample_id": donor,
                    "age": age,
                    "sex": "female" if idx % 2 else "male",
                    "chemistry": "flex_v1",
                    "collection_method": "brush",
                    "site": "",
                    "usable_for_ora_training": True,
                    "passes_strict_ora_training_rule": idx < 32,
                }
            )

        associations = pd.DataFrame(
            {
                "feature": ["clr__old_state", "clr__young_state"],
                "beta_per_10_years": [0.4, -0.4],
                "fdr": [0.01, 0.02],
                "direction": ["positive", "negative"],
            }
        )

        result = run_compositional_age_model(
            pd.DataFrame(count_rows),
            pd.DataFrame(manifest_rows),
            age_associations=associations,
            min_scenario_donors=20,
            min_nonzero_donors=3,
        )

        old_state = result.summary[result.summary["cell_state"].eq("Old State")].iloc[0]
        self.assertEqual(old_state["age_direction"], "positive")
        self.assertTrue(bool(old_state["direction_concordant_with_ora_age_association"]))
        self.assertTrue(bool(old_state["directionally_stable_in_sensitivity"]))
        self.assertIn("primary_all_healthy", set(result.sensitivity["scenario"]))
        self.assertIn("strict_threshold", set(result.sensitivity["scenario"]))
        self.assertIn("single_flex_v1_brush", set(result.sensitivity["scenario"]))


if __name__ == "__main__":
    unittest.main()
