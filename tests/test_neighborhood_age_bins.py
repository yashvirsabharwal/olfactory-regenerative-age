import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ora.neighborhood_age_bins import AgeBinConfig, summarize_neighborhood_age_bins


class NeighborhoodAgeBinTests(unittest.TestCase):
    def test_summarize_neighborhood_age_bins_checks_sign_direction(self):
        donors = pd.DataFrame(
            {
                "donor_id": ["young1", "young2", "old1", "old2"],
                "age": [35, 40, 78, 82],
            }
        )
        memberships = pd.DataFrame(
            {
                "neighborhood_id": [0] * 8 + [1] * 8,
                "donor_id": ["young1"] * 3 + ["young2"] * 3 + ["old1"] + ["old2"] + ["young1"] + ["young2"] + ["old1"] * 3 + ["old2"] * 3,
            }
        )
        da = pd.DataFrame(
            {
                "neighborhood_id": [0, 1],
                "top_fine_celltype": ["Early_iOSN", "Cycling_HBC"],
                "top_coarse_celltype": ["Olf_iOSNs", "Resp_HBC"],
                "age_coef": [-0.8, 0.7],
                "age_fdr": [0.04, 0.05],
                "status": ["tested", "tested"],
            }
        )
        config = AgeBinConfig(
            run_name="toy",
            bins=(("young", 0, 50), ("old", 50, float("inf"))),
        )

        neighborhoods, summary = summarize_neighborhood_age_bins(memberships, donors, da_table=da, config=config)

        self.assertEqual(neighborhoods.shape[0], 2)
        self.assertTrue(neighborhoods["bin_agrees_with_regression"].all())
        self.assertEqual(
            int(summary.loc[summary["metric"].eq("negative_sig_bin_agreement"), "value"].iloc[0]),
            1,
        )
        self.assertEqual(
            int(summary.loc[summary["metric"].eq("positive_sig_bin_agreement"), "value"].iloc[0]),
            1,
        )


if __name__ == "__main__":
    unittest.main()
