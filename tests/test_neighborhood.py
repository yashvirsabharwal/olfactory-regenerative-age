import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from ora.neighborhood import NeighborhoodConfig, run_neighborhood_da, summarize_neighborhood_da


class NeighborhoodTests(unittest.TestCase):
    def test_run_neighborhood_da_returns_age_statistics(self):
        rng = np.random.default_rng(1)
        donors = [f"d{idx:02d}" for idx in range(30)]
        donor_metadata = pd.DataFrame(
            {
                "donor_id": donors,
                "age": np.linspace(25, 80, len(donors)),
                "sex": ["female", "male"] * 15,
                "chemistry": ["flex_v1"] * 15 + ["flex_v2"] * 15,
                "collection_method": ["device"] * 30,
            }
        )
        cell_rows = []
        embedding_rows = []
        for donor_idx, donor in enumerate(donors):
            age = donor_metadata.loc[donor_idx, "age"]
            n_old_state = 6 if age > 55 else 1
            n_base_state = 6
            for _ in range(n_old_state):
                cell_rows.append({"donor_id": donor, "fine_celltype": "old_state", "coarse_celltype": "epithelial"})
                embedding_rows.append(rng.normal(loc=[3.0, 0.0], scale=0.1, size=2))
            for _ in range(n_base_state):
                cell_rows.append({"donor_id": donor, "fine_celltype": "base_state", "coarse_celltype": "epithelial"})
                embedding_rows.append(rng.normal(loc=[0.0, 0.0], scale=0.1, size=2))
        neighborhoods, summary = run_neighborhood_da(
            np.vstack(embedding_rows),
            pd.DataFrame(cell_rows),
            donor_metadata,
            config=NeighborhoodConfig(n_neighborhoods=12, n_neighbors=25, min_donors=8, seed=2),
        )

        self.assertEqual(neighborhoods.shape[0], 12)
        self.assertIn("age_fdr", neighborhoods.columns)
        self.assertGreaterEqual(int(summary.loc[summary["metric"].eq("neighborhoods_tested"), "value"].iloc[0]), 1)

    def test_summarize_neighborhood_da_handles_empty_input(self):
        summary = summarize_neighborhood_da(pd.DataFrame())
        self.assertEqual(summary.loc[0, "metric"], "neighborhoods")
        self.assertEqual(summary.loc[0, "value"], 0)


if __name__ == "__main__":
    unittest.main()
