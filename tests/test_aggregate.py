from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.aggregate import aggregate_cell_counts, build_cell_state_features


def config():
    return {
        "columns": {
            "donor_id": ["donor_id"],
            "sample_id": ["sample_id"],
            "age": ["age"],
            "disease": ["disease"],
            "coarse_cell_type": ["coarse_cell_type"],
            "fine_cell_type": ["fine_cell_type"],
        },
        "lineage_cell_types": {
            "quiescent_hbc": ["Quiescent HBC"],
            "activated_hbc": ["Activated HBC"],
            "early_inp": ["Early INP"],
            "late_inp": ["Late INP"],
            "early_iosn": ["Early iOSN"],
            "late_iosn": ["Late iOSN"],
            "early_mature_mosn": ["Early mature mOSN"],
            "fully_mature_mosn": ["Fully mature mOSN"],
            "stressed_mosn": ["Stressed mOSN"],
        },
    }


class AggregateTests(unittest.TestCase):
    def test_counts_and_features(self):
        obs = pd.DataFrame(
            [
                ["d1", "s1", 40, "healthy", "HBC", "Quiescent HBC"],
                ["d1", "s1", 40, "healthy", "Neuron", "Early iOSN"],
                ["d1", "s1", 40, "healthy", "Neuron", "Fully mature mOSN"],
                ["d2", "s2", 70, "healthy", "HBC", "Activated HBC"],
                ["d2", "s2", 70, "healthy", "Neuron", "Stressed mOSN"],
            ],
            columns=["donor_id", "sample_id", "age", "disease", "coarse_cell_type", "fine_cell_type"],
        )

        counts = aggregate_cell_counts(obs, config())
        features = build_cell_state_features(counts, config())

        self.assertEqual(int(counts["cell_count"].sum()), 5)
        self.assertIn("prop__fully_mature_mosn", features.columns)
        self.assertIn("clr__stressed_mosn", features.columns)
        self.assertIn("ratio__immature_to_mature", features.columns)
        self.assertTrue(np.isfinite(features.filter(regex="^(prop|clr|ratio)__").to_numpy()).all())


if __name__ == "__main__":
    unittest.main()

