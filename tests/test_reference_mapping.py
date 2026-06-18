from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.reference_mapping import mapped_label_donor_features, mapping_qc_by_sample, normalized_entropy


class ReferenceMappingTests(unittest.TestCase):
    def test_mapped_label_donor_features_emits_prop_and_clr_columns(self):
        obs = pd.DataFrame(
            {
                "dataset_id": ["x", "x", "x", "x"],
                "sample_id": ["s1", "s1", "s1", "s1"],
                "donor_id": ["d1", "d1", "d1", "d1"],
                "age": [70, 70, 70, 70],
                "disease_group": ["healthy", "healthy", "healthy", "healthy"],
                "scanvi_predicted_label": ["Mature OSN", "Mature OSN", "HBC", "HBC"],
                "scanvi_label_confidence": [0.9, 0.8, 0.7, 0.6],
            }
        )

        features = mapped_label_donor_features(
            obs,
            label_column="scanvi_predicted_label",
            confidence_column="scanvi_label_confidence",
        )

        self.assertEqual(features.shape[0], 1)
        self.assertAlmostEqual(float(features.loc[0, "prop__mature_osn"]), 0.5)
        self.assertIn("clr__hbc", features.columns)
        self.assertAlmostEqual(float(features.loc[0, "mean_label_confidence"]), 0.75)

    def test_mapping_qc_by_sample_marks_confidence_status(self):
        obs = pd.DataFrame(
            {
                "dataset_id": ["x", "x"],
                "sample_id": ["s1", "s1"],
                "donor_id": ["d1", "d1"],
                "age": [70, 70],
                "disease_group": ["healthy", "healthy"],
                "scanvi_predicted_label": ["a", "b"],
                "scanvi_label_confidence": [0.8, 0.9],
                "scanvi_label_entropy": [0.1, 0.2],
            }
        )

        qc = mapping_qc_by_sample(
            obs,
            label_column="scanvi_predicted_label",
            confidence_column="scanvi_label_confidence",
            entropy_column="scanvi_label_entropy",
        )

        self.assertEqual(qc.loc[0, "status"], "high_confidence")
        self.assertEqual(int(qc.loc[0, "n_labels"]), 2)

    def test_normalized_entropy_bounds(self):
        entropy = normalized_entropy(np.array([[1.0, 0.0], [0.5, 0.5]]))

        self.assertAlmostEqual(float(entropy[0]), 0.0, places=5)
        self.assertAlmostEqual(float(entropy[1]), 1.0, places=5)


if __name__ == "__main__":
    unittest.main()
