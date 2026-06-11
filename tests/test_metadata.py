from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.metadata import (
    build_manifest,
    collection_method_group,
    disease_group,
    parse_age_series,
    resolve_columns,
)


def config():
    return {
        "columns": {
            "donor_id": ["donor_id"],
            "sample_id": ["sample_id"],
            "age": ["age"],
            "sex": ["sex"],
            "race_ethnicity": ["race_ethnicity"],
            "disease": ["disease_condition"],
            "chemistry": ["chemistry"],
            "collection_method": ["device_usage"],
            "coarse_cell_type": ["coarse_cell_type"],
            "fine_cell_type": ["fine_cell_type"],
        },
        "healthy_values": ["healthy", "cognitively normal"],
        "ndd_values": {
            "ad": ["Alzheimer's disease"],
            "pd": ["Parkinson's disease"],
        },
        "lineage_cell_types": {
            "quiescent_hbc": ["Quiescent HBC"],
            "activated_hbc": ["Activated HBC"],
            "fully_mature_mosn": ["Fully mature mOSN"],
            "stressed_mosn": ["Stressed mOSN"],
        },
        "collection_method_values": {
            "device": ["T", "device"],
            "brush": ["F", "brush"],
        },
    }


class MetadataTests(unittest.TestCase):
    def test_resolve_columns_from_aliases(self):
        obs = pd.DataFrame(
            {
                "donor_id": [],
                "sample_id": [],
                "age": [],
                "disease_condition": [],
                "coarse_cell_type": [],
                "fine_cell_type": [],
            }
        )
        resolved = resolve_columns(list(obs.columns), config())

        self.assertEqual(resolved.donor_id, "donor_id")
        self.assertEqual(resolved.disease, "disease_condition")

    def test_build_manifest_excludes_ndd_from_training(self):
        obs = pd.DataFrame(
            [
                ["d1", "s1", 44, "healthy", "F", "White", "v2", "device", "HBC", "Quiescent HBC"],
                ["d1", "s1", 44, "healthy", "F", "White", "v2", "device", "Neuron", "Fully mature mOSN"],
                ["d2", "s2", 75, "Alzheimer's disease", "M", "White", "v2", "device", "Neuron", "Stressed mOSN"],
                ["d3", "s3", None, "healthy", "F", "Asian", "v1", "brush", "HBC", "Activated HBC"],
            ],
            columns=[
                "donor_id",
                "sample_id",
                "age",
                "disease_condition",
                "sex",
                "race_ethnicity",
                "chemistry",
                "device_usage",
                "coarse_cell_type",
                "fine_cell_type",
            ],
        )

        manifest = build_manifest(obs, config())
        trainable = set(manifest.loc[manifest["usable_for_ora_training"], "donor_id"])

        self.assertEqual(trainable, {"d1"})
        self.assertEqual(manifest.loc[manifest["donor_id"].eq("d2"), "disease_group"].iloc[0], "ad")
        self.assertEqual(int(manifest.loc[manifest["donor_id"].eq("d1"), "mature_neurons"].iloc[0]), 1)

    def test_disease_group_normalization(self):
        self.assertEqual(disease_group("Cognitively Normal", config()), "healthy")
        self.assertEqual(disease_group("Parkinson disease", config()), "pd")

    def test_parse_age_from_development_stage(self):
        parsed = parse_age_series(pd.Series(["65-year-old stage", "unknown", "20-year-old stage"]))

        self.assertEqual(parsed.iloc[0], 65)
        self.assertTrue(pd.isna(parsed.iloc[1]))
        self.assertEqual(parsed.iloc[2], 20)

    def test_collection_method_normalization(self):
        self.assertEqual(collection_method_group("T", config()), "device")
        self.assertEqual(collection_method_group("F", config()), "brush")


if __name__ == "__main__":
    unittest.main()
