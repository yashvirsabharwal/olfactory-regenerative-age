from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config


class ConfigTests(unittest.TestCase):
    def test_default_gateway_config_loads(self):
        config = load_config(Path(__file__).resolve().parents[1] / "configs" / "gateway.yaml")

        self.assertEqual(config["source"]["doi"], "10.64898/2026.06.10.731272")
        self.assertEqual(config["source"]["dataset_id"], "16c13603-a36e-467a-a8e8-3d118f2eef45")
        self.assertEqual(config["source"]["dataset_version_id"], "060f346b-c744-4845-9191-65595dd893b9")
        self.assertEqual(config["source"]["h5ad_filesize_bytes"], 27651770343)
        self.assertEqual(config["paper_defaults"]["donors_total"], 202)
        self.assertIn("donor_id", config["columns"])
        self.assertIn("fully_mature_mosn", config["lineage_cell_types"])


if __name__ == "__main__":
    unittest.main()
