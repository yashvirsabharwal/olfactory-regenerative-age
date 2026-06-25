from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.aggregate import compute_lineage_ratios


def config():
    return {
        "lineage_cell_types": {
            "quiescent_hbc": ["Quiescent HBC"],
            "activated_hbc": ["Activated HBC"],
            "cycling_hbc": ["Cycling HBC"],
            "early_inp": ["Early INP"],
            "late_inp": ["Late INP"],
            "early_iosn": ["Early iOSN"],
            "late_iosn": ["Late iOSN"],
            "early_mature_mosn": ["Early mature mOSN"],
            "fully_mature_mosn": ["Fully mature mOSN"],
            "stressed_mosn": ["Stressed mOSN"],
        }
    }


class LineageRatioTests(unittest.TestCase):
    def test_all_lineage_ratio_formulas(self):
        donor_counts = pd.DataFrame(
            {
                "Quiescent HBC": [4],
                "Activated HBC": [2],
                "Cycling HBC": [1],
                "Early INP": [3],
                "Late INP": [1],
                "Early iOSN": [5],
                "Late iOSN": [1],
                "Early mature mOSN": [6],
                "Fully mature mOSN": [2],
                "Stressed mOSN": [1],
                "Goblet": [9],
            },
            index=pd.Index(["d1"], name="donor_id"),
        )

        ratios = compute_lineage_ratios(donor_counts, config(), pseudocount=0.5).set_index("donor_id")
        row = ratios.loc["d1"]

        self.assertNotIn("ratio__hbc_to_inp", ratios.columns)
        self.assertNotIn("ratio__iosn_to_mosn", ratios.columns)
        self.assertAlmostEqual(row["ratio__neuronal_fraction"], 19 / 35.5)
        self.assertAlmostEqual(row["ratio__mature_neuron_fraction"], 9 / 35.5)
        self.assertAlmostEqual(row["ratio__immature_to_mature"], 6 / 9.5)
        self.assertAlmostEqual(row["ratio__progenitor_to_neuron"], 4 / 15.5)
        self.assertAlmostEqual(row["ratio__activated_to_quiescent_hbc"], 2 / 4.5)
        self.assertAlmostEqual(row["ratio__inp_to_activated_hbc"], 4 / 2.5)
        self.assertAlmostEqual(row["ratio__inp_to_iosn"], 6 / 4.5)
        self.assertAlmostEqual(row["ratio__mature_mosn_to_iosn"], 8 / 6.5)
        self.assertAlmostEqual(row["ratio__stressed_to_mature_mosn"], 1 / 8.5)
        self.assertAlmostEqual(row["ratio__lineage_fraction"], 26 / 35.5)


if __name__ == "__main__":
    unittest.main()
