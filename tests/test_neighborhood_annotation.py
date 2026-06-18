import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ora.neighborhood_annotation import annotate_neighborhood_table, build_neighborhood_annotation


class NeighborhoodAnnotationTests(unittest.TestCase):
    def test_annotate_neighborhood_table_adds_theme_and_claim_gate(self):
        neighborhoods = pd.DataFrame(
            {
                "neighborhood_id": [1, 2],
                "top_fine_celltype": ["Early_iOSN", "Mucous_gland"],
                "top_coarse_celltype": ["Olf_iOSNs", "Resp_Secretory"],
                "age_coef": [-0.8, 0.4],
                "age_fdr": [0.04, 0.20],
                "status": ["tested", "tested"],
            }
        )

        annotated = annotate_neighborhood_table(neighborhoods, run_name="lineage_matched")

        self.assertEqual(annotated.loc[0, "fine_theme"], "neuronal lineage maturation")
        self.assertTrue(bool(annotated.loc[0, "is_age_associated_fdr_0_10"]))
        self.assertEqual(annotated.loc[0, "claim_gate"], "matched_regenerative_neuronal_support")
        self.assertEqual(annotated.loc[1, "claim_gate"], "not_significant")

    def test_build_neighborhood_annotation_summarizes_significant_themes(self):
        neighborhoods = pd.DataFrame(
            {
                "neighborhood_id": [1, 2, 3],
                "top_fine_celltype": ["Early_iOSN", "Late_iOSN", "Naive_CD8"],
                "top_coarse_celltype": ["Olf_iOSNs", "Olf_iOSNs", "Tcell"],
                "age_coef": [-0.8, -0.6, -0.7],
                "age_fdr": [0.04, 0.03, 0.02],
                "status": ["tested", "tested", "tested"],
            }
        )

        top, summary = build_neighborhood_annotation({"all_matched": neighborhoods}, top_n=2)

        self.assertEqual(top.shape[0], 2)
        self.assertIn("neuronal lineage maturation", set(summary["fine_theme"]))
        self.assertIn("immune/inflammatory compartment", set(summary["fine_theme"]))


if __name__ == "__main__":
    unittest.main()
