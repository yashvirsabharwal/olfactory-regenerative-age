import unittest

import pandas as pd

from ora.foundation_benchmark import (
    build_donor_split_table,
    build_gene_manifest,
    select_foundation_subset_indices,
)


def _config():
    return {
        "columns": {
            "donor_id": ["donor_id"],
            "sample_id": ["sample_id"],
            "age": ["development_stage"],
            "disease": ["condition"],
            "coarse_cell_type": ["coarse_celltype"],
            "fine_cell_type": ["fine_celltype"],
            "sex": ["sex"],
            "chemistry": ["flex_version"],
            "collection_method": ["device_guided"],
        },
        "healthy_values": ["healthy", "normal"],
    }


class FoundationBenchmarkTests(unittest.TestCase):
    def test_select_lineage_subset_uses_olfactory_lineage_labels(self):
        obs = pd.DataFrame(
            {
                "donor_id": ["D1", "D1", "D2", "D2"],
                "sample_id": ["S1", "S1", "S2", "S2"],
                "coarse_celltype": ["Olf_mOSNs", "Resp_Secretory", "Tcell", "Olf_Sus"],
                "fine_celltype": ["Early_iOSN", "Goblet", "NaiveB", "Fully_mature_mOSN"],
                "_ora_row_position": [0, 1, 2, 3],
            }
        )

        selected = select_foundation_subset_indices(obs, _config(), "lineage", 10, seed=1)

        self.assertEqual(selected.tolist(), [0, 3])

    def test_stratified_sampling_is_deterministic_and_bounded(self):
        obs = pd.DataFrame(
            {
                "donor_id": [f"D{i % 4}" for i in range(100)],
                "sample_id": [f"S{i % 4}" for i in range(100)],
                "coarse_celltype": ["Resp_Ciliated"] * 100,
                "fine_celltype": ["Multiciliated" if i % 2 else "Club" for i in range(100)],
                "_ora_row_position": list(range(100)),
            }
        )

        first = select_foundation_subset_indices(obs, _config(), "epithelial", 17, seed=5)
        second = select_foundation_subset_indices(obs, _config(), "epithelial", 17, seed=5)

        self.assertEqual(first.tolist(), second.tolist())
        self.assertEqual(first.size, 17)
        self.assertTrue((first[:-1] <= first[1:]).all())

    def test_build_donor_split_table_exports_train_and_test_rows(self):
        manifest = pd.DataFrame(
            {
                "donor_id": [f"D{i}" for i in range(12)],
                "sample_id": [f"S{i}" for i in range(12)],
                "age": [30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85],
                "sex": ["female", "male"] * 6,
                "chemistry": ["flex_v2"] * 12,
                "collection_method": ["device"] * 12,
                "disease_group": ["healthy"] * 12,
                "usable_for_ora_training": [True] * 12,
            }
        )

        splits = build_donor_split_table(manifest, _config(), {"outer_cv_folds": 3, "random_seed": 13})

        self.assertEqual(set(splits["split"]), {"train", "test"})
        self.assertEqual(splits["fold"].nunique(), 3)
        for fold, frame in splits.groupby("fold"):
            test_donors = set(frame.loc[frame["split"].eq("test"), "donor_id"])
            train_donors = set(frame.loc[frame["split"].eq("train"), "donor_id"])
            self.assertTrue(test_donors)
            self.assertFalse(test_donors & train_donors, msg=f"overlap in fold {fold}")

    def test_gene_manifest_prefers_feature_name_and_marks_protein_coding(self):
        var = pd.DataFrame(
            {
                "feature_name": ["GENE1", "GENE2"],
                "feature_type": ["protein_coding", "lncRNA"],
                "feature_biotype": ["gene", "gene"],
            },
            index=["ENSG1", "ENSG2"],
        )

        genes = build_gene_manifest(var, var.index)

        self.assertEqual(genes["gene_symbol"].tolist(), ["GENE1", "GENE2"])
        self.assertEqual(genes["is_protein_coding"].tolist(), [True, False])


if __name__ == "__main__":
    unittest.main()
