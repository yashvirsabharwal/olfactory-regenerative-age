from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd

from ora.niche_signaling import (
    build_niche_age_associations,
    parse_niche_interactions,
    score_niche_interactions,
)


class TestNicheSignaling(unittest.TestCase):
    def test_scores_sender_receiver_interaction_from_pseudobulk(self):
        config = {
            "interactions": {
                "tnf_test": {
                    "family": "TNF",
                    "ligand_genes": ["TNF"],
                    "receptor_genes": ["TNFRSF1A"],
                    "sender_groups": ["immune"],
                    "receiver_groups": ["hbc"],
                    "expected_age_direction": "positive",
                }
            }
        }
        interactions = parse_niche_interactions(config)
        metadata = pd.DataFrame(
            {
                "pseudobulk_id": ["PB1", "PB2", "PB3", "PB4"],
                "donor_id": ["D1", "D1", "D2", "D2"],
                "sample_id": ["D1", "D1", "D2", "D2"],
                "disease_group": ["healthy", "healthy", "healthy", "healthy"],
                "coarse_cell_type": ["Dendritic", "Resp_HBC", "Dendritic", "Resp_HBC"],
                "fine_cell_type": ["Macrophage", "Quiescent_HBC", "Macrophage", "Quiescent_HBC"],
                "n_cells": [10, 20, 10, 20],
                "sum_n_counts": [100.0, 100.0, 100.0, 100.0],
            }
        )
        counts = pd.DataFrame(
            {
                "gene_id": ["g1", "g2"],
                "gene_symbol": ["TNF", "TNFRSF1A"],
                "gene_index": [0, 1],
                "PB1": [10, 0],
                "PB2": [0, 25],
                "PB3": [20, 0],
                "PB4": [0, 30],
            }
        )
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "counts.tsv"
            counts.to_csv(path, sep="\t", index=False)
            donor_scores, coverage = score_niche_interactions(
                counts_path=path,
                metadata=metadata,
                interactions=interactions,
                chunksize=1,
            )

        self.assertEqual(coverage.loc[0, "n_present"], 2)
        self.assertEqual(donor_scores.shape[0], 2)
        self.assertTrue(np.isfinite(donor_scores["interaction_score"]).all())
        self.assertGreater(donor_scores.loc[donor_scores["donor_id"].eq("D2"), "interaction_score"].iloc[0], 0)

    def test_age_association_returns_primary_row(self):
        donor_scores = pd.DataFrame(
            {
                "interaction_id": ["x"] * 8,
                "family": ["TNF"] * 8,
                "sender_group": ["immune"] * 8,
                "receiver_group": ["hbc"] * 8,
                "donor_id": [f"D{i}" for i in range(8)],
                "interaction_score": np.arange(8, dtype=float),
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": [f"D{i}" for i in range(8)],
                "age": np.arange(40, 48, dtype=float),
                "usable_for_ora_training": [True] * 8,
                "passes_strict_ora_training_rule": [False] * 8,
            }
        )
        interactions = pd.DataFrame(
            {
                "interaction_id": ["x"],
                "family": ["TNF"],
                "expected_age_direction": ["positive"],
                "ligand_genes": [("TNF",)],
                "receptor_genes": [("TNFRSF1A",)],
            }
        )
        result = build_niche_age_associations(
            donor_scores=donor_scores,
            manifest=manifest,
            interactions=interactions,
            covariates=(),
        )
        primary = result[result["analysis_set"].eq("primary")]
        self.assertEqual(primary.shape[0], 1)
        self.assertEqual(primary.iloc[0]["direction"], "positive")


if __name__ == "__main__":
    unittest.main()
