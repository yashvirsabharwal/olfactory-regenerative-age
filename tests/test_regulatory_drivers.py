import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from ora.regulatory_drivers import (
    build_driver_age_associations,
    build_driver_ora_correlations,
    build_regulatory_driver_map,
    parse_driver_metadata,
    score_regulatory_driver_activity,
)


class RegulatoryDriverTests(unittest.TestCase):
    def test_driver_scoring_and_rank_tables(self):
        config = {
            "gene_sets": {
                "ascl1_neurogenic": {
                    "driver": "ASCL1",
                    "driver_class": "transcription_factor",
                    "target_theme": "neural_progenitor",
                    "expected_age_direction": "negative",
                    "source": "literature",
                    "citation": "PMID:1",
                    "description": "ASCL1 targets",
                    "genes": ["ASCL1", "NEUROD1"],
                }
            }
        }
        driver_metadata = parse_driver_metadata(config)
        donors = [f"d{i}" for i in range(10)]
        metadata = pd.DataFrame(
            {
                "pseudobulk_id": [f"PB{i:03d}" for i in range(10)],
                "donor_id": donors,
                "sample_id": donors,
                "disease_group": ["healthy"] * 10,
                "coarse_cell_type": ["Olf_INPs"] * 10,
                "fine_cell_type": ["Early_INP"] * 10,
                "n_cells": [100] * 10,
                "sum_n_counts": [1000] * 10,
                "include_for_de": [True] * 10,
            }
        )
        decreasing = np.linspace(100, 10, 10).astype(int)
        counts = pd.DataFrame(
            {
                "gene_id": ["g1", "g2", "g3"],
                "gene_symbol": ["ASCL1", "NEUROD1", "OTHER"],
                "gene_index": [0, 1, 2],
                **{f"PB{i:03d}": [decreasing[i], decreasing[i], 1] for i in range(10)},
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "age": np.linspace(30, 80, 10),
                "usable_for_ora_training": [True] * 10,
                "passes_strict_ora_training_rule": [True] * 10,
            }
        )
        ora_scores = pd.DataFrame(
            {
                "donor_id": donors,
                "model": ["ridge"] * 10,
                "ora": np.linspace(80, 30, 10),
                "oraa": np.linspace(1, -1, 10),
            }
        )

        with tempfile.TemporaryDirectory() as tmp:
            counts_path = Path(tmp) / "counts.tsv.gz"
            counts.to_csv(counts_path, sep="\t", index=False)
            activity, donor_activity, coverage = score_regulatory_driver_activity(
                counts_path=counts_path,
                metadata=metadata,
                driver_metadata=driver_metadata,
                chunksize=2,
            )

        age = build_driver_age_associations(
            donor_activity=donor_activity,
            manifest=manifest,
            driver_metadata=driver_metadata,
        )
        corr = build_driver_ora_correlations(
            donor_activity=donor_activity,
            ora_scores=ora_scores,
            driver_metadata=driver_metadata,
        )
        driver_map = build_regulatory_driver_map(
            driver_metadata=driver_metadata,
            coverage=coverage,
            age_associations=age,
            ora_correlations=corr,
        )

        self.assertEqual(float(coverage.iloc[0]["coverage_fraction"]), 1.0)
        self.assertIn("gbc_inp", set(activity["lineage_group"]))
        primary = age[
            age["analysis_set"].eq("primary")
            & age["lineage_group"].eq("gbc_inp")
        ].iloc[0]
        self.assertEqual(primary["direction"], "negative")
        self.assertEqual(primary["observed_vs_expected"], "aligned")
        self.assertEqual(corr[corr["score_metric"].eq("oraa")].iloc[0]["direction"], "positive")
        self.assertEqual(driver_map.iloc[0]["driver"], "ASCL1")


if __name__ == "__main__":
    unittest.main()
