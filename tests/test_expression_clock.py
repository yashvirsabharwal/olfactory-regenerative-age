import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from ora.expression_clock import load_donor_logcpm_expression, run_expression_clock_baseline


class ExpressionClockTests(unittest.TestCase):
    def test_load_donor_logcpm_expression_aggregates_pseudobulk_groups(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            counts = pd.DataFrame(
                {
                    "gene_id": ["g1", "g2"],
                    "gene_symbol": ["G1", "G2"],
                    "gene_index": [0, 1],
                    "PB1": [10, 0],
                    "PB2": [5, 2],
                    "PB3": [0, 8],
                }
            )
            metadata = pd.DataFrame(
                {
                    "pseudobulk_id": ["PB1", "PB2", "PB3"],
                    "donor_id": ["D1", "D1", "D2"],
                }
            )
            counts_path = root / "counts.tsv.gz"
            metadata_path = root / "metadata.tsv"
            counts.to_csv(counts_path, sep="\t", index=False)
            metadata.to_csv(metadata_path, sep="\t", index=False)

            expression, gene_qc = load_donor_logcpm_expression(
                counts_path,
                metadata_path,
                min_detection_donors=1,
                chunksize=1,
            )

            self.assertEqual(expression.shape, (2, 2))
            self.assertEqual(expression.index.tolist(), ["D1", "D2"])
            self.assertEqual(gene_qc["gene_id"].tolist(), ["g1", "g2"])

    def test_run_expression_clock_baseline_returns_summary_and_feasibility(self):
        rng = np.random.default_rng(7)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            donors = [f"D{i:02d}" for i in range(14)]
            genes = [f"g{i:02d}" for i in range(8)]
            rows = []
            for gene_idx, gene in enumerate(genes):
                row = {"gene_id": gene, "gene_symbol": gene.upper(), "gene_index": gene_idx}
                for donor_idx, donor in enumerate(donors):
                    age_signal = donor_idx * (1 if gene_idx < 3 else 0)
                    row[f"PB{donor_idx:02d}"] = int(max(0, rng.poisson(10 + age_signal)))
                rows.append(row)
            counts = pd.DataFrame(rows)
            metadata = pd.DataFrame(
                {
                    "pseudobulk_id": [f"PB{i:02d}" for i in range(len(donors))],
                    "donor_id": donors,
                }
            )
            manifest = pd.DataFrame(
                {
                    "donor_id": donors,
                    "sample_id": donors,
                    "age": np.linspace(30, 82, len(donors)),
                    "usable_for_ora_training": [True] * len(donors),
                    "disease_group": ["healthy"] * len(donors),
                }
            )
            counts_path = root / "counts.tsv.gz"
            metadata_path = root / "metadata.tsv"
            counts.to_csv(counts_path, sep="\t", index=False)
            metadata.to_csv(metadata_path, sep="\t", index=False)

            result = run_expression_clock_baseline(
                counts_path=counts_path,
                metadata_path=metadata_path,
                manifest=manifest,
                model_config={"outer_cv_folds": 3, "random_seed": 3, "age_bins": {"all": [0, 120]}},
                models=["ridge"],
                repeats=2,
                n_pcs=3,
                top_variable_genes=5,
                min_detection_donors=1,
                chunksize=2,
            )

            self.assertIn("mae_mean", result.performance_summary.columns)
            self.assertEqual(result.performance_summary["model"].tolist(), ["ridge"])
            self.assertEqual(set(result.feasibility["baseline"]), {"public_scageclock_direct", "fold_internal_pseudobulk_expression_pcs"})
            self.assertFalse(result.predictions["ora"].isna().any())


if __name__ == "__main__":
    unittest.main()
