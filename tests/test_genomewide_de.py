from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.genomewide_de import audit_genomewide_de, summarize_genomewide_de


class GenomewideDeTests(unittest.TestCase):
    def test_summarize_genomewide_de_counts_significant_and_sex_linked_hits(self):
        de = pd.DataFrame(
            {
                "contrast": ["ad_vs_healthy", "ad_vs_healthy", "pd_vs_healthy"],
                "fine_cell_type": ["qHBC", "qHBC", "mOSN"],
                "gene_symbol": ["MAFF", "USP9Y", "SNCA"],
                "log2fc": [1.0, 2.0, -1.5],
                "p_value": [0.002, 0.001, 0.2],
                "fdr": [0.02, 0.01, 0.3],
            }
        )
        run_summary = pd.DataFrame(
            {
                "contrast": ["ad_vs_healthy", "pd_vs_healthy", "pd_vs_healthy"],
                "fine_cell_type": ["qHBC", "mOSN", "qHBC"],
                "status": ["ok", "ok", "too_few_donors"],
            }
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            de_path = Path(tmpdir) / "de.tsv"
            summary_path = Path(tmpdir) / "run_summary.tsv"
            de.to_csv(de_path, sep="\t", index=False)
            run_summary.to_csv(summary_path, sep="\t", index=False)

            summary, top_hits = summarize_genomewide_de(de_path, summary_path, top_n=2)

        ad = summary[summary["contrast"].eq("ad_vs_healthy")].iloc[0]
        self.assertEqual(ad["tested_rows"], 2)
        self.assertEqual(ad["significant_rows"], 2)
        self.assertEqual(ad["significant_genes"], 2)
        self.assertEqual(ad["sex_linked_significant_rows"], 1)
        self.assertEqual(ad["ok_cell_state_models"], 1)
        self.assertEqual(top_hits["gene_symbol"].tolist(), ["USP9Y", "MAFF"])

    def test_audit_genomewide_de_summarizes_flags_and_matching(self):
        de = pd.DataFrame(
            {
                "contrast": ["pd_vs_healthy", "pd_vs_healthy", "ad_vs_healthy"],
                "fine_cell_type": ["qHBC", "qHBC", "mOSN"],
                "gene_symbol": ["USP9Y", "RPL3", "HBA1"],
                "log2fc": [5.0, 1.0, 2.0],
                "p_value": [0.001, 0.002, 0.003],
                "fdr": [0.01, 0.02, 0.03],
            }
        )
        run_summary = pd.DataFrame(
            {
                "contrast": ["pd_vs_healthy", "ad_vs_healthy"],
                "fine_cell_type": ["qHBC", "mOSN"],
                "n_case": [2, 5],
                "n_control": [9, 20],
                "n_genes_tested": [100, 100],
                "status": ["ok", "ok"],
            }
        )
        metadata = pd.DataFrame(
            {
                "pseudobulk_id": ["h1", "h2", "a1", "p1"],
                "donor_id": ["h1", "h2", "a1", "p1"],
                "disease_group": ["healthy", "healthy", "ad", "pd"],
                "fine_cell_type": ["qHBC", "qHBC", "qHBC", "qHBC"],
                "chemistry": ["flex_v2"] * 4,
                "collection_method": ["device"] * 4,
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            de_path = root / "de.tsv"
            summary_path = root / "summary.tsv"
            metadata_path = root / "metadata.tsv"
            de.to_csv(de_path, sep="\t", index=False)
            run_summary.to_csv(summary_path, sep="\t", index=False)
            metadata.to_csv(metadata_path, sep="\t", index=False)

            audit, balance, matched = audit_genomewide_de(
                de_path,
                summary_path,
                metadata_path,
                min_case_donors=1,
                min_control_donors=2,
            )

        pd_row = audit[audit["contrast"].eq("pd_vs_healthy")].iloc[0]
        self.assertEqual(int(pd_row["is_sex_linked_initial_significant_rows"]), 1)
        self.assertEqual(int(pd_row["is_ribosomal_significant_rows"]), 1)
        self.assertIn("balance_status", balance.columns)
        self.assertTrue(bool(matched["ready_for_matched_de"].any()))


if __name__ == "__main__":
    unittest.main()
