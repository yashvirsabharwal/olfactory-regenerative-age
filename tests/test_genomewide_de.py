from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.genomewide_de import summarize_genomewide_de


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


if __name__ == "__main__":
    unittest.main()
