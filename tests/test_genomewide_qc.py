from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.genomewide_qc import summarize_genomewide_pseudobulk


class GenomewideQCTests(unittest.TestCase):
    def test_summarize_genomewide_pseudobulk_tracks_alignment_and_detection(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            counts = root / "counts.tsv"
            metadata = root / "metadata.tsv"
            genes = root / "genes.tsv"
            counts.write_text(
                "\t".join(["gene_id", "gene_symbol", "gene_index", "PB1", "PB2", "PB3"])
                + "\n"
                + "ENSG1\tTP63\t0\t10\t0\t5\n"
                + "ENSG2\tOMP\t1\t0\t0\t20\n"
                + "ENSG3\tSNCA\t2\t2\t2\t2\n",
                encoding="utf-8",
            )
            metadata.write_text(
                "pseudobulk_id\tdonor_id\tdisease_group\tfine_cell_type\tn_cells\tsum_n_counts\n"
                "PB1\td1\thealthy\tqHBC\t10\t12\n"
                "PB2\td2\thealthy\tqHBC\t8\t2\n"
                "PB3\td3\tad\tqHBC\t12\t27\n",
                encoding="utf-8",
            )
            genes.write_text(
                "gene_index\tgene_id\tgene_symbol\n0\tENSG1\tTP63\n1\tENSG2\tOMP\n2\tENSG3\tSNCA\n",
                encoding="utf-8",
            )

            result = summarize_genomewide_pseudobulk(counts, metadata, genes, chunksize=2)

        self.assertTrue(bool(result.summary.loc[0, "matrix_columns_match_metadata"]))
        self.assertEqual(int(result.summary.loc[0, "n_genes"]), 3)
        self.assertEqual(int(result.summary.loc[0, "n_groups"]), 3)
        self.assertEqual(result.gene_qc.iloc[0]["gene_symbol"], "OMP")
        self.assertEqual(int(result.gene_qc.loc[result.gene_qc["gene_symbol"].eq("TP63"), "detected_groups"].iloc[0]), 2)
        self.assertEqual(result.group_qc["matrix_total_count"].tolist(), [12.0, 2.0, 27.0])
        self.assertEqual(set(result.disease_summary["disease_group"]), {"healthy", "ad"})


if __name__ == "__main__":
    unittest.main()
