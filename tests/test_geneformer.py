import unittest

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ora.geneformer import (
    geneformer_gene_arrays,
    geneformer_gene_coverage,
    select_geneformer_cells,
    tokenize_geneformer_matrix,
)


class GeneformerHelperTests(unittest.TestCase):
    def test_gene_arrays_and_coverage_require_token_and_median_overlap(self):
        genes = np.array(["ENSG1", "ENSG2", "ENSG3"])
        token_dict = {"<pad>": 0, "ENSG1": 11, "ENSG2": 22}
        median_dict = {"ENSG1": 2.0, "ENSG3": 3.0}

        indices, tokens, medians = geneformer_gene_arrays(genes, token_dict, median_dict)
        coverage = geneformer_gene_coverage(
            genes,
            token_dict,
            median_dict,
            model_family="geneformer",
            checkpoint="test",
        )

        self.assertEqual(indices.tolist(), [0])
        self.assertEqual(tokens.tolist(), [11])
        self.assertEqual(medians.tolist(), [2.0])
        self.assertEqual(int(coverage.loc[0, "usable_genes"]), 1)
        self.assertAlmostEqual(float(coverage.loc[0, "usable_gene_fraction"]), 1 / 3)
        self.assertEqual(coverage.loc[0, "status"], "low_gene_coverage")

    def test_select_geneformer_cells_is_deterministic_and_bounded(self):
        obs = pd.DataFrame(
            {
                "donor_id": [f"D{i % 3}" for i in range(30)],
                "fine_celltype": ["HBC" if i % 2 else "iOSN" for i in range(30)],
            }
        )

        first = select_geneformer_cells(obs, max_cells=11, seed=19)
        second = select_geneformer_cells(obs, max_cells=11, seed=19)

        self.assertEqual(first.tolist(), second.tolist())
        self.assertEqual(first.size, 11)
        self.assertTrue(np.all(first[:-1] <= first[1:]))

    def test_tokenize_geneformer_matrix_ranks_normalized_nonzero_genes(self):
        matrix = sp.csr_matrix(
            [
                [5.0, 4.0, 1.0],
                [0.0, 0.0, 0.0],
                [2.0, 8.0, 1.0],
            ]
        )
        gene_tokens = np.array([10, 20, 30])
        gene_medians = np.array([1.0, 2.0, 1.0])

        sequences, valid_rows = tokenize_geneformer_matrix(
            matrix,
            n_counts=np.array([10.0, 1.0, 11.0]),
            gene_tokens=gene_tokens,
            gene_medians=gene_medians,
            target_sum=10.0,
            max_length=2,
        )

        self.assertEqual(valid_rows, [0, 2])
        self.assertEqual([seq.tolist() for seq in sequences], [[10, 20], [20, 10]])


if __name__ == "__main__":
    unittest.main()
