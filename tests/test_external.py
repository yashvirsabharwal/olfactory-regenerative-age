from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.external import (
    external_dataset_summary,
    feature_matrix_contract_summary,
    parse_published_gene_lists,
    published_gene_list_coverage,
)


class ExternalValidationTests(unittest.TestCase):
    def test_external_dataset_summary_tracks_file_readiness(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            feature_matrix = root / "features.tsv"
            feature_matrix.write_text("donor_id\tage\nD1\t40\n", encoding="utf-8")
            config = {
                "datasets": {
                    "toy": {
                        "title": "Toy validation",
                        "status": "configured",
                        "validation_use": "unit test",
                        "species": "human",
                        "tissue": "olfactory",
                        "disease_context": ["healthy"],
                        "expected_level": "donor_feature_matrix",
                        "required_files": {"feature_matrix": "features.tsv", "expression": None, "metadata": None},
                    }
                }
            }

            summary = external_dataset_summary(config, base_dir=root)

        self.assertEqual(summary.loc[0, "dataset_id"], "toy")
        self.assertTrue(bool(summary.loc[0, "ready_for_feature_validation"]))
        self.assertFalse(bool(summary.loc[0, "ready_for_raw_adapter"]))
        self.assertIn("metadata", summary.loc[0, "files_missing"])

    def test_published_gene_list_coverage_uses_feature_name(self):
        config = {
            "published_gene_lists": {
                "regen": {"description": "Regeneration", "genes": ["TP63", "OMP", "MISSING"]}
            }
        }
        var = pd.DataFrame({"feature_name": ["TP63", "OMP"]}, index=["ENSG1", "ENSG2"])

        coverage = published_gene_list_coverage(config, var, var.index, ["feature_name"])

        self.assertEqual(coverage.loc[0, "gene_list"], "regen")
        self.assertEqual(int(coverage.loc[0, "n_present"]), 2)
        self.assertEqual(coverage.loc[0, "missing_genes"], "MISSING")

    def test_parse_gene_lists_and_contract_summary(self):
        config = {
            "published_gene_lists": {"aging": {"genes": ["CDKN1A", "CDKN2A"]}},
            "feature_matrix_contract": {
                "required_columns": ["donor_id", "age"],
                "optional_covariates": ["sex"],
                "accepted_feature_prefixes": ["prop__"],
            },
        }

        gene_lists = parse_published_gene_lists(config)
        contract = feature_matrix_contract_summary(config)

        self.assertEqual(gene_lists[0].genes, ("CDKN1A", "CDKN2A"))
        self.assertEqual(set(contract["kind"]), {"required_column", "optional_covariate", "accepted_feature_prefix"})


if __name__ == "__main__":
    unittest.main()
