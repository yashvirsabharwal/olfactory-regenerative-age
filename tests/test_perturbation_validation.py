import unittest

from ora.perturbation_validation import (
    build_minimum_experiment_table,
    build_perturbation_candidate_matrix,
    build_perturbation_search_log,
)


class PerturbationValidationTests(unittest.TestCase):
    def test_candidate_matrix_ranks_and_flags_directness(self):
        config = {
            "candidates": [
                {
                    "accession": "GSE2",
                    "source_url": "https://example.org/2",
                    "title": "Context",
                    "tissue_model": "nasal ALI",
                    "assay": "bulk RNA-seq",
                    "perturbations": "IFN",
                    "mechanism_match": "IFN",
                    "data_access": "public",
                    "adapter_decision": "context",
                    "priority": 2,
                    "directness": "nasal_ali_mechanism_context",
                    "limitations": "not olfactory or aging",
                    "recommended_next_step": "score modules",
                },
                {
                    "accession": "GSE1",
                    "source_url": "https://example.org/1",
                    "title": "Organoid",
                    "tissue_model": "olfactory organoid",
                    "assay": "scRNA-seq",
                    "perturbations": "aging cytokine model",
                    "mechanism_match": "aging",
                    "data_access": "public",
                    "adapter_decision": "adapter",
                    "priority": 1,
                    "directness": "olfactory_relevant_organoid_context",
                    "limitations": "none",
                    "recommended_next_step": "adapt",
                },
            ]
        }

        matrix = build_perturbation_candidate_matrix(config)

        self.assertEqual(matrix["accession"].tolist(), ["GSE1", "GSE2"])
        self.assertTrue(bool(matrix.loc[0, "usable_for_direct_ora_mechanism"]))
        self.assertFalse(bool(matrix.loc[1, "usable_for_direct_ora_mechanism"]))

    def test_search_log_and_experiment_table_flatten(self):
        config = {
            "search_log": [
                {
                    "search_date": "2026-06-25",
                    "database_or_resource": "NCBI GEO",
                    "query_or_filter": "olfactory organoid",
                    "result_summary": "candidate found",
                }
            ],
            "minimum_experiment": [
                {
                    "experiment_id": "adult_oe",
                    "model": "adult OE organoids",
                    "perturbations": "TNF",
                    "timepoints": "24h",
                    "readout": "scRNA-seq",
                    "target_n": "6 donors",
                }
            ],
        }

        log = build_perturbation_search_log(config)
        experiment = build_minimum_experiment_table(config)

        self.assertEqual(log.loc[0, "database_or_resource"], "NCBI GEO")
        self.assertEqual(experiment.loc[0, "experiment_id"], "adult_oe")
        self.assertIn("Causal support", experiment.loc[0, "claim_boundary"])


if __name__ == "__main__":
    unittest.main()
