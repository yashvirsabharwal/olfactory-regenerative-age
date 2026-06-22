import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.publication_tables import build_publication_tables, render_publication_table_index


class PublicationTablesTests(unittest.TestCase):
    def test_build_publication_tables_compacts_core_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "cohort_summary.tsv").write_text(
                "cohort\tdonors\tcells\tsamples\nhealthy\t2\t100\t2\nad\t1\t20\t1\n", encoding="utf-8"
            )
            (root / "ora_model_card.tsv").write_text(
                "model\tfeature_set\trole\tn\trepeats\tmae_mean\tmae_ci_low\tmae_ci_high\tspearman_r_mean\tcalibration_slope\tpermutation_p_mae\tlimitations\n"
                "xgboost\tcomposition\tpreferred_benchmark\t2\t1\t1.2\t1.0\t1.4\t0.5\t0.2\t0.01\ttoy\n",
                encoding="utf-8",
            )
            (root / "external_validation_evidence.tsv").write_text(
                "dataset_id\taccession\tevidence_type\tfeature_level\treadiness_class\tvalidation_strength\tn_samples\tn_donors\tsupports_primary_claim\tlimitation\tnext_action\n"
                "oliva\tGSE\tconfigured\tsingle_cell\tready\tsmall_n\t6\t6\tno\ttoy\tmap\n",
                encoding="utf-8",
            )
            (root / "scvi_embedding_claim_gates.tsv").write_text(
                "model\trole\tcells\tfine_label_purity\tlatent_dimensions\tclaim_gate\n"
                "full_4m_reduced\tprimary\t10\t0.9\t10\tprimary\n",
                encoding="utf-8",
            )
            (root / "milo_full_4m_lineage_summary.tsv").write_text(
                "metric\tvalue\tdetail\nage_fdr_lt_0_10\t3\tx\nneighborhoods_tested\t10\tx\n",
                encoding="utf-8",
            )
            (root / "pseudobulk_genomewide_de_audit.tsv").write_text(
                "contrast\ttested_rows\tsignificant_rows\tis_sex_linked_initial_significant_rows\tis_hemoglobin_significant_rows\tis_immunoglobulin_significant_rows\n"
                "ad_vs_healthy\t10\t2\t0\t0\t1\n",
                encoding="utf-8",
            )
            (root / "ndd_ora_projection_summary.tsv").write_text(
                "model\tdisease_group\tdonors\tmean_age\tmean_ora\tmean_oraa\tsd_oraa\n"
                "xgboost\tad\t1\t70\t60\t-10\t1\n",
                encoding="utf-8",
            )

            tables = build_publication_tables(root)
            index = render_publication_table_index(tables)

        self.assertIn("manuscript_table_model_card", tables)
        self.assertEqual(tables["manuscript_table_cohort"].shape[0], 2)
        self.assertIn("python_lineage_full", tables["manuscript_table_latent_neighborhood_gates"]["analysis"].tolist())
        self.assertIn("Publication Table Bundle", index)


if __name__ == "__main__":
    unittest.main()
