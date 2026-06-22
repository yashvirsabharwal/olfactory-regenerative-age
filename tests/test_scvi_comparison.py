import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.scvi_comparison import compare_scvi_validation_runs


class ScviComparisonTests(unittest.TestCase):
    def test_compare_scvi_validation_runs_builds_claim_gates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            full = tmp / "full.tsv"
            seed = tmp / "seed.tsv"
            full.write_text(
                "\n".join(
                    [
                        "check\tstatus\tdetail\trecommendation",
                        "pilot_h5ad\tok\t4028275 cells x 3000 HVGs\tok",
                        "embedding_dimensions\tok\tX_scvi:(4028275, 10)\tok",
                        "neighbor_label_purity__fine_celltype\tok\tmean_same_label=0.88;k=16\tok",
                        "neighbor_label_purity__coarse_celltype\tok\tmean_same_label=0.98;k=16\tok",
                        "neighbor_mixing_entropy__flex_version\tok\tnormalized_entropy=0.66;levels=2\tok",
                        "neighbor_mixing_entropy__device_guided\tok\tnormalized_entropy=0.51;levels=2\tok",
                        "neighbor_mixing_entropy__condition\tok\tnormalized_entropy=0.31;levels=3\tok",
                        "neighbor_mixing_entropy__sex\tok\tnormalized_entropy=0.54;levels=3\tok",
                        "marker_continuity__immature_osn\tok\tpresent_genes=3;top_label=Fully_mature_mOSN;top_decile_enrichment=7.10\tok",
                        "marker_continuity__immune\tlimited\tpresent_genes=3;top_label=Multiciliated;top_decile_enrichment=1.00\treview",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            seed.write_text(
                "\n".join(
                    [
                        "check\tstatus\tdetail\trecommendation",
                        "pilot_h5ad\tok\t250000 cells x 3003 HVGs\tok",
                        "embedding_dimensions\tok\tX_scvi:(250000, 10)\tok",
                        "neighbor_label_purity__fine_celltype\tok\tmean_same_label=0.72;k=16\tok",
                        "neighbor_label_purity__coarse_celltype\tok\tmean_same_label=0.94;k=16\tok",
                        "neighbor_mixing_entropy__flex_version\tok\tnormalized_entropy=0.76;levels=2\tok",
                        "neighbor_mixing_entropy__device_guided\tok\tnormalized_entropy=0.59;levels=2\tok",
                        "neighbor_mixing_entropy__condition\tok\tnormalized_entropy=0.51;levels=3\tok",
                        "neighbor_mixing_entropy__sex\tok\tnormalized_entropy=0.61;levels=3\tok",
                        "marker_continuity__immature_osn\tok\tpresent_genes=3;top_label=Early_iOSN;top_decile_enrichment=5.09\tok",
                        "marker_continuity__immune\tok\tpresent_genes=3;top_label=matureDC;top_decile_enrichment=5.11\tok",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            summary, markers, note = compare_scvi_validation_runs(
                {"full_4m_reduced": full, "stratified_250k_seed13": seed}
            )

        primary = summary[summary["model"].eq("full_4m_reduced")].iloc[0]
        immune = markers[markers["marker"].eq("immune")].iloc[0]
        immature = markers[markers["marker"].eq("immature_osn")].iloc[0]
        self.assertEqual(primary["claim_gate"], "primary_with_technical_caveat")
        self.assertEqual(immune["claim_gate"], "guarded")
        self.assertEqual(immature["claim_gate"], "guarded")
        self.assertIn("full 4M reduced scVI model represents 4,028,275 cells", note)


if __name__ == "__main__":
    unittest.main()
