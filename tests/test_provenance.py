from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.provenance import command_manifest_table, output_provenance_table


class ProvenanceTests(unittest.TestCase):
    def test_output_provenance_marks_existing_and_missing_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "out.txt").write_text("hello\n", encoding="utf-8")
            (root / "model").mkdir()
            (root / "model" / "model.pt").write_text("weights\n", encoding="utf-8")
            manifest = {
                "commands": {
                    "toy": {
                        "description": "Toy command",
                        "command": "make toy",
                        "inputs": ["in.txt"],
                        "outputs": ["out.txt", "model", "missing.txt"],
                    }
                }
            }

            commands = command_manifest_table(manifest)
            provenance = output_provenance_table(manifest, base_dir=root)

        self.assertEqual(commands.loc[0, "stage"], "toy")
        self.assertEqual(int(provenance["exists"].sum()), 2)
        self.assertIn("sha256", set(provenance["checksum_status"]))
        self.assertIn("directory", set(provenance["checksum_status"]))
        self.assertIn("missing", set(provenance["checksum_status"]))


if __name__ == "__main__":
    unittest.main()
