from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.release import build_release_manifest, render_release_manifest_markdown


class ReleaseManifestTests(unittest.TestCase):
    def test_build_release_manifest_records_status_and_overrides(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "results" / "table.tsv"
            output.parent.mkdir(parents=True)
            output.write_text("a\tb\n1\t2\n")
            command_manifest = {
                "commands": {
                    "toy": {
                        "command": "make toy",
                        "inputs": ["input.tsv"],
                        "outputs": ["results/table.tsv", "missing.tsv"],
                    }
                }
            }

            manifest = build_release_manifest(
                command_manifest,
                base_dir=root,
                extra_artifacts=[
                    {
                        "path": "results/table.tsv",
                        "category": "locally_generated",
                        "required_for_review": True,
                        "notes": "Toy output.",
                    },
                    {
                        "path": "source/raw.h5ad",
                        "category": "source_data",
                        "required_for_review": True,
                    },
                ],
            )

            table = manifest[manifest["path"].eq("results/table.tsv")].iloc[0]
            missing = manifest[manifest["path"].eq("missing.tsv")].iloc[0]
            source = manifest[manifest["path"].eq("source/raw.h5ad")].iloc[0]
            self.assertEqual(table["artifact_status"], "present")
            self.assertEqual(table["category"], "locally_generated")
            self.assertEqual(table["generating_command"], "make toy")
            self.assertEqual(missing["artifact_status"], "missing")
            self.assertEqual(source["category"], "source_data")
            self.assertEqual(source["artifact_status"], "missing")

            markdown = render_release_manifest_markdown(manifest)
            self.assertIn("Required Items Needing Attention", markdown)
            self.assertIn("missing.tsv", markdown)

    def test_external_archives_use_larger_checksum_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            external = root / "data" / "external" / "raw.tar"
            source = root / "data" / "raw" / "gateway.h5ad"
            external.parent.mkdir(parents=True)
            source.parent.mkdir(parents=True)
            external.write_bytes(b"external")
            source.write_bytes(b"source")
            manifest = build_release_manifest(
                {"commands": {}},
                base_dir=root,
                checksum_max_bytes=1,
                extra_artifacts=[
                    {
                        "path": "data/external/raw.tar",
                        "category": "external_archive",
                        "required_for_review": True,
                    },
                    {
                        "path": "data/raw/gateway.h5ad",
                        "category": "source_data",
                        "required_for_review": True,
                    },
                ],
            )

        external_row = manifest[manifest["path"].eq("data/external/raw.tar")].iloc[0]
        source_row = manifest[manifest["path"].eq("data/raw/gateway.h5ad")].iloc[0]
        self.assertEqual(external_row["checksum_status"], "sha256")
        self.assertEqual(source_row["checksum_status"], "skipped_large_file")


if __name__ == "__main__":
    unittest.main()
