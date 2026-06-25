from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.archive import build_archive_review_package, render_archive_review_markdown


class ArchiveReviewPackageTests(unittest.TestCase):
    def test_archive_review_flags_deferred_and_large_artifacts(self):
        manifest = pd.DataFrame(
            [
                {
                    "path": "README.md",
                    "category": "git_tracked",
                    "artifact_status": "present",
                    "required_for_review": True,
                    "size_bytes": 100,
                    "sha256": "abc",
                    "checksum_status": "sha256",
                    "archive_uri": "",
                    "notes": "Code readme.",
                },
                {
                    "path": "data/processed/full.h5ad",
                    "category": "large_artifact",
                    "artifact_status": "deferred",
                    "required_for_review": True,
                    "size_bytes": 10,
                    "sha256": "",
                    "checksum_status": "skipped_large_file",
                    "archive_uri": "",
                    "notes": "Large latent artifact.",
                },
                {
                    "path": "results/table.tsv",
                    "category": "locally_generated",
                    "artifact_status": "present",
                    "required_for_review": True,
                    "size_bytes": 10,
                    "sha256": "def",
                    "checksum_status": "sha256",
                    "archive_uri": "https://example.org/archive",
                    "notes": "Archived table.",
                },
            ]
        )

        package = build_archive_review_package(manifest)
        large = package[package["path"].eq("data/processed/full.h5ad")].iloc[0]
        archived = package[package["path"].eq("results/table.tsv")].iloc[0]
        self.assertIn("deferred", large["blocking_issue"])
        self.assertEqual(archived["blocking_issue"], "")
        markdown = render_archive_review_markdown(package)
        self.assertIn("Blocking Items", markdown)
        self.assertIn("data/processed/full.h5ad", markdown)


if __name__ == "__main__":
    unittest.main()
