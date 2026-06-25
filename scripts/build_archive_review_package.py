#!/usr/bin/env python3
"""Build reviewer-access archive action tables from the release manifest."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.archive import build_archive_review_package, render_archive_review_markdown
from ora.config import load_config
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--release-manifest", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--markdown-out", default=None)
    args = parser.parse_args()

    outputs = load_config(args.config).get("outputs", {})
    release_manifest = args.release_manifest or outputs.get(
        "release_artifact_manifest_tsv",
        "results/reports/release_artifact_manifest.tsv",
    )
    out_path = args.out or outputs.get(
        "archive_review_package_tsv",
        "results/reports/archive_review_package.tsv",
    )
    markdown_path = args.markdown_out or outputs.get(
        "archive_review_package_md",
        "results/reports/archive_review_package.md",
    )
    package = build_archive_review_package(pd.read_csv(release_manifest, sep="\t"))
    package.to_csv(ensure_parent(out_path), sep="\t", index=False)
    ensure_parent(markdown_path).write_text(render_archive_review_markdown(package), encoding="utf-8")
    blockers = int(package["blocking_issue"].astype(str).ne("").sum()) if "blocking_issue" in package else 0
    print(f"Wrote archive review package: {out_path} ({package.shape[0]} rows; {blockers} blockers)")
    print(f"Wrote archive review summary: {markdown_path}")


if __name__ == "__main__":
    main()
