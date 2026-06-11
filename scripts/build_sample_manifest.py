#!/usr/bin/env python3
"""Build a donor/sample manifest from Gateway H5AD metadata."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config, project_path
from ora.io import load_obs
from ora.metadata import build_manifest, resolve_columns, summarize_cohort
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--summary-out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    h5ad_path = project_path(args.h5ad or config["source"]["h5ad_path"])
    obs = load_obs(h5ad_path)
    columns = resolve_columns(list(obs.columns), config)
    manifest = build_manifest(obs, config, columns)
    summary = summarize_cohort(manifest)

    outputs = config.get("outputs", {})
    manifest_path = args.out or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv")
    summary_path = args.summary_out or outputs.get("cohort_summary_tsv", "results/tables/cohort_summary.tsv")
    manifest.to_csv(ensure_parent(manifest_path), sep="\t", index=False)
    summary.to_csv(ensure_parent(summary_path), sep="\t", index=False)
    print(f"Wrote manifest: {manifest_path} ({manifest.shape[0]} donor/sample rows)")
    print(f"Wrote cohort summary: {summary_path}")


if __name__ == "__main__":
    main()
