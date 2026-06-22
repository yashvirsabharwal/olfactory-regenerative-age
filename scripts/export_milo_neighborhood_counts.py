#!/usr/bin/env python3
"""Export neighborhood-by-donor count inputs for edgeR parity testing."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.neighborhood_parity import export_neighborhood_count_inputs
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memberships", required=True)
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--counts-out", required=True)
    parser.add_argument("--design-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--donor-query", default=None, help="Optional pandas query applied to donor metadata after healthy/age filtering.")
    parser.add_argument("--include-disease", action="store_true", help="Do not restrict donor metadata to healthy donors.")
    args = parser.parse_args()

    memberships = pd.read_csv(args.memberships, sep="\t", usecols=["neighborhood_id", "donor_id"])
    manifest = pd.read_csv(args.manifest, sep="\t")
    if not args.include_disease and "is_healthy" in manifest:
        manifest = manifest[manifest["is_healthy"].astype(str).str.lower().eq("true")].copy()
    manifest = manifest[manifest["age"].notna()].copy()
    if args.donor_query:
        try:
            manifest = manifest.query(args.donor_query).copy()
        except Exception as exc:
            raise SystemExit(f"Invalid --donor-query `{args.donor_query}`: {exc}") from exc
    counts, design, summary = export_neighborhood_count_inputs(memberships, manifest)
    counts.to_csv(ensure_parent(args.counts_out), sep="\t", index=False)
    design.to_csv(ensure_parent(args.design_out), sep="\t", index=False)
    summary.to_csv(ensure_parent(args.summary_out), sep="\t", index=False)
    print(f"Wrote neighborhood counts: {args.counts_out} ({counts.shape[0]} neighborhoods x {counts.shape[1] - 1} donors)")
    print(f"Wrote neighborhood design: {args.design_out} ({design.shape[0]} donors)")
    print(f"Wrote neighborhood export summary: {args.summary_out}")


if __name__ == "__main__":
    main()
