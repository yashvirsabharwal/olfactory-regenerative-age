#!/usr/bin/env python3
"""Summarize age-bin robustness for Milo-style neighborhood memberships."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.neighborhood_age_bins import AgeBinConfig, summarize_neighborhood_age_bins
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memberships", required=True)
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--da-table", required=True)
    parser.add_argument("--run-name", default="neighborhood_run")
    parser.add_argument("--neighborhoods-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--donor-query", default=None, help="Optional pandas query applied to donor metadata after healthy/age filtering.")
    parser.add_argument("--include-disease", action="store_true", help="Do not restrict donor metadata to healthy donors.")
    parser.add_argument(
        "--age-bins",
        default="lt45:0:45,45_59:45:60,60_74:60:75,75_plus:75:inf",
        help="Comma-separated label:lower:upper bins. Upper bounds are exclusive; use inf for open-ended.",
    )
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
    da_table = pd.read_csv(args.da_table, sep="\t")
    config = AgeBinConfig(run_name=args.run_name, bins=_parse_bins(args.age_bins))
    neighborhoods, summary = summarize_neighborhood_age_bins(
        memberships,
        manifest,
        da_table=da_table,
        config=config,
    )
    neighborhoods.to_csv(ensure_parent(args.neighborhoods_out), sep="\t", index=False)
    summary.to_csv(ensure_parent(args.summary_out), sep="\t", index=False)
    print(f"Wrote age-bin neighborhood robustness: {args.neighborhoods_out} ({neighborhoods.shape[0]} rows)")
    print(f"Wrote age-bin summary: {args.summary_out} ({summary.shape[0]} rows)")


def _parse_bins(spec: str) -> tuple[tuple[str, float, float], ...]:
    bins = []
    for item in spec.split(","):
        if not item.strip():
            continue
        parts = item.split(":")
        if len(parts) != 3:
            raise SystemExit(f"Invalid age-bin spec `{item}`; expected label:lower:upper")
        label, lower, upper = parts
        upper_value = float("inf") if upper.lower() in {"inf", "infinity"} else float(upper)
        bins.append((label, float(lower), upper_value))
    if not bins:
        raise SystemExit("At least one age bin is required.")
    return tuple(bins)


if __name__ == "__main__":
    main()
