#!/usr/bin/env python3
"""Audit genome-wide pseudobulk DE for donor balance and sentinel gene classes."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.genomewide_de import audit_genomewide_de
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--de", default=None)
    parser.add_argument("--run-summary", default=None)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--min-case-donors", type=int, default=3)
    parser.add_argument("--min-control-donors", type=int, default=10)
    parser.add_argument("--audit-out", default=None)
    parser.add_argument("--donor-balance-out", default=None)
    parser.add_argument("--matched-feasibility-out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    de_path = args.de or outputs.get("pseudobulk_genomewide_edger_tsv", "results/tables/pseudobulk_genomewide_edger.tsv.gz")
    run_summary_path = args.run_summary or outputs.get(
        "pseudobulk_genomewide_edger_summary_tsv",
        "results/tables/pseudobulk_genomewide_edger_summary.tsv",
    )
    metadata_path = args.metadata or outputs.get(
        "pseudobulk_genomewide_metadata_tsv",
        "data/processed/pseudobulk_genomewide_metadata.tsv",
    )
    manifest_path = args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv")
    audit_out = args.audit_out or outputs.get(
        "pseudobulk_genomewide_de_audit_tsv",
        "results/tables/pseudobulk_genomewide_de_audit.tsv",
    )
    donor_balance_out = args.donor_balance_out or outputs.get(
        "pseudobulk_genomewide_donor_balance_tsv",
        "results/tables/pseudobulk_genomewide_donor_balance.tsv",
    )
    matched_out = args.matched_feasibility_out or outputs.get(
        "pseudobulk_genomewide_matched_feasibility_tsv",
        "results/tables/pseudobulk_genomewide_matched_feasibility.tsv",
    )

    audit, donor_balance, matched = audit_genomewide_de(
        de_path,
        run_summary_path,
        metadata_path,
        manifest_path,
        fdr_threshold=args.fdr_threshold,
        min_case_donors=args.min_case_donors,
        min_control_donors=args.min_control_donors,
    )
    audit.to_csv(ensure_parent(audit_out), sep="\t", index=False)
    donor_balance.to_csv(ensure_parent(donor_balance_out), sep="\t", index=False)
    matched.to_csv(ensure_parent(matched_out), sep="\t", index=False)
    print(f"Wrote genome-wide DE audit: {audit_out} ({audit.shape[0]} rows)")
    print(f"Wrote genome-wide donor-balance audit: {donor_balance_out} ({donor_balance.shape[0]} rows)")
    print(f"Wrote matched DE feasibility: {matched_out} ({matched.shape[0]} rows)")


if __name__ == "__main__":
    main()
