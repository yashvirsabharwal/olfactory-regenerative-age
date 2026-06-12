#!/usr/bin/env python3
"""Run covariate-adjusted donor-level tests on targeted pseudobulk counts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.pseudobulk import DEFAULT_CONTRASTS, DEFAULT_COVARIATES, parse_contrasts, run_covariate_pseudobulk_de
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--counts", default=None)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--coverage", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--contrasts", nargs="+", default=list(DEFAULT_CONTRASTS))
    parser.add_argument("--covariates", nargs="+", default=list(DEFAULT_COVARIATES))
    parser.add_argument("--min-donors", type=int, default=3)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    counts_path = args.counts or outputs.get("pseudobulk_counts_tsv", "data/processed/pseudobulk_counts.tsv.gz")
    metadata_path = args.metadata or outputs.get("pseudobulk_metadata_tsv", "data/processed/pseudobulk_metadata.tsv")
    manifest_path = args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv")
    coverage_path = args.coverage or outputs.get("pseudobulk_gene_coverage_tsv", "results/tables/pseudobulk_gene_coverage.tsv")
    out = args.out or outputs.get("pseudobulk_covariate_de_tsv", "results/tables/pseudobulk_covariate_de.tsv")

    counts = pd.read_csv(counts_path, sep="\t")
    metadata = pd.read_csv(metadata_path, sep="\t")
    manifest = pd.read_csv(manifest_path, sep="\t")
    genes = _genes_from_coverage(coverage_path) or sorted(counts["gene"].dropna().astype(str).unique().tolist())
    result = run_covariate_pseudobulk_de(
        counts,
        metadata,
        manifest,
        genes=genes,
        contrasts=parse_contrasts(args.contrasts),
        covariates=args.covariates,
        min_donors=args.min_donors,
    )
    result.to_csv(ensure_parent(out), sep="\t", index=False)
    ok_tests = int(result["status"].eq("ok").sum()) if "status" in result else 0
    print(f"Wrote covariate-adjusted pseudobulk DE: {out} ({ok_tests} ok tests)")


def _genes_from_coverage(path: str | Path | None) -> list[str]:
    if not path:
        return []
    candidate = Path(path)
    if not candidate.exists():
        return []
    coverage = pd.read_csv(candidate, sep="\t")
    if coverage.empty or "present_genes" not in coverage:
        return []
    genes: list[str] = []
    for value in coverage["present_genes"].dropna().astype(str):
        genes.extend(gene for gene in value.split(",") if gene)
    return list(dict.fromkeys(genes))


if __name__ == "__main__":
    main()
