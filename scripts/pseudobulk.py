#!/usr/bin/env python3
"""Pseudobulk export, QC, DE summary, and audit command group."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.genomewide_de import audit_genomewide_de, summarize_genomewide_de
from ora.genomewide_qc import summarize_genomewide_pseudobulk
from ora.pseudobulk import (
    DEFAULT_CONTRASTS,
    DEFAULT_COVARIATES,
    DEFAULT_GROUPBY,
    aggregate_targeted_pseudobulk_h5ad,
    export_genomewide_pseudobulk_h5ad,
    genes_from_gene_sets,
    parse_contrasts,
    run_covariate_pseudobulk_de,
)
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_targeted(subparsers)
    _add_export_genomewide(subparsers)
    _add_qc(subparsers)
    _add_de_summary(subparsers)
    _add_de_audit(subparsers)
    _add_covariate_de(subparsers)
    args = parser.parse_args()
    args.func(args)


def _add_targeted(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("targeted")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--gene-sets", default="configs/gene_sets.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--genes", nargs="*", default=None)
    parser.add_argument("--groupby", nargs="+", default=list(DEFAULT_GROUPBY))
    parser.add_argument("--contrasts", nargs="+", default=list(DEFAULT_CONTRASTS))
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--min-donors", type=int, default=3)
    parser.add_argument("--apply-qc", action="store_true")
    parser.add_argument("--counts-out", default=None)
    parser.add_argument("--metadata-out", default=None)
    parser.add_argument("--coverage-out", default=None)
    parser.add_argument("--de-out", default=None)
    parser.set_defaults(func=_targeted)


def _targeted(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    outputs = config.get("outputs", {})
    result = aggregate_targeted_pseudobulk_h5ad(
        args.h5ad or config.get("source", {}).get("h5ad_path", "data/raw/gateway.h5ad"),
        config,
        args.genes or genes_from_gene_sets(load_config(args.gene_sets)),
        groupby=args.groupby,
        chunk_size=args.chunk_size,
        apply_qc=args.apply_qc,
        contrasts=parse_contrasts(args.contrasts),
        min_donors=args.min_donors,
    )
    counts_out = args.counts_out or outputs.get("pseudobulk_counts_tsv", "data/processed/pseudobulk_counts.tsv.gz")
    metadata_out = args.metadata_out or outputs.get("pseudobulk_metadata_tsv", "data/processed/pseudobulk_metadata.tsv")
    coverage_out = args.coverage_out or outputs.get("pseudobulk_gene_coverage_tsv", "results/tables/pseudobulk_gene_coverage.tsv")
    de_out = args.de_out or outputs.get("pseudobulk_de_tsv", "results/tables/pseudobulk_de.tsv")
    result.counts.to_csv(ensure_parent(counts_out), sep="\t", index=False)
    result.metadata.to_csv(ensure_parent(metadata_out), sep="\t", index=False)
    result.coverage.to_csv(ensure_parent(coverage_out), sep="\t", index=False)
    result.de.to_csv(ensure_parent(de_out), sep="\t", index=False)
    print(f"Wrote pseudobulk DE: {de_out} ({result.de.shape[0]} rows)")


def _add_export_genomewide(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("export-genomewide")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--groupby", nargs="+", default=list(DEFAULT_GROUPBY))
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--gene-chunk-size", type=int, default=500)
    parser.add_argument("--min-cells-per-group", type=int, default=10)
    parser.add_argument("--min-donors-per-cell-state", type=int, default=3)
    parser.add_argument("--apply-qc", action="store_true")
    parser.add_argument("--counts-out", default=None)
    parser.add_argument("--metadata-out", default=None)
    parser.add_argument("--genes-out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.set_defaults(func=_export_genomewide)


def _export_genomewide(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    outputs = config.get("outputs", {})
    counts_out = args.counts_out or outputs.get("pseudobulk_genomewide_counts_tsv", "data/processed/pseudobulk_genomewide_counts.tsv.gz")
    metadata_out = args.metadata_out or outputs.get("pseudobulk_genomewide_metadata_tsv", "data/processed/pseudobulk_genomewide_metadata.tsv")
    genes_out = args.genes_out or outputs.get("pseudobulk_genomewide_genes_tsv", "data/processed/pseudobulk_genomewide_genes.tsv")
    summary_out = args.summary_out or outputs.get("pseudobulk_genomewide_summary_tsv", "results/tables/pseudobulk_genomewide_summary.tsv")
    result = export_genomewide_pseudobulk_h5ad(
        args.h5ad or config.get("source", {}).get("h5ad_path", "data/raw/gateway.h5ad"),
        config,
        counts_out=counts_out,
        metadata_out=metadata_out,
        genes_out=genes_out,
        groupby=args.groupby,
        chunk_size=args.chunk_size,
        gene_chunk_size=args.gene_chunk_size,
        apply_qc=args.apply_qc,
        min_cells_per_group=args.min_cells_per_group,
        min_donors_per_cell_state=args.min_donors_per_cell_state,
    )
    result.summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    print(f"Wrote genome-wide pseudobulk export summary: {summary_out}")


def _add_qc(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("qc")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--counts", default=None)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--genes", default=None)
    parser.add_argument("--chunksize", type=int, default=500)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--gene-qc-out", default=None)
    parser.add_argument("--group-qc-out", default=None)
    parser.add_argument("--disease-summary-out", default=None)
    parser.add_argument("--cell-state-summary-out", default=None)
    parser.set_defaults(func=_qc)


def _qc(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    outputs = config.get("outputs", {})
    result = summarize_genomewide_pseudobulk(
        args.counts or outputs.get("pseudobulk_genomewide_counts_tsv", "data/processed/pseudobulk_genomewide_counts.tsv.gz"),
        args.metadata or outputs.get("pseudobulk_genomewide_metadata_tsv", "data/processed/pseudobulk_genomewide_metadata.tsv"),
        args.genes or outputs.get("pseudobulk_genomewide_genes_tsv", "data/processed/pseudobulk_genomewide_genes.tsv"),
        chunksize=args.chunksize,
    )
    result.summary.to_csv(ensure_parent(args.summary_out or outputs.get("pseudobulk_genomewide_qc_summary_tsv", "results/tables/pseudobulk_genomewide_qc_summary.tsv")), sep="\t", index=False)
    result.gene_qc.to_csv(ensure_parent(args.gene_qc_out or outputs.get("pseudobulk_genomewide_gene_qc_tsv", "results/tables/pseudobulk_genomewide_gene_qc.tsv")), sep="\t", index=False)
    result.group_qc.to_csv(ensure_parent(args.group_qc_out or outputs.get("pseudobulk_genomewide_group_qc_tsv", "results/tables/pseudobulk_genomewide_group_qc.tsv")), sep="\t", index=False)
    result.disease_summary.to_csv(ensure_parent(args.disease_summary_out or outputs.get("pseudobulk_genomewide_disease_summary_tsv", "results/tables/pseudobulk_genomewide_disease_summary.tsv")), sep="\t", index=False)
    result.cell_state_summary.to_csv(ensure_parent(args.cell_state_summary_out or outputs.get("pseudobulk_genomewide_cell_state_summary_tsv", "results/tables/pseudobulk_genomewide_cell_state_summary.tsv")), sep="\t", index=False)
    print("Wrote genome-wide pseudobulk QC outputs")


def _add_de_summary(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("de-summary")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--de", default=None)
    parser.add_argument("--run-summary", default=None)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--top-hits-out", default=None)
    parser.set_defaults(func=_de_summary)


def _de_summary(args: argparse.Namespace) -> None:
    outputs = load_config(args.config).get("outputs", {})
    summary, top_hits = summarize_genomewide_de(
        args.de or outputs.get("pseudobulk_genomewide_edger_tsv", "results/tables/pseudobulk_genomewide_edger.tsv.gz"),
        args.run_summary or outputs.get("pseudobulk_genomewide_edger_summary_tsv", "results/tables/pseudobulk_genomewide_edger_summary.tsv"),
        fdr_threshold=args.fdr_threshold,
        top_n=args.top_n,
    )
    summary_out = args.summary_out or outputs.get("pseudobulk_genomewide_de_summary_tsv", "results/tables/pseudobulk_genomewide_de_summary.tsv")
    top_hits_out = args.top_hits_out or outputs.get("pseudobulk_genomewide_de_top_hits_tsv", "results/tables/pseudobulk_genomewide_de_top_hits.tsv")
    summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    top_hits.to_csv(ensure_parent(top_hits_out), sep="\t", index=False)
    print(f"Wrote genome-wide DE summary: {summary_out}")


def _add_de_audit(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("de-audit")
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
    parser.set_defaults(func=_de_audit)


def _de_audit(args: argparse.Namespace) -> None:
    outputs = load_config(args.config).get("outputs", {})
    audit, donor_balance, matched = audit_genomewide_de(
        args.de or outputs.get("pseudobulk_genomewide_edger_tsv", "results/tables/pseudobulk_genomewide_edger.tsv.gz"),
        args.run_summary or outputs.get("pseudobulk_genomewide_edger_summary_tsv", "results/tables/pseudobulk_genomewide_edger_summary.tsv"),
        args.metadata or outputs.get("pseudobulk_genomewide_metadata_tsv", "data/processed/pseudobulk_genomewide_metadata.tsv"),
        args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv"),
        fdr_threshold=args.fdr_threshold,
        min_case_donors=args.min_case_donors,
        min_control_donors=args.min_control_donors,
    )
    audit.to_csv(ensure_parent(args.audit_out or outputs.get("pseudobulk_genomewide_de_audit_tsv", "results/tables/pseudobulk_genomewide_de_audit.tsv")), sep="\t", index=False)
    donor_balance.to_csv(ensure_parent(args.donor_balance_out or outputs.get("pseudobulk_genomewide_donor_balance_tsv", "results/tables/pseudobulk_genomewide_donor_balance.tsv")), sep="\t", index=False)
    matched.to_csv(ensure_parent(args.matched_feasibility_out or outputs.get("pseudobulk_genomewide_matched_feasibility_tsv", "results/tables/pseudobulk_genomewide_matched_feasibility.tsv")), sep="\t", index=False)
    print(f"Wrote genome-wide DE audit ({audit.shape[0]} rows)")


def _add_covariate_de(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("covariate-de")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--counts", default=None)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--coverage", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--contrasts", nargs="+", default=list(DEFAULT_CONTRASTS))
    parser.add_argument("--covariates", nargs="+", default=list(DEFAULT_COVARIATES))
    parser.add_argument("--min-donors", type=int, default=3)
    parser.set_defaults(func=_covariate_de)


def _covariate_de(args: argparse.Namespace) -> None:
    outputs = load_config(args.config).get("outputs", {})
    counts_path = args.counts or outputs.get("pseudobulk_counts_tsv", "data/processed/pseudobulk_counts.tsv.gz")
    metadata_path = args.metadata or outputs.get("pseudobulk_metadata_tsv", "data/processed/pseudobulk_metadata.tsv")
    manifest_path = args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv")
    coverage_path = args.coverage or outputs.get("pseudobulk_gene_coverage_tsv", "results/tables/pseudobulk_gene_coverage.tsv")
    out = args.out or outputs.get("pseudobulk_covariate_de_tsv", "results/tables/pseudobulk_covariate_de.tsv")
    counts = pd.read_csv(counts_path, sep="\t")
    result = run_covariate_pseudobulk_de(
        counts,
        pd.read_csv(metadata_path, sep="\t"),
        pd.read_csv(manifest_path, sep="\t"),
        genes=_genes_from_coverage(coverage_path) or sorted(counts["gene"].dropna().astype(str).unique().tolist()),
        contrasts=parse_contrasts(args.contrasts),
        covariates=args.covariates,
        min_donors=args.min_donors,
    )
    result.to_csv(ensure_parent(out), sep="\t", index=False)
    print(f"Wrote covariate-adjusted pseudobulk DE: {out} ({result.shape[0]} rows)")


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
