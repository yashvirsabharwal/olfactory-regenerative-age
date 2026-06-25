#!/usr/bin/env python3
"""Run curated ligand-receptor niche signaling analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.niche_signaling import (
    build_niche_age_associations,
    build_niche_ora_associations,
    build_niche_priority_table,
    parse_niche_interactions,
    render_niche_signaling_report,
    score_niche_interactions,
    write_niche_signaling_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interaction-config", default="configs/niche_ligand_receptor_pairs.yaml")
    parser.add_argument("--counts", default="data/processed/pseudobulk_genomewide_counts.tsv.gz")
    parser.add_argument("--metadata", default="data/processed/pseudobulk_genomewide_metadata.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--ora-scores", default="results/tables/ora_augmented_candidate_repeated_cv_scores.tsv")
    parser.add_argument("--chunksize", type=int, default=1000)
    parser.add_argument("--donor-scores-out", default="results/tables/niche_ligand_receptor_donor_scores.tsv")
    parser.add_argument("--coverage-out", default="results/tables/niche_ligand_receptor_gene_coverage.tsv")
    parser.add_argument("--age-out", default="results/tables/niche_ligand_receptor_age_associations.tsv")
    parser.add_argument("--ora-out", default="results/tables/niche_ligand_receptor_ora_associations.tsv")
    parser.add_argument("--priority-out", default="results/tables/niche_driver_priority_table.tsv")
    parser.add_argument("--report-out", default="docs/niche_signaling_feasibility.md")
    parser.add_argument("--figure-pdf", default="results/figures/extended_data_niche_signaling.pdf")
    parser.add_argument("--figure-png", default="results/figures/extended_data_niche_signaling.png")
    args = parser.parse_args()

    interactions = parse_niche_interactions(load_config(args.interaction_config))
    metadata = pd.read_csv(args.metadata, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    ora_scores = pd.read_csv(args.ora_scores, sep="\t")

    donor_scores, coverage = score_niche_interactions(
        counts_path=args.counts,
        metadata=metadata,
        interactions=interactions,
        chunksize=args.chunksize,
    )
    age = build_niche_age_associations(
        donor_scores=donor_scores,
        manifest=manifest,
        interactions=interactions,
    )
    ora = build_niche_ora_associations(donor_scores=donor_scores, ora_scores=ora_scores)
    priority = build_niche_priority_table(
        interactions=interactions,
        coverage=coverage,
        age_associations=age,
        ora_associations=ora,
    )
    report = render_niche_signaling_report(
        interactions=interactions,
        coverage=coverage,
        priority=priority,
    )
    write_niche_signaling_outputs(
        donor_scores=donor_scores,
        coverage=coverage,
        age_associations=age,
        ora_associations=ora,
        priority=priority,
        report=report,
        donor_scores_out=args.donor_scores_out,
        coverage_out=args.coverage_out,
        age_out=args.age_out,
        ora_out=args.ora_out,
        priority_out=args.priority_out,
        report_out=args.report_out,
        figure_pdf=args.figure_pdf,
        figure_png=args.figure_png,
    )
    print(
        "Wrote niche signaling analysis: "
        f"{args.priority_out} ({priority.shape[0]} sender-receiver hypotheses), "
        f"{args.age_out} ({age.shape[0]} age rows)"
    )


if __name__ == "__main__":
    main()
