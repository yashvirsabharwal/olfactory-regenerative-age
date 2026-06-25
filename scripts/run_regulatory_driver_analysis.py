#!/usr/bin/env python3
"""Run first-pass regulatory-driver pseudobulk target-program analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.regulatory_drivers import (
    build_driver_age_associations,
    build_driver_ora_correlations,
    build_regulatory_driver_map,
    parse_driver_metadata,
    render_regulatory_driver_feasibility,
    score_regulatory_driver_activity,
    write_regulatory_driver_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver-gene-sets", default="configs/regulatory_driver_gene_sets.yaml")
    parser.add_argument("--counts", default="data/processed/pseudobulk_genomewide_counts.tsv.gz")
    parser.add_argument("--metadata", default="data/processed/pseudobulk_genomewide_metadata.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--ora-scores", default="results/tables/ora_augmented_candidate_repeated_cv_scores.tsv")
    parser.add_argument("--chunksize", type=int, default=1000)
    parser.add_argument("--activity-out", default="results/tables/regulatory_driver_activity.tsv")
    parser.add_argument("--donor-activity-out", default="results/tables/regulatory_driver_donor_activity.tsv")
    parser.add_argument("--coverage-out", default="results/tables/regulatory_driver_gene_coverage.tsv")
    parser.add_argument("--driver-map-out", default="results/tables/regulatory_driver_map.tsv")
    parser.add_argument("--age-out", default="results/tables/regulatory_driver_age_associations.tsv")
    parser.add_argument("--ora-out", default="results/tables/regulatory_driver_ora_correlations.tsv")
    parser.add_argument("--feasibility-out", default="docs/regulatory_driver_feasibility.md")
    parser.add_argument("--figure-pdf", default="results/figures/extended_data_regulatory_drivers.pdf")
    parser.add_argument("--figure-png", default="results/figures/extended_data_regulatory_drivers.png")
    args = parser.parse_args()

    driver_metadata = parse_driver_metadata(load_config(args.driver_gene_sets))
    metadata = pd.read_csv(args.metadata, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    ora_scores = pd.read_csv(args.ora_scores, sep="\t")

    activity, donor_activity, coverage = score_regulatory_driver_activity(
        counts_path=args.counts,
        metadata=metadata,
        driver_metadata=driver_metadata,
        chunksize=args.chunksize,
    )
    age = build_driver_age_associations(
        donor_activity=donor_activity,
        manifest=manifest,
        driver_metadata=driver_metadata,
    )
    ora = build_driver_ora_correlations(
        donor_activity=donor_activity,
        ora_scores=ora_scores,
        driver_metadata=driver_metadata,
    )
    driver_map = build_regulatory_driver_map(
        driver_metadata=driver_metadata,
        coverage=coverage,
        age_associations=age,
        ora_correlations=ora,
    )
    feasibility = render_regulatory_driver_feasibility(driver_map=driver_map, coverage=coverage)
    write_regulatory_driver_outputs(
        activity=activity,
        donor_activity=donor_activity,
        coverage=coverage,
        driver_map=driver_map,
        age_associations=age,
        ora_correlations=ora,
        feasibility_note=feasibility,
        activity_out=args.activity_out,
        donor_activity_out=args.donor_activity_out,
        coverage_out=args.coverage_out,
        driver_map_out=args.driver_map_out,
        age_out=args.age_out,
        ora_out=args.ora_out,
        feasibility_out=args.feasibility_out,
        figure_pdf=args.figure_pdf,
        figure_png=args.figure_png,
    )
    print(
        "Wrote regulatory driver analysis: "
        f"{args.driver_map_out} ({driver_map.shape[0]} drivers), "
        f"{args.age_out} ({age.shape[0]} age rows), "
        f"{args.ora_out} ({ora.shape[0]} ORA rows)"
    )


if __name__ == "__main__":
    main()
