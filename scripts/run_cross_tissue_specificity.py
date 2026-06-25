#!/usr/bin/env python3
"""Run first-pass ORA cross-tissue specificity classification."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.cross_tissue import (
    build_cross_tissue_candidate_matrix,
    build_cross_tissue_specificity_summary,
    build_ora_cross_tissue_feature_classification,
    render_cross_tissue_specificity_plan,
    write_cross_tissue_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", default="data/processed/ora_augmented_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--feature-stability", default="results/tables/ora_augmented_candidate_repeated_cv_feature_stability.tsv")
    parser.add_argument("--feature-interpretation", default="results/tables/ora_feature_interpretation.tsv")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--candidate-out", default="results/tables/cross_tissue_candidate_matrix.tsv")
    parser.add_argument("--classification-out", default="results/tables/ora_cross_tissue_feature_classification.tsv")
    parser.add_argument("--summary-out", default="results/tables/ora_cross_tissue_specificity.tsv")
    parser.add_argument("--plan-out", default="docs/cross_tissue_specificity_plan.md")
    parser.add_argument("--figure-pdf", default="results/figures/extended_data_cross_tissue_specificity.pdf")
    parser.add_argument("--figure-png", default="results/figures/extended_data_cross_tissue_specificity.png")
    args = parser.parse_args()

    external_config = load_config(args.external_config)
    feature_matrix = pd.read_csv(args.features, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    feature_stability = _read_optional_table(args.feature_stability)
    feature_interpretation = _read_optional_table(args.feature_interpretation)
    candidate_matrix = build_cross_tissue_candidate_matrix(external_config)
    classification = build_ora_cross_tissue_feature_classification(
        feature_matrix=feature_matrix,
        manifest=manifest,
        feature_stability=feature_stability,
        feature_interpretation=feature_interpretation,
        comparator_matrix=candidate_matrix,
    )
    summary = build_cross_tissue_specificity_summary(classification)
    plan = render_cross_tissue_specificity_plan(candidate_matrix, summary)
    write_cross_tissue_outputs(
        candidate_matrix=candidate_matrix,
        classification=classification,
        summary=summary,
        plan_md=plan,
        candidate_out=args.candidate_out,
        classification_out=args.classification_out,
        summary_out=args.summary_out,
        plan_out=args.plan_out,
        figure_pdf=args.figure_pdf,
        figure_png=args.figure_png,
    )
    print(f"Wrote cross-tissue candidate matrix: {args.candidate_out} ({candidate_matrix.shape[0]} rows)")
    print(f"Wrote ORA cross-tissue classification: {args.classification_out} ({classification.shape[0]} features)")
    print(f"Wrote ORA cross-tissue summary: {args.summary_out}")
    print(f"Wrote cross-tissue specificity plan: {args.plan_out}")


def _read_optional_table(path: str) -> pd.DataFrame | None:
    candidate = Path(path)
    if not candidate.exists():
        return None
    return pd.read_csv(candidate, sep="\t")


if __name__ == "__main__":
    main()
