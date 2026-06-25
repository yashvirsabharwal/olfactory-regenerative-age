#!/usr/bin/env python3
"""Run leakage-safe aging-clock-style expression baselines."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.expression_clock import run_expression_clock_baseline, write_expression_clock_outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", default="data/processed/pseudobulk_genomewide_counts.tsv.gz")
    parser.add_argument("--metadata", default="data/processed/pseudobulk_genomewide_metadata.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument(
        "--models",
        nargs="+",
        default=["ridge", "hist_gradient_boosting", "xgboost", "catboost", "boosted_ensemble"],
    )
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--n-pcs", type=int, default=20)
    parser.add_argument("--top-variable-genes", type=int, default=5000)
    parser.add_argument("--min-detection-donors", type=int, default=20)
    parser.add_argument("--chunksize", type=int, default=750)
    parser.add_argument("--summary-out", default="results/tables/aging_clock_baseline_performance.tsv")
    parser.add_argument("--repeat-performance-out", default="results/tables/aging_clock_baseline_repeat_performance.tsv")
    parser.add_argument("--scores-out", default="results/tables/aging_clock_baseline_scores.tsv")
    parser.add_argument("--feature-stability-out", default="results/tables/aging_clock_baseline_feature_stability.tsv")
    parser.add_argument("--feasibility-out", default="results/tables/aging_clock_baseline_feasibility.tsv")
    parser.add_argument("--comparison-out", default="results/tables/aging_clock_baseline_model_comparison.tsv")
    args = parser.parse_args()

    result = run_expression_clock_baseline(
        counts_path=args.counts,
        metadata_path=args.metadata,
        manifest=pd.read_csv(args.manifest, sep="\t"),
        model_config=load_config(args.model_config),
        models=args.models,
        repeats=args.repeats,
        n_pcs=args.n_pcs,
        top_variable_genes=args.top_variable_genes,
        min_detection_donors=args.min_detection_donors,
        chunksize=args.chunksize,
        comparison_summaries=[
            ("feature_family_best", "results/tables/ora_feature_family_ablation_summary.tsv"),
            ("ora_scvi_hybrid", "results/tables/ora_scvi_hybrid_age_model_summary.tsv"),
            ("scvi_donor_embedding", "results/tables/scvi_donor_embedding_age_model_summary.tsv"),
            ("composition_plus_modules", "results/tables/ora_augmented_candidate_repeated_cv_summary.tsv"),
        ],
    )
    write_expression_clock_outputs(
        result,
        summary_out=args.summary_out,
        repeat_performance_out=args.repeat_performance_out,
        scores_out=args.scores_out,
        feature_stability_out=args.feature_stability_out,
        feasibility_out=args.feasibility_out,
        comparison_out=args.comparison_out,
    )
    best = result.performance_summary.sort_values("mae_mean").iloc[0]
    print(
        "Wrote aging-clock-style expression baseline: "
        f"{args.summary_out}; best={best['model']} MAE={best['mae_mean']:.3f}"
    )


if __name__ == "__main__":
    main()
