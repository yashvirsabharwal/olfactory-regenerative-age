#!/usr/bin/env python3
"""Run matched feature-family ablation benchmarks."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config  # noqa: E402
from ora.feature_ablation import run_feature_family_ablation, write_ablation_figure  # noqa: E402
from ora.utils import ensure_parent  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", default="data/processed/ora_scvi_hybrid_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--models", nargs="*", default=["hist_gradient_boosting", "xgboost", "catboost", "boosted_ensemble"])
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--n-permutations", type=int, default=20)
    parser.add_argument("--permutation-repeats", type=int, default=1)
    parser.add_argument("--random-seed", type=int, default=20260624)
    parser.add_argument("--matrix-dir", default=None)
    parser.add_argument("--pseudobulk-counts", default=None)
    parser.add_argument("--pseudobulk-metadata", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--deltas-out", default=None)
    parser.add_argument("--feasibility-out", default=None)
    parser.add_argument("--repeat-performance-out", default=None)
    parser.add_argument("--scores-out", default=None)
    parser.add_argument("--figure-pdf", default=None)
    parser.add_argument("--figure-png", default=None)
    args = parser.parse_args()

    model_config = load_config(args.model_config)
    outputs = load_config(args.gateway_config).get("outputs", {})
    result = run_feature_family_ablation(
        feature_matrix=pd.read_csv(args.features, sep="\t"),
        manifest=pd.read_csv(args.manifest, sep="\t"),
        model_config=model_config,
        output_dir=args.matrix_dir or outputs.get("feature_family_ablation_matrix_dir", "data/processed/feature_family_ablation"),
        pseudobulk_counts_path=args.pseudobulk_counts or outputs.get("pseudobulk_genomewide_counts_tsv", "data/processed/pseudobulk_genomewide_counts.tsv.gz"),
        pseudobulk_metadata_path=args.pseudobulk_metadata or outputs.get("pseudobulk_genomewide_metadata_tsv", "data/processed/pseudobulk_genomewide_metadata.tsv"),
        models=args.models,
        repeats=args.repeats,
        n_permutations=args.n_permutations,
        permutation_repeats=args.permutation_repeats,
        random_seed=args.random_seed,
    )
    summary_out = args.summary_out or outputs.get("ora_feature_family_ablation_summary_tsv", "results/tables/ora_feature_family_ablation_summary.tsv")
    deltas_out = args.deltas_out or outputs.get("ora_feature_family_ablation_deltas_tsv", "results/tables/ora_feature_family_ablation_deltas.tsv")
    feasibility_out = args.feasibility_out or outputs.get("ora_feature_family_ablation_feasibility_tsv", "results/tables/ora_feature_family_ablation_feasibility.tsv")
    repeat_out = args.repeat_performance_out or outputs.get("ora_feature_family_ablation_repeat_performance_tsv", "results/tables/ora_feature_family_ablation_repeat_performance.tsv")
    scores_out = args.scores_out or outputs.get("ora_feature_family_ablation_scores_tsv", "results/tables/ora_feature_family_ablation_scores.tsv")
    figure_pdf = args.figure_pdf or outputs.get("feature_family_ablation_pdf", "results/figures/manuscript_figure_feature_ablation.pdf")
    figure_png = args.figure_png or outputs.get("feature_family_ablation_png", "results/figures/manuscript_figure_feature_ablation.png")

    result.summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.deltas.to_csv(ensure_parent(deltas_out), sep="\t", index=False)
    result.feasibility.to_csv(ensure_parent(feasibility_out), sep="\t", index=False)
    result.repeat_performance.to_csv(ensure_parent(repeat_out), sep="\t", index=False)
    result.scores.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    write_ablation_figure(result.summary, figure_pdf, figure_png)
    print(f"Wrote feature-family ablation summary: {summary_out} ({result.summary.shape[0]} rows)")
    print(f"Wrote feature-family ablation deltas: {deltas_out}")
    print(f"Wrote feature-family ablation feasibility: {feasibility_out}")
    print(f"Wrote feature-family ablation figure: {figure_pdf}")


if __name__ == "__main__":
    main()
