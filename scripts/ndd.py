#!/usr/bin/env python3
"""NDD projection and guardrail command group."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.age_model import project_ora_models
from ora.config import load_config
from ora.ndd import (
    compare_ndd_feature_sets,
    donor_projection_appendix,
    ndd_label_permutation,
    ndd_projection_diagnostics,
    summarize_ndd_projection_uncertainty,
)
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_project(subparsers)
    _add_feature_sensitivity(subparsers)
    _add_uncertainty(subparsers)
    _add_diagnostics(subparsers)
    _add_label_permutation(subparsers)
    args = parser.parse_args()
    args.func(args)


def _add_project(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("project")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--allow-fallback", action="store_true")
    parser.add_argument("--features", default=None)
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--scores-out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.set_defaults(func=_project)


def _project(args: argparse.Namespace) -> None:
    gateway_config = load_config(args.gateway_config)
    model_config = load_config(args.model_config)
    model_config["allow_fallback"] = bool(args.allow_fallback)
    outputs = gateway_config.get("outputs", {})
    features_path = args.features or outputs.get("ora_augmented_feature_matrix_tsv", "data/processed/ora_augmented_feature_matrix.tsv")
    manifest_path = args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv")
    scores_out = args.scores_out or outputs.get("ndd_ora_projection_tsv", "results/tables/ndd_ora_projection.tsv")
    summary_out = args.summary_out or outputs.get("ndd_ora_projection_summary_tsv", "results/tables/ndd_ora_projection_summary.tsv")
    result = project_ora_models(pd.read_csv(features_path, sep="\t"), pd.read_csv(manifest_path, sep="\t"), model_config)
    result.predictions.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    result.summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    print(f"Wrote frozen ORA projections: {scores_out} ({result.predictions['donor_id'].nunique()} donors)")


def _add_feature_sensitivity(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("feature-sensitivity")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--allow-fallback", action="store_true")
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--composition-features", default=None)
    parser.add_argument("--augmented-features", default=None)
    parser.add_argument("--composition-projection-out", default=None)
    parser.add_argument("--composition-summary-out", default=None)
    parser.add_argument("--augmented-projection-out", default=None)
    parser.add_argument("--augmented-summary-out", default=None)
    parser.add_argument("--comparison-out", default=None)
    parser.add_argument("--donor-appendix-out", default=None)
    parser.set_defaults(func=_feature_sensitivity)


def _feature_sensitivity(args: argparse.Namespace) -> None:
    gateway_config = load_config(args.gateway_config)
    model_config = load_config(args.model_config)
    model_config["allow_fallback"] = bool(args.allow_fallback)
    outputs = gateway_config.get("outputs", {})
    manifest = pd.read_csv(args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv"), sep="\t")
    composition = _project_feature_set(
        args.composition_features or outputs.get("ora_feature_matrix_tsv", "data/processed/ora_feature_matrix.tsv"),
        manifest,
        model_config,
        "composition",
    )
    augmented = _project_feature_set(
        args.augmented_features or outputs.get("ora_augmented_feature_matrix_tsv", "data/processed/ora_augmented_feature_matrix.tsv"),
        manifest,
        model_config,
        "augmented",
    )
    comparison = compare_ndd_feature_sets({"composition": composition.predictions, "augmented": augmented.predictions})
    appendix = pd.concat(
        [donor_projection_appendix(composition.predictions), donor_projection_appendix(augmented.predictions)],
        ignore_index=True,
    )
    paths = {
        "composition_projection": args.composition_projection_out or outputs.get("ndd_ora_projection_composition_tsv", "results/tables/ndd_ora_projection_composition.tsv"),
        "composition_summary": args.composition_summary_out or outputs.get("ndd_ora_projection_composition_summary_tsv", "results/tables/ndd_ora_projection_composition_summary.tsv"),
        "augmented_projection": args.augmented_projection_out or outputs.get("ndd_ora_projection_tsv", "results/tables/ndd_ora_projection.tsv"),
        "augmented_summary": args.augmented_summary_out or outputs.get("ndd_ora_projection_summary_tsv", "results/tables/ndd_ora_projection_summary.tsv"),
        "comparison": args.comparison_out or outputs.get("ndd_ora_projection_feature_comparison_tsv", "results/tables/ndd_ora_projection_feature_comparison.tsv"),
        "appendix": args.donor_appendix_out or outputs.get("ndd_ora_projection_donor_appendix_tsv", "results/tables/ndd_ora_projection_donor_appendix.tsv"),
    }
    composition.predictions.to_csv(ensure_parent(paths["composition_projection"]), sep="\t", index=False)
    composition.summary.to_csv(ensure_parent(paths["composition_summary"]), sep="\t", index=False)
    augmented.predictions.to_csv(ensure_parent(paths["augmented_projection"]), sep="\t", index=False)
    augmented.summary.to_csv(ensure_parent(paths["augmented_summary"]), sep="\t", index=False)
    comparison.to_csv(ensure_parent(paths["comparison"]), sep="\t", index=False)
    appendix.to_csv(ensure_parent(paths["appendix"]), sep="\t", index=False)
    print(f"Wrote NDD feature-set comparison: {paths['comparison']} ({comparison.shape[0]} rows)")


def _add_uncertainty(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("uncertainty")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--projection", default=None)
    parser.add_argument("--n-bootstrap", type=int, default=5000)
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--uncertainty-out", default=None)
    parser.add_argument("--context-out", default=None)
    parser.set_defaults(func=_uncertainty)


def _uncertainty(args: argparse.Namespace) -> None:
    outputs = load_config(args.gateway_config).get("outputs", {})
    projection_path = args.projection or outputs.get("ndd_ora_projection_tsv", "results/tables/ndd_ora_projection.tsv")
    result = summarize_ndd_projection_uncertainty(
        pd.read_csv(projection_path, sep="\t"),
        n_bootstrap=args.n_bootstrap,
        random_seed=args.random_seed,
    )
    uncertainty_out = args.uncertainty_out or outputs.get("ndd_ora_projection_uncertainty_tsv", "results/tables/ndd_ora_projection_uncertainty.tsv")
    context_out = args.context_out or outputs.get("ndd_ora_projection_context_tsv", "results/tables/ndd_ora_projection_context.tsv")
    result.uncertainty.to_csv(ensure_parent(uncertainty_out), sep="\t", index=False)
    result.context.to_csv(ensure_parent(context_out), sep="\t", index=False)
    print(f"Wrote NDD projection uncertainty: {uncertainty_out} ({result.uncertainty.shape[0]} rows)")


def _add_diagnostics(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("diagnostics")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--projection", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--min-donors-ok", type=int, default=2)
    parser.set_defaults(func=_diagnostics)


def _diagnostics(args: argparse.Namespace) -> None:
    outputs = load_config(args.gateway_config).get("outputs", {})
    projection_path = args.projection or outputs.get("ndd_ora_projection_tsv", "results/tables/ndd_ora_projection.tsv")
    out_path = args.out or outputs.get("ndd_ora_projection_diagnostics_tsv", "results/tables/ndd_ora_projection_diagnostics.tsv")
    diagnostics = ndd_projection_diagnostics(pd.read_csv(projection_path, sep="\t"), min_donors_ok=args.min_donors_ok)
    diagnostics.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(f"Wrote NDD projection diagnostics: {out_path} ({diagnostics.shape[0]} rows)")


def _add_label_permutation(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("label-permutation")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--projection", default=None)
    parser.add_argument("--n-permutations", type=int, default=5000)
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--out", default=None)
    parser.set_defaults(func=_label_permutation)


def _label_permutation(args: argparse.Namespace) -> None:
    outputs = load_config(args.gateway_config).get("outputs", {})
    projection_path = args.projection or outputs.get("ndd_ora_projection_tsv", "results/tables/ndd_ora_projection.tsv")
    out_path = args.out or outputs.get("ndd_label_permutation_tsv", "results/tables/ndd_label_permutation.tsv")
    result = ndd_label_permutation(
        pd.read_csv(projection_path, sep="\t"),
        n_permutations=args.n_permutations,
        random_seed=args.random_seed,
    )
    result.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(f"Wrote NDD label permutation: {out_path} ({result.shape[0]} rows)")


def _project_feature_set(feature_path: str, manifest: pd.DataFrame, model_config: dict, feature_set: str):
    result = project_ora_models(pd.read_csv(feature_path, sep="\t"), manifest, model_config)
    result.predictions.insert(0, "feature_set", feature_set)
    result.summary.insert(0, "feature_set", feature_set)
    return result


if __name__ == "__main__":
    main()
