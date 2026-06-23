#!/usr/bin/env python3
"""ORA modeling and interpretation command group."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.age_model import train_ora_models, train_ora_models_repeated
from ora.config import load_config
from ora.diagnostics import summarize_ora_diagnostics
from ora.interpretation import build_feature_interpretation
from ora.model_compare import compare_feature_set_deltas, rank_feature_set_summaries
from ora.permutation import run_permutation_null
from ora.sensitivity import run_ora_sensitivity
from ora.stacking import DEFAULT_BASE_MODELS, run_stacked_ora
from ora.stats import run_age_associations
from ora.tuning import run_nested_tuning
from ora.utils import ensure_parent


DEFAULT_PERMUTATION_MODELS = ["random_forest", "xgboost", "catboost", "boosted_ensemble"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_age_associations(subparsers)
    _add_train(subparsers)
    _add_repeated(subparsers)
    _add_diagnostics(subparsers)
    _add_feature_set_comparison(subparsers)
    _add_permutation_null(subparsers)
    _add_nested_tuning(subparsers)
    _add_stacking(subparsers)
    _add_interpretation(subparsers)
    _add_sensitivity(subparsers)
    args = parser.parse_args()
    args.func(args)


def _add_age_associations(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("age-associations")
    parser.add_argument("--features", default="data/processed/donor_cell_state_features.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--out", default="results/tables/age_cell_state_associations.tsv")
    parser.add_argument("--all-donors", action="store_true")
    parser.set_defaults(func=_age_associations)


def _age_associations(args: argparse.Namespace) -> None:
    load_config(args.config)
    features = pd.read_csv(args.features, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    if not args.all_donors and "usable_for_ora_training" in manifest.columns:
        manifest = manifest[manifest["usable_for_ora_training"].astype(bool)].copy()
    covariates = [
        col
        for col in ["sex", "race_ethnicity", "chemistry", "collection_method", "site", "total_cells"]
        if col in manifest.columns
    ]
    feature_columns = [
        col
        for col in features.columns
        if col.startswith("prop__") or col.startswith("clr__") or col.startswith("ratio__")
    ]
    results = run_age_associations(features, manifest, feature_columns=feature_columns, covariates=covariates)
    results.to_csv(ensure_parent(args.out), sep="\t", index=False)
    print(f"Wrote age associations: {args.out} ({results.shape[0]} tests)")


def _add_train(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("train")
    parser.add_argument("--features", default="data/processed/ora_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--config", default="configs/models.yaml")
    parser.add_argument("--performance-out", default="results/tables/ora_model_performance.tsv")
    parser.add_argument("--scores-out", default="results/tables/donor_ora_scores.tsv")
    parser.add_argument("--importance-out", default="results/tables/ora_feature_importance.tsv")
    parser.set_defaults(func=_train)


def _train(args: argparse.Namespace) -> None:
    result = train_ora_models(
        pd.read_csv(args.features, sep="\t"),
        pd.read_csv(args.manifest, sep="\t"),
        load_config(args.config),
    )
    result.performance.to_csv(ensure_parent(args.performance_out), sep="\t", index=False)
    result.predictions.to_csv(ensure_parent(args.scores_out), sep="\t", index=False)
    result.feature_importance.to_csv(ensure_parent(args.importance_out), sep="\t", index=False)
    print(f"Wrote ORA performance: {args.performance_out}")
    print(f"Wrote donor ORA scores: {args.scores_out}")
    print(f"Wrote feature importance: {args.importance_out}")


def _add_repeated(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("repeated")
    parser.add_argument("--features", default="data/processed/ora_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--repeats", type=int, default=None)
    parser.add_argument("--models", nargs="*", default=None)
    parser.add_argument("--repeat-performance-out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--scores-out", default=None)
    parser.add_argument("--feature-stability-out", default=None)
    parser.set_defaults(func=_repeated)


def _repeated(args: argparse.Namespace) -> None:
    model_config = load_config(args.model_config)
    if args.models:
        model_config["model_names"] = args.models
    outputs = load_config(args.gateway_config).get("outputs", {})
    result = train_ora_models_repeated(
        pd.read_csv(args.features, sep="\t"),
        pd.read_csv(args.manifest, sep="\t"),
        model_config,
        repeats=args.repeats,
    )
    repeat_performance_out = args.repeat_performance_out or outputs.get("ora_repeated_cv_performance_tsv", "results/tables/ora_repeated_cv_performance.tsv")
    summary_out = args.summary_out or outputs.get("ora_repeated_cv_summary_tsv", "results/tables/ora_repeated_cv_summary.tsv")
    scores_out = args.scores_out or outputs.get("ora_repeated_cv_scores_tsv", "results/tables/ora_repeated_cv_scores.tsv")
    feature_stability_out = args.feature_stability_out or outputs.get("ora_repeated_cv_feature_stability_tsv", "results/tables/ora_repeated_cv_feature_stability.tsv")
    result.repeat_performance.to_csv(ensure_parent(repeat_performance_out), sep="\t", index=False)
    result.performance_summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.predictions.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    result.feature_stability.to_csv(ensure_parent(feature_stability_out), sep="\t", index=False)
    print(f"Wrote repeated-CV summary: {summary_out}")


def _add_diagnostics(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("diagnostics")
    parser.add_argument("--scores", default=None)
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--calibration-out", default=None)
    parser.add_argument("--age-bin-out", default=None)
    parser.add_argument("--residuals-out", default=None)
    parser.add_argument("--calibrated-scores-out", default=None)
    parser.set_defaults(func=_diagnostics)


def _diagnostics(args: argparse.Namespace) -> None:
    model_config = load_config(args.model_config)
    outputs = load_config(args.gateway_config).get("outputs", {})
    scores_path = args.scores or outputs.get("donor_ora_scores_tsv", "results/tables/donor_ora_scores.tsv")
    manifest_path = args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv")
    result = summarize_ora_diagnostics(
        pd.read_csv(scores_path, sep="\t"),
        model_config=model_config,
        manifest=pd.read_csv(manifest_path, sep="\t"),
    )
    calibration_out = args.calibration_out or outputs.get("ora_calibration_tsv", "results/tables/ora_calibration.tsv")
    age_bin_out = args.age_bin_out or outputs.get("ora_age_bin_errors_tsv", "results/tables/ora_age_bin_errors.tsv")
    residuals_out = args.residuals_out or outputs.get("ora_residual_diagnostics_tsv", "results/tables/ora_residual_diagnostics.tsv")
    calibrated_scores_out = args.calibrated_scores_out or outputs.get("ora_calibrated_scores_tsv", "results/tables/ora_calibrated_scores.tsv")
    result.calibration.to_csv(ensure_parent(calibration_out), sep="\t", index=False)
    result.age_bin_errors.to_csv(ensure_parent(age_bin_out), sep="\t", index=False)
    result.residual_diagnostics.to_csv(ensure_parent(residuals_out), sep="\t", index=False)
    result.calibrated_scores.to_csv(ensure_parent(calibrated_scores_out), sep="\t", index=False)
    print(f"Wrote ORA calibration: {calibration_out} ({result.calibration.shape[0]} models)")


def _add_feature_set_comparison(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("feature-set-comparison")
    parser.add_argument("--base-summary", default="results/tables/ora_repeated_cv_summary.tsv")
    parser.add_argument("--augmented-summary", default="results/tables/ora_augmented_repeated_cv_summary.tsv")
    parser.add_argument("--base-label", default="composition")
    parser.add_argument("--augmented-label", default="composition_plus_modules")
    parser.add_argument("--combined-out", default="results/tables/ora_feature_set_model_comparison.tsv")
    parser.add_argument("--delta-out", default="results/tables/ora_feature_set_model_deltas.tsv")
    parser.set_defaults(func=_feature_set_comparison)


def _feature_set_comparison(args: argparse.Namespace) -> None:
    base = pd.read_csv(args.base_summary, sep="\t")
    augmented = pd.read_csv(args.augmented_summary, sep="\t")
    combined = rank_feature_set_summaries({args.base_label: base, args.augmented_label: augmented})
    deltas = compare_feature_set_deltas(base, augmented, base_label=args.base_label, augmented_label=args.augmented_label)
    combined.to_csv(ensure_parent(args.combined_out), sep="\t", index=False)
    deltas.to_csv(ensure_parent(args.delta_out), sep="\t", index=False)
    print(f"Wrote combined model comparison: {args.combined_out} ({combined.shape[0]} rows)")


def _add_permutation_null(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("permutation-null")
    parser.add_argument("--features", default="data/processed/ora_augmented_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--observed-summary", default=None)
    parser.add_argument("--models", nargs="+", default=DEFAULT_PERMUTATION_MODELS)
    parser.add_argument("--n-permutations", type=int, default=100)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--random-seed", type=int, default=20260615)
    parser.add_argument("--repeat-performance-out", default=None)
    parser.add_argument("--permutation-summary-out", default=None)
    parser.add_argument("--empirical-summary-out", default=None)
    parser.set_defaults(func=_permutation_null)


def _permutation_null(args: argparse.Namespace) -> None:
    model_config = load_config(args.model_config)
    model_config["model_names"] = args.models
    outputs = load_config(args.gateway_config).get("outputs", {})
    result = run_permutation_null(
        pd.read_csv(args.features, sep="\t"),
        pd.read_csv(args.manifest, sep="\t"),
        model_config,
        n_permutations=args.n_permutations,
        repeats=args.repeats,
        random_seed=args.random_seed,
        observed_summary=pd.read_csv(args.observed_summary, sep="\t") if args.observed_summary else None,
    )
    repeat_performance_out = args.repeat_performance_out or outputs.get("ora_permutation_repeat_performance_tsv", "results/tables/ora_permutation_repeat_performance.tsv")
    permutation_summary_out = args.permutation_summary_out or outputs.get("ora_permutation_summary_tsv", "results/tables/ora_permutation_summary.tsv")
    empirical_summary_out = args.empirical_summary_out or outputs.get("ora_permutation_empirical_tsv", "results/tables/ora_permutation_empirical.tsv")
    result.repeat_performance.to_csv(ensure_parent(repeat_performance_out), sep="\t", index=False)
    result.permutation_summary.to_csv(ensure_parent(permutation_summary_out), sep="\t", index=False)
    result.empirical_summary.to_csv(ensure_parent(empirical_summary_out), sep="\t", index=False)
    print(f"Wrote empirical null summary: {empirical_summary_out} ({result.empirical_summary.shape[0]} rows)")


def _add_nested_tuning(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("nested-tuning")
    parser.add_argument("--features", default="data/processed/ora_augmented_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--tuning-config", default=None)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--models", nargs="+", default=["xgboost", "catboost"])
    parser.add_argument("--allow-fallback", action="store_true")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--inner-folds", type=int, default=3)
    parser.add_argument("--max-candidates", type=int, default=None)
    parser.add_argument("--performance-out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--scores-out", default=None)
    parser.add_argument("--trace-out", default=None)
    parser.add_argument("--selected-params-out", default=None)
    parser.set_defaults(func=_nested_tuning)


def _nested_tuning(args: argparse.Namespace) -> None:
    if not args.allow_fallback:
        _require_native_backends(args.models)
    model_config = load_config(args.model_config)
    model_config["model_names"] = args.models
    outputs = load_config(args.gateway_config).get("outputs", {})
    result = run_nested_tuning(
        pd.read_csv(args.features, sep="\t"),
        pd.read_csv(args.manifest, sep="\t"),
        model_config,
        load_config(args.tuning_config) if args.tuning_config else {},
        repeats=args.repeats,
        inner_folds=args.inner_folds,
        max_candidates=args.max_candidates,
    )
    _write_nested_outputs(args, outputs, result)


def _write_nested_outputs(args: argparse.Namespace, outputs: dict, result) -> None:
    performance_out = args.performance_out or outputs.get("ora_nested_tuning_performance_tsv", "results/tables/ora_nested_tuning_performance.tsv")
    summary_out = args.summary_out or outputs.get("ora_nested_tuning_summary_tsv", "results/tables/ora_nested_tuning_summary.tsv")
    scores_out = args.scores_out or outputs.get("ora_nested_tuning_scores_tsv", "results/tables/ora_nested_tuning_scores.tsv")
    trace_out = args.trace_out or outputs.get("ora_nested_tuning_trace_tsv", "results/tables/ora_nested_tuning_trace.tsv")
    selected_params_out = args.selected_params_out or outputs.get("ora_nested_tuning_selected_params_tsv", "results/tables/ora_nested_tuning_selected_params.tsv")
    result.performance.to_csv(ensure_parent(performance_out), sep="\t", index=False)
    result.performance_summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.predictions.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    result.tuning_trace.to_csv(ensure_parent(trace_out), sep="\t", index=False)
    result.selected_params.to_csv(ensure_parent(selected_params_out), sep="\t", index=False)
    print(f"Wrote nested tuning summary: {summary_out} ({result.performance_summary.shape[0]} rows)")


def _add_stacking(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("stacking")
    parser.add_argument("--features", default="data/processed/ora_augmented_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--base-models", nargs="+", default=DEFAULT_BASE_MODELS)
    parser.add_argument("--allow-fallback", action="store_true")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--inner-folds", type=int, default=3)
    parser.add_argument("--performance-out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--scores-out", default=None)
    parser.add_argument("--weights-out", default=None)
    parser.set_defaults(func=_stacking)


def _stacking(args: argparse.Namespace) -> None:
    if not args.allow_fallback:
        _require_native_backends(args.base_models)
    outputs = load_config(args.gateway_config).get("outputs", {})
    result = run_stacked_ora(
        pd.read_csv(args.features, sep="\t"),
        pd.read_csv(args.manifest, sep="\t"),
        load_config(args.model_config),
        base_models=args.base_models,
        repeats=args.repeats,
        inner_folds=args.inner_folds,
    )
    performance_out = args.performance_out or outputs.get("ora_stacking_performance_tsv", "results/tables/ora_stacking_performance.tsv")
    summary_out = args.summary_out or outputs.get("ora_stacking_summary_tsv", "results/tables/ora_stacking_summary.tsv")
    scores_out = args.scores_out or outputs.get("ora_stacking_scores_tsv", "results/tables/ora_stacking_scores.tsv")
    weights_out = args.weights_out or outputs.get("ora_stacking_weights_tsv", "results/tables/ora_stacking_weights.tsv")
    result.performance.to_csv(ensure_parent(performance_out), sep="\t", index=False)
    result.performance_summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.predictions.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    result.meta_weights.to_csv(ensure_parent(weights_out), sep="\t", index=False)
    print(f"Wrote stacking summary: {summary_out} ({result.performance_summary.shape[0]} rows)")


def _add_interpretation(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("interpretation")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--feature-stability", default=None)
    parser.add_argument("--associations", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--top-per-model", type=int, default=12)
    parser.add_argument("--top-n", type=int, default=30)
    parser.set_defaults(func=_interpretation)


def _interpretation(args: argparse.Namespace) -> None:
    outputs = load_config(args.gateway_config).get("outputs", {})
    stability_path = args.feature_stability or outputs.get("ora_augmented_repeated_cv_feature_stability_tsv", "results/tables/ora_augmented_repeated_cv_feature_stability.tsv")
    associations_path = args.associations or outputs.get("age_associations_tsv", "results/tables/age_cell_state_associations.tsv")
    out_path = args.out or outputs.get("ora_feature_interpretation_tsv", "results/tables/ora_feature_interpretation.tsv")
    interpretation = build_feature_interpretation(
        pd.read_csv(stability_path, sep="\t"),
        _read_optional_tsv(associations_path),
        top_per_model=args.top_per_model,
        top_n=args.top_n,
    )
    interpretation.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(f"Wrote ORA feature interpretation: {out_path} ({interpretation.shape[0]} rows)")


def _add_sensitivity(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("sensitivity")
    parser.add_argument("--features", default="data/processed/ora_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--min-cell-thresholds", nargs="*", type=int, default=[500, 1000, 5000, 10000])
    parser.add_argument("--min-train-donors", type=int, default=20)
    parser.add_argument("--models", nargs="*", default=None)
    parser.add_argument("--scenarios-out", default=None)
    parser.add_argument("--performance-out", default=None)
    parser.add_argument("--scores-out", default=None)
    parser.set_defaults(func=_sensitivity)


def _sensitivity(args: argparse.Namespace) -> None:
    model_config = load_config(args.model_config)
    if args.models:
        model_config["model_names"] = args.models
    outputs = load_config(args.gateway_config).get("outputs", {})
    result = run_ora_sensitivity(
        pd.read_csv(args.features, sep="\t"),
        pd.read_csv(args.manifest, sep="\t"),
        model_config,
        min_cell_thresholds=args.min_cell_thresholds,
        min_train_donors=args.min_train_donors,
    )
    scenarios_out = args.scenarios_out or outputs.get("ora_sensitivity_scenarios_tsv", "results/tables/ora_sensitivity_scenarios.tsv")
    performance_out = args.performance_out or outputs.get("ora_sensitivity_performance_tsv", "results/tables/ora_sensitivity_performance.tsv")
    scores_out = args.scores_out or outputs.get("ora_sensitivity_scores_tsv", "results/tables/ora_sensitivity_scores.tsv")
    result.scenarios.to_csv(ensure_parent(scenarios_out), sep="\t", index=False)
    result.performance.to_csv(ensure_parent(performance_out), sep="\t", index=False)
    result.scores.to_csv(ensure_parent(scores_out), sep="\t", index=False)
    print(f"Wrote ORA sensitivity performance: {performance_out}")


def _require_native_backends(model_names: list[str]) -> None:
    required = {"xgboost": "xgboost", "lightgbm": "lightgbm", "catboost": "catboost"}
    missing = []
    for model_name in model_names:
        package = required.get(model_name)
        if package is None:
            continue
        try:
            __import__(package)
        except Exception as exc:  # pragma: no cover
            missing.append(f"{model_name}: {exc}")
    if missing:
        joined = "\n".join(missing)
        raise SystemExit(
            "Native booster backend check failed. Fix optional runtime dependencies or pass "
            f"--allow-fallback to use sklearn fallback models.\n{joined}"
        )


def _read_optional_tsv(path: str | None) -> pd.DataFrame | None:
    if not path:
        return None
    candidate = Path(path)
    if not candidate.exists():
        return None
    return pd.read_csv(candidate, sep="\t")


if __name__ == "__main__":
    main()
