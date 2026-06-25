#!/usr/bin/env python3
"""Build and summarize donor-level scVI embedding age baselines."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.scvi_donor import (  # noqa: E402
    build_scvi_donor_embedding_features,
    compare_scvi_donor_baseline,
    summarize_scvi_state_importance,
)
from ora.features import (  # noqa: E402
    feature_kind_counts,
    merge_donor_feature_matrices,
    summarize_feature_family_stability,
)
from ora.utils import ensure_parent  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_features(subparsers)
    _add_hybrid_features(subparsers)
    _add_state_importance(subparsers)
    _add_family_importance(subparsers)
    _add_compare(subparsers)
    args = parser.parse_args()
    args.func(args)


def _add_features(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("features")
    parser.add_argument("--h5ad", default="data/processed/gateway_scvi_full_4m_reduced.h5ad")
    parser.add_argument("--embedding-key", default="X_scvi")
    parser.add_argument("--donor-col", default="donor_id")
    parser.add_argument("--cell-state-col", default="fine_celltype")
    parser.add_argument("--top-cell-states", type=int, default=12)
    parser.add_argument("--min-cells-per-donor", type=int, default=20)
    parser.add_argument("--min-cells-per-state", type=int, default=20)
    parser.add_argument("--chunk-size", type=int, default=250_000)
    parser.add_argument("--no-state-features", action="store_true")
    parser.add_argument("--features-out", default="data/processed/scvi_donor_embedding_features.tsv")
    parser.add_argument("--qc-out", default="results/tables/scvi_donor_embedding_feature_qc.tsv")
    parser.set_defaults(func=_features)


def _features(args: argparse.Namespace) -> None:
    features, qc = build_scvi_donor_embedding_features(
        args.h5ad,
        embedding_key=args.embedding_key,
        donor_col=args.donor_col,
        cell_state_col=args.cell_state_col,
        top_cell_states=args.top_cell_states,
        min_cells_per_donor=args.min_cells_per_donor,
        min_cells_per_state=args.min_cells_per_state,
        chunk_size=args.chunk_size,
        include_state_features=not args.no_state_features,
    )
    features.to_csv(ensure_parent(args.features_out), sep="\t", index=False)
    qc.to_csv(ensure_parent(args.qc_out), sep="\t", index=False)
    print(f"Wrote scVI donor features: {args.features_out} ({features.shape[0]} donors x {features.shape[1] - 1} features)")
    print(f"Wrote scVI donor feature QC: {args.qc_out} ({qc.shape[0]} rows)")


def _add_hybrid_features(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("hybrid-features")
    parser.add_argument("--ora-features", default="data/processed/ora_augmented_feature_matrix.tsv")
    parser.add_argument("--scvi-features", default="data/processed/scvi_donor_embedding_features.tsv")
    parser.add_argument("--out", default="data/processed/ora_scvi_hybrid_feature_matrix.tsv")
    parser.set_defaults(func=_hybrid_features)


def _hybrid_features(args: argparse.Namespace) -> None:
    hybrid = merge_donor_feature_matrices(
        [
            pd.read_csv(args.ora_features, sep="\t"),
            pd.read_csv(args.scvi_features, sep="\t"),
        ]
    )
    counts = feature_kind_counts(hybrid)
    hybrid.to_csv(ensure_parent(args.out), sep="\t", index=False)
    print(
        "Wrote ORA+scVI hybrid features: "
        f"{args.out} ({hybrid.shape[0]} donors; "
        f"{counts['composition']} composition, {counts['module']} module, "
        f"{counts['scvi_global']} scVI global, {counts['scvi_cell_state']} scVI cell-state)"
    )


def _add_state_importance(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("state-importance")
    parser.add_argument("--feature-stability", default="results/tables/scvi_donor_embedding_feature_stability.tsv")
    parser.add_argument("--out", default="results/tables/scvi_donor_embedding_state_importance.tsv")
    parser.set_defaults(func=_state_importance)


def _state_importance(args: argparse.Namespace) -> None:
    summary = summarize_scvi_state_importance(pd.read_csv(args.feature_stability, sep="\t"))
    summary.to_csv(ensure_parent(args.out), sep="\t", index=False)
    print(f"Wrote scVI donor state importance: {args.out} ({summary.shape[0]} rows)")


def _add_family_importance(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("family-importance")
    parser.add_argument("--feature-stability", default="results/tables/ora_scvi_hybrid_feature_stability.tsv")
    parser.add_argument("--out", default="results/tables/ora_scvi_hybrid_feature_family_importance.tsv")
    parser.set_defaults(func=_family_importance)


def _family_importance(args: argparse.Namespace) -> None:
    summary = summarize_feature_family_stability(pd.read_csv(args.feature_stability, sep="\t"))
    summary.to_csv(ensure_parent(args.out), sep="\t", index=False)
    print(f"Wrote hybrid feature-family importance: {args.out} ({summary.shape[0]} rows)")


def _add_compare(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("compare")
    parser.add_argument(
        "--summary",
        nargs=2,
        action="append",
        metavar=("LABEL", "PATH"),
        required=True,
        help="Feature-set label and repeated-CV summary TSV. May be repeated.",
    )
    parser.add_argument("--out", default="results/tables/scvi_donor_embedding_model_comparison.tsv")
    parser.set_defaults(func=_compare)


def _compare(args: argparse.Namespace) -> None:
    comparison = compare_scvi_donor_baseline(args.summary)
    comparison.to_csv(ensure_parent(args.out), sep="\t", index=False)
    print(f"Wrote scVI donor model comparison: {args.out} ({comparison.shape[0]} rows)")


if __name__ == "__main__":
    main()
