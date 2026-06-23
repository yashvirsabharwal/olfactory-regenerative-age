#!/usr/bin/env python3
"""Data download, inspection, aggregation, feature, and module commands."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.aggregate import aggregate_cell_counts, build_cell_state_features
from ora.config import load_config, project_path
from ora.download import download_cellxgene_dataset, download_direct_url, infer_download_mode
from ora.features import build_ora_feature_matrix, feature_kind_counts
from ora.io import inspect_h5ad, load_obs
from ora.metadata import build_manifest, resolve_columns, summarize_cohort
from ora.modules import DEFAULT_GROUPBY, score_gene_sets_h5ad
from ora.utils import ensure_parent, write_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_download(subparsers)
    _add_inspect(subparsers)
    _add_manifest(subparsers)
    _add_aggregate(subparsers)
    _add_features(subparsers)
    _add_modules(subparsers)
    args = parser.parse_args()
    args.func(args)


def _add_download(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("download")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--dataset-id", default=None)
    parser.add_argument("--url", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.set_defaults(func=_download)


def _download(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    source = config.get("source", {})
    dataset_id = args.dataset_id or source.get("dataset_id")
    output_path = project_path(args.out or source.get("h5ad_path", "data/raw/gateway.h5ad"))
    direct_url = args.url or source.get("direct_h5ad_url")
    mode = infer_download_mode(direct_url, dataset_id)
    if mode == "missing":
        raise SystemExit("No dataset ID or direct H5AD URL was provided.")
    if args.dry_run:
        print(f"mode: {mode}")
        print(f"dataset_id: {dataset_id}")
        print(f"url: {direct_url}")
        print(f"output: {output_path}")
        if source.get("h5ad_filesize_bytes"):
            print(f"expected_bytes: {source['h5ad_filesize_bytes']}")
        return
    if output_path.exists() and not args.force:
        raise SystemExit(f"Output already exists: {output_path}. Pass --force to overwrite.")
    path = download_direct_url(str(direct_url), output_path) if mode == "url" else download_cellxgene_dataset(str(dataset_id), output_path)
    print(f"Wrote H5AD: {path}")


def _add_inspect(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("inspect")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--schema-out", default=None)
    parser.add_argument("--obs-columns-out", default=None)
    parser.add_argument("--var-columns-out", default=None)
    parser.set_defaults(func=_inspect)


def _inspect(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    h5ad_path = args.h5ad or config["source"]["h5ad_path"]
    schema, obs_table, var_table = inspect_h5ad(project_path(h5ad_path))
    outputs = config.get("outputs", {})
    schema_path = args.schema_out or outputs.get("schema_json", "results/reports/h5ad_schema.json")
    obs_columns_path = args.obs_columns_out or outputs.get("obs_columns_tsv", "data/metadata/gateway_obs_columns.tsv")
    var_columns_path = args.var_columns_out or outputs.get("var_columns_tsv", "data/metadata/gateway_var_columns.tsv")
    write_json(schema, schema_path)
    obs_table.to_csv(ensure_parent(obs_columns_path), sep="\t", index=False)
    var_table.to_csv(ensure_parent(var_columns_path), sep="\t", index=False)
    print(f"Inspected {h5ad_path}: {schema['n_obs']} cells x {schema['n_vars']} genes")


def _add_manifest(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("manifest")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.set_defaults(func=_manifest)


def _manifest(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    obs = load_obs(project_path(args.h5ad or config["source"]["h5ad_path"]))
    manifest = build_manifest(obs, config, resolve_columns(list(obs.columns), config))
    summary = summarize_cohort(manifest)
    outputs = config.get("outputs", {})
    manifest_path = args.out or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv")
    summary_path = args.summary_out or outputs.get("cohort_summary_tsv", "results/tables/cohort_summary.tsv")
    manifest.to_csv(ensure_parent(manifest_path), sep="\t", index=False)
    summary.to_csv(ensure_parent(summary_path), sep="\t", index=False)
    print(f"Wrote manifest: {manifest_path} ({manifest.shape[0]} donor/sample rows)")


def _add_aggregate(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("aggregate")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--counts-out", default=None)
    parser.add_argument("--features-out", default=None)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    parser.set_defaults(func=_aggregate)


def _aggregate(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    obs = load_obs(project_path(args.h5ad or config["source"]["h5ad_path"]))
    counts = aggregate_cell_counts(obs, config, resolve_columns(list(obs.columns), config))
    features = build_cell_state_features(counts, config, pseudocount=args.pseudocount)
    outputs = config.get("outputs", {})
    counts_path = args.counts_out or outputs.get("cell_counts_tsv", "data/processed/donor_cell_state_counts.tsv")
    features_path = args.features_out or outputs.get("cell_features_tsv", "data/processed/donor_cell_state_features.tsv")
    counts.to_csv(ensure_parent(counts_path), sep="\t", index=False)
    features.to_csv(ensure_parent(features_path), sep="\t", index=False)
    print(f"Wrote cell counts: {counts_path} ({counts.shape[0]} rows)")
    print(f"Wrote cell-state features: {features_path} ({features.shape[1] - 1} features)")


def _add_features(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("features")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--features", default=None)
    parser.add_argument("--module-features", default=None)
    parser.add_argument("--include-modules", action="store_true")
    parser.add_argument("--out", default=None)
    parser.set_defaults(func=_features)


def _features(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    outputs = config.get("outputs", {})
    feature_path = args.features or outputs.get("cell_features_tsv", "data/processed/donor_cell_state_features.tsv")
    include_modules = args.include_modules or bool(args.module_features)
    if include_modules:
        out_path = args.out or outputs.get("ora_augmented_feature_matrix_tsv", "data/processed/ora_augmented_feature_matrix.tsv")
        module_path = args.module_features or outputs.get("donor_module_features_tsv", "data/processed/donor_module_features.tsv")
    else:
        out_path = args.out or outputs.get("ora_feature_matrix_tsv", "data/processed/ora_feature_matrix.tsv")
        module_path = None
    matrix = build_ora_feature_matrix(
        pd.read_csv(feature_path, sep="\t"),
        pd.read_csv(module_path, sep="\t") if module_path else None,
    )
    matrix.to_csv(ensure_parent(out_path), sep="\t", index=False)
    counts = feature_kind_counts(matrix)
    print(f"Wrote ORA feature matrix: {out_path} ({counts['composition']} composition, {counts['module']} module)")


def _add_modules(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("modules")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--gene-sets", default="configs/gene_sets.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--coverage-out", default=None)
    parser.add_argument("--donor-features-out", default=None)
    parser.add_argument("--chunk-size", type=int, default=None)
    parser.add_argument("--layer", default=None)
    parser.add_argument("--groupby", nargs="+", default=list(DEFAULT_GROUPBY))
    parser.add_argument("--apply-qc", action="store_true")
    parser.add_argument("--no-log1p", action="store_true")
    parser.set_defaults(func=_modules)


def _modules(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    outputs = config.get("outputs", {})
    result = score_gene_sets_h5ad(
        args.h5ad or config.get("source", {}).get("h5ad_path", "data/raw/gateway.h5ad"),
        config,
        load_config(args.gene_sets),
        groupby=args.groupby,
        chunk_size=args.chunk_size,
        layer=args.layer,
        log1p=not args.no_log1p,
        apply_qc=args.apply_qc,
    )
    summary_out = args.summary_out or outputs.get("module_score_summary_tsv", "results/tables/module_score_summary.tsv")
    coverage_out = args.coverage_out or outputs.get("module_gene_coverage_tsv", "results/tables/module_gene_coverage.tsv")
    donor_features_out = args.donor_features_out or outputs.get("donor_module_features_tsv", "data/processed/donor_module_features.tsv")
    result.summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.coverage.to_csv(ensure_parent(coverage_out), sep="\t", index=False)
    result.donor_features.to_csv(ensure_parent(donor_features_out), sep="\t", index=False)
    print(f"Wrote donor module features: {donor_features_out} ({result.donor_features.shape[1] - 1} features)")


if __name__ == "__main__":
    main()
