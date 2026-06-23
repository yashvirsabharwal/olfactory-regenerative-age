#!/usr/bin/env python3
"""External validation command group."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.external import (
    build_external_10x_marker_mapped_anndata,
    external_candidate_matrix,
    external_dataset_summary,
    external_mapped_feature_concordance,
    external_marker_age_concordance,
    external_validation_evidence_summary,
    feature_matrix_contract_summary,
    inspect_external_archive,
    parse_geo_series_matrix_metadata,
    published_gene_list_coverage,
    score_external_10x_marker_composition,
    score_external_10x_modules,
    validate_external_feature_matrix,
)
from ora.io import read_h5ad_backed
from ora.modules import DEFAULT_SYMBOL_COLUMNS
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_inspect_archive(subparsers)
    _add_score_modules(subparsers)
    _add_score_markers(subparsers)
    _add_build_mapped(subparsers)
    _add_compare_mapped(subparsers)
    _add_compare_marker_age(subparsers)
    _add_summarize_validation(subparsers)
    _add_candidate_matrix(subparsers)
    _add_evidence(subparsers)
    _add_validate_feature_matrix(subparsers)
    args = parser.parse_args()
    args.func(args)


def _add_inspect_archive(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("inspect-archive")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--archive", required=True)
    parser.add_argument("--dataset-id", default="external")
    parser.add_argument("--out", default=None)
    parser.set_defaults(func=_inspect_archive)


def _inspect_archive(args: argparse.Namespace) -> None:
    config = load_config(args.gateway_config)
    out_path = args.out or config.get("outputs", {}).get(
        "external_raw_inventory_tsv",
        "results/tables/external_raw_inventory.tsv",
    )
    inventory = inspect_external_archive(args.archive, dataset_id=args.dataset_id)
    inventory.to_csv(ensure_parent(out_path), sep="\t", index=False)
    roles = ",".join(sorted(inventory["role"].dropna().astype(str).unique())) if not inventory.empty else "none"
    print(f"Wrote external raw inventory: {out_path} ({inventory.shape[0]} files; roles: {roles})")


def _add_score_modules(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("score-modules")
    parser.add_argument("--archive", default="data/external/GSE184117_RAW.tar")
    parser.add_argument("--metadata", default="data/external/GSE184117_series_matrix.txt.gz")
    parser.add_argument("--dataset-id", default="oliva_2022")
    parser.add_argument("--gene-set-config", default="configs/gene_sets.yaml")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--sample-metadata-out", default=None)
    parser.add_argument("--sample-qc-out", default=None)
    parser.add_argument("--module-scores-out", default=None)
    parser.add_argument("--module-contrasts-out", default=None)
    parser.set_defaults(func=_score_modules)


def _score_modules(args: argparse.Namespace) -> None:
    external_config = load_config(args.external_config)
    gene_set_config = load_config(args.gene_set_config)
    gene_set_config["published_gene_lists"] = external_config.get("published_gene_lists", {})
    outputs = external_config.get("outputs", {})
    sample_metadata_out = args.sample_metadata_out or outputs.get(
        "external_sample_metadata_tsv",
        "results/tables/external_sample_metadata.tsv",
    )
    sample_qc_out = args.sample_qc_out or outputs.get(
        "external_10x_sample_qc_tsv",
        "results/tables/external_10x_sample_qc.tsv",
    )
    module_scores_out = args.module_scores_out or outputs.get(
        "external_10x_module_scores_tsv",
        "results/tables/external_10x_module_scores.tsv",
    )
    module_contrasts_out = args.module_contrasts_out or outputs.get(
        "external_10x_module_contrasts_tsv",
        "results/tables/external_10x_module_contrasts.tsv",
    )
    metadata = parse_geo_series_matrix_metadata(args.metadata, dataset_id=args.dataset_id)
    qc, scores, contrasts = score_external_10x_modules(
        args.archive,
        metadata,
        gene_set_config,
        dataset_id=args.dataset_id,
    )
    metadata.to_csv(ensure_parent(sample_metadata_out), sep="\t", index=False)
    qc.to_csv(ensure_parent(sample_qc_out), sep="\t", index=False)
    scores.to_csv(ensure_parent(module_scores_out), sep="\t", index=False)
    contrasts.to_csv(ensure_parent(module_contrasts_out), sep="\t", index=False)
    usable = int(metadata["usable_for_external_validation"].sum()) if "usable_for_external_validation" in metadata else 0
    print(f"Wrote external sample metadata: {sample_metadata_out} ({usable} usable biopsy samples)")
    print(f"Wrote external 10x sample QC: {sample_qc_out} ({qc.shape[0]} samples)")
    print(f"Wrote external 10x module scores: {module_scores_out} ({scores.shape[0]} rows)")
    print(f"Wrote external 10x module contrasts: {module_contrasts_out} ({contrasts.shape[0]} modules)")


def _add_score_markers(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("score-markers")
    parser.add_argument("--archive", default="data/external/GSE184117_RAW.tar")
    parser.add_argument("--metadata", default="data/external/GSE184117_series_matrix.txt.gz")
    parser.add_argument("--dataset-id", default="oliva_2022")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--sample-metadata-out", default=None)
    parser.add_argument("--marker-coverage-out", default=None)
    parser.add_argument("--marker-composition-out", default=None)
    parser.add_argument("--marker-contrasts-out", default=None)
    parser.set_defaults(func=_score_markers)


def _score_markers(args: argparse.Namespace) -> None:
    external_config = load_config(args.external_config)
    outputs = external_config.get("outputs", {})
    sample_metadata_out = args.sample_metadata_out or outputs.get(
        "external_sample_metadata_tsv",
        "results/tables/external_sample_metadata.tsv",
    )
    marker_coverage_out = args.marker_coverage_out or outputs.get(
        "external_10x_marker_coverage_tsv",
        "results/tables/external_10x_marker_coverage.tsv",
    )
    marker_composition_out = args.marker_composition_out or outputs.get(
        "external_10x_marker_composition_tsv",
        "results/tables/external_10x_marker_composition.tsv",
    )
    marker_contrasts_out = args.marker_contrasts_out or outputs.get(
        "external_10x_marker_contrasts_tsv",
        "results/tables/external_10x_marker_contrasts.tsv",
    )
    metadata = parse_geo_series_matrix_metadata(args.metadata, dataset_id=args.dataset_id)
    coverage, composition, contrasts = score_external_10x_marker_composition(
        args.archive,
        metadata,
        external_config.get("marker_panels"),
        dataset_id=args.dataset_id,
    )
    metadata.to_csv(ensure_parent(sample_metadata_out), sep="\t", index=False)
    coverage.to_csv(ensure_parent(marker_coverage_out), sep="\t", index=False)
    composition.to_csv(ensure_parent(marker_composition_out), sep="\t", index=False)
    contrasts.to_csv(ensure_parent(marker_contrasts_out), sep="\t", index=False)
    print(f"Wrote external marker coverage: {marker_coverage_out} ({coverage.shape[0]} rows)")
    print(f"Wrote external marker composition: {marker_composition_out} ({composition.shape[0]} rows)")
    print(f"Wrote external marker contrasts: {marker_contrasts_out} ({contrasts.shape[0]} marker panels)")


def _add_build_mapped(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("build-mapped")
    parser.add_argument("--archive", default="data/external/GSE184117_RAW.tar")
    parser.add_argument("--metadata", default="data/external/GSE184117_series_matrix.txt.gz")
    parser.add_argument("--dataset-id", default="oliva_2022")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--out", default=None)
    parser.add_argument("--mapping-qc-out", default=None)
    parser.add_argument("--donor-features-out", default=None)
    parser.set_defaults(func=_build_mapped)


def _build_mapped(args: argparse.Namespace) -> None:
    external_config = load_config(args.external_config)
    outputs = external_config.get("outputs", {})
    out = args.out or outputs.get("external_10x_mapped_h5ad", "data/processed/gse184117_marker_mapped.h5ad")
    mapping_qc_out = args.mapping_qc_out or outputs.get(
        "external_10x_mapping_qc_tsv",
        "results/tables/external_10x_mapping_qc.tsv",
    )
    donor_features_out = args.donor_features_out or outputs.get(
        "external_10x_mapped_donor_features_tsv",
        "data/processed/gse184117_mapped_donor_features.tsv",
    )
    metadata = parse_geo_series_matrix_metadata(args.metadata, dataset_id=args.dataset_id)
    adata, mapping_qc, donor_features = build_external_10x_marker_mapped_anndata(
        args.archive,
        metadata,
        external_config.get("marker_panels"),
        dataset_id=args.dataset_id,
    )
    ensure_parent(out)
    adata.write_h5ad(out)
    mapping_qc.to_csv(ensure_parent(mapping_qc_out), sep="\t", index=False)
    donor_features.to_csv(ensure_parent(donor_features_out), sep="\t", index=False)
    print(f"Wrote external marker-mapped AnnData: {out} ({adata.n_obs} cells x {adata.n_vars} genes)")


def _add_compare_mapped(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("compare-mapped")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--mapped-features", default=None)
    parser.add_argument("--age-associations", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--direct-feature-map", action="store_true")
    parser.set_defaults(func=_compare_mapped)


def _compare_mapped(args: argparse.Namespace) -> None:
    external_config = load_config(args.external_config)
    gateway_config = load_config(args.gateway_config)
    mapped_path = args.mapped_features or external_config.get("outputs", {}).get(
        "external_10x_mapped_donor_features_tsv",
        "data/processed/gse184117_mapped_donor_features.tsv",
    )
    age_path = args.age_associations or gateway_config.get("outputs", {}).get(
        "age_associations_tsv",
        "results/tables/age_cell_state_associations.tsv",
    )
    out_path = args.out or external_config.get("outputs", {}).get(
        "external_mapped_feature_concordance_tsv",
        "results/tables/external_mapped_feature_concordance.tsv",
    )
    mapped = pd.read_csv(mapped_path, sep="\t")
    age = pd.read_csv(age_path, sep="\t")
    feature_map = _direct_feature_map(mapped, age) if args.direct_feature_map else None
    concordance = external_mapped_feature_concordance(mapped, age, marker_to_gateway=feature_map)
    concordance.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(f"Wrote external mapped-feature concordance: {out_path} ({concordance.shape[0]} rows)")


def _add_compare_marker_age(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("compare-marker-age")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--marker-contrasts", default=None)
    parser.add_argument("--age-associations", default=None)
    parser.add_argument("--out", default=None)
    parser.set_defaults(func=_compare_marker_age)


def _compare_marker_age(args: argparse.Namespace) -> None:
    config = load_config(args.gateway_config)
    outputs = config.get("outputs", {})
    marker_path = args.marker_contrasts or outputs.get(
        "external_10x_marker_contrasts_tsv",
        "results/tables/external_10x_marker_contrasts.tsv",
    )
    age_path = args.age_associations or outputs.get(
        "age_associations_tsv",
        "results/tables/age_cell_state_associations.tsv",
    )
    out_path = args.out or outputs.get(
        "external_marker_age_concordance_tsv",
        "results/tables/external_marker_age_concordance.tsv",
    )
    concordance = external_marker_age_concordance(
        pd.read_csv(marker_path, sep="\t"),
        pd.read_csv(age_path, sep="\t"),
    )
    concordance.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(f"Wrote external marker-age concordance: {out_path} ({concordance.shape[0]} rows)")


def _add_summarize_validation(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("summarize-validation")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--gene-table", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--coverage-out", default=None)
    parser.add_argument("--contract-out", default=None)
    parser.set_defaults(func=_summarize_validation)


def _summarize_validation(args: argparse.Namespace) -> None:
    external_config = load_config(args.external_config)
    gateway_config = load_config(args.gateway_config)
    outputs = external_config.get("outputs", {})
    summary_out = args.summary_out or outputs.get("validation_summary_tsv", "results/tables/external_validation_summary.tsv")
    coverage_out = args.coverage_out or outputs.get("gene_list_coverage_tsv", "results/tables/external_gene_list_coverage.tsv")
    contract_out = args.contract_out or outputs.get("feature_contract_tsv", "results/tables/external_feature_contract.tsv")
    dataset_summary = external_dataset_summary(external_config)
    dataset_summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    feature_matrix_contract_summary(external_config).to_csv(ensure_parent(contract_out), sep="\t", index=False)
    var, var_names = _load_var(args.h5ad or gateway_config.get("source", {}).get("h5ad_path"), args.gene_table)
    coverage = published_gene_list_coverage(
        external_config,
        var,
        var_names,
        external_config.get("score", {}).get("var_symbol_columns", DEFAULT_SYMBOL_COLUMNS),
    )
    coverage.to_csv(ensure_parent(coverage_out), sep="\t", index=False)
    ready = int(dataset_summary["ready_for_feature_validation"].sum()) if "ready_for_feature_validation" in dataset_summary else 0
    print(f"Wrote external validation summary: {summary_out} ({dataset_summary.shape[0]} datasets; {ready} feature-ready)")


def _add_candidate_matrix(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("candidate-matrix")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--out", default=None)
    parser.set_defaults(func=_candidate_matrix)


def _candidate_matrix(args: argparse.Namespace) -> None:
    config = load_config(args.external_config)
    out = args.out or config.get("outputs", {}).get(
        "external_candidate_matrix_tsv",
        "results/tables/external_candidate_matrix.tsv",
    )
    matrix = external_candidate_matrix(config)
    matrix.to_csv(ensure_parent(out), sep="\t", index=False)
    print(f"Wrote external candidate matrix: {out} ({matrix.shape[0]} datasets)")


def _add_evidence(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("evidence")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--validation-summary", default=None)
    parser.add_argument("--sample-metadata", default=None)
    parser.add_argument("--module-contrasts", default=None)
    parser.add_argument("--marker-contrasts", default=None)
    parser.add_argument("--mapped-features", default=None)
    parser.add_argument("--scanvi-features", default=None)
    parser.add_argument("--out", default=None)
    parser.set_defaults(func=_evidence)


def _evidence(args: argparse.Namespace) -> None:
    config = load_config(args.external_config)
    outputs = config.get("outputs", {})
    summary_path = args.validation_summary or outputs.get(
        "validation_summary_tsv",
        "results/tables/external_validation_summary.tsv",
    )
    evidence_out = args.out or outputs.get(
        "external_validation_evidence_tsv",
        "results/tables/external_validation_evidence.tsv",
    )
    dataset_summary = _read_optional_tsv(summary_path)
    if dataset_summary is None:
        dataset_summary = external_dataset_summary(config)
    evidence = external_validation_evidence_summary(
        config,
        dataset_summary,
        sample_metadata=_read_optional_tsv(args.sample_metadata or outputs.get("external_sample_metadata_tsv", "results/tables/external_sample_metadata.tsv")),
        module_contrasts=_read_optional_tsv(args.module_contrasts or outputs.get("external_10x_module_contrasts_tsv", "results/tables/external_10x_module_contrasts.tsv")),
        marker_contrasts=_read_optional_tsv(args.marker_contrasts or outputs.get("external_10x_marker_contrasts_tsv", "results/tables/external_10x_marker_contrasts.tsv")),
        mapped_features=_read_optional_tsv(args.mapped_features or outputs.get("external_10x_mapped_donor_features_tsv", "data/processed/gse184117_mapped_donor_features.tsv")),
        scanvi_features=_read_optional_tsv(args.scanvi_features or outputs.get("external_scanvi_donor_features_tsv", "data/processed/gse184117_scanvi_donor_features.tsv")),
    )
    evidence.to_csv(ensure_parent(evidence_out), sep="\t", index=False)
    print(f"Wrote external validation evidence ledger: {evidence_out} ({evidence.shape[0]} rows)")


def _add_validate_feature_matrix(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("validate-feature-matrix")
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--feature-matrix", default=None)
    parser.add_argument("--gateway-features", default=None)
    parser.add_argument("--dataset-id", default="external_feature_matrix")
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--harmonization-out", default=None)
    parser.set_defaults(func=_validate_feature_matrix)


def _validate_feature_matrix(args: argparse.Namespace) -> None:
    external_config = load_config(args.external_config)
    gateway_config = load_config(args.gateway_config)
    gateway_features_path = args.gateway_features or gateway_config.get("outputs", {}).get(
        "ora_augmented_feature_matrix_tsv",
        "data/processed/ora_augmented_feature_matrix.tsv",
    )
    gateway_features = pd.read_csv(gateway_features_path, sep="\t") if Path(gateway_features_path).exists() else None
    summary, harmonization = validate_external_feature_matrix(
        args.feature_matrix,
        external_config,
        gateway_features=gateway_features,
        dataset_id=args.dataset_id,
    )
    outputs = gateway_config.get("outputs", {})
    summary_out = args.summary_out or outputs.get("external_feature_validation_tsv", "results/tables/external_feature_validation.tsv")
    harmonization_out = args.harmonization_out or outputs.get("external_feature_harmonization_tsv", "results/tables/external_feature_harmonization.tsv")
    summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    harmonization.to_csv(ensure_parent(harmonization_out), sep="\t", index=False)
    status = summary.loc[0, "status"] if not summary.empty else "empty"
    print(f"Wrote external feature validation summary: {summary_out} ({status})")
    print(f"Wrote external feature harmonization: {harmonization_out} ({harmonization.shape[0]} rows)")


def _direct_feature_map(mapped_features: pd.DataFrame, age_associations: pd.DataFrame) -> dict[str, tuple[str, ...]]:
    if "feature" not in age_associations:
        return {}
    gateway_features = set(age_associations["feature"].astype(str))
    mapping: dict[str, list[str]] = {}
    for feature in mapped_features.columns:
        if not str(feature).startswith(("prop__", "clr__")) or str(feature) not in gateway_features:
            continue
        panel = str(feature).split("__", 1)[1]
        mapping.setdefault(panel, []).append(str(feature))
    return {panel: tuple(features) for panel, features in mapping.items()}


def _load_var(h5ad_path: str | None, gene_table_path: str | None) -> tuple[pd.DataFrame, pd.Index]:
    if h5ad_path:
        adata = read_h5ad_backed(h5ad_path)
        try:
            return adata.var.copy(), pd.Index(adata.var_names.astype(str))
        finally:
            close = getattr(adata, "file", None)
            if close is not None:
                close.close()
    if gene_table_path:
        sep = "," if str(gene_table_path).endswith(".csv") else "\t"
        table = pd.read_csv(gene_table_path, sep=sep)
        if "gene_id" in table:
            var = table.set_index("gene_id", drop=False)
            return var, pd.Index(var.index.astype(str))
        if "gene_symbol" in table:
            return table, pd.Index(table["gene_symbol"].astype(str))
        if "feature_name" in table:
            return table, pd.Index(table["feature_name"].astype(str))
        raise ValueError("--gene-table must include gene_id, gene_symbol, or feature_name.")
    raise ValueError("Provide --h5ad or --gene-table to resolve published gene-list coverage.")


def _read_optional_tsv(path: str | None) -> pd.DataFrame | None:
    if not path:
        return None
    candidate = Path(path)
    if not candidate.exists():
        return None
    return pd.read_csv(candidate, sep="\t")


if __name__ == "__main__":
    main()
