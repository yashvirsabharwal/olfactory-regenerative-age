"""External validation registry and gene-list coverage helpers."""

from __future__ import annotations

import csv
import gzip
import re
import tarfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse
from scipy import stats

from .modules import GeneSet, parse_gene_sets, resolve_gene_sets


DEFAULT_EXTERNAL_MARKER_PANELS: dict[str, tuple[str, ...]] = {
    "quiescent_hbc": ("TP63", "KRT5", "KRT14", "CXCL14"),
    "activated_hbc": ("KRT6A", "KRT6B", "KRT16", "SERPINB3"),
    "cycling": ("MKI67", "TOP2A", "UBE2C", "HMGB2"),
    "neural_progenitor": ("ASCL1", "NEUROG1", "NEUROD1", "STMN2"),
    "immature_osn": ("GAP43", "DCX", "TUBB3", "STMN2"),
    "mature_osn": ("OMP", "ADCY3", "GNAL", "RTP1"),
    "sustentacular": ("CYP2A13", "CYP2J2", "MUC1", "SOX2"),
    "bowman_gland": ("LTF", "BPIFB1", "MUC5B", "BPIFA1"),
    "respiratory_ciliated": ("FOXJ1", "PIFO", "TPPP3", "DNAH5"),
    "goblet_secretory": ("MUC5AC", "MUC5B", "SPDEF", "AGR2"),
    "immune": ("PTPRC", "LST1", "TYROBP", "CD3D", "MS4A1"),
}

DEFAULT_EXTERNAL_MARKER_GATEWAY_MAP: dict[str, tuple[str, ...]] = {
    "quiescent_hbc": ("prop__quiescent_hbc", "clr__quiescent_hbc"),
    "cycling": ("prop__cycling_hbc", "clr__cycling_hbc"),
    "neural_progenitor": ("prop__early_inp", "clr__early_inp", "prop__late_inp", "clr__late_inp"),
    "immature_osn": ("prop__early_iosn", "clr__early_iosn", "prop__late_iosn", "clr__late_iosn"),
    "mature_osn": ("prop__early_mature_mosn", "clr__early_mature_mosn", "prop__fully_mature_mosn", "clr__fully_mature_mosn"),
    "sustentacular": ("prop__olf_sus", "clr__olf_sus"),
    "bowman_gland": ("prop__bowman_gland", "clr__bowman_gland"),
    "respiratory_ciliated": ("prop__multiciliated", "clr__multiciliated", "prop__deuterosomal", "clr__deuterosomal"),
    "goblet_secretory": ("prop__goblet", "clr__goblet", "prop__proliferating_secretory", "clr__proliferating_secretory"),
    "immune": ("prop__inflammatory", "clr__inflammatory", "prop__antigen_presenting", "clr__antigen_presenting"),
}


def external_dataset_summary(config: dict[str, Any], base_dir: str | Path = ".") -> pd.DataFrame:
    """Summarize configured external datasets and whether required files exist."""

    rows = []
    root = Path(base_dir)
    for dataset_id, spec in config.get("datasets", {}).items():
        required_files = spec.get("required_files", {}) if isinstance(spec, dict) else {}
        resolved_files = {
            name: _resolve_optional_path(path, root)
            for name, path in required_files.items()
        }
        missing = [name for name, path in resolved_files.items() if path is None or not path.exists()]
        provided = [name for name, path in resolved_files.items() if path is not None and path.exists()]
        rows.append(
            {
                "dataset_id": dataset_id,
                "title": spec.get("title", dataset_id),
                "status": spec.get("status", "unknown"),
                "validation_use": spec.get("validation_use", ""),
                "species": spec.get("species", ""),
                "tissue": spec.get("tissue", ""),
                "disease_context": ",".join(str(item) for item in spec.get("disease_context", [])),
                "expected_level": spec.get("expected_level", ""),
                "files_provided": ",".join(provided),
                "files_missing": ",".join(missing),
                "ready_for_feature_validation": "feature_matrix" in provided,
                "ready_for_raw_adapter": bool({"expression", "metadata"}.issubset(set(provided))),
                "readiness_class": _readiness_class(spec, provided, missing),
                "notes": spec.get("notes", ""),
            }
        )
    return pd.DataFrame(rows)


def inspect_external_archive(archive_path: str | Path, dataset_id: str = "external") -> pd.DataFrame:
    """Inventory a raw external TAR archive without extracting it."""

    path = Path(archive_path)
    if not path.exists():
        return pd.DataFrame(
            [
                {
                    "dataset_id": dataset_id,
                    "archive": str(path),
                    "member": "",
                    "size_bytes": 0,
                    "role": "missing_archive",
                    "sample_guess": "",
                    "status": "missing_archive",
                }
            ]
        )

    rows = []
    with tarfile.open(path) as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            role = _archive_member_role(member.name)
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "archive": str(path),
                    "member": member.name,
                    "size_bytes": int(member.size),
                    "role": role,
                    "sample_guess": _sample_guess(member.name, role),
                    "status": "ok",
                }
            )
    if not rows:
        rows.append(
            {
                "dataset_id": dataset_id,
                "archive": str(path),
                "member": "",
                "size_bytes": 0,
                "role": "empty_archive",
                "sample_guess": "",
                "status": "empty_archive",
            }
        )
    return pd.DataFrame(rows).sort_values(["sample_guess", "role", "member"]).reset_index(drop=True)


def parse_geo_series_matrix_metadata(
    series_matrix_path: str | Path,
    *,
    dataset_id: str = "external",
) -> pd.DataFrame:
    """Parse GEO series-matrix sample metadata into a donor/sample manifest."""

    path = Path(series_matrix_path)
    if not path.exists():
        return pd.DataFrame(
            [
                {
                    "dataset_id": dataset_id,
                    "sample_id": "",
                    "donor_id": "",
                    "age": np.nan,
                    "disease_state": "",
                    "disease_group": "",
                    "sample_class": "",
                    "status": "missing_metadata",
                }
            ]
        )

    sample_rows = _read_geo_sample_rows(path)
    accessions = _first_geo_row(sample_rows, "Sample_geo_accession")
    if not accessions:
        return pd.DataFrame(columns=_geo_metadata_columns())

    rows = []
    for idx, sample_id in enumerate(accessions):
        title = _geo_value(sample_rows, "Sample_title", idx)
        description = _geo_join(sample_rows, "Sample_description", idx)
        characteristics = _geo_characteristics(sample_rows, idx)
        supplementary = [
            _geo_value(sample_rows, f"Sample_supplementary_file_{pos}", idx)
            for pos in range(1, 5)
            if _geo_value(sample_rows, f"Sample_supplementary_file_{pos}", idx)
        ]
        sample_prefix = _sample_prefix_from_urls(supplementary)
        subject_id = _subject_id_from_prefix(sample_prefix)
        disease_state = characteristics.get("disease state", "")
        sample_class = _external_sample_class(title, description, sample_prefix)
        age = pd.to_numeric(characteristics.get("age", ""), errors="coerce")
        disease_group = _external_disease_group(disease_state, sample_class)
        rows.append(
            {
                "dataset_id": dataset_id,
                "sample_id": str(sample_id),
                "donor_id": subject_id or str(sample_id),
                "subject_id": subject_id,
                "sample_title": title,
                "sample_description": description,
                "age": age,
                "disease_state": disease_state,
                "disease_group": disease_group,
                "sample_class": sample_class,
                "is_presbyosmic": disease_group == "presbyosmia",
                "usable_for_external_validation": bool(sample_class != "culture" and pd.notna(age) and disease_group in {"healthy", "presbyosmia"}),
                "source_name": _geo_value(sample_rows, "Sample_source_name_ch1", idx),
                "organism": _geo_value(sample_rows, "Sample_organism_ch1", idx),
                "instrument_model": _geo_value(sample_rows, "Sample_instrument_model", idx),
                "library_strategy": _geo_value(sample_rows, "Sample_library_strategy", idx),
                "chemistry": "10x_3p_v3_1",
                "collection_method": "biopsy_dissociation",
                "sample_prefix": sample_prefix,
                "supplementary_files": ",".join(supplementary),
                "status": "ok",
            }
        )
    return pd.DataFrame(rows, columns=_geo_metadata_columns())


def score_external_10x_modules(
    archive_path: str | Path,
    sample_metadata: pd.DataFrame,
    gene_set_config: dict[str, Any],
    *,
    dataset_id: str = "external",
    pseudocount_scale: float = 10_000.0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score configured gene modules from raw external 10x matrices."""

    archive_path = Path(archive_path)
    gene_sets = parse_gene_sets(gene_set_config)
    if not gene_sets:
        raise ValueError("No gene sets configured for external 10x scoring.")
    if not archive_path.exists():
        empty_qc = pd.DataFrame(columns=_external_10x_qc_columns())
        empty_scores = pd.DataFrame(columns=_external_10x_score_columns())
        empty_contrasts = pd.DataFrame(columns=_external_10x_contrast_columns())
        return empty_qc, empty_scores, empty_contrasts

    metadata = sample_metadata.copy()
    if "sample_id" not in metadata:
        metadata["sample_id"] = ""
    metadata_by_accession = metadata.set_index("sample_id", drop=False)
    metadata_by_prefix = (
        metadata[metadata.get("sample_prefix", pd.Series("", index=metadata.index)).astype(str).ne("")]
        .drop_duplicates("sample_prefix")
        .set_index("sample_prefix", drop=False)
    )
    inventory = inspect_external_archive(archive_path, dataset_id=dataset_id)
    sample_roles = _archive_sample_roles(inventory)
    qc_rows = []
    score_rows = []
    with tarfile.open(archive_path) as archive:
        for sample_prefix, roles in sorted(sample_roles.items()):
            if not {"matrix", "features", "barcodes"}.issubset(roles):
                continue
            meta = _metadata_for_external_prefix(sample_prefix, metadata_by_accession, metadata_by_prefix)
            matrix = _read_10x_matrix_from_tar(archive, roles["matrix"])
            var = _read_10x_features_from_tar(archive, roles["features"])
            n_barcodes = _count_gzip_lines_from_tar(archive, roles["barcodes"])
            gene_means = _mean_log1p_cpm_by_gene(matrix, scale=pseudocount_scale)
            resolved, coverage = resolve_gene_sets(var, pd.Index(var["gene_id"].astype(str)), gene_sets, ["feature_name", "gene_id"])
            sample_id = str(meta.get("sample_id", _accession_from_prefix(sample_prefix)))
            qc_rows.append(
                {
                    "dataset_id": dataset_id,
                    "sample_id": sample_id,
                    "donor_id": meta.get("donor_id", ""),
                    "sample_prefix": sample_prefix,
                    "age": meta.get("age", np.nan),
                    "disease_state": meta.get("disease_state", ""),
                    "disease_group": meta.get("disease_group", ""),
                    "sample_class": meta.get("sample_class", ""),
                    "n_cells": int(matrix.shape[1]),
                    "n_barcodes": int(n_barcodes),
                    "n_genes": int(matrix.shape[0]),
                    "detected_genes": int((np.asarray(matrix.sum(axis=1)).ravel() > 0).sum()),
                    "total_counts": float(matrix.sum()),
                    "median_counts_per_cell": float(np.median(np.asarray(matrix.sum(axis=0)).ravel())),
                    "status": "ok",
                }
            )
            for _, cov in coverage.iterrows():
                module = str(cov["module"])
                indices = resolved.get(module, [])
                score = float(np.nanmean(gene_means[indices])) if indices else np.nan
                score_rows.append(
                    {
                        "dataset_id": dataset_id,
                        "sample_id": sample_id,
                        "donor_id": meta.get("donor_id", ""),
                        "sample_prefix": sample_prefix,
                        "age": meta.get("age", np.nan),
                        "disease_state": meta.get("disease_state", ""),
                        "disease_group": meta.get("disease_group", ""),
                        "sample_class": meta.get("sample_class", ""),
                        "module": module,
                        "description": cov.get("description", ""),
                        "n_requested": int(cov["n_requested"]),
                        "n_present": int(cov["n_present"]),
                        "coverage_fraction": float(cov["coverage_fraction"]),
                        "missing_genes": cov.get("missing_genes", ""),
                        "mean_log1p_cpm": score,
                        "status": "ok" if indices else "no_genes_present",
                    }
                )
    qc = pd.DataFrame(qc_rows, columns=_external_10x_qc_columns())
    scores = pd.DataFrame(score_rows, columns=_external_10x_score_columns())
    contrasts = external_module_contrasts(scores)
    return qc, scores, contrasts


def score_external_10x_marker_composition(
    archive_path: str | Path,
    sample_metadata: pd.DataFrame,
    marker_panels: dict[str, Any] | None = None,
    *,
    dataset_id: str = "external",
    pseudocount_scale: float = 10_000.0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Estimate coarse marker-panel composition from raw external 10x matrices."""

    archive_path = Path(archive_path)
    panels = _normalize_marker_panels(marker_panels)
    if not panels:
        raise ValueError("No marker panels configured for external 10x marker composition.")
    if not archive_path.exists():
        empty_coverage = pd.DataFrame(columns=_external_marker_coverage_columns())
        empty_composition = pd.DataFrame(columns=_external_marker_composition_columns())
        empty_contrasts = pd.DataFrame(columns=_external_marker_contrast_columns())
        return empty_coverage, empty_composition, empty_contrasts

    metadata = sample_metadata.copy()
    if "sample_id" not in metadata:
        metadata["sample_id"] = ""
    metadata_by_accession = metadata.set_index("sample_id", drop=False)
    metadata_by_prefix = (
        metadata[metadata.get("sample_prefix", pd.Series("", index=metadata.index)).astype(str).ne("")]
        .drop_duplicates("sample_prefix")
        .set_index("sample_prefix", drop=False)
    )
    inventory = inspect_external_archive(archive_path, dataset_id=dataset_id)
    sample_roles = _archive_sample_roles(inventory)
    coverage_rows = []
    composition_rows = []
    with tarfile.open(archive_path) as archive:
        for sample_prefix, roles in sorted(sample_roles.items()):
            if not {"matrix", "features", "barcodes"}.issubset(roles):
                continue
            meta = _metadata_for_external_prefix(sample_prefix, metadata_by_accession, metadata_by_prefix)
            matrix = _read_10x_matrix_from_tar(archive, roles["matrix"])
            var = _read_10x_features_from_tar(archive, roles["features"])
            sample_id = str(meta.get("sample_id", _accession_from_prefix(sample_prefix)))
            panel_scores, coverage = _score_marker_panels_by_cell(
                matrix,
                var,
                panels,
                scale=pseudocount_scale,
            )
            total_cells = int(matrix.shape[1])
            assignments, max_scores, margins = _assign_marker_panels(panel_scores)
            for panel_name, coverage_row in coverage.items():
                coverage_rows.append(
                    {
                        "dataset_id": dataset_id,
                        "sample_id": sample_id,
                        "sample_prefix": sample_prefix,
                        "marker_panel": panel_name,
                        **coverage_row,
                    }
                )
            for panel_name in [*panels.keys(), "unassigned_marker_low"]:
                if panel_name == "unassigned_marker_low":
                    mask = assignments == panel_name
                    panel_score = np.full(total_cells, np.nan)
                    n_present = 0
                else:
                    mask = assignments == panel_name
                    panel_score = panel_scores.get(panel_name, np.full(total_cells, np.nan))
                    n_present = int(coverage[panel_name]["n_present"])
                n_assigned = int(mask.sum())
                composition_rows.append(
                    {
                        "dataset_id": dataset_id,
                        "sample_id": sample_id,
                        "donor_id": meta.get("donor_id", ""),
                        "sample_prefix": sample_prefix,
                        "age": meta.get("age", np.nan),
                        "disease_state": meta.get("disease_state", ""),
                        "disease_group": meta.get("disease_group", ""),
                        "sample_class": meta.get("sample_class", ""),
                        "marker_panel": panel_name,
                        "n_cells_assigned": n_assigned,
                        "total_cells": total_cells,
                        "fraction_cells": float(n_assigned / total_cells) if total_cells else np.nan,
                        "mean_marker_score": float(np.nanmean(panel_score[mask])) if n_assigned and panel_name in panels else np.nan,
                        "median_marker_score": float(np.nanmedian(panel_score[mask])) if n_assigned and panel_name in panels else np.nan,
                        "mean_assignment_margin": float(np.nanmean(margins[mask])) if n_assigned else np.nan,
                        "mean_top_marker_score": float(np.nanmean(max_scores[mask])) if n_assigned else np.nan,
                        "n_markers_present": n_present,
                        "status": "marker_only",
                    }
                )
    coverage = pd.DataFrame(coverage_rows, columns=_external_marker_coverage_columns())
    composition = pd.DataFrame(composition_rows, columns=_external_marker_composition_columns())
    contrasts = external_marker_composition_contrasts(composition)
    return coverage, composition, contrasts


def build_external_10x_marker_mapped_anndata(
    archive_path: str | Path,
    sample_metadata: pd.DataFrame,
    marker_panels: dict[str, Any] | None = None,
    *,
    dataset_id: str = "external",
    pseudocount_scale: float = 10_000.0,
) -> tuple[Any, pd.DataFrame, pd.DataFrame]:
    """Build a harmonized external AnnData with marker-reference cell-state labels."""

    try:
        import anndata as ad  # type: ignore
    except ModuleNotFoundError as exc:
        raise SystemExit("anndata is required to build external 10x AnnData.") from exc

    archive_path = Path(archive_path)
    panels = _normalize_marker_panels(marker_panels)
    metadata = sample_metadata.copy()
    metadata_by_accession = metadata.set_index("sample_id", drop=False) if "sample_id" in metadata else pd.DataFrame()
    metadata_by_prefix = (
        metadata[metadata.get("sample_prefix", pd.Series("", index=metadata.index)).astype(str).ne("")]
        .drop_duplicates("sample_prefix")
        .set_index("sample_prefix", drop=False)
        if not metadata.empty
        else pd.DataFrame()
    )
    inventory = inspect_external_archive(archive_path, dataset_id=dataset_id)
    sample_roles = _archive_sample_roles(inventory)
    matrices = []
    obs_frames = []
    reference_var: pd.DataFrame | None = None
    qc_rows = []
    with tarfile.open(archive_path) as archive:
        for sample_prefix, roles in sorted(sample_roles.items()):
            if not {"matrix", "features", "barcodes"}.issubset(roles):
                continue
            meta = _metadata_for_external_prefix(sample_prefix, metadata_by_accession, metadata_by_prefix)
            if bool(meta.get("usable_for_external_validation", False)) is False:
                continue
            matrix = _read_10x_matrix_from_tar(archive, roles["matrix"])
            var = _read_10x_features_from_tar(archive, roles["features"])
            barcodes = _read_10x_barcodes_from_tar(archive, roles["barcodes"])
            if reference_var is None:
                reference_var = var.copy()
            elif not _features_compatible(reference_var, var):
                raise ValueError(f"10x feature table for {sample_prefix} does not match the first sample.")
            panel_scores, coverage = _score_marker_panels_by_cell(matrix, var, panels, scale=pseudocount_scale)
            assignments, max_scores, margins = _assign_marker_panels(panel_scores)
            sample_id = str(meta.get("sample_id", _accession_from_prefix(sample_prefix)))
            obs_frames.append(
                pd.DataFrame(
                    {
                        "dataset_id": dataset_id,
                        "sample_id": sample_id,
                        "donor_id": meta.get("donor_id", ""),
                        "sample_prefix": sample_prefix,
                        "barcode": barcodes,
                        "age": meta.get("age", np.nan),
                        "disease_state": meta.get("disease_state", ""),
                        "disease_group": meta.get("disease_group", ""),
                        "sample_class": meta.get("sample_class", ""),
                        "chemistry": meta.get("chemistry", ""),
                        "collection_method": meta.get("collection_method", ""),
                        "mapped_cell_state": assignments,
                        "mapping_confidence": max_scores,
                        "mapping_margin": margins,
                    }
                )
            )
            matrices.append(matrix.T.tocsr())
            present_markers = sorted(
                {
                    gene
                    for row in coverage.values()
                    for gene in str(row.get("present_genes", "")).split(",")
                    if gene
                }
            )
            qc_rows.append(
                {
                    "dataset_id": dataset_id,
                    "sample_id": sample_id,
                    "donor_id": meta.get("donor_id", ""),
                    "sample_prefix": sample_prefix,
                    "age": meta.get("age", np.nan),
                    "disease_group": meta.get("disease_group", ""),
                    "sample_class": meta.get("sample_class", ""),
                    "n_cells": int(matrix.shape[1]),
                    "n_genes": int(matrix.shape[0]),
                    "mapped_fraction": float(np.mean(assignments != "unassigned_marker_low")),
                    "mean_mapping_confidence": float(np.mean(max_scores)),
                    "mean_mapping_margin": float(np.mean(margins)),
                    "n_marker_genes_present": int(len(present_markers)),
                    "status": "marker_reference_mapped",
                }
            )
    if not matrices or reference_var is None:
        empty = ad.AnnData(X=sparse.csr_matrix((0, 0)))
        return empty, pd.DataFrame(columns=_external_mapping_qc_columns()), pd.DataFrame()
    obs = pd.concat(obs_frames, ignore_index=True)
    obs.index = _unique_obs_names(obs["sample_id"], obs["barcode"])
    var = reference_var.copy()
    var.index = var["gene_id"].astype(str)
    adata = ad.AnnData(X=sparse.vstack(matrices, format="csr"), obs=obs, var=var)
    adata.uns["dataset_id"] = dataset_id
    adata.uns["mapping_method"] = "marker_reference"
    qc = pd.DataFrame(qc_rows, columns=_external_mapping_qc_columns())
    donor_features = external_mapped_donor_features(obs)
    return adata, qc, donor_features


def external_mapped_donor_features(obs: pd.DataFrame) -> pd.DataFrame:
    """Convert mapped external cell labels into ORA-compatible donor/sample features."""

    if obs.empty or "mapped_cell_state" not in obs:
        return pd.DataFrame()
    rows = []
    for keys, group in obs.groupby(["dataset_id", "sample_id", "donor_id", "age", "disease_group"], dropna=False, observed=True):
        dataset_id, sample_id, donor_id, age, disease_group = keys
        counts = group["mapped_cell_state"].astype(str).value_counts()
        total = float(counts.sum())
        fractions = counts / total if total else counts.astype(float)
        feature_row: dict[str, Any] = {
            "dataset_id": dataset_id,
            "sample_id": sample_id,
            "donor_id": donor_id,
            "age": age,
            "disease_group": disease_group,
            "n_cells": int(total),
        }
        panels = sorted(label for label in fractions.index if label != "unassigned_marker_low")
        if panels:
            geometric = float(np.exp(np.mean(np.log([float(fractions.get(panel, 0.0)) + 1e-6 for panel in panels]))))
        else:
            geometric = np.nan
        for panel in panels:
            fraction = float(fractions.get(panel, 0.0))
            feature_row[f"prop__{panel}"] = fraction
            feature_row[f"clr__{panel}"] = float(np.log((fraction + 1e-6) / geometric)) if geometric > 0 else np.nan
        feature_row["ratio__mature_osn_to_quiescent_hbc"] = _safe_ratio(fractions, "mature_osn", "quiescent_hbc")
        feature_row["ratio__immature_osn_to_quiescent_hbc"] = _safe_ratio(fractions, "immature_osn", "quiescent_hbc")
        feature_row["ratio__mature_osn_to_sustentacular"] = _safe_ratio(fractions, "mature_osn", "sustentacular")
        rows.append(feature_row)
    return pd.DataFrame(rows).sort_values(["dataset_id", "sample_id"]).reset_index(drop=True)


def external_mapped_feature_concordance(
    mapped_features: pd.DataFrame,
    age_associations: pd.DataFrame,
    marker_to_gateway: dict[str, tuple[str, ...]] | None = None,
) -> pd.DataFrame:
    """Compare marker-reference mapped external features with Gateway age-association directions."""

    columns = _external_mapped_feature_concordance_columns()
    if mapped_features.empty or age_associations.empty:
        return pd.DataFrame(columns=columns)
    mapping = marker_to_gateway or DEFAULT_EXTERNAL_MARKER_GATEWAY_MAP
    assoc = age_associations.copy()
    if "feature" not in assoc:
        return pd.DataFrame(columns=columns)
    assoc = assoc.set_index("feature", drop=False)
    rows = []
    feature_cols = [col for col in mapped_features.columns if str(col).startswith(("prop__", "clr__"))]
    disease = mapped_features.get("disease_group", pd.Series("", index=mapped_features.index)).astype(str)
    for external_feature in feature_cols:
        panel = str(external_feature).split("__", 1)[1]
        healthy = pd.to_numeric(mapped_features.loc[disease.eq("healthy"), external_feature], errors="coerce").dropna().to_numpy()
        case = pd.to_numeric(mapped_features.loc[disease.eq("presbyosmia"), external_feature], errors="coerce").dropna().to_numpy()
        delta = float(np.mean(case) - np.mean(healthy)) if healthy.size and case.size else np.nan
        external_direction = "positive" if delta > 0 else "negative" if delta < 0 else "flat"
        p_value = _mann_whitney_pvalue(case, healthy)
        candidate_gateway_features = [feature for feature in mapping.get(panel, ()) if feature.startswith(external_feature.split("__", 1)[0])]
        if not candidate_gateway_features:
            candidate_gateway_features = list(mapping.get(panel, ()))
        for gateway_feature in candidate_gateway_features:
            if gateway_feature not in assoc.index:
                rows.append(
                    _empty_mapped_feature_concordance_row(
                        external_feature=external_feature,
                        gateway_feature=gateway_feature,
                        external_direction=external_direction,
                        delta=delta,
                        p_value=p_value,
                        n_healthy=healthy.size,
                        n_case=case.size,
                        status="missing_gateway_age_association",
                    )
                )
                continue
            age = assoc.loc[gateway_feature]
            gateway_direction = str(age.get("direction", ""))
            concordance = _direction_concordance(external_direction, gateway_direction)
            rows.append(
                {
                    "external_feature": external_feature,
                    "gateway_feature": gateway_feature,
                    "external_direction": external_direction,
                    "gateway_age_direction": gateway_direction,
                    "concordance": concordance,
                    "external_delta": delta,
                    "external_p_value": p_value,
                    "gateway_beta_per_10_years": pd.to_numeric(age.get("beta_per_10_years"), errors="coerce"),
                    "gateway_p_value": pd.to_numeric(age.get("p_value"), errors="coerce"),
                    "gateway_fdr": pd.to_numeric(age.get("fdr"), errors="coerce"),
                    "n_healthy": int(healthy.size),
                    "n_presbyosmia": int(case.size),
                    "status": _mapped_feature_concordance_status(concordance, healthy.size, case.size),
                    "interpretation": _mapped_feature_concordance_interpretation(
                        external_feature,
                        gateway_feature,
                        concordance,
                    ),
                }
            )
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["concordance", "gateway_fdr", "external_p_value", "external_feature", "gateway_feature"],
        ascending=[True, True, True, True, True],
        na_position="last",
    ).reset_index(drop=True)


def external_marker_composition_contrasts(composition: pd.DataFrame) -> pd.DataFrame:
    """Summarize descriptive marker-composition contrasts for external 10x samples."""

    columns = _external_marker_contrast_columns()
    if composition.empty:
        return pd.DataFrame(columns=columns)
    frame = composition.copy()
    frame = frame[frame.get("sample_class", "").astype(str).eq("biopsy")]
    frame = frame[frame["disease_group"].isin(["healthy", "presbyosmia"])]
    frame["fraction_cells"] = pd.to_numeric(frame["fraction_cells"], errors="coerce")
    rows = []
    for panel, group in frame.groupby("marker_panel", observed=True):
        healthy = group[group["disease_group"].eq("healthy")]["fraction_cells"].dropna().to_numpy(dtype=float)
        presby = group[group["disease_group"].eq("presbyosmia")]["fraction_cells"].dropna().to_numpy(dtype=float)
        diff = float(np.mean(presby) - np.mean(healthy)) if healthy.size and presby.size else np.nan
        p_value = _mann_whitney_pvalue(presby, healthy)
        rows.append(
            {
                "marker_panel": panel,
                "n_healthy": int(healthy.size),
                "n_presbyosmia": int(presby.size),
                "mean_fraction_healthy": float(np.mean(healthy)) if healthy.size else np.nan,
                "mean_fraction_presbyosmia": float(np.mean(presby)) if presby.size else np.nan,
                "presbyosmia_minus_healthy": diff,
                "p_value": p_value,
                "direction": "higher_in_presbyosmia" if diff > 0 else "lower_in_presbyosmia" if diff < 0 else "flat",
                "status": "marker_only_small_n" if min(healthy.size, presby.size) < 5 else "marker_only",
            }
        )
    return pd.DataFrame(rows, columns=columns).sort_values(["p_value", "marker_panel"], na_position="last").reset_index(drop=True)


def external_marker_age_concordance(
    marker_contrasts: pd.DataFrame,
    age_associations: pd.DataFrame,
    marker_to_gateway: dict[str, tuple[str, ...]] | None = None,
) -> pd.DataFrame:
    """Compare external marker-panel shifts with Gateway donor-level age-association directions."""

    columns = _external_marker_age_concordance_columns()
    if marker_contrasts.empty or age_associations.empty:
        return pd.DataFrame(columns=columns)
    mapping = marker_to_gateway or DEFAULT_EXTERNAL_MARKER_GATEWAY_MAP
    assoc = age_associations.copy()
    if "feature" not in assoc:
        return pd.DataFrame(columns=columns)
    assoc = assoc.set_index("feature", drop=False)
    rows = []
    for _, contrast in marker_contrasts.iterrows():
        panel = str(contrast.get("marker_panel", ""))
        external_direction = _external_direction_to_age_direction(str(contrast.get("direction", "")))
        for gateway_feature in mapping.get(panel, ()):
            if gateway_feature not in assoc.index:
                rows.append(
                    _empty_marker_age_concordance_row(
                        contrast,
                        panel=panel,
                        gateway_feature=gateway_feature,
                        external_direction=external_direction,
                        status="missing_gateway_age_association",
                    )
                )
                continue
            age = assoc.loc[gateway_feature]
            gateway_direction = str(age.get("direction", ""))
            concordance = _direction_concordance(external_direction, gateway_direction)
            rows.append(
                {
                    "marker_panel": panel,
                    "gateway_feature": gateway_feature,
                    "external_direction": external_direction,
                    "gateway_age_direction": gateway_direction,
                    "concordance": concordance,
                    "external_delta": pd.to_numeric(contrast.get("presbyosmia_minus_healthy"), errors="coerce"),
                    "external_p_value": pd.to_numeric(contrast.get("p_value"), errors="coerce"),
                    "gateway_beta_per_10_years": pd.to_numeric(age.get("beta_per_10_years"), errors="coerce"),
                    "gateway_p_value": pd.to_numeric(age.get("p_value"), errors="coerce"),
                    "gateway_fdr": pd.to_numeric(age.get("fdr"), errors="coerce"),
                    "n_healthy": _safe_int(contrast.get("n_healthy")),
                    "n_presbyosmia": _safe_int(contrast.get("n_presbyosmia")),
                    "status": _marker_age_concordance_status(concordance, contrast.get("status", "")),
                    "interpretation": _marker_age_concordance_interpretation(panel, gateway_feature, concordance),
                }
            )
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["concordance", "gateway_fdr", "external_p_value", "marker_panel", "gateway_feature"],
        ascending=[True, True, True, True, True],
        na_position="last",
    ).reset_index(drop=True)


def external_module_contrasts(scores: pd.DataFrame) -> pd.DataFrame:
    """Summarize descriptive module contrasts between presbyosmic and healthy samples."""

    columns = _external_10x_contrast_columns()
    if scores.empty:
        return pd.DataFrame(columns=columns)
    frame = scores.copy()
    frame = frame[frame.get("sample_class", "").astype(str).eq("biopsy")]
    frame = frame[frame["disease_group"].isin(["healthy", "presbyosmia"])]
    frame["mean_log1p_cpm"] = pd.to_numeric(frame["mean_log1p_cpm"], errors="coerce")
    rows = []
    for module, group in frame.groupby("module", observed=True):
        healthy = group[group["disease_group"].eq("healthy")]["mean_log1p_cpm"].dropna().to_numpy(dtype=float)
        presby = group[group["disease_group"].eq("presbyosmia")]["mean_log1p_cpm"].dropna().to_numpy(dtype=float)
        diff = float(np.mean(presby) - np.mean(healthy)) if healthy.size and presby.size else np.nan
        p_value = _mann_whitney_pvalue(presby, healthy)
        rows.append(
            {
                "module": module,
                "n_healthy": int(healthy.size),
                "n_presbyosmia": int(presby.size),
                "mean_healthy": float(np.mean(healthy)) if healthy.size else np.nan,
                "mean_presbyosmia": float(np.mean(presby)) if presby.size else np.nan,
                "presbyosmia_minus_healthy": diff,
                "p_value": p_value,
                "direction": "higher_in_presbyosmia" if diff > 0 else "lower_in_presbyosmia" if diff < 0 else "flat",
                "status": "descriptive_small_n" if min(healthy.size, presby.size) < 5 else "ok",
            }
        )
    return pd.DataFrame(rows, columns=columns).sort_values(["p_value", "module"], na_position="last").reset_index(drop=True)


def validate_external_feature_matrix(
    feature_matrix_path: str | Path | None,
    config: dict[str, Any],
    *,
    gateway_features: pd.DataFrame | None = None,
    dataset_id: str = "external",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate a donor-level external feature matrix and compare feature coverage to Gateway."""

    contract = config.get("feature_matrix_contract", {})
    required = [str(col) for col in contract.get("required_columns", [])]
    optional = [str(col) for col in contract.get("optional_covariates", [])]
    prefixes = [str(prefix) for prefix in contract.get("accepted_feature_prefixes", [])]
    path = Path(feature_matrix_path) if feature_matrix_path not in {None, ""} else None
    if path is None or not path.exists():
        summary = pd.DataFrame(
            [
                {
                    "dataset_id": dataset_id,
                    "path": "" if path is None else str(path),
                    "rows": 0,
                    "columns": 0,
                    "donors": 0,
                    "feature_columns": 0,
                    "missing_required_columns": ",".join(required),
                    "optional_covariates_present": "",
                    "status": "missing_feature_matrix",
                }
            ]
        )
        return summary, pd.DataFrame(columns=_harmonization_columns())

    frame = _read_table(path)
    columns = set(frame.columns.astype(str))
    missing_required = [col for col in required if col not in columns]
    feature_cols = _feature_columns(frame, prefixes)
    gateway_cols = set(_feature_columns(gateway_features, prefixes)) if gateway_features is not None else set()
    harmonization = _feature_harmonization(feature_cols, gateway_cols, dataset_id=dataset_id)
    summary = pd.DataFrame(
        [
            {
                "dataset_id": dataset_id,
                "path": str(path),
                "rows": int(frame.shape[0]),
                "columns": int(frame.shape[1]),
                "donors": int(frame["donor_id"].nunique()) if "donor_id" in frame else 0,
                "feature_columns": int(len(feature_cols)),
                "missing_required_columns": ",".join(missing_required),
                "optional_covariates_present": ",".join(col for col in optional if col in columns),
                "status": "ok" if not missing_required and feature_cols else "invalid_contract",
            }
        ]
    )
    return summary, harmonization


def parse_published_gene_lists(config: dict[str, Any]) -> list[GeneSet]:
    """Parse published validation gene lists from external dataset config."""

    gene_sets = []
    for name, spec in config.get("published_gene_lists", {}).items():
        genes = spec.get("genes", []) if isinstance(spec, dict) else spec
        description = spec.get("description", "") if isinstance(spec, dict) else ""
        cleaned = tuple(str(gene).strip() for gene in genes if str(gene).strip())
        if cleaned:
            gene_sets.append(GeneSet(name=str(name), genes=cleaned, description=str(description)))
    return gene_sets


def published_gene_list_coverage(
    config: dict[str, Any],
    var: pd.DataFrame,
    var_names: pd.Index,
    symbol_columns: list[str] | tuple[str, ...],
) -> pd.DataFrame:
    """Resolve published validation gene lists against a matrix gene index."""

    gene_sets = parse_published_gene_lists(config)
    if not gene_sets:
        return pd.DataFrame(
            columns=["module", "description", "n_requested", "n_present", "coverage_fraction", "present_genes", "missing_genes"]
        )
    _, coverage = resolve_gene_sets(var, var_names, gene_sets, symbol_columns)
    return coverage.rename(columns={"module": "gene_list"})


def feature_matrix_contract_summary(config: dict[str, Any]) -> pd.DataFrame:
    """Return the configured external feature matrix contract as a table."""

    contract = config.get("feature_matrix_contract", {})
    rows = []
    for column in contract.get("required_columns", []):
        rows.append({"field": column, "kind": "required_column"})
    for column in contract.get("optional_covariates", []):
        rows.append({"field": column, "kind": "optional_covariate"})
    for prefix in contract.get("accepted_feature_prefixes", []):
        rows.append({"field": prefix, "kind": "accepted_feature_prefix"})
    return pd.DataFrame(rows)


def external_validation_evidence_summary(
    config: dict[str, Any],
    dataset_summary: pd.DataFrame,
    *,
    sample_metadata: pd.DataFrame | None = None,
    module_contrasts: pd.DataFrame | None = None,
    marker_contrasts: pd.DataFrame | None = None,
    mapped_features: pd.DataFrame | None = None,
    scanvi_features: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build a claim-gated evidence ledger for external validation datasets."""

    rows: list[dict[str, Any]] = []
    summary_by_dataset = (
        dataset_summary.set_index("dataset_id", drop=False)
        if dataset_summary is not None and not dataset_summary.empty and "dataset_id" in dataset_summary
        else pd.DataFrame()
    )
    for dataset_id, spec in config.get("datasets", {}).items():
        summary_row = (
            summary_by_dataset.loc[dataset_id].to_dict()
            if not summary_by_dataset.empty and dataset_id in summary_by_dataset.index
            else {}
        )
        rows.append(_configured_external_evidence_row(dataset_id, spec, summary_row))

    if module_contrasts is not None and not module_contrasts.empty:
        rows.append(
            _external_result_evidence_row(
                config,
                sample_metadata,
                module_contrasts,
                evidence_type="raw_10x_sample_module_scores",
                feature_column="module",
                status_column="status",
                feature_level="sample_module",
                validation_strength="descriptive_sanity_only",
                limitation="Sample-level 10x module scores; n=3 versus n=3 and public cell labels are unavailable.",
                next_action="Resolve cell labels or reference-map cells before donor-level ORA feature replication.",
            )
        )
    if marker_contrasts is not None and not marker_contrasts.empty:
        rows.append(
            _external_result_evidence_row(
                config,
                sample_metadata,
                marker_contrasts,
                evidence_type="raw_10x_marker_only_composition",
                feature_column="marker_panel",
                status_column="status",
                feature_level="marker_panel_fraction",
                validation_strength="marker_only_sanity",
                limitation="Marker-only coarse composition; useful for direction checks, not a substitute for cell labels.",
                next_action="Use reference mapping or manual annotation to convert raw cells into Gateway-compatible states.",
            )
        )
    if mapped_features is not None and not mapped_features.empty:
        rows.append(_mapped_feature_evidence_row(config, sample_metadata, mapped_features, method="marker_reference"))
    if scanvi_features is not None and not scanvi_features.empty:
        rows.append(_mapped_feature_evidence_row(config, sample_metadata, scanvi_features, method="scanvi_scarches"))

    columns = [
        "dataset_id",
        "title",
        "accession",
        "validation_use",
        "evidence_type",
        "feature_level",
        "readiness_class",
        "status",
        "validation_strength",
        "n_samples",
        "n_donors",
        "n_healthy",
        "n_case",
        "case_group",
        "n_features",
        "supports_primary_claim",
        "supports_ndd_claim",
        "limitation",
        "next_action",
        "source_url",
    ]
    return pd.DataFrame(rows, columns=columns)


def _mapped_feature_evidence_row(
    config: dict[str, Any],
    sample_metadata: pd.DataFrame | None,
    mapped_features: pd.DataFrame,
    *,
    method: str,
) -> dict[str, Any]:
    dataset_id = str(mapped_features["dataset_id"].dropna().iloc[0]) if "dataset_id" in mapped_features else "oliva_2022"
    spec = config.get("datasets", {}).get(dataset_id, {})
    disease = mapped_features.get("disease_group", pd.Series(dtype=str)).astype(str)
    n_healthy = int(disease.eq("healthy").sum())
    n_case = int(disease.eq("presbyosmia").sum())
    prefixes = config.get("feature_matrix_contract", {}).get("accepted_feature_prefixes", ["prop__", "clr__", "ratio__"])
    n_features = len(_feature_columns(mapped_features, prefixes))
    is_scanvi = method == "scanvi_scarches"
    return {
        "dataset_id": dataset_id,
        "title": spec.get("title", dataset_id),
        "accession": spec.get("accession", ""),
        "validation_use": spec.get("validation_use", ""),
        "evidence_type": "raw_10x_scanvi_scarches_mapped_features" if is_scanvi else "raw_10x_marker_reference_mapped_features",
        "feature_level": "donor_sample_composition_features",
        "readiness_class": "feature_ready_scanvi_scarches" if is_scanvi else "feature_ready_marker_reference",
        "status": "small_n_mapped_features" if min(n_healthy, n_case) < 5 else "mapped_features",
        "validation_strength": "scanvi_mapped_feature_candidate" if is_scanvi else "mapped_feature_candidate",
        "n_samples": _count_unique(mapped_features, "sample_id"),
        "n_donors": int(mapped_features["donor_id"].nunique()) if "donor_id" in mapped_features else 0,
        "n_healthy": n_healthy,
        "n_case": n_case,
        "case_group": "presbyosmia",
        "n_features": n_features,
        "supports_primary_claim": "candidate_after_replication_test",
        "supports_ndd_claim": "no",
        "limitation": (
            "scANVI/scArches-mapped cell states from 6 biopsy samples; true reference transfer is now available, "
            "but n=3 versus n=3 keeps claims small-n."
            if is_scanvi
            else "Marker-reference mapped cell states from 6 biopsy samples; superseded by scANVI/scArches mapping for label transfer."
        ),
        "next_action": (
            "Compare mapped-feature directions with Gateway aging and expand independent donor-level validation."
            if is_scanvi
            else "Use as a conservative baseline against scANVI/scArches mapping."
        ),
        "source_url": spec.get("source_url", ""),
    }


def _configured_external_evidence_row(
    dataset_id: str,
    spec: dict[str, Any],
    summary_row: dict[str, Any],
) -> dict[str, Any]:
    readiness_class = str(summary_row.get("readiness_class", "not_evaluated"))
    expected_level = str(spec.get("expected_level", ""))
    status = str(spec.get("status", "unknown"))
    ready_feature = _as_bool(summary_row.get("ready_for_feature_validation", False))
    ready_raw = _as_bool(summary_row.get("ready_for_raw_adapter", False))
    evidence_type = "configured_dataset"
    if ready_feature:
        validation_strength = "candidate_direct_validation"
        supports_primary = "candidate_after_replication_test"
        limitation = "Donor-level feature matrix is available but replication statistics still need to be run."
        next_action = "Run feature harmonization and ORA signature replication tests."
    elif ready_raw:
        validation_strength = "candidate_after_cell_mapping"
        supports_primary = "candidate_after_adapter"
        limitation = "Raw expression and metadata are available, but Gateway-compatible donor features are not built yet."
        next_action = "Resolve cell labels or reference-map raw cells, then emit donor composition/module features."
    elif "bulk" in expected_level:
        validation_strength = "marker_context_only"
        supports_primary = "context_only"
        limitation = "Bulk tissue data cannot validate donor-level single-cell ORA features."
        next_action = "Use for olfactory/respiratory marker sanity checks only."
    else:
        validation_strength = "blocked"
        supports_primary = "no"
        limitation = "Required expression, metadata, or donor-feature files are not available locally."
        next_action = "Add source files, access notes, or a donor-feature matrix before validation."

    ndd_contexts = {str(item).lower() for item in spec.get("disease_context", [])}
    supports_ndd = "candidate_after_adapter" if ndd_contexts & {"ad", "pd", "neurodegeneration"} and ready_raw else "no"
    if "bulk" in expected_level:
        supports_ndd = "context_only"
    return {
        "dataset_id": dataset_id,
        "title": spec.get("title", dataset_id),
        "accession": spec.get("accession", ""),
        "validation_use": spec.get("validation_use", ""),
        "evidence_type": evidence_type,
        "feature_level": expected_level,
        "readiness_class": readiness_class,
        "status": status,
        "validation_strength": validation_strength,
        "n_samples": "",
        "n_donors": "",
        "n_healthy": "",
        "n_case": "",
        "case_group": _case_group_from_context(spec.get("disease_context", [])),
        "n_features": "",
        "supports_primary_claim": supports_primary,
        "supports_ndd_claim": supports_ndd,
        "limitation": limitation,
        "next_action": next_action,
        "source_url": spec.get("source_url", ""),
    }


def _external_result_evidence_row(
    config: dict[str, Any],
    sample_metadata: pd.DataFrame | None,
    contrasts: pd.DataFrame,
    *,
    evidence_type: str,
    feature_column: str,
    status_column: str,
    feature_level: str,
    validation_strength: str,
    limitation: str,
    next_action: str,
) -> dict[str, Any]:
    dataset_id = _external_result_dataset_id(sample_metadata, default="oliva_2022")
    spec = config.get("datasets", {}).get(dataset_id, {})
    usable = _usable_external_samples(sample_metadata)
    status_values = _join_unique(contrasts.get(status_column, pd.Series(dtype=str)))
    return {
        "dataset_id": dataset_id,
        "title": spec.get("title", dataset_id),
        "accession": spec.get("accession", ""),
        "validation_use": spec.get("validation_use", ""),
        "evidence_type": evidence_type,
        "feature_level": feature_level,
        "readiness_class": "sanity_check_generated",
        "status": status_values,
        "validation_strength": validation_strength,
        "n_samples": _count_unique(usable, "sample_id"),
        "n_donors": _count_unique(usable, "donor_id"),
        "n_healthy": _max_numeric(contrasts, "n_healthy"),
        "n_case": _max_numeric(contrasts, "n_presbyosmia"),
        "case_group": "presbyosmia",
        "n_features": _count_unique(contrasts, feature_column),
        "supports_primary_claim": "sanity_only",
        "supports_ndd_claim": "no",
        "limitation": limitation,
        "next_action": next_action,
        "source_url": spec.get("source_url", ""),
    }


def _external_result_dataset_id(sample_metadata: pd.DataFrame | None, *, default: str) -> str:
    if sample_metadata is None or sample_metadata.empty or "dataset_id" not in sample_metadata:
        return default
    values = sample_metadata["dataset_id"].dropna().astype(str).unique().tolist()
    return values[0] if len(values) == 1 else default


def _usable_external_samples(sample_metadata: pd.DataFrame | None) -> pd.DataFrame:
    if sample_metadata is None or sample_metadata.empty:
        return pd.DataFrame()
    frame = sample_metadata.copy()
    if "usable_for_external_validation" in frame:
        usable = frame["usable_for_external_validation"].astype(str).str.lower().isin(["true", "1", "yes"])
        frame = frame[usable]
    if "sample_class" in frame:
        frame = frame[frame["sample_class"].astype(str).ne("culture")]
    return frame


def _count_unique(frame: pd.DataFrame, column: str) -> int | str:
    if frame.empty or column not in frame:
        return ""
    return int(frame[column].dropna().astype(str).nunique())


def _max_numeric(frame: pd.DataFrame, column: str) -> int | str:
    if frame.empty or column not in frame:
        return ""
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    return int(values.max()) if not values.empty else ""


def _join_unique(series: pd.Series) -> str:
    values = [str(value) for value in series.dropna().astype(str).unique().tolist() if str(value)]
    return ",".join(sorted(values))


def _as_bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes", "y"}
    return bool(value)


def _case_group_from_context(contexts: Any) -> str:
    values = [str(item).lower() for item in contexts] if isinstance(contexts, list) else [str(contexts).lower()]
    for candidate in ["ad", "pd", "presbyosmia", "aging"]:
        if candidate in values:
            return candidate
    return ""


def _external_direction_to_age_direction(direction: str) -> str:
    if direction == "higher_in_presbyosmia":
        return "positive"
    if direction == "lower_in_presbyosmia":
        return "negative"
    return "flat"


def _direction_concordance(external_direction: str, gateway_direction: str) -> str:
    if external_direction not in {"positive", "negative"} or gateway_direction not in {"positive", "negative"}:
        return "not_evaluable"
    return "concordant" if external_direction == gateway_direction else "discordant"


def _marker_age_concordance_status(concordance: str, contrast_status: Any) -> str:
    suffix = "small_n_marker_only" if str(contrast_status) == "marker_only_small_n" else "marker_only"
    if concordance == "concordant":
        return f"concordant_{suffix}"
    if concordance == "discordant":
        return f"discordant_{suffix}"
    return f"not_evaluable_{suffix}"


def _marker_age_concordance_interpretation(panel: str, gateway_feature: str, concordance: str) -> str:
    if concordance == "concordant":
        return f"{panel} marker shift in GSE184117 has the same direction as Gateway aging for {gateway_feature}."
    if concordance == "discordant":
        return f"{panel} marker shift in GSE184117 has the opposite direction from Gateway aging for {gateway_feature}."
    return f"{panel} marker shift cannot be directionally compared with Gateway aging for {gateway_feature}."


def _mapped_feature_concordance_status(concordance: str, n_healthy: int, n_case: int) -> str:
    suffix = "small_n_mapped" if min(n_healthy, n_case) < 5 else "mapped"
    if concordance == "concordant":
        return f"concordant_{suffix}"
    if concordance == "discordant":
        return f"discordant_{suffix}"
    return f"not_evaluable_{suffix}"


def _mapped_feature_concordance_interpretation(
    external_feature: str,
    gateway_feature: str,
    concordance: str,
) -> str:
    if concordance == "concordant":
        return f"{external_feature} shift in GSE184117 has the same direction as Gateway aging for {gateway_feature}."
    if concordance == "discordant":
        return f"{external_feature} shift in GSE184117 has the opposite direction from Gateway aging for {gateway_feature}."
    return f"{external_feature} shift cannot be directionally compared with Gateway aging for {gateway_feature}."


def _empty_mapped_feature_concordance_row(
    *,
    external_feature: str,
    gateway_feature: str,
    external_direction: str,
    delta: float,
    p_value: float,
    n_healthy: int,
    n_case: int,
    status: str,
) -> dict[str, Any]:
    return {
        "external_feature": external_feature,
        "gateway_feature": gateway_feature,
        "external_direction": external_direction,
        "gateway_age_direction": "",
        "concordance": "not_evaluable",
        "external_delta": delta,
        "external_p_value": p_value,
        "gateway_beta_per_10_years": np.nan,
        "gateway_p_value": np.nan,
        "gateway_fdr": np.nan,
        "n_healthy": int(n_healthy),
        "n_presbyosmia": int(n_case),
        "status": status,
        "interpretation": f"Gateway age association is missing for {gateway_feature}.",
    }


def _empty_marker_age_concordance_row(
    contrast: pd.Series,
    *,
    panel: str,
    gateway_feature: str,
    external_direction: str,
    status: str,
) -> dict[str, Any]:
    return {
        "marker_panel": panel,
        "gateway_feature": gateway_feature,
        "external_direction": external_direction,
        "gateway_age_direction": "",
        "concordance": "not_evaluable",
        "external_delta": pd.to_numeric(contrast.get("presbyosmia_minus_healthy"), errors="coerce"),
        "external_p_value": pd.to_numeric(contrast.get("p_value"), errors="coerce"),
        "gateway_beta_per_10_years": np.nan,
        "gateway_p_value": np.nan,
        "gateway_fdr": np.nan,
        "n_healthy": _safe_int(contrast.get("n_healthy")),
        "n_presbyosmia": _safe_int(contrast.get("n_presbyosmia")),
        "status": status,
        "interpretation": f"No Gateway age-association row was available for mapped feature {gateway_feature}.",
    }


def _safe_int(value: Any) -> int:
    numeric = pd.to_numeric(value, errors="coerce")
    return 0 if pd.isna(numeric) else int(numeric)


def _read_geo_sample_rows(path: Path) -> dict[str, list[list[str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    rows: dict[str, list[list[str]]] = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("!series_matrix_table_begin"):
                break
            if not raw.startswith("!Sample_"):
                continue
            cells = next(csv.reader([raw.rstrip("\n")], delimiter="\t", quotechar='"'))
            if not cells:
                continue
            key = cells[0].lstrip("!")
            rows.setdefault(key, []).append(cells[1:])
    return rows


def _first_geo_row(rows: dict[str, list[list[str]]], key: str) -> list[str]:
    values = rows.get(key, [])
    return values[0] if values else []


def _geo_value(rows: dict[str, list[list[str]]], key: str, idx: int) -> str:
    values = rows.get(key, [])
    if not values or idx >= len(values[0]):
        return ""
    return str(values[0][idx])


def _geo_join(rows: dict[str, list[list[str]]], key: str, idx: int) -> str:
    values = []
    for row in rows.get(key, []):
        if idx < len(row) and str(row[idx]).strip():
            values.append(str(row[idx]).strip())
    return " | ".join(values)


def _geo_characteristics(rows: dict[str, list[list[str]]], idx: int) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for row in rows.get("Sample_characteristics_ch1", []):
        if idx >= len(row):
            continue
        value = str(row[idx]).strip()
        if ":" not in value:
            continue
        key, item = value.split(":", 1)
        parsed[key.strip().lower()] = item.strip()
    return parsed


def _sample_prefix_from_urls(urls: list[str]) -> str:
    for url in urls:
        name = Path(url).name
        for suffix in [
            ".matrix.mtx.gz",
            ".barcodes.tsv.gz",
            ".features.tsv.gz",
            "_matrix.mtx.gz",
            "_barcodes.tsv.gz",
            "_features.tsv.gz",
        ]:
            if name.endswith(suffix):
                return name[: -len(suffix)]
    return ""


def _subject_id_from_prefix(prefix: str) -> str:
    match = re.search(r"_(H\d{4}_\d+)_", prefix)
    return match.group(1) if match else ""


def _accession_from_prefix(prefix: str) -> str:
    match = re.match(r"(GSM\d+)", prefix)
    return match.group(1) if match else ""


def _contains_token(value: str, token: str) -> bool:
    return token.lower() in str(value).lower()


def _external_disease_group(disease_state: str, sample_class: str) -> str:
    if sample_class == "culture":
        return "culture"
    value = disease_state.lower()
    if "normosmia" in value:
        return "healthy"
    if "hyposmia" in value or "anosmia" in value or "presby" in value:
        return "presbyosmia"
    return "unknown"


def _external_sample_class(title: str, description: str, sample_prefix: str) -> str:
    first_description = str(description).split("|", 1)[0].strip().lower()
    if _contains_token(title, "culture") or first_description == "culture" or sample_prefix.lower().endswith("_culture"):
        return "culture"
    return "biopsy"


def _archive_sample_roles(inventory: pd.DataFrame) -> dict[str, dict[str, str]]:
    roles: dict[str, dict[str, str]] = {}
    if inventory.empty:
        return roles
    for _, row in inventory.iterrows():
        sample = str(row.get("sample_guess", ""))
        role = str(row.get("role", ""))
        member = str(row.get("member", ""))
        if not sample or role not in {"matrix", "features", "barcodes"}:
            continue
        roles.setdefault(sample, {})[role] = member
    return roles


def _metadata_for_external_prefix(
    sample_prefix: str,
    metadata_by_accession: pd.DataFrame,
    metadata_by_prefix: pd.DataFrame,
) -> pd.Series:
    if sample_prefix in metadata_by_prefix.index:
        return metadata_by_prefix.loc[sample_prefix]
    accession = _accession_from_prefix(sample_prefix)
    if accession and accession in metadata_by_accession.index:
        return metadata_by_accession.loc[accession]
    return pd.Series({"sample_id": accession, "donor_id": "", "sample_prefix": sample_prefix})


def _read_10x_matrix_from_tar(archive: tarfile.TarFile, member_name: str):
    extracted = archive.extractfile(member_name)
    if extracted is None:
        raise FileNotFoundError(member_name)
    with gzip.GzipFile(fileobj=extracted) as handle:
        return mmread(handle).tocsr()


def _read_10x_features_from_tar(archive: tarfile.TarFile, member_name: str) -> pd.DataFrame:
    extracted = archive.extractfile(member_name)
    if extracted is None:
        raise FileNotFoundError(member_name)
    with gzip.GzipFile(fileobj=extracted) as handle:
        frame = pd.read_csv(handle, sep="\t", header=None, names=["gene_id", "feature_name", "feature_type"])
    return frame


def _read_10x_barcodes_from_tar(archive: tarfile.TarFile, member_name: str) -> list[str]:
    extracted = archive.extractfile(member_name)
    if extracted is None:
        raise FileNotFoundError(member_name)
    with gzip.GzipFile(fileobj=extracted) as handle:
        return [line.decode("utf-8").strip() for line in handle if line.strip()]


def _features_compatible(left: pd.DataFrame, right: pd.DataFrame) -> bool:
    left_ids = left["gene_id"].astype(str).tolist() if "gene_id" in left else []
    right_ids = right["gene_id"].astype(str).tolist() if "gene_id" in right else []
    left_names = left["feature_name"].astype(str).tolist() if "feature_name" in left else []
    right_names = right["feature_name"].astype(str).tolist() if "feature_name" in right else []
    return left_ids == right_ids and left_names == right_names


def _count_gzip_lines_from_tar(archive: tarfile.TarFile, member_name: str) -> int:
    extracted = archive.extractfile(member_name)
    if extracted is None:
        raise FileNotFoundError(member_name)
    with gzip.GzipFile(fileobj=extracted) as handle:
        return sum(1 for _ in handle)


def _mean_log1p_cpm_by_gene(matrix, *, scale: float = 10_000.0) -> np.ndarray:
    cell_totals = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    multipliers = np.zeros_like(cell_totals, dtype=float)
    valid = cell_totals > 0
    multipliers[valid] = scale / cell_totals[valid]
    normalized = matrix.multiply(multipliers)
    normalized.data = np.log1p(normalized.data)
    return np.asarray(normalized.mean(axis=1)).ravel()


def _normalize_marker_panels(marker_panels: dict[str, Any] | None) -> dict[str, tuple[str, ...]]:
    source = marker_panels or DEFAULT_EXTERNAL_MARKER_PANELS
    panels: dict[str, tuple[str, ...]] = {}
    for panel_name, spec in source.items():
        genes = spec.get("genes", []) if isinstance(spec, dict) else spec
        cleaned = tuple(str(gene).strip().upper() for gene in genes if str(gene).strip())
        if cleaned:
            panels[str(panel_name)] = cleaned
    return panels


def _score_marker_panels_by_cell(
    matrix,
    var: pd.DataFrame,
    panels: dict[str, tuple[str, ...]],
    *,
    scale: float,
) -> tuple[dict[str, np.ndarray], dict[str, dict[str, object]]]:
    symbol_to_indices = _feature_symbol_indices(var)
    cell_totals = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    multipliers = np.zeros_like(cell_totals, dtype=float)
    valid = cell_totals > 0
    multipliers[valid] = scale / cell_totals[valid]
    panel_scores: dict[str, np.ndarray] = {}
    coverage: dict[str, dict[str, object]] = {}
    n_cells = int(matrix.shape[1])
    for panel_name, genes in panels.items():
        indices: list[int] = []
        present: list[str] = []
        for gene in genes:
            hits = symbol_to_indices.get(gene.upper(), [])
            if hits:
                indices.extend(hits)
                present.append(gene)
        missing = [gene for gene in genes if gene not in set(present)]
        if indices:
            subset = matrix[indices, :].multiply(multipliers)
            subset.data = np.log1p(subset.data)
            score = np.asarray(subset.mean(axis=0)).ravel()
        else:
            score = np.zeros(n_cells, dtype=float)
        panel_scores[panel_name] = score
        coverage[panel_name] = {
            "n_requested": int(len(genes)),
            "n_present": int(len(set(present))),
            "coverage_fraction": float(len(set(present)) / len(genes)) if genes else np.nan,
            "present_genes": ",".join(sorted(set(present))),
            "missing_genes": ",".join(missing),
            "status": "ok" if present else "no_markers_present",
        }
    return panel_scores, coverage


def _feature_symbol_indices(var: pd.DataFrame) -> dict[str, list[int]]:
    columns = [col for col in ["feature_name", "gene_name", "gene_symbol", "gene_id"] if col in var.columns]
    symbol_to_indices: dict[str, list[int]] = {}
    for idx, row in var.reset_index(drop=True).iterrows():
        for column in columns:
            value = str(row.get(column, "")).strip().upper()
            if value:
                symbol_to_indices.setdefault(value, []).append(int(idx))
    return symbol_to_indices


def _assign_marker_panels(panel_scores: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    panel_names = list(panel_scores)
    if not panel_names:
        empty = np.array([], dtype=object)
        return empty, np.array([], dtype=float), np.array([], dtype=float)
    matrix = np.vstack([panel_scores[name] for name in panel_names])
    best_idx = np.argmax(matrix, axis=0)
    max_scores = matrix[best_idx, np.arange(matrix.shape[1])]
    if matrix.shape[0] > 1:
        second_scores = np.partition(matrix, -2, axis=0)[-2, :]
    else:
        second_scores = np.zeros_like(max_scores)
    margins = max_scores - second_scores
    assignments = np.array([panel_names[idx] for idx in best_idx], dtype=object)
    assignments[max_scores <= 0] = "unassigned_marker_low"
    return assignments, max_scores, margins


def _unique_obs_names(sample_ids: pd.Series, barcodes: pd.Series) -> list[str]:
    return [f"{sample_id}:{barcode}" for sample_id, barcode in zip(sample_ids.astype(str), barcodes.astype(str), strict=False)]


def _safe_ratio(fractions: pd.Series, numerator: str, denominator: str) -> float:
    return float((float(fractions.get(numerator, 0.0)) + 1e-6) / (float(fractions.get(denominator, 0.0)) + 1e-6))


def _mann_whitney_pvalue(case: np.ndarray, control: np.ndarray) -> float:
    if case.size == 0 or control.size == 0:
        return np.nan
    try:
        return float(stats.mannwhitneyu(case, control, alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def _geo_metadata_columns() -> list[str]:
    return [
        "dataset_id",
        "sample_id",
        "donor_id",
        "subject_id",
        "sample_title",
        "sample_description",
        "age",
        "disease_state",
        "disease_group",
        "sample_class",
        "is_presbyosmic",
        "usable_for_external_validation",
        "source_name",
        "organism",
        "instrument_model",
        "library_strategy",
        "chemistry",
        "collection_method",
        "sample_prefix",
        "supplementary_files",
        "status",
    ]


def _external_10x_qc_columns() -> list[str]:
    return [
        "dataset_id",
        "sample_id",
        "donor_id",
        "sample_prefix",
        "age",
        "disease_state",
        "disease_group",
        "sample_class",
        "n_cells",
        "n_barcodes",
        "n_genes",
        "detected_genes",
        "total_counts",
        "median_counts_per_cell",
        "status",
    ]


def _external_10x_score_columns() -> list[str]:
    return [
        "dataset_id",
        "sample_id",
        "donor_id",
        "sample_prefix",
        "age",
        "disease_state",
        "disease_group",
        "sample_class",
        "module",
        "description",
        "n_requested",
        "n_present",
        "coverage_fraction",
        "missing_genes",
        "mean_log1p_cpm",
        "status",
    ]


def _external_10x_contrast_columns() -> list[str]:
    return [
        "module",
        "n_healthy",
        "n_presbyosmia",
        "mean_healthy",
        "mean_presbyosmia",
        "presbyosmia_minus_healthy",
        "p_value",
        "direction",
        "status",
    ]


def _external_marker_coverage_columns() -> list[str]:
    return [
        "dataset_id",
        "sample_id",
        "sample_prefix",
        "marker_panel",
        "n_requested",
        "n_present",
        "coverage_fraction",
        "present_genes",
        "missing_genes",
        "status",
    ]


def _external_marker_composition_columns() -> list[str]:
    return [
        "dataset_id",
        "sample_id",
        "donor_id",
        "sample_prefix",
        "age",
        "disease_state",
        "disease_group",
        "sample_class",
        "marker_panel",
        "n_cells_assigned",
        "total_cells",
        "fraction_cells",
        "mean_marker_score",
        "median_marker_score",
        "mean_assignment_margin",
        "mean_top_marker_score",
        "n_markers_present",
        "status",
    ]


def _external_marker_contrast_columns() -> list[str]:
    return [
        "marker_panel",
        "n_healthy",
        "n_presbyosmia",
        "mean_fraction_healthy",
        "mean_fraction_presbyosmia",
        "presbyosmia_minus_healthy",
        "p_value",
        "direction",
        "status",
    ]


def _external_mapping_qc_columns() -> list[str]:
    return [
        "dataset_id",
        "sample_id",
        "donor_id",
        "sample_prefix",
        "age",
        "disease_group",
        "sample_class",
        "n_cells",
        "n_genes",
        "mapped_fraction",
        "mean_mapping_confidence",
        "mean_mapping_margin",
        "n_marker_genes_present",
        "status",
    ]


def _external_marker_age_concordance_columns() -> list[str]:
    return [
        "marker_panel",
        "gateway_feature",
        "external_direction",
        "gateway_age_direction",
        "concordance",
        "external_delta",
        "external_p_value",
        "gateway_beta_per_10_years",
        "gateway_p_value",
        "gateway_fdr",
        "n_healthy",
        "n_presbyosmia",
        "status",
        "interpretation",
    ]


def _external_mapped_feature_concordance_columns() -> list[str]:
    return [
        "external_feature",
        "gateway_feature",
        "external_direction",
        "gateway_age_direction",
        "concordance",
        "external_delta",
        "external_p_value",
        "gateway_beta_per_10_years",
        "gateway_p_value",
        "gateway_fdr",
        "n_healthy",
        "n_presbyosmia",
        "status",
        "interpretation",
    ]


def _readiness_class(spec: dict[str, Any], provided: list[str], missing: list[str]) -> str:
    expected = str(spec.get("expected_level", "")).lower()
    status = str(spec.get("status", "")).lower()
    provided_set = set(provided)
    if "feature_matrix" in provided_set:
        return "ready_feature_matrix"
    if {"expression", "metadata"}.issubset(provided_set):
        return "ready_raw_adapter"
    if "expression" in provided_set and "metadata" in missing:
        return "blocked_metadata"
    if "bulk" in expected:
        return "marker_only"
    if "download_ready" in status:
        return "blocked_file_pending"
    if "metadata_pending" in status:
        return "blocked_metadata"
    if missing:
        return "blocked_missing_files"
    return "unknown"


def _archive_member_role(name: str) -> str:
    lower = Path(name).name.lower()
    if "matrix.mtx" in lower:
        return "matrix"
    if "barcodes.tsv" in lower:
        return "barcodes"
    if "features.tsv" in lower or "genes.tsv" in lower:
        return "features"
    if lower.endswith((".csv", ".tsv", ".txt")) and any(token in lower for token in ["meta", "annot", "cell", "sample"]):
        return "metadata"
    if lower.endswith((".h5ad", ".h5")):
        return "matrix_h5"
    return "other"


def _sample_guess(name: str, role: str) -> str:
    if role in {"other", "empty_archive", "missing_archive"}:
        return ""
    stem = Path(name).name
    for suffix in [
        "_matrix.mtx.gz",
        "_matrix.mtx",
        "_barcodes.tsv.gz",
        "_barcodes.tsv",
        "_features.tsv.gz",
        "_features.tsv",
        "_genes.tsv.gz",
        "_genes.tsv",
        ".matrix.mtx.gz",
        ".matrix.mtx",
        ".barcodes.tsv.gz",
        ".barcodes.tsv",
        ".features.tsv.gz",
        ".features.tsv",
        ".genes.tsv.gz",
        ".genes.tsv",
    ]:
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return Path(name).parent.name if Path(name).parent.name != "." else Path(name).stem


def _read_table(path: Path) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    sep = "," if ".csv" in suffixes else "\t"
    return pd.read_csv(path, sep=sep)


def _feature_columns(frame: pd.DataFrame | None, prefixes: list[str]) -> list[str]:
    if frame is None:
        return []
    return [str(col) for col in frame.columns if any(str(col).startswith(prefix) for prefix in prefixes)]


def _feature_harmonization(feature_cols: list[str], gateway_cols: set[str], *, dataset_id: str) -> pd.DataFrame:
    if not gateway_cols:
        rows = [
            {
                "dataset_id": dataset_id,
                "feature": feature,
                "status": "external_only",
                "feature_kind": _feature_kind(feature),
            }
            for feature in sorted(feature_cols)
        ]
        return pd.DataFrame(rows, columns=_harmonization_columns())
    external = set(feature_cols)
    rows = []
    for feature in sorted(gateway_cols | external):
        if feature in gateway_cols and feature in external:
            status = "matched"
        elif feature in gateway_cols:
            status = "missing_from_external"
        else:
            status = "external_only"
        rows.append(
            {
                "dataset_id": dataset_id,
                "feature": feature,
                "status": status,
                "feature_kind": _feature_kind(feature),
            }
        )
    return pd.DataFrame(rows, columns=_harmonization_columns())


def _feature_kind(feature: str) -> str:
    if "__" in feature:
        return feature.split("__", 1)[0]
    return "metadata_or_unknown"


def _harmonization_columns() -> list[str]:
    return ["dataset_id", "feature", "status", "feature_kind"]


def _resolve_optional_path(path: object, base_dir: Path) -> Path | None:
    if path in {None, ""}:
        return None
    candidate = Path(str(path))
    if candidate.is_absolute():
        return candidate
    return base_dir / candidate
