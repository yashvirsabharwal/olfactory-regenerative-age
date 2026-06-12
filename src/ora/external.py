"""External validation registry and gene-list coverage helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .modules import GeneSet, resolve_gene_sets


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
                "notes": spec.get("notes", ""),
            }
        )
    return pd.DataFrame(rows)


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


def _resolve_optional_path(path: object, base_dir: Path) -> Path | None:
    if path in {None, ""}:
        return None
    candidate = Path(str(path))
    if candidate.is_absolute():
        return candidate
    return base_dir / candidate
