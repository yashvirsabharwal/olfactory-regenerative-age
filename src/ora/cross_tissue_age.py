"""Cross-tissue CELLxGENE donor-level age-effect estimates."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .modules import parse_gene_sets, resolve_gene_sets
from .utils import ensure_parent, normalize_token


DEFAULT_MIN_DONORS = 4
DEFAULT_MIN_CELLS_PER_DONOR = 20
DEFAULT_ADULT_MIN_AGE = 18.0
PRIMARY_SCOPE = "adult_only"
CONTEXT_SCOPE = "all_stages_context"


@dataclass(frozen=True)
class CrossTissueDatasetSpec:
    dataset_id: str
    collection_id: str
    title: str
    tissue_class: str
    h5ad_path: str
    source_url: str = ""


def parse_cellxgene_age(value: object) -> float:
    """Parse CELLxGENE development-stage labels into approximate years."""

    text = str(value or "").strip().lower()
    if not text or text in {"nan", "none", "unknown"}:
        return np.nan
    year_match = re.search(r"(\d+(?:\.\d+)?)\s*-\s*year\s*-\s*old", text)
    if year_match:
        return float(year_match.group(1))
    direct_year_match = re.search(r"(\d+(?:\.\d+)?)\s*year", text)
    if direct_year_match:
        return float(direct_year_match.group(1))
    week_match = re.search(r"(\d+(?:\.\d+)?)\s*(?:st|nd|rd|th)?\s*week", text)
    if week_match:
        return float(week_match.group(1)) / 52.0
    month_match = re.search(r"(\d+(?:\.\d+)?)\s*month", text)
    if month_match:
        return float(month_match.group(1)) / 12.0
    numeric_match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*", text)
    if numeric_match:
        return float(numeric_match.group(1))
    return np.nan


def classify_age_context(value: object, age_years: float, adult_min_age: float = DEFAULT_ADULT_MIN_AGE) -> str:
    """Classify whether an age label is usable for adult-aging inference."""

    text = str(value or "").lower()
    if "post-fertilization" in text or "embry" in text or "fetal" in text:
        return "prenatal"
    if not np.isfinite(age_years):
        return "unknown"
    if age_years < adult_min_age:
        return "child_or_developmental"
    return "adult"


def harmonize_cell_group(label: object) -> str:
    """Map public CELLxGENE cell labels to coarse ORA-relevant comparator groups."""

    text = normalize_token(label)
    if not text:
        return "unknown"
    if "goblet" in text:
        return "goblet"
    if "club" in text:
        return "club"
    if "ionocyte" in text:
        return "ionocyte"
    if "multiciliated" in text or "ciliated" in text or "deuterosomal" in text:
        return "ciliated"
    if "suprabasal" in text:
        return "suprabasal"
    if "basal" in text:
        return "basal"
    if "secretory" in text or "glandular" in text or "mucous" in text or "serous" in text:
        return "secretory"
    if "respiratory tract epithelial" in text or "respiratory epithelial" in text:
        return "respiratory_epithelial"
    if "epithelial" in text or "pneumocyte" in text or "alveolar type" in text:
        return "other_epithelial"
    if "macrophage" in text or "monocyte" in text or "myeloid" in text or "neutrophil" in text:
        return "myeloid"
    if "antigen presenting" in text or "dendritic" in text:
        return "antigen_presenting"
    if "t cell" in text:
        return "t_cell"
    if "natural killer" in text or "nk cell" in text:
        return "nk_cell"
    if "b cell" in text or "plasma" in text:
        return "b_cell"
    if "mast" in text:
        return "mast_cell"
    if "immune" in text:
        return "immune_other"
    if "fibroblast" in text or "mesenchymal" in text:
        return "mesenchymal"
    if "endothelial" in text:
        return "endothelial"
    return "other"


def parse_dataset_specs(config: dict[str, Any]) -> list[CrossTissueDatasetSpec]:
    """Parse selected comparator datasets from a YAML config dictionary."""

    specs = []
    for row in config.get("datasets", []):
        if not isinstance(row, dict):
            continue
        specs.append(
            CrossTissueDatasetSpec(
                dataset_id=str(row.get("dataset_id", "")).strip(),
                collection_id=str(row.get("collection_id", "")).strip(),
                title=str(row.get("title", "")).strip(),
                tissue_class=str(row.get("tissue_class", "")).strip(),
                h5ad_path=str(row.get("h5ad_path", "")).strip(),
                source_url=str(row.get("source_url", "")).strip(),
            )
        )
    return [spec for spec in specs if spec.dataset_id and spec.h5ad_path]


def build_cellxgene_asset_inventory(dataset_specs: list[CrossTissueDatasetSpec]) -> pd.DataFrame:
    """Inspect local comparator H5AD files and summarize metadata readiness."""

    import anndata as ad

    rows: list[dict[str, object]] = []
    for spec in dataset_specs:
        path = Path(spec.h5ad_path)
        if not path.exists():
            rows.append(_missing_inventory_row(spec, "missing_h5ad"))
            continue
        adata = ad.read_h5ad(path, backed="r")
        try:
            obs = adata.obs
            age = pd.to_numeric(
                obs.get("development_stage", pd.Series(index=obs.index, dtype=object)).map(parse_cellxgene_age),
                errors="coerce",
            )
            age_context = [
                classify_age_context(stage, age_year)
                for stage, age_year in zip(
                    obs.get("development_stage", pd.Series(index=obs.index, dtype=object)),
                    age,
                    strict=False,
                )
            ]
            donor_col = _first_present(obs, ["donor_id", "sample_id", "subject_id"])
            label_col = _first_present(obs, ["cell_type", "author_cell_type", "author_cell_type_final"])
            donors = obs[donor_col].astype(str) if donor_col else pd.Series(dtype=object)
            rows.append(
                {
                    "dataset_id": spec.dataset_id,
                    "collection_id": spec.collection_id,
                    "title": spec.title,
                    "tissue_class": spec.tissue_class,
                    "h5ad_path": str(path),
                    "source_url": spec.source_url,
                    "status": "ok",
                    "n_cells": int(adata.n_obs),
                    "n_genes": int(adata.n_vars),
                    "n_donors": int(donors.nunique()) if donor_col else 0,
                    "n_adult_donors": _n_donors_with_context(obs, donor_col, age_context, "adult"),
                    "age_min_years": float(np.nanmin(age)) if np.isfinite(age).any() else np.nan,
                    "age_max_years": float(np.nanmax(age)) if np.isfinite(age).any() else np.nan,
                    "age_contexts": ";".join(sorted(set(str(item) for item in age_context if str(item)))),
                    "donor_column": donor_col,
                    "age_column": "development_stage" if "development_stage" in obs.columns else "",
                    "cell_label_column": label_col,
                    "expression_scale": "cellxgene_X_as_provided",
                    "primary_age_scope": PRIMARY_SCOPE,
                    "notes": _inventory_notes(age_context, obs, donor_col),
                }
            )
        finally:
            close = getattr(adata, "file", None)
            if close is not None:
                close.close()
    return pd.DataFrame(rows)


def build_cellxgene_donor_feature_matrix(
    dataset_specs: list[CrossTissueDatasetSpec],
    gene_set_config: dict[str, Any],
    *,
    chunk_size: int = 20_000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build donor-level composition and module features from selected public H5ADs."""

    import anndata as ad

    feature_tables: list[pd.DataFrame] = []
    coverage_tables: list[pd.DataFrame] = []
    for spec in dataset_specs:
        path = Path(spec.h5ad_path)
        if not path.exists():
            continue
        adata = ad.read_h5ad(path, backed="r")
        try:
            metadata = _cellxgene_cell_metadata(adata.obs, spec)
            composition = _donor_composition_features(metadata)
            modules, coverage = _donor_module_features(
                adata=adata,
                metadata=metadata,
                gene_set_config=gene_set_config,
                dataset_id=spec.dataset_id,
                chunk_size=chunk_size,
            )
            donor_features = composition.merge(modules, on="donor_id", how="left")
            donor_features.insert(0, "dataset_id", spec.dataset_id)
            donor_features.insert(1, "dataset_title", spec.title)
            donor_features.insert(2, "tissue_class", spec.tissue_class)
            feature_tables.append(donor_features)
            coverage.insert(0, "dataset_id", spec.dataset_id)
            coverage.insert(1, "dataset_title", spec.title)
            coverage.insert(2, "tissue_class", spec.tissue_class)
            coverage_tables.append(coverage)
        finally:
            close = getattr(adata, "file", None)
            if close is not None:
                close.close()
    donor_matrix = pd.concat(feature_tables, ignore_index=True) if feature_tables else pd.DataFrame()
    if not donor_matrix.empty:
        composition_cols = [col for col in donor_matrix.columns if col.startswith("cell_group__")]
        donor_matrix[composition_cols] = donor_matrix[composition_cols].fillna(0.0)
    coverage_matrix = pd.concat(coverage_tables, ignore_index=True) if coverage_tables else pd.DataFrame()
    return donor_matrix, coverage_matrix


def estimate_cross_tissue_age_effects(
    donor_features: pd.DataFrame,
    *,
    min_donors: int = DEFAULT_MIN_DONORS,
    min_cells_per_donor: int = DEFAULT_MIN_CELLS_PER_DONOR,
    adult_min_age: float = DEFAULT_ADULT_MIN_AGE,
) -> pd.DataFrame:
    """Estimate per-dataset donor-level age effects for cross-tissue features."""

    if donor_features.empty:
        return pd.DataFrame()
    feature_cols = [
        col
        for col in donor_features.columns
        if col.startswith("cell_group__") or col.startswith("module_score__")
    ]
    rows: list[dict[str, object]] = []
    for dataset_id, dataset_frame in donor_features.groupby("dataset_id", sort=False):
        frame = dataset_frame.copy()
        frame["age_years"] = pd.to_numeric(frame["age_years"], errors="coerce")
        frame["n_cells"] = pd.to_numeric(frame["n_cells"], errors="coerce").fillna(0)
        scopes = [
            (PRIMARY_SCOPE, frame[frame["age_context"].eq("adult")].copy()),
            (CONTEXT_SCOPE, frame.copy()),
        ]
        for scope, scoped in scopes:
            scoped = scoped[(scoped["n_cells"] >= min_cells_per_donor) & scoped["age_years"].notna()].copy()
            scope_note = _scope_note(scope, scoped, adult_min_age)
            for feature in feature_cols:
                values = pd.to_numeric(scoped.get(feature, pd.Series(dtype=float)), errors="coerce")
                mask = values.notna() & scoped["age_years"].notna()
                n = int(mask.sum())
                base = _effect_base_row(dataset_id, scoped, feature, scope, n, scope_note)
                if n < min_donors:
                    rows.append({**base, **_empty_effect("insufficient_donors")})
                    continue
                y = values[mask].to_numpy(dtype=float)
                if float(np.nanstd(y)) <= 0:
                    rows.append({**base, **_empty_effect("constant_feature")})
                    continue
                x = (scoped.loc[mask, "age_years"].to_numpy(dtype=float) - float(scoped.loc[mask, "age_years"].mean())) / 10.0
                result = stats.linregress(x, y)
                direction = "positive" if result.slope > 0 else "negative" if result.slope < 0 else "flat"
                rows.append(
                    {
                        **base,
                        "beta_per_10_years": float(result.slope),
                        "p_value": float(result.pvalue),
                        "r_value": float(result.rvalue),
                        "direction": direction,
                        "status": "ok" if scope == PRIMARY_SCOPE else "context_only",
                    }
                )
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    table["fdr_within_dataset_scope"] = np.nan
    for _, idx in table.groupby(["dataset_id", "analysis_scope"]).groups.items():
        table.loc[list(idx), "fdr_within_dataset_scope"] = _benjamini_hochberg(table.loc[list(idx), "p_value"])
    return table.sort_values(["dataset_id", "analysis_scope", "feature"]).reset_index(drop=True)


def summarize_ora_cross_tissue_age_effects(
    classification: pd.DataFrame,
    effects: pd.DataFrame,
) -> pd.DataFrame:
    """Join measured comparator age effects back to ORA features."""

    if classification.empty:
        return pd.DataFrame()
    ok_effects = (
        effects[
            effects["analysis_scope"].eq(PRIMARY_SCOPE)
            & effects["status"].eq("ok")
            & effects["feature"].notna()
        ].copy()
        if not effects.empty
        else pd.DataFrame()
    )
    rows: list[dict[str, object]] = []
    for _, row in classification.iterrows():
        feature = str(row.get("feature", ""))
        mapped = map_ora_feature_to_cross_tissue_feature(feature)
        matched = ok_effects[ok_effects["feature"].isin(mapped)] if mapped and not ok_effects.empty else pd.DataFrame()
        gateway_direction = str(row.get("gateway_age_direction", "not_tested"))
        directions = sorted(set(matched["direction"].astype(str))) if not matched.empty else []
        rows.append(
            {
                "feature": feature,
                "specificity_class": row.get("specificity_class", ""),
                "gateway_age_direction": gateway_direction,
                "gateway_age_beta_per_10_years": row.get("gateway_age_beta_per_10_years", np.nan),
                "mapped_cross_tissue_features": ";".join(mapped),
                "n_adult_comparator_effects": int(matched.shape[0]),
                "comparator_dataset_ids": ";".join(sorted(set(matched["dataset_id"].astype(str)))) if not matched.empty else "",
                "comparator_directions": ";".join(directions),
                "n_positive_comparators": int(matched["direction"].eq("positive").sum()) if not matched.empty else 0,
                "n_negative_comparators": int(matched["direction"].eq("negative").sum()) if not matched.empty else 0,
                "n_gateway_direction_concordant": _n_concordant(matched, gateway_direction),
                "min_comparator_p_value": float(matched["p_value"].min()) if not matched.empty else np.nan,
                "max_abs_comparator_beta": float(matched["beta_per_10_years"].abs().max()) if not matched.empty else np.nan,
                "external_age_effect_status": _summary_status(feature, mapped, matched, row),
                "top_comparator_effect": _top_effect_label(matched),
            }
        )
    return pd.DataFrame(rows)


def render_cross_tissue_age_effect_report(
    *,
    inventory: pd.DataFrame,
    effects: pd.DataFrame,
    ora_summary: pd.DataFrame,
) -> str:
    """Render a short research note for the cross-tissue age-effect stage."""

    lines = [
        "# Cross-Tissue Age-Effect Estimates",
        "",
        "Updated: 2026-06-25",
        "",
        "## Scope",
        "",
        "This analysis estimates donor-level age effects in selected public CELLxGENE nasal, bronchial, and lung comparator datasets. The primary claim gate uses adult donors only; fetal or pediatric lung stages are retained only as context.",
        "",
        "## Comparator Assets",
        "",
        "| Dataset | Tissue | Cells | Donors | Adult donors | Age range | Status | Notes |",
        "| --- | --- | ---: | ---: | ---: | --- | --- | --- |",
    ]
    for _, row in inventory.iterrows():
        age_range = _format_age_range(row.get("age_min_years", np.nan), row.get("age_max_years", np.nan))
        lines.append(
            "| "
            f"{row.get('title', row.get('dataset_id', ''))} | {row.get('tissue_class', '')} | "
            f"{int(row.get('n_cells', 0) or 0)} | {int(row.get('n_donors', 0) or 0)} | "
            f"{int(row.get('n_adult_donors', 0) or 0)} | {age_range} | {row.get('status', '')} | "
            f"{row.get('notes', '')} |"
        )
    lines.extend(["", "## Measured Adult Effects", ""])
    if effects.empty:
        lines.append("No effect estimates were generated.")
    else:
        adult_ok = effects[effects["analysis_scope"].eq(PRIMARY_SCOPE) & effects["status"].eq("ok")]
        n_nominal = int((adult_ok["p_value"] < 0.05).sum()) if "p_value" in adult_ok else 0
        n_fdr = (
            int((adult_ok["fdr_within_dataset_scope"] < 0.05).sum())
            if "fdr_within_dataset_scope" in adult_ok
            else 0
        )
        by_dataset = (
            adult_ok.groupby(["dataset_id", "tissue_class"])["feature"]
            .count()
            .reset_index(name="n_tested_features")
            .sort_values(["dataset_id", "tissue_class"])
        )
        lines.extend(["| Dataset | Tissue | Tested adult effects |", "| --- | --- | ---: |"])
        for _, row in by_dataset.iterrows():
            lines.append(f"| {row['dataset_id']} | {row['tissue_class']} | {int(row['n_tested_features'])} |")
        if by_dataset.empty:
            lines.append("| none | | 0 |")
        lines.extend(
            [
                "",
                f"Nominal adult comparator effects with p<0.05: {n_nominal}.",
                f"Adult comparator effects with within-dataset/scope FDR<0.05: {n_fdr}.",
            ]
        )
    lines.extend(["", "## ORA Feature-Level Status", ""])
    if ora_summary.empty:
        lines.append("No ORA feature summary was generated.")
    else:
        status_counts = (
            ora_summary["external_age_effect_status"].value_counts().reset_index(name="n_features")
        )
        status_counts.columns = ["status", "n_features"]
        lines.extend(["| Status | ORA features |", "| --- | ---: |"])
        for _, row in status_counts.iterrows():
            lines.append(f"| {row['status']} | {int(row['n_features'])} |")
    lines.extend(
        [
            "",
            "## Interpretation Guardrails",
            "",
            "- Adult nasal and bronchial estimates are donor-level comparator evidence for shared airway epithelial and immune programs.",
            "- The LungMAP asset contains prenatal/child/adult stages, so its all-stage estimates are context-only unless an adult-only subset reaches the donor threshold.",
            "- These estimates test cross-tissue sharing of feature families; they do not turn ORA into an absolute age clock or prove causal regeneration mechanisms.",
            "",
        ]
    )
    return "\n".join(lines)


def write_cross_tissue_age_outputs(
    *,
    inventory: pd.DataFrame,
    donor_features: pd.DataFrame,
    module_coverage: pd.DataFrame,
    effects: pd.DataFrame,
    ora_summary: pd.DataFrame,
    report_md: str,
    inventory_out: str | Path,
    donor_features_out: str | Path,
    module_coverage_out: str | Path,
    effects_out: str | Path,
    ora_summary_out: str | Path,
    report_out: str | Path,
) -> None:
    """Write all cross-tissue age-effect outputs."""

    inventory.to_csv(ensure_parent(inventory_out), sep="\t", index=False)
    donor_features.to_csv(ensure_parent(donor_features_out), sep="\t", index=False)
    module_coverage.to_csv(ensure_parent(module_coverage_out), sep="\t", index=False)
    effects.to_csv(ensure_parent(effects_out), sep="\t", index=False)
    ora_summary.to_csv(ensure_parent(ora_summary_out), sep="\t", index=False)
    ensure_parent(report_out).write_text(report_md, encoding="utf-8")


def map_ora_feature_to_cross_tissue_feature(feature: str) -> list[str]:
    """Map an ORA feature name to one or more measured cross-tissue feature names."""

    if feature.startswith("module_score__"):
        label = feature.removeprefix("module_score__")
        return [f"module_score__{label}"]
    if "__" not in feature:
        return []
    kind, label = feature.split("__", 1)
    if kind not in {"prop", "clr", "ratio"}:
        return []
    mapped = _ORA_CELL_LABEL_MAP.get(label.lower())
    if mapped is None:
        mapped = harmonize_cell_group(label)
        if mapped in {"other", "unknown"}:
            return []
    if isinstance(mapped, str):
        return [f"cell_group__{mapped}"]
    return [f"cell_group__{item}" for item in mapped]


_ORA_CELL_LABEL_MAP: dict[str, str | list[str]] = {
    "goblet": "goblet",
    "club": "club",
    "ionocyte": "ionocyte",
    "mv_ionocyte": "ionocyte",
    "multiciliated": "ciliated",
    "deuterosomal": "ciliated",
    "mucous_gland": "secretory",
    "bowman_gland": "secretory",
    "serous": "secretory",
    "proliferating_secretory": "secretory",
    "suprabasal": "suprabasal",
    "quiescent_hbc": "basal",
    "cycling_hbc": "basal",
    "activated_hbc": "basal",
    "classical": "myeloid",
    "nonclassical": "myeloid",
    "macrophage": "myeloid",
    "proliferating_mac": "myeloid",
    "cdc1": "antigen_presenting",
    "cdc2": "antigen_presenting",
    "maturedc": "antigen_presenting",
    "pdc": "antigen_presenting",
    "antigen_presenting": "antigen_presenting",
    "cd56bright": "nk_cell",
    "proliferating_nk": "nk_cell",
    "naiveb": "b_cell",
    "gcb": "b_cell",
    "plasmacell": "b_cell",
    "naive_cd4": "t_cell",
    "naive_cd8": "t_cell",
    "cm_cd4": "t_cell",
    "em_cd8": "t_cell",
    "cytotoxic_cd4": "t_cell",
    "cytotoxic_cd8": "t_cell",
    "proliferating_tcell": "t_cell",
    "inflammatory": ["myeloid", "antigen_presenting"],
    "complement_high": "myeloid",
    "lipid_associated": "myeloid",
}


def _cellxgene_cell_metadata(obs: pd.DataFrame, spec: CrossTissueDatasetSpec) -> pd.DataFrame:
    donor_col = _first_present(obs, ["donor_id", "sample_id", "subject_id"])
    if not donor_col:
        raise KeyError(f"No donor/sample column found for {spec.dataset_id}.")
    label_col = _first_present(obs, ["cell_type", "author_cell_type", "author_cell_type_final"])
    if not label_col:
        raise KeyError(f"No cell label column found for {spec.dataset_id}.")
    stage_col = "development_stage" if "development_stage" in obs.columns else ""
    ages = (
        pd.to_numeric(obs[stage_col].map(parse_cellxgene_age), errors="coerce")
        if stage_col
        else pd.Series(np.nan, index=obs.index)
    )
    metadata = pd.DataFrame(
        {
            "donor_id": obs[donor_col].astype(str),
            "sample_id": obs["sample_id"].astype(str) if "sample_id" in obs.columns else obs[donor_col].astype(str),
            "age_years": ages.astype(float),
            "development_stage": obs[stage_col].astype(str) if stage_col else "",
            "age_context": [
                classify_age_context(stage, age)
                for stage, age in zip(obs[stage_col] if stage_col else [""] * obs.shape[0], ages, strict=False)
            ],
            "sex": obs["sex"].astype(str) if "sex" in obs.columns else "unknown",
            "smoking_status": obs["smoking_status"].astype(str) if "smoking_status" in obs.columns else "unknown",
            "tissue": obs["tissue"].astype(str) if "tissue" in obs.columns else spec.tissue_class,
            "cell_label": obs[label_col].astype(str),
            "cell_group": obs[label_col].map(harmonize_cell_group),
        },
        index=obs.index,
    )
    return metadata


def _donor_composition_features(metadata: pd.DataFrame) -> pd.DataFrame:
    donor_meta = (
        metadata.groupby("donor_id", sort=True)
        .agg(
            n_cells=("donor_id", "size"),
            age_years=("age_years", "median"),
            development_stage=("development_stage", _mode_join),
            age_context=("age_context", _mode_join),
            sex=("sex", _mode_join),
            smoking_status=("smoking_status", _mode_join),
            tissue=("tissue", _mode_join),
        )
        .reset_index()
    )
    counts = (
        metadata.groupby(["donor_id", "cell_group"], sort=True, observed=False)
        .size()
        .reset_index(name="n")
    )
    wide = counts.pivot(index="donor_id", columns="cell_group", values="n").fillna(0)
    totals = wide.sum(axis=1).replace(0, np.nan)
    proportions = wide.div(totals, axis=0)
    proportions.columns = [f"cell_group__{col}" for col in proportions.columns]
    proportions = proportions.reset_index()
    return donor_meta.merge(proportions, on="donor_id", how="left").fillna(0)


def _donor_module_features(
    *,
    adata: Any,
    metadata: pd.DataFrame,
    gene_set_config: dict[str, Any],
    dataset_id: str,
    chunk_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_sets = parse_gene_sets(gene_set_config)
    symbol_columns = gene_set_config.get("score", {}).get(
        "var_symbol_columns",
        ["feature_name", "gene_symbol", "gene_name", "symbol"],
    )
    resolved, coverage = resolve_gene_sets(adata.var, adata.var_names, gene_sets, symbol_columns)
    module_indices = {module: idxs for module, idxs in resolved.items() if idxs}
    coverage["status"] = np.where(coverage["n_present"] > 0, "ok", "no_genes_present")
    if not module_indices:
        return pd.DataFrame({"donor_id": sorted(metadata["donor_id"].unique())}), coverage
    union_indices = sorted({idx for idxs in module_indices.values() for idx in idxs})
    position = {idx: pos for pos, idx in enumerate(union_indices)}
    module_positions = {
        module: [position[idx] for idx in idxs]
        for module, idxs in module_indices.items()
    }
    partials = []
    matrix = adata.X
    n_obs = int(adata.n_obs)
    chunk_size = max(1, int(chunk_size))
    for start in range(0, n_obs, chunk_size):
        stop = min(start + chunk_size, n_obs)
        x = _as_dense(matrix[start:stop, union_indices])
        if x.size == 0:
            continue
        scores = {
            f"module_score__{module}": x[:, positions].mean(axis=1)
            for module, positions in module_positions.items()
            if positions
        }
        score_frame = pd.DataFrame(scores)
        score_frame["donor_id"] = metadata.iloc[start:stop]["donor_id"].to_numpy()
        partials.append(score_frame)
    if not partials:
        return pd.DataFrame({"donor_id": sorted(metadata["donor_id"].unique())}), coverage
    scored = pd.concat(partials, ignore_index=True)
    module_cols = [col for col in scored.columns if col.startswith("module_score__")]
    donor_features = scored.groupby("donor_id", sort=True)[module_cols].mean().reset_index()
    coverage["dataset_id_check"] = dataset_id
    return donor_features, coverage


def _effect_base_row(
    dataset_id: str,
    frame: pd.DataFrame,
    feature: str,
    scope: str,
    n: int,
    scope_note: str,
) -> dict[str, object]:
    return {
        "dataset_id": dataset_id,
        "dataset_title": _single_value(frame.get("dataset_title", pd.Series(dtype=object))),
        "tissue_class": _single_value(frame.get("tissue_class", pd.Series(dtype=object))),
        "analysis_scope": scope,
        "feature": feature,
        "feature_family": "module_score" if feature.startswith("module_score__") else "cell_composition",
        "n_donors": n,
        "age_min_years": float(frame["age_years"].min()) if n else np.nan,
        "age_max_years": float(frame["age_years"].max()) if n else np.nan,
        "mean_cells_per_donor": float(frame["n_cells"].mean()) if n else np.nan,
        "scope_note": scope_note,
    }


def _empty_effect(status: str) -> dict[str, object]:
    return {
        "beta_per_10_years": np.nan,
        "p_value": np.nan,
        "r_value": np.nan,
        "direction": "not_tested",
        "status": status,
    }


def _scope_note(scope: str, scoped: pd.DataFrame, adult_min_age: float) -> str:
    if scope == PRIMARY_SCOPE:
        return f"adult donors with age >= {adult_min_age:g} years"
    contexts = ";".join(sorted(set(scoped.get("age_context", pd.Series(dtype=object)).astype(str))))
    return f"context-only all available stages: {contexts}"


def _benjamini_hochberg(p_values: pd.Series) -> np.ndarray:
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    q = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return q
    valid_values = values[valid]
    order = np.argsort(valid_values)
    ranked = valid_values[order]
    n = ranked.size
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    valid_indices = np.where(valid)[0]
    q[valid_indices[order]] = adjusted
    return q


def _n_concordant(matched: pd.DataFrame, gateway_direction: str) -> int:
    if matched.empty or gateway_direction not in {"positive", "negative"}:
        return 0
    return int(matched["direction"].eq(gateway_direction).sum())


def _summary_status(feature: str, mapped: list[str], matched: pd.DataFrame, row: pd.Series) -> str:
    specificity = str(row.get("specificity_class", ""))
    if not mapped:
        if specificity == "olfactory_specific":
            return "expected_no_cross_tissue_cell_state_mapping"
        return "no_cross_tissue_feature_mapping"
    if matched.empty:
        return "mapped_but_no_adult_comparator_effect"
    if specificity == "olfactory_specific":
        return "measured_airway_context_for_olfactory_feature"
    return "measured_adult_comparator_age_effect"


def _top_effect_label(matched: pd.DataFrame) -> str:
    if matched.empty:
        return ""
    ranked = matched.assign(abs_beta=matched["beta_per_10_years"].abs()).sort_values(
        ["abs_beta", "feature"],
        ascending=[False, True],
    )
    row = ranked.iloc[0]
    return (
        f"{row.get('dataset_id', '')}:{row.get('feature', '')}:"
        f"{row.get('direction', '')}:beta={row.get('beta_per_10_years', np.nan):.4g}"
    )


def _format_age_range(age_min: object, age_max: object) -> str:
    if not np.isfinite(pd.to_numeric(pd.Series([age_min]), errors="coerce").iloc[0]):
        return "unknown"
    min_value = float(age_min)
    max_value = float(age_max)
    return f"{min_value:.1f}-{max_value:.1f}"


def _inventory_notes(age_context: list[str], obs: pd.DataFrame, donor_col: str) -> str:
    contexts = set(age_context)
    notes = []
    if contexts - {"adult"}:
        notes.append("non-adult stages present; adult-only estimates are primary")
    if donor_col:
        adult_donors = _n_donors_with_context(obs, donor_col, age_context, "adult")
        if adult_donors < DEFAULT_MIN_DONORS:
            notes.append("adult donor count below primary age-effect threshold")
    return "; ".join(notes) or "adult donor-level age metadata available"


def _missing_inventory_row(spec: CrossTissueDatasetSpec, status: str) -> dict[str, object]:
    return {
        "dataset_id": spec.dataset_id,
        "collection_id": spec.collection_id,
        "title": spec.title,
        "tissue_class": spec.tissue_class,
        "h5ad_path": spec.h5ad_path,
        "source_url": spec.source_url,
        "status": status,
        "n_cells": 0,
        "n_genes": 0,
        "n_donors": 0,
        "n_adult_donors": 0,
        "age_min_years": np.nan,
        "age_max_years": np.nan,
        "age_contexts": "",
        "donor_column": "",
        "age_column": "",
        "cell_label_column": "",
        "expression_scale": "",
        "primary_age_scope": PRIMARY_SCOPE,
        "notes": "Download the public H5AD before running estimates.",
    }


def _n_donors_with_context(obs: pd.DataFrame, donor_col: str, age_context: list[str], context: str) -> int:
    if not donor_col:
        return 0
    frame = pd.DataFrame({"donor_id": obs[donor_col].astype(str), "age_context": age_context})
    return int(frame[frame["age_context"].eq(context)]["donor_id"].nunique())


def _first_present(frame: pd.DataFrame, candidates: list[str]) -> str:
    for col in candidates:
        if col in frame.columns:
            return col
    return ""


def _mode_join(values: pd.Series) -> str:
    clean = values.dropna().astype(str)
    if clean.empty:
        return ""
    counts = clean.value_counts()
    top_count = counts.iloc[0]
    return ";".join(sorted(counts[counts.eq(top_count)].index.tolist()))


def _single_value(values: pd.Series) -> str:
    clean = values.dropna().astype(str).unique().tolist()
    return clean[0] if clean else ""


def _as_dense(matrix: Any) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        return np.asarray(matrix.toarray(), dtype=float)
    return np.asarray(matrix, dtype=float)
