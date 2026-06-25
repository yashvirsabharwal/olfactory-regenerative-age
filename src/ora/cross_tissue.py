"""Cross-tissue specificity planning and ORA feature classification."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .utils import ensure_parent


SPECIFICITY_CLASSES = (
    "olfactory_specific",
    "airway_nasal_shared",
    "pan_epithelial_regenerative",
    "immune_inflammatory_shared",
    "technical_yield_associated",
    "not_comparable",
)


@dataclass(frozen=True)
class FeatureClassification:
    specificity_class: str
    biology_theme: str
    confidence: str
    comparable_tissue_classes: str
    evidence_source: str
    rationale: str


def build_cross_tissue_candidate_matrix(config: dict[str, Any]) -> pd.DataFrame:
    """Build the M6.1 cross-tissue comparator candidate matrix."""

    rows: list[dict[str, object]] = []
    for row in config.get("public_data_exhaustion", {}).get("candidates", []):
        if not isinstance(row, dict) or str(row.get("species", "")).lower() != "human":
            continue
        decision = str(row.get("inclusion_decision", "")).lower()
        tissue = str(row.get("tissue", "")).lower()
        if not _is_cross_tissue_candidate(decision, tissue):
            continue
        tissue_class = _candidate_tissue_class(row)
        rows.append(
            {
                "dataset_id": _dataset_id(row.get("accession_or_dataset", "")),
                "accession_or_dataset": row.get("accession_or_dataset", ""),
                "source_url": row.get("source_url", ""),
                "tissue_class": tissue_class,
                "tissue": row.get("tissue", ""),
                "assay": row.get("assay", ""),
                "species": row.get("species", ""),
                "donor_or_sample_count": row.get("donor_or_sample_count", ""),
                "age_availability": row.get("age_availability", ""),
                "counts_availability": row.get("counts_availability", ""),
                "labels_availability": row.get("labels_availability", ""),
                "inclusion_decision": row.get("inclusion_decision", ""),
                "comparator_role": _comparator_role(tissue_class, row),
                "feature_harmonization_strategy": _harmonization_strategy(tissue_class, row),
                "local_status": _candidate_local_status(row),
                "notes": row.get("notes", ""),
            }
        )
    rows.extend(_placeholder_comparator_rows())
    columns = [
        "dataset_id",
        "accession_or_dataset",
        "source_url",
        "tissue_class",
        "tissue",
        "assay",
        "species",
        "donor_or_sample_count",
        "age_availability",
        "counts_availability",
        "labels_availability",
        "inclusion_decision",
        "comparator_role",
        "feature_harmonization_strategy",
        "local_status",
        "notes",
    ]
    return pd.DataFrame(rows, columns=columns)


def build_ora_cross_tissue_feature_classification(
    *,
    feature_matrix: pd.DataFrame,
    manifest: pd.DataFrame,
    feature_stability: pd.DataFrame | None = None,
    feature_interpretation: pd.DataFrame | None = None,
    comparator_matrix: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Classify ORA features by likely tissue specificity and comparator readiness."""

    feature_cols = [
        col
        for col in feature_matrix.columns
        if col != "donor_id" and pd.api.types.is_numeric_dtype(feature_matrix[col])
    ]
    age_assoc = _feature_age_associations(feature_matrix[["donor_id", *feature_cols]], manifest)
    stability = _summarize_feature_stability(feature_stability)
    interpretation = (
        feature_interpretation.set_index("feature", drop=False)
        if feature_interpretation is not None and "feature" in feature_interpretation
        else pd.DataFrame()
    )
    comparator_lookup = _comparator_lookup(comparator_matrix)
    rows: list[dict[str, object]] = []
    for feature in feature_cols:
        parsed = parse_feature_name(feature)
        existing_theme = ""
        existing_caution = ""
        if not interpretation.empty and feature in interpretation.index:
            existing_theme = str(interpretation.loc[feature].get("biology_theme", ""))
            existing_caution = str(interpretation.loc[feature].get("caution", ""))
        classification = classify_feature(feature, existing_theme=existing_theme)
        class_comparators = comparator_lookup.get(classification.specificity_class, [])
        assoc_row = age_assoc.loc[feature].to_dict() if feature in age_assoc.index else {}
        stability_row = stability.loc[feature].to_dict() if feature in stability.index else {}
        rows.append(
            {
                "feature": feature,
                "feature_kind": parsed["feature_kind"],
                "feature_label": parsed["feature_label"],
                "biology_theme": existing_theme or classification.biology_theme,
                "specificity_class": classification.specificity_class,
                "classification_confidence": classification.confidence,
                "comparable_tissue_classes": classification.comparable_tissue_classes,
                "candidate_comparator_ids": ";".join(class_comparators),
                "external_age_effect_status": _external_age_effect_status(
                    classification.specificity_class,
                    class_comparators,
                ),
                "gateway_age_beta_per_10_years": assoc_row.get("beta_per_10_years", np.nan),
                "gateway_age_p_value": assoc_row.get("p_value", np.nan),
                "gateway_age_fdr": assoc_row.get("fdr", np.nan),
                "gateway_age_direction": assoc_row.get("direction", "not_tested"),
                "gateway_age_status": assoc_row.get("status", "not_tested"),
                "max_abs_importance": stability_row.get("max_abs_importance", 0.0),
                "top_model": stability_row.get("top_model", ""),
                "max_selection_fraction": stability_row.get("max_selection_fraction", 0.0),
                "evidence_source": classification.evidence_source,
                "existing_interpretation_caution": existing_caution,
                "classification_rationale": classification.rationale,
            }
        )
    table = pd.DataFrame(rows)
    table["importance_rank"] = (
        table["max_abs_importance"].rank(method="min", ascending=False).astype(int)
    )
    return table.sort_values(["importance_rank", "feature"]).reset_index(drop=True)


def build_cross_tissue_specificity_summary(classification: pd.DataFrame) -> pd.DataFrame:
    """Summarize feature specificity classes for reviewer-facing reporting."""

    rows: list[dict[str, object]] = []
    for specificity_class in SPECIFICITY_CLASSES:
        frame = classification[classification["specificity_class"].eq(specificity_class)].copy()
        if frame.empty:
            rows.append(
                {
                    "specificity_class": specificity_class,
                    "n_features": 0,
                    "n_gateway_age_fdr_lt_0_05": 0,
                    "n_gateway_positive": 0,
                    "n_gateway_negative": 0,
                    "max_abs_importance": 0.0,
                    "top_features": "",
                    "external_age_effect_statuses": "",
                    "interpretation": _class_interpretation(specificity_class),
                }
            )
            continue
        top = frame.sort_values(["max_abs_importance", "feature"], ascending=[False, True]).head(6)
        rows.append(
            {
                "specificity_class": specificity_class,
                "n_features": int(frame.shape[0]),
                "n_gateway_age_fdr_lt_0_05": int((frame["gateway_age_fdr"] < 0.05).sum()),
                "n_gateway_positive": int(frame["gateway_age_direction"].eq("positive").sum()),
                "n_gateway_negative": int(frame["gateway_age_direction"].eq("negative").sum()),
                "max_abs_importance": float(frame["max_abs_importance"].max()),
                "top_features": ";".join(top["feature"].astype(str).tolist()),
                "external_age_effect_statuses": ";".join(
                    sorted(set(frame["external_age_effect_status"].astype(str)))
                ),
                "interpretation": _class_interpretation(specificity_class),
            }
        )
    return pd.DataFrame(rows)


def write_cross_tissue_specificity_figure(
    summary: pd.DataFrame,
    *,
    pdf_out: str | Path,
    png_out: str | Path | None = None,
) -> None:
    """Write a compact bar figure summarizing feature specificity classes."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402

    ordered = summary.copy()
    colors = {
        "olfactory_specific": "#3b6ea8",
        "airway_nasal_shared": "#55a868",
        "pan_epithelial_regenerative": "#c44e52",
        "immune_inflammatory_shared": "#8172b2",
        "technical_yield_associated": "#8c8c8c",
        "not_comparable": "#ccb974",
    }
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.barh(
        ordered["specificity_class"],
        ordered["n_features"],
        color=[colors.get(item, "#4c72b0") for item in ordered["specificity_class"]],
    )
    ax.set_xlabel("ORA feature count")
    ax.set_ylabel("Cross-tissue specificity class")
    ax.set_title("First-pass ORA cross-tissue specificity classification")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(ensure_parent(pdf_out))
    if png_out is not None:
        fig.savefig(ensure_parent(png_out), dpi=200)
    plt.close(fig)


def render_cross_tissue_specificity_plan(
    candidate_matrix: pd.DataFrame,
    summary: pd.DataFrame | None = None,
) -> str:
    """Render the M6.1 cross-tissue specificity plan."""

    lines = [
        "# Cross-Tissue Specificity Plan",
        "",
        "Updated: 2026-06-25",
        "",
        "## Goal",
        "",
        "Classify ORA signals as olfactory-specific, airway/nasal shared, pan-epithelial regenerative, immune/inflammatory shared, technical/yield-associated, or not comparable. The first pass below uses curated public-resource metadata and Gateway feature evidence; measured nasal/bronchial CELLxGENE comparator age effects are generated by `make cross-tissue-age-effects`.",
        "",
        "## Data Access Strategy",
        "",
        "- Use the curated public-data registry in `configs/external_datasets.yaml` as the source of truth for selected airway/nasal/lung comparators.",
        "- Use CELLxGENE Census for scalable future pulls of curated single-cell objects where age and cell labels are present.",
        "- Keep bulk/spatial datasets as marker/context checks, not direct donor-level single-cell replication.",
        "- Add explicit query-required rows for skin, gut, and PBMC/blood so the specificity claim remains bounded.",
        "",
        "## Comparator Panel",
        "",
        "| Dataset | Tissue class | Assay | Age metadata | Role | Status |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in candidate_matrix.iterrows():
        dataset = str(row.get("accession_or_dataset", ""))
        url = str(row.get("source_url", ""))
        label = f"[{dataset}]({url})" if url else dataset
        lines.append(
            "| "
            f"{label} | {row.get('tissue_class', '')} | {row.get('assay', '')} | "
            f"{row.get('age_availability', '')} | {row.get('comparator_role', '')} | "
            f"{row.get('local_status', '')} |"
        )
    if summary is not None and not summary.empty:
        lines.extend(
            [
                "",
                "## Current First-Pass Classification",
                "",
                "| Specificity class | Features | Age-associated features | Top features |",
                "| --- | ---: | ---: | --- |",
            ]
        )
        for _, row in summary.iterrows():
            lines.append(
                "| "
                f"{row.get('specificity_class', '')} | {row.get('n_features', 0)} | "
                f"{row.get('n_gateway_age_fdr_lt_0_05', 0)} | {row.get('top_features', '')} |"
            )
    lines.extend(
        [
            "",
            "## Claim Rules",
            "",
            "- `olfactory_specific`: allowed only for OSN, INP/neurogenic, olfactory-sustentacular, and olfactory-transduction signals unless external airway/lung age effects show the same feature.",
            "- `airway_nasal_shared`: ciliated, goblet, club, ionocyte, tuft, and respiratory/secretory signals should be described as shared airway/nasal epithelial context.",
            "- `pan_epithelial_regenerative`: basal activation, cycling, stress/senescence, and injury-repair programs should be described as regenerative epithelial themes, not olfactory-exclusive mechanisms.",
            "- `immune_inflammatory_shared`: immune cell-state and inflammatory module signals should be separated from epithelial regeneration claims.",
            "- `technical_yield_associated`: cell-yield and technical features remain controls and cannot support biological specificity.",
            "- `not_comparable`: disease-risk modules or ambiguous labels require manual review before claim use.",
            "",
            "## Companion Age-Effect Stage",
            "",
            "See `docs/cross_tissue_age_effects.md`, `results/tables/cross_tissue_age_effects.tsv`, and `results/tables/ora_cross_tissue_age_effect_summary.tsv` for donor-level public CELLxGENE nasal/bronchial comparator estimates. Remaining specificity expansion is mainly non-airway skin/gut/blood comparator querying and larger adult lung resources.",
            "",
        ]
    )
    return "\n".join(lines)


def write_cross_tissue_outputs(
    *,
    candidate_matrix: pd.DataFrame,
    classification: pd.DataFrame,
    summary: pd.DataFrame,
    plan_md: str,
    candidate_out: str | Path,
    classification_out: str | Path,
    summary_out: str | Path,
    plan_out: str | Path,
    figure_pdf: str | Path,
    figure_png: str | Path | None = None,
) -> None:
    """Write all M6 first-pass specificity artifacts."""

    candidate_matrix.to_csv(ensure_parent(candidate_out), sep="\t", index=False)
    classification.to_csv(ensure_parent(classification_out), sep="\t", index=False)
    summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    ensure_parent(plan_out).write_text(plan_md, encoding="utf-8")
    write_cross_tissue_specificity_figure(summary, pdf_out=figure_pdf, png_out=figure_png)


def parse_feature_name(feature: str) -> dict[str, str]:
    """Parse a project feature name into kind and label."""

    if feature.startswith("module_score__"):
        return {"feature_kind": "module_score", "feature_label": feature.removeprefix("module_score__")}
    if "__" in feature:
        kind, label = feature.split("__", 1)
        return {"feature_kind": kind, "feature_label": label}
    return {"feature_kind": "unknown", "feature_label": feature}


def classify_feature(feature: str, *, existing_theme: str = "") -> FeatureClassification:
    """Classify one ORA feature by likely cross-tissue specificity."""

    parsed = parse_feature_name(feature)
    label = parsed["feature_label"].lower()
    kind = parsed["feature_kind"]
    text = f"{kind} {label} {existing_theme.lower()}"
    if _contains_any(text, _TECHNICAL_TOKENS):
        return FeatureClassification(
            "technical_yield_associated",
            "technical/yield",
            "high",
            "Gateway technical covariates",
            "feature_name",
            "Feature maps to cell yield, chemistry, collection, or other technical covariates.",
        )
    if _contains_any(text, _IMMUNE_TOKENS):
        return FeatureClassification(
            "immune_inflammatory_shared",
            "immune/inflammatory",
            "high",
            "blood/PBMC; lung/nasal immune compartments; pan-tissue immune atlases",
            "feature_name;existing_interpretation",
            "Immune and inflammatory features are expected to be shared across tissues.",
        )
    if _contains_any(text, _OLFACTORY_TOKENS):
        return FeatureClassification(
            "olfactory_specific",
            "olfactory sensory/regenerative lineage",
            "high",
            "olfactory epithelium; airway/lung as negative/context comparators",
            "feature_name;module_definition",
            "Feature is tied to OSN, INP/neurogenic, sustentacular, or olfactory-transduction biology.",
        )
    if _contains_any(text, _AIRWAY_TOKENS):
        return FeatureClassification(
            "airway_nasal_shared",
            "airway/nasal epithelial",
            "high",
            "nasal respiratory; bronchial; lung airway",
            "feature_name;cell_state_label",
            "Feature maps to respiratory, ciliated, goblet, club, ionocyte, tuft, or secretory airway states.",
        )
    if _contains_any(text, _PAN_EPITHELIAL_TOKENS):
        return FeatureClassification(
            "pan_epithelial_regenerative",
            "pan-epithelial regeneration/stress",
            "medium",
            "nasal/airway; lung; skin; gut epithelium",
            "feature_name;module_definition",
            "Basal, cycling, injury-repair, stress, or senescence features are plausible shared epithelial programs.",
        )
    return FeatureClassification(
        "not_comparable",
        existing_theme or "uncertain",
        "low",
        "manual review required",
        "feature_name",
        "Feature label is ambiguous or disease-context-specific for cross-tissue specificity.",
    )


_IMMUNE_TOKENS = {
    "immune",
    "inflammatory",
    "antigen",
    "cdc1",
    "cdc2",
    "maturedc",
    "pdc",
    "macrophage",
    "monocyte",
    "classical",
    "nonclassical",
    "tcell",
    "cd4",
    "cd8",
    "nk",
    "bcell",
    "naiveb",
    "gcb",
    "plasmacell",
    "cytotoxic",
    "complement",
    "neuroinflammation",
    "lipid_associated",
}

_OLFACTORY_TOKENS = {
    "osn",
    "mosn",
    "iosn",
    "olf_sus",
    "olfactory",
    "omp",
    "adcy3",
    "neuron",
    "neuronal",
    "immature_neuron",
    "mature_olfactory_neuron",
    "cilia_olfactory_transduction",
    "neuroblast",
    "inp",
    "immature_to_mature",
    "lineage_fraction",
    "pre_dysfunctional",
    "dysfunctional",
    "stressed_mosn",
}

_AIRWAY_TOKENS = {
    "goblet",
    "club",
    "ciliated",
    "multiciliated",
    "deuterosomal",
    "respiratory",
    "secretory",
    "mucous",
    "serous",
    "bowman",
    "gland",
    "ionocyte",
    "tuft",
    "mv_type",
    "mv_ionocyte",
    "chemosensory",
}

_PAN_EPITHELIAL_TOKENS = {
    "hbc",
    "basal",
    "suprabasal",
    "cycling",
    "proliferating",
    "proliferation",
    "activation",
    "injury",
    "stress",
    "senescence",
    "sasp",
    "repair",
    "metaplasia",
}

_TECHNICAL_TOKENS = {
    "technical",
    "yield",
    "total_cells",
    "lineage_cells",
    "mature_neurons",
    "chemistry",
    "collection",
    "site",
}


def _feature_age_associations(feature_matrix: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    merged = feature_matrix.merge(manifest[["donor_id", "age", "usable_for_ora_training"]], on="donor_id", how="inner")
    if "usable_for_ora_training" in merged:
        merged = merged[merged["usable_for_ora_training"].astype(bool)]
    merged = merged[pd.to_numeric(merged["age"], errors="coerce").notna()].copy()
    merged["age"] = pd.to_numeric(merged["age"], errors="coerce")
    rows: list[dict[str, object]] = []
    for feature in [col for col in feature_matrix.columns if col != "donor_id"]:
        values = pd.to_numeric(merged[feature], errors="coerce")
        mask = values.notna() & merged["age"].notna()
        n = int(mask.sum())
        if n < 6:
            rows.append(_age_assoc_row(feature, n, status="insufficient"))
            continue
        y = values[mask].to_numpy(dtype=float)
        if np.nanstd(y) <= 0:
            rows.append(_age_assoc_row(feature, n, status="constant"))
            continue
        x = (merged.loc[mask, "age"].to_numpy(dtype=float) - float(merged.loc[mask, "age"].mean())) / 10.0
        result = stats.linregress(x, y)
        direction = "positive" if result.slope > 0 else "negative" if result.slope < 0 else "flat"
        rows.append(
            {
                "feature": feature,
                "n": n,
                "beta_per_10_years": float(result.slope),
                "p_value": float(result.pvalue),
                "direction": direction,
                "status": "ok",
            }
        )
    table = pd.DataFrame(rows).set_index("feature", drop=False)
    table["fdr"] = _benjamini_hochberg(table["p_value"])
    return table


def _age_assoc_row(feature: str, n: int, *, status: str) -> dict[str, object]:
    return {
        "feature": feature,
        "n": n,
        "beta_per_10_years": np.nan,
        "p_value": np.nan,
        "direction": "not_tested",
        "status": status,
    }


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


def _summarize_feature_stability(feature_stability: pd.DataFrame | None) -> pd.DataFrame:
    if feature_stability is None or feature_stability.empty:
        return pd.DataFrame()
    rows = []
    for feature, frame in feature_stability.groupby("feature", sort=False):
        ranked = frame.sort_values(["abs_mean_importance", "selection_fraction"], ascending=False)
        top = ranked.iloc[0]
        rows.append(
            {
                "feature": feature,
                "max_abs_importance": float(ranked["abs_mean_importance"].max()),
                "top_model": top.get("model", ""),
                "max_selection_fraction": float(ranked["selection_fraction"].max()),
            }
        )
    return pd.DataFrame(rows).set_index("feature", drop=False)


def _comparator_lookup(comparator_matrix: pd.DataFrame | None) -> dict[str, list[str]]:
    lookup: dict[str, list[str]] = {key: [] for key in SPECIFICITY_CLASSES}
    if comparator_matrix is None or comparator_matrix.empty:
        return lookup
    for _, row in comparator_matrix.iterrows():
        dataset_id = str(row.get("dataset_id", ""))
        tissue_class = str(row.get("tissue_class", ""))
        status = str(row.get("local_status", ""))
        if not dataset_id or status == "query_required":
            continue
        if tissue_class in {"nasal", "airway", "lung", "airway_lung", "olfactory_respiratory"}:
            lookup["airway_nasal_shared"].append(dataset_id)
        if tissue_class in {"nasal", "airway", "lung", "airway_lung"}:
            lookup["pan_epithelial_regenerative"].append(dataset_id)
            lookup["immune_inflammatory_shared"].append(dataset_id)
        if tissue_class in {"olfactory_respiratory"}:
            lookup["olfactory_specific"].append(dataset_id)
    return {key: sorted(set(values)) for key, values in lookup.items()}


def _external_age_effect_status(specificity_class: str, comparator_ids: list[str]) -> str:
    if specificity_class == "technical_yield_associated":
        return "internal_technical_control_only"
    if specificity_class == "not_comparable":
        return "not_comparable_manual_review"
    if comparator_ids:
        return "comparator_selected_external_age_effect_pending"
    return "comparator_query_required"


def _contains_any(text: str, tokens: set[str]) -> bool:
    return any(token in text for token in tokens)


def _is_cross_tissue_candidate(decision: str, tissue: str) -> bool:
    if "olfactory" in tissue and "respiratory" not in tissue:
        return False
    terms = {
        "specificity",
        "comparator",
        "lung aging",
        "nasal",
        "airway",
        "respiratory",
        "oe-vs-re",
        "sampling-method",
    }
    tissue_terms = {"lung", "nasal", "respiratory", "bronchial", "airway", "sinus"}
    return any(term in decision for term in terms) or any(term in tissue for term in tissue_terms)


def _candidate_tissue_class(row: dict[str, Any]) -> str:
    tissue = str(row.get("tissue", "")).lower()
    if "olfactory" in tissue and "respiratory" in tissue:
        return "olfactory_respiratory"
    if "lung" in tissue and ("airway" in tissue or "nose" in tissue or "respiratory" in tissue):
        return "airway_lung"
    if "lung" in tissue:
        return "lung"
    if "bronchial" in tissue or "airway" in tissue:
        return "airway"
    if "nasal" in tissue or "sinus" in tissue:
        return "nasal"
    return "other_context"


def _dataset_id(accession_or_dataset: object) -> str:
    text = str(accession_or_dataset).strip()
    if not text:
        return ""
    return (
        text.lower()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("+", "plus")
        .replace(".", "")
        .replace("-", "_")
    )


def _comparator_role(tissue_class: str, row: dict[str, Any]) -> str:
    assay = str(row.get("assay", "")).lower()
    if "bulk" in assay:
        return "marker_context_only"
    if "spatial" in assay or "geomx" in assay:
        return "spatial_context_only"
    if tissue_class == "olfactory_respiratory":
        return "olfactory_vs_respiratory_marker_control"
    if tissue_class in {"nasal", "airway", "airway_lung", "lung"}:
        return "single_cell_cross_tissue_age_comparator"
    return "context_only"


def _harmonization_strategy(tissue_class: str, row: dict[str, Any]) -> str:
    assay = str(row.get("assay", "")).lower()
    if "bulk" in assay or "spatial" in assay or "geomx" in assay:
        return "score marker/module panels only; do not treat as donor-level single-cell replication"
    if tissue_class in {"nasal", "airway", "airway_lung", "lung"}:
        return "map coarse epithelial/immune labels and score ORA module markers by donor age"
    return "manual harmonization required"


def _candidate_local_status(row: dict[str, Any]) -> str:
    accession = str(row.get("accession_or_dataset", "")).lower()
    source_url = str(row.get("source_url", "")).lower()
    if "af500bec-a5d6-4569-9f50-4314f2e7a011" in source_url or "xu et al" in accession:
        return "adult_age_effects_computed"
    if "625f6bf4-2f33-4942-962e-35243d284837" in source_url or "lungmap broad-age" in accession:
        return "context_only_insufficient_adult_donors"
    return "candidate_selected_external_effects_pending"


def _placeholder_comparator_rows() -> list[dict[str, object]]:
    placeholders = [
        ("skin_epithelium_query_required", "Skin epithelium aging CELLxGENE query", "skin"),
        ("gut_epithelium_query_required", "Gut epithelium aging CELLxGENE query", "gut"),
        ("pbmc_blood_query_required", "PBMC/blood immune aging query", "blood_immune"),
    ]
    rows = []
    for dataset_id, label, tissue_class in placeholders:
        rows.append(
            {
                "dataset_id": dataset_id,
                "accession_or_dataset": label,
                "source_url": "https://cellxgene.cziscience.com/",
                "tissue_class": tissue_class,
                "tissue": tissue_class.replace("_", " "),
                "assay": "CELLxGENE query",
                "species": "human",
                "donor_or_sample_count": "pending query",
                "age_availability": "query required",
                "counts_availability": "query required",
                "labels_availability": "query required",
                "inclusion_decision": "query required for non-airway specificity",
                "comparator_role": "missing_non_airway_specificity_panel",
                "feature_harmonization_strategy": "select adult human datasets with donor age, labels, and counts",
                "local_status": "query_required",
                "notes": "Placeholder documents that full specificity versus non-airway tissues is not yet complete.",
            }
        )
    return rows


def _class_interpretation(specificity_class: str) -> str:
    return {
        "olfactory_specific": "Olfactory-lineage signals; strongest candidates for tissue-specific ORA biology.",
        "airway_nasal_shared": "Likely shared nasal/airway epithelial composition or barrier/secretory programs.",
        "pan_epithelial_regenerative": "Basal, cycling, stress, and repair programs likely shared across epithelia.",
        "immune_inflammatory_shared": "Immune and inflammatory signals likely shared across tissue immune compartments.",
        "technical_yield_associated": "Technical/yield signals; control only, not biological specificity evidence.",
        "not_comparable": "Ambiguous or disease-context features needing manual review.",
    }.get(specificity_class, "")
