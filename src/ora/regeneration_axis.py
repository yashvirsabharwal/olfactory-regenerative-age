"""Regeneration-axis curation and evidence joins for ORA features."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .utils import ensure_parent


THEMES = (
    "basal_quiescence",
    "basal_activation",
    "cycling_proliferation",
    "neural_progenitor",
    "immature_osn",
    "mature_osn",
    "sustentacular_barrier_detox",
    "bowman_gland_secretory",
    "respiratory_metaplasia_ciliated_goblet",
    "immune_inflammatory",
    "stress_senescence",
    "ecm_remodeling",
    "technical_yield",
)

THEME_LABELS = {
    "basal_quiescence": "Basal quiescence",
    "basal_activation": "Basal activation",
    "cycling_proliferation": "Cycling/proliferation",
    "neural_progenitor": "Neural progenitor",
    "immature_osn": "Immature OSN",
    "mature_osn": "Mature OSN",
    "sustentacular_barrier_detox": "Sustentacular/barrier/detox",
    "bowman_gland_secretory": "Bowman gland/secretory",
    "respiratory_metaplasia_ciliated_goblet": "Respiratory metaplasia/ciliated/goblet",
    "immune_inflammatory": "Immune/inflammatory",
    "stress_senescence": "Stress/senescence",
    "ecm_remodeling": "ECM/remodeling",
    "technical_yield": "Technical/yield",
}

THEME_EXPECTED_DIRECTIONS = {
    "basal_quiescence": "positive",
    "basal_activation": "context_dependent",
    "cycling_proliferation": "context_dependent",
    "neural_progenitor": "negative",
    "immature_osn": "negative",
    "mature_osn": "negative",
    "sustentacular_barrier_detox": "context_dependent",
    "bowman_gland_secretory": "context_dependent",
    "respiratory_metaplasia_ciliated_goblet": "positive",
    "immune_inflammatory": "positive",
    "stress_senescence": "positive",
    "ecm_remodeling": "positive",
    "technical_yield": "unknown",
}

THEME_INTERPRETATIONS = {
    "basal_quiescence": "Reserve basal-cell identity/quiescence features.",
    "basal_activation": "Basal activation and injury-response features.",
    "cycling_proliferation": "Cell-cycle or proliferating-state features.",
    "neural_progenitor": "Globose-basal/intermediate neural progenitor features.",
    "immature_osn": "Immature olfactory sensory neuron features.",
    "mature_osn": "Mature olfactory sensory neuron and transduction features.",
    "sustentacular_barrier_detox": "Support/barrier/detox epithelial features.",
    "bowman_gland_secretory": "Bowman gland and secretory epithelial features.",
    "respiratory_metaplasia_ciliated_goblet": "Respiratory-like ciliated/goblet airway features.",
    "immune_inflammatory": "Immune-cell and inflammatory-state features.",
    "stress_senescence": "Cell stress, dysfunction, senescence, or disease-prior features.",
    "ecm_remodeling": "Extracellular-matrix and remodeling features.",
    "technical_yield": "Yield, assay, or sampling-sensitive features.",
}

THEME_CITATIONS = {
    "basal_quiescence": "PMID:21677159;PMID:27560601;PMID:37260223;PMID:42228282",
    "basal_activation": "PMID:21677159;PMID:37260223;PMID:42228282",
    "cycling_proliferation": "PMID:28506465;PMID:29934351;PMID:42228282",
    "neural_progenitor": "PMID:28506465;PMID:24920630;PMID:21486944;PMID:42228282",
    "immature_osn": "PMID:28506465;PMID:29934351;PMID:38903957",
    "mature_osn": "PMID:28506465;PMID:29934351;PMID:38903957",
    "sustentacular_barrier_detox": "PMID:38903957;PMID:42228282",
    "bowman_gland_secretory": "PMID:38903957;PMID:42228282",
    "respiratory_metaplasia_ciliated_goblet": "PMID:41461651;PMID:38903957",
    "immune_inflammatory": "PMID:41461651;PMID:29934351;PMID:42228282",
    "stress_senescence": "PMID:29934351;PMID:41461651;PMID:42228282",
    "ecm_remodeling": "PMID:38903957;PMID:42228282",
    "technical_yield": "project_provenance",
}


@dataclass(frozen=True)
class RegenerationFeatureCuration:
    primary_theme: str
    secondary_theme: str
    expected_aging_direction: str
    expected_direction_basis: str
    confidence: str
    evidence_source: str
    evidence_citations: str
    biological_role: str
    interpretation_caution: str


def parse_feature_name(feature: str) -> dict[str, str]:
    """Split ORA feature names into kind and human-readable labels."""

    if "__" in feature:
        kind, label = feature.split("__", 1)
    else:
        kind, label = "feature", feature
    return {
        "feature_kind": kind,
        "feature_label": label,
        "feature_display": label.replace("_", " "),
    }


def classify_regeneration_feature(feature: str) -> RegenerationFeatureCuration:
    """Map one ORA feature to the controlled regeneration-axis vocabulary."""

    parsed = parse_feature_name(feature)
    label = parsed["feature_label"].lower()
    kind = parsed["feature_kind"]
    theme = _primary_theme(label, kind)
    secondary = _secondary_theme(label, theme)
    expected = _expected_direction(label, theme)
    confidence = _confidence(label, kind, theme)
    evidence = _evidence_source(kind, label)
    citations = _citations(theme, label)
    role = _biological_role(label, theme)
    caution = _caution(label, kind, theme, expected)
    return RegenerationFeatureCuration(
        primary_theme=theme,
        secondary_theme=secondary,
        expected_aging_direction=expected,
        expected_direction_basis=_expected_basis(theme, expected),
        confidence=confidence,
        evidence_source=evidence,
        evidence_citations=citations,
        biological_role=role,
        interpretation_caution=caution,
    )


def build_regeneration_feature_resource_map(feature_matrix: pd.DataFrame) -> pd.DataFrame:
    """Build the static feature-to-regeneration-theme map for every ORA feature."""

    rows = []
    for feature in _numeric_feature_columns(feature_matrix):
        parsed = parse_feature_name(feature)
        curation = classify_regeneration_feature(feature)
        rows.append(
            {
                "feature": feature,
                **parsed,
                "primary_theme": curation.primary_theme,
                "primary_theme_label": THEME_LABELS[curation.primary_theme],
                "secondary_theme": curation.secondary_theme,
                "expected_aging_direction": curation.expected_aging_direction,
                "expected_direction_basis": curation.expected_direction_basis,
                "confidence": curation.confidence,
                "evidence_source": curation.evidence_source,
                "evidence_citations": curation.evidence_citations,
                "biological_role": curation.biological_role,
                "interpretation_caution": curation.interpretation_caution,
            }
        )
    return pd.DataFrame(rows, columns=_resource_columns())


def build_regeneration_axis_feature_map(
    *,
    feature_matrix: pd.DataFrame,
    manifest: pd.DataFrame,
    feature_stability: pd.DataFrame | None = None,
    cross_tissue_classification: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Join regeneration curation to Gateway age, importance, and specificity evidence."""

    resource = build_regeneration_feature_resource_map(feature_matrix)
    feature_cols = resource["feature"].astype(str).tolist()
    age = compute_feature_age_associations(feature_matrix[["donor_id", *feature_cols]], manifest)
    stability = summarize_feature_stability(feature_stability)
    cross_tissue = _cross_tissue_lookup(cross_tissue_classification)
    rows = []
    for _, row in resource.iterrows():
        feature = str(row["feature"])
        age_row = age.loc[feature].to_dict() if feature in age.index else {}
        stability_row = stability.loc[feature].to_dict() if feature in stability.index else {}
        cross_row = cross_tissue.loc[feature].to_dict() if feature in cross_tissue.index else {}
        observed_direction = age_row.get("direction", "not_tested")
        fdr = age_row.get("fdr", np.nan)
        rows.append(
            {
                **row.to_dict(),
                "gateway_age_beta_per_10_years": age_row.get("beta_per_10_years", np.nan),
                "gateway_age_p_value": age_row.get("p_value", np.nan),
                "gateway_age_fdr": fdr,
                "gateway_age_direction": observed_direction,
                "gateway_age_status": age_row.get("status", "not_tested"),
                "observed_vs_expected": _observed_vs_expected(
                    str(row["expected_aging_direction"]),
                    str(observed_direction),
                    fdr,
                ),
                "max_abs_importance": stability_row.get("max_abs_importance", 0.0),
                "top_model": stability_row.get("top_model", ""),
                "max_selection_fraction": stability_row.get("max_selection_fraction", 0.0),
                "mean_selection_fraction": stability_row.get("mean_selection_fraction", 0.0),
                "supporting_models": stability_row.get("supporting_models", ""),
                "specificity_class": cross_row.get("specificity_class", "not_classified"),
                "cross_tissue_confidence": cross_row.get("classification_confidence", ""),
                "external_age_effect_status": cross_row.get("external_age_effect_status", ""),
            }
        )
    table = pd.DataFrame(rows, columns=_result_columns())
    table["importance_rank"] = (
        table["max_abs_importance"].rank(method="min", ascending=False).astype(int)
    )
    return table.sort_values(["importance_rank", "feature"]).reset_index(drop=True)


def build_regeneration_axis_theme_summary(feature_map: pd.DataFrame) -> pd.DataFrame:
    """Summarize feature-map evidence by controlled regeneration theme."""

    rows = []
    for theme in THEMES:
        frame = feature_map[feature_map["primary_theme"].eq(theme)].copy()
        if frame.empty:
            rows.append(_empty_theme_row(theme))
            continue
        top = frame.sort_values(["max_abs_importance", "feature"], ascending=[False, True]).head(6)
        observed = frame["gateway_age_direction"].astype(str)
        expected = frame["expected_aging_direction"].astype(str)
        fdr = pd.to_numeric(frame["gateway_age_fdr"], errors="coerce")
        rows.append(
            {
                "primary_theme": theme,
                "primary_theme_label": THEME_LABELS[theme],
                "n_features": int(frame.shape[0]),
                "n_age_tested": int(frame["gateway_age_status"].eq("ok").sum()),
                "n_gateway_age_fdr_lt_0_05": int((fdr < 0.05).sum()),
                "n_observed_positive": int(observed.eq("positive").sum()),
                "n_observed_negative": int(observed.eq("negative").sum()),
                "n_expected_positive": int(expected.eq("positive").sum()),
                "n_expected_negative": int(expected.eq("negative").sum()),
                "n_direction_aligned": int(frame["observed_vs_expected"].eq("aligned").sum()),
                "n_direction_opposite": int(frame["observed_vs_expected"].eq("opposite").sum()),
                "max_abs_importance": float(frame["max_abs_importance"].max()),
                "median_abs_age_beta_per_10_years": float(
                    pd.to_numeric(frame["gateway_age_beta_per_10_years"], errors="coerce").abs().median()
                ),
                "top_features": ";".join(top["feature"].astype(str).tolist()),
                "theme_interpretation": THEME_INTERPRETATIONS[theme],
                "theme_citations": THEME_CITATIONS[theme],
            }
        )
    return pd.DataFrame(rows, columns=_summary_columns())


def compute_feature_age_associations(feature_matrix: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    """Run simple donor-level feature~age associations for all numeric ORA features."""

    if "donor_id" not in feature_matrix or "donor_id" not in manifest:
        raise KeyError("feature_matrix and manifest must both contain donor_id.")
    if "age" not in manifest:
        raise KeyError("manifest must contain age.")
    manifest_cols = ["donor_id", "age"]
    if "usable_for_ora_training" in manifest:
        manifest_cols.append("usable_for_ora_training")
    frame = feature_matrix.merge(manifest[manifest_cols], on="donor_id", how="left")
    if "usable_for_ora_training" in frame:
        frame = frame[frame["usable_for_ora_training"].fillna(False).astype(bool)].copy()
    frame["age"] = pd.to_numeric(frame["age"], errors="coerce")
    rows = []
    for feature in _numeric_feature_columns(feature_matrix):
        values = pd.to_numeric(frame[feature], errors="coerce")
        mask = frame["age"].notna() & values.notna()
        n = int(mask.sum())
        if n < 8 or values.loc[mask].nunique(dropna=True) < 2:
            rows.append(_age_row(feature, n, np.nan, np.nan, np.nan, np.nan, "not_tested"))
            continue
        slope, intercept, r_value, p_value, stderr = stats.linregress(frame.loc[mask, "age"], values.loc[mask])
        del intercept, r_value
        beta = float(slope * 10.0)
        direction = "positive" if beta > 0 else "negative" if beta < 0 else "flat"
        rows.append(_age_row(feature, n, beta, float(stderr * 10.0), np.nan, float(p_value), direction))
    result = pd.DataFrame(rows)
    if result.empty:
        return result.set_index("feature")
    ok = result["p_value"].notna()
    result["fdr"] = np.nan
    if ok.any():
        result.loc[ok, "fdr"] = _benjamini_hochberg(result.loc[ok, "p_value"].to_numpy(dtype=float))
    return result.set_index("feature", drop=False)


def summarize_feature_stability(feature_stability: pd.DataFrame | None) -> pd.DataFrame:
    """Summarize repeated-CV model support per feature."""

    columns = [
        "feature",
        "max_abs_importance",
        "top_model",
        "max_selection_fraction",
        "mean_selection_fraction",
        "supporting_models",
    ]
    if feature_stability is None or feature_stability.empty or "feature" not in feature_stability:
        return pd.DataFrame(columns=columns).set_index("feature")
    frame = feature_stability.copy()
    frame["abs_importance"] = pd.to_numeric(
        frame.get("abs_mean_importance", frame.get("mean_importance", 0.0)),
        errors="coerce",
    ).abs()
    frame["selection_fraction"] = pd.to_numeric(
        frame.get("selection_fraction", 0.0),
        errors="coerce",
    ).fillna(0.0)
    rows = []
    for feature, group in frame.groupby("feature", observed=True):
        group = group.sort_values(["abs_importance", "selection_fraction"], ascending=[False, False])
        top = group.iloc[0]
        supported = group[group["selection_fraction"].gt(0)]
        models = sorted(supported.get("model", pd.Series(dtype=str)).astype(str).unique().tolist())
        rows.append(
            {
                "feature": feature,
                "max_abs_importance": float(group["abs_importance"].max()),
                "top_model": str(top.get("model", "")),
                "max_selection_fraction": float(group["selection_fraction"].max()),
                "mean_selection_fraction": float(group["selection_fraction"].mean()),
                "supporting_models": ";".join(models),
            }
        )
    return pd.DataFrame(rows, columns=columns).set_index("feature", drop=False)


def write_regeneration_axis_figure(
    summary: pd.DataFrame,
    *,
    pdf_out: str | Path,
    png_out: str | Path | None = None,
) -> None:
    """Write a compact feature-count/importance figure by regeneration theme."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402

    plot = summary[summary["n_features"].gt(0)].copy()
    plot = plot.sort_values(["n_features", "max_abs_importance"], ascending=[True, True])
    colors = _theme_colors(plot["primary_theme"].tolist())
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), sharey=True)
    axes[0].barh(plot["primary_theme_label"], plot["n_features"], color=colors)
    axes[0].set_xlabel("ORA feature count")
    axes[0].set_ylabel("Regeneration-axis theme")
    axes[0].set_title("Feature map coverage")
    axes[0].grid(axis="x", alpha=0.25)
    axes[1].barh(plot["primary_theme_label"], plot["max_abs_importance"], color=colors)
    axes[1].set_xlabel("Maximum repeated-CV importance")
    axes[1].set_title("Strongest model evidence")
    axes[1].grid(axis="x", alpha=0.25)
    fig.suptitle("ORA regeneration-axis feature map", y=0.98)
    fig.tight_layout()
    fig.savefig(ensure_parent(pdf_out))
    if png_out is not None:
        fig.savefig(ensure_parent(png_out), dpi=200)
    plt.close(fig)


def write_regeneration_axis_outputs(
    *,
    resource_map: pd.DataFrame,
    feature_map: pd.DataFrame,
    summary: pd.DataFrame,
    resource_out: str | Path,
    feature_map_out: str | Path,
    summary_out: str | Path,
    figure_pdf: str | Path,
    figure_png: str | Path | None = None,
) -> None:
    """Write all M5.1 regeneration-axis artifacts."""

    resource_map.to_csv(ensure_parent(resource_out), sep="\t", index=False)
    feature_map.to_csv(ensure_parent(feature_map_out), sep="\t", index=False)
    summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    write_regeneration_axis_figure(summary, pdf_out=figure_pdf, png_out=figure_png)


def _primary_theme(label: str, kind: str) -> str:
    if any(term in label for term in ["quiescent_hbc", "hbc_identity"]):
        return "basal_quiescence"
    if any(term in label for term in ["activation", "activated_hbc", "suprabasal"]):
        return "basal_activation"
    if any(term in label for term in ["cycling", "proliferating"]):
        return "cycling_proliferation"
    if any(term in label for term in ["early_inp", "late_inp", "progenitor", "neuroblast"]):
        return "neural_progenitor"
    if any(term in label for term in ["early_iosn", "late_iosn", "immature"]):
        return "immature_osn"
    if any(
        term in label
        for term in [
            "early_mature_mosn",
            "fully_mature_mosn",
            "mature_",
            "mature_mosn",
            "olfactory_transduction",
            "neuronal_fraction",
            "lineage_fraction",
            "olfactory_neuron",
        ]
    ):
        return "mature_osn"
    if any(term in label for term in ["stressed", "pre_dysfunctional", "senescence", "sasp"]):
        return "stress_senescence"
    if any(term in label for term in ["ad_amyloid", "pd_synuclein"]):
        return "stress_senescence"
    if _is_immune(label):
        return "immune_inflammatory"
    if any(
        term in label
        for term in [
            "goblet",
            "multiciliated",
            "deuterosomal",
            "ciliated",
            "club",
            "ionocyte",
            "tuft",
            "mv_",
        ]
    ):
        return "respiratory_metaplasia_ciliated_goblet"
    if any(term in label for term in ["bowman", "mucous_gland", "serous", "secretory"]):
        return "bowman_gland_secretory"
    if any(term in label for term in ["olf_sus", "sustentacular", "detox"]):
        return "sustentacular_barrier_detox"
    if any(term in label for term in ["ecm", "collagen", "matrix", "remodel"]):
        return "ecm_remodeling"
    if kind in {"ratio", "prop", "clr"}:
        return "technical_yield" if "total" in label or "n_cells" in label else "stress_senescence"
    return "technical_yield"


def _secondary_theme(label: str, primary_theme: str) -> str:
    secondary = []
    if _is_immune(label) and primary_theme != "immune_inflammatory":
        secondary.append("immune_inflammatory")
    if any(term in label for term in ["cycling", "proliferating"]) and primary_theme != "cycling_proliferation":
        secondary.append("cycling_proliferation")
    if any(term in label for term in ["stress", "senescence", "dysfunctional", "sasp"]):
        if primary_theme != "stress_senescence":
            secondary.append("stress_senescence")
    if "ratio" in label and primary_theme not in {"technical_yield"}:
        secondary.append("technical_yield")
    return ";".join(dict.fromkeys(secondary))


def _expected_direction(label: str, theme: str) -> str:
    if "activated_to_quiescent_hbc" in label:
        return "context_dependent"
    if "immature_to_mature" in label:
        return "positive"
    if "mature_mosn_to_iosn" in label:
        return "negative"
    if "stressed_to_mature_mosn" in label:
        return "positive"
    if "inp_to_activated_hbc" in label:
        return "negative"
    if "inp_to_iosn" in label:
        return "context_dependent"
    if "progenitor_to_neuron" in label:
        return "context_dependent"
    if "lineage_fraction" in label:
        return "negative"
    if "mature_neuron_fraction" in label or "neuronal_fraction" in label:
        return "negative"
    if any(term in label for term in ["ad_amyloid", "pd_synuclein"]):
        return "unknown"
    return THEME_EXPECTED_DIRECTIONS[theme]


def _expected_basis(theme: str, expected: str) -> str:
    if expected == "unknown":
        return "No prespecified direction; retained as exploratory or technical evidence."
    if expected == "context_dependent":
        return "Direction depends on injury, sampling context, or numerator/denominator balance."
    if theme in {"neural_progenitor", "immature_osn", "mature_osn"}:
        return "Aging is expected to reduce regenerative neurogenesis and OSN maturation."
    if theme == "respiratory_metaplasia_ciliated_goblet":
        return "Aging/inflammation is expected to promote respiratory-like metaplasia."
    if theme in {"immune_inflammatory", "stress_senescence", "ecm_remodeling"}:
        return "Aging is expected to increase inflammatory, stress, and remodeling signals."
    if theme == "basal_quiescence":
        return "Reserve basal-cell/quiescent-state features can increase with impaired turnover."
    return f"Expected direction follows curated {THEME_LABELS[theme].lower()} literature."


def _confidence(label: str, kind: str, theme: str) -> str:
    if any(term in label for term in ["ad_amyloid", "pd_synuclein"]):
        return "low"
    if kind == "module_score":
        return "medium"
    if kind == "ratio":
        return "medium"
    if theme in {
        "basal_quiescence",
        "basal_activation",
        "neural_progenitor",
        "immature_osn",
        "mature_osn",
        "immune_inflammatory",
    }:
        return "high"
    if theme in {"respiratory_metaplasia_ciliated_goblet", "stress_senescence"}:
        return "medium"
    return "medium"


def _evidence_source(kind: str, label: str) -> str:
    sources = []
    if kind == "module_score":
        sources.extend(["module_definition", "literature"])
    elif kind == "ratio":
        sources.extend(["feature_engineering", "cell_state_annotation", "literature"])
    else:
        sources.extend(["cell_state_annotation", "literature"])
    if _is_immune(label):
        sources.append("cross_tissue_prior")
    return ";".join(dict.fromkeys(sources))


def _citations(theme: str, label: str) -> str:
    citations = [THEME_CITATIONS[theme]]
    if any(term in label for term in ["notch", "lgr5"]):
        citations.append("PMID:29739871;PMID:33390928")
    if any(term in label for term in ["wnt", "lgr5"]):
        citations.append("PMID:24920630;PMID:21486944")
    if any(term in label for term in ["yap", "tead"]):
        citations.append("PMID:35148842")
    return ";".join(dict.fromkeys(";".join(citations).split(";")))


def _biological_role(label: str, theme: str) -> str:
    display = label.replace("_", " ")
    if theme == "stress_senescence" and any(term in label for term in ["ad_amyloid", "pd_synuclein"]):
        return f"{display} disease-prior module; use only as exploratory stress-context evidence."
    return f"{display} mapped to {THEME_LABELS[theme].lower()}."


def _caution(label: str, kind: str, theme: str, expected: str) -> str:
    if any(term in label for term in ["ad_amyloid", "pd_synuclein"]):
        return "Disease-prior module; not evidence of disease prediction, diagnosis, or causality."
    if kind == "ratio":
        return "Derived ratio; validate numerator and denominator state stability separately."
    if kind == "module_score":
        return "Module score supports biological interpretation but is not independent validation."
    if theme == "immune_inflammatory":
        return "Separate inflammatory/immune composition from olfactory-regeneration-specific claims."
    if theme == "technical_yield" or expected == "context_dependent":
        return "Interpret with donor yield, chemistry, collection-method, and context sensitivity checks."
    return "Associational donor-level feature; do not interpret as direct lineage flux."


def _is_immune(label: str) -> bool:
    immune_terms = [
        "antigen_presenting",
        "macrophage",
        "proliferating_mac",
        "maturedc",
        "proliferating_dc",
        "lipid_associated",
        "classical",
        "nonclassical",
        "cdc",
        "pdc",
        "dendritic",
        "cd56",
        "nk",
        "immune",
        "inflammation",
        "inflammatory",
        "complement",
        "cytotoxic",
        "naive",
        "tcell",
        "cd4",
        "cd8",
        "plasma",
        "gcb",
        "neuroinflammation",
    ]
    return any(term in label for term in immune_terms)


def _numeric_feature_columns(feature_matrix: pd.DataFrame) -> list[str]:
    return [
        col
        for col in feature_matrix.columns
        if col != "donor_id" and pd.api.types.is_numeric_dtype(feature_matrix[col])
    ]


def _age_row(
    feature: str,
    n: int,
    beta: float,
    standard_error: float,
    t_value: float,
    p_value: float,
    direction: str,
) -> dict[str, object]:
    return {
        "feature": feature,
        "n": n,
        "beta_per_10_years": beta,
        "standard_error": standard_error,
        "t_value": t_value,
        "p_value": p_value,
        "direction": direction,
        "status": "ok" if direction in {"positive", "negative", "flat"} else "not_tested",
    }


def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    n = len(ranked)
    adjusted = np.empty(n, dtype=float)
    cumulative = 1.0
    for idx in range(n - 1, -1, -1):
        rank = idx + 1
        cumulative = min(cumulative, ranked[idx] * n / rank)
        adjusted[idx] = cumulative
    result = np.empty(n, dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def _cross_tissue_lookup(classification: pd.DataFrame | None) -> pd.DataFrame:
    columns = [
        "feature",
        "specificity_class",
        "classification_confidence",
        "external_age_effect_status",
    ]
    if classification is None or classification.empty or "feature" not in classification:
        return pd.DataFrame(columns=columns).set_index("feature")
    present = [col for col in columns if col in classification.columns]
    return classification[present].drop_duplicates("feature").set_index("feature", drop=False)


def _observed_vs_expected(expected: str, observed: str, fdr: object) -> str:
    if observed not in {"positive", "negative"}:
        return "not_tested"
    fdr_value = pd.to_numeric(fdr, errors="coerce")
    if pd.isna(fdr_value) or fdr_value >= 0.05:
        return "observed_not_fdr_significant"
    if expected not in {"positive", "negative"}:
        return "no_directional_prior"
    return "aligned" if expected == observed else "opposite"


def _empty_theme_row(theme: str) -> dict[str, object]:
    return {
        "primary_theme": theme,
        "primary_theme_label": THEME_LABELS[theme],
        "n_features": 0,
        "n_age_tested": 0,
        "n_gateway_age_fdr_lt_0_05": 0,
        "n_observed_positive": 0,
        "n_observed_negative": 0,
        "n_expected_positive": 0,
        "n_expected_negative": 0,
        "n_direction_aligned": 0,
        "n_direction_opposite": 0,
        "max_abs_importance": 0.0,
        "median_abs_age_beta_per_10_years": np.nan,
        "top_features": "",
        "theme_interpretation": THEME_INTERPRETATIONS[theme],
        "theme_citations": THEME_CITATIONS[theme],
    }


def _theme_colors(themes: list[str]) -> list[str]:
    palette = {
        "basal_quiescence": "#4c78a8",
        "basal_activation": "#f58518",
        "cycling_proliferation": "#e45756",
        "neural_progenitor": "#72b7b2",
        "immature_osn": "#54a24b",
        "mature_osn": "#b279a2",
        "sustentacular_barrier_detox": "#ff9da6",
        "bowman_gland_secretory": "#9d755d",
        "respiratory_metaplasia_ciliated_goblet": "#bab0ac",
        "immune_inflammatory": "#7f7f7f",
        "stress_senescence": "#d4a72c",
        "ecm_remodeling": "#5f9ed1",
        "technical_yield": "#8cd17d",
    }
    return [palette.get(theme, "#4c78a8") for theme in themes]


def _resource_columns() -> list[str]:
    return [
        "feature",
        "feature_kind",
        "feature_label",
        "feature_display",
        "primary_theme",
        "primary_theme_label",
        "secondary_theme",
        "expected_aging_direction",
        "expected_direction_basis",
        "confidence",
        "evidence_source",
        "evidence_citations",
        "biological_role",
        "interpretation_caution",
    ]


def _result_columns() -> list[str]:
    return [
        *_resource_columns(),
        "gateway_age_beta_per_10_years",
        "gateway_age_p_value",
        "gateway_age_fdr",
        "gateway_age_direction",
        "gateway_age_status",
        "observed_vs_expected",
        "max_abs_importance",
        "top_model",
        "max_selection_fraction",
        "mean_selection_fraction",
        "supporting_models",
        "specificity_class",
        "cross_tissue_confidence",
        "external_age_effect_status",
        "importance_rank",
    ]


def _summary_columns() -> list[str]:
    return [
        "primary_theme",
        "primary_theme_label",
        "n_features",
        "n_age_tested",
        "n_gateway_age_fdr_lt_0_05",
        "n_observed_positive",
        "n_observed_negative",
        "n_expected_positive",
        "n_expected_negative",
        "n_direction_aligned",
        "n_direction_opposite",
        "max_abs_importance",
        "median_abs_age_beta_per_10_years",
        "top_features",
        "theme_interpretation",
        "theme_citations",
    ]
