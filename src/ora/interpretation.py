"""Biological interpretation helpers for ORA feature summaries."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


DEFAULT_MODEL_ORDER = [
    "catboost",
    "xgboost",
    "hist_gradient_boosting",
    "boosted_ensemble",
    "random_forest",
    "ridge",
    "lasso",
    "elastic_net",
]


def build_feature_interpretation(
    feature_stability: pd.DataFrame,
    associations: pd.DataFrame | None = None,
    *,
    top_per_model: int = 12,
    top_n: int = 30,
    model_order: list[str] | tuple[str, ...] = DEFAULT_MODEL_ORDER,
) -> pd.DataFrame:
    """Map stable ORA features to broad biological themes for manuscript triage."""

    columns = _interpretation_columns()
    if feature_stability.empty or not {"model", "feature", "selection_fraction"}.issubset(feature_stability.columns):
        return pd.DataFrame(columns=columns)

    frame = feature_stability.copy()
    frame["selection_fraction"] = pd.to_numeric(frame["selection_fraction"], errors="coerce").fillna(0.0)
    importance_col = "abs_mean_importance" if "abs_mean_importance" in frame else "mean_importance"
    frame["abs_importance"] = pd.to_numeric(frame.get(importance_col), errors="coerce").abs().fillna(0.0)
    frame = frame[frame["selection_fraction"].gt(0) & frame["abs_importance"].gt(0)].copy()
    if frame.empty:
        return pd.DataFrame(columns=columns)

    selected = []
    model_rank = {model: idx for idx, model in enumerate(model_order)}
    for model, group in frame.groupby("model", observed=True):
        if str(model) == "null_model":
            continue
        ranked = group.sort_values(["selection_fraction", "abs_importance", "feature"], ascending=[False, False, True])
        selected.append(ranked.head(top_per_model))
    if not selected:
        return pd.DataFrame(columns=columns)

    top = pd.concat(selected, ignore_index=True)
    association_lookup = _association_lookup(associations)
    rows = []
    for feature, group in top.groupby("feature", observed=True):
        best = group.sort_values(["abs_importance", "selection_fraction"], ascending=[False, False]).iloc[0]
        models = sorted(group["model"].astype(str).unique().tolist(), key=lambda item: model_rank.get(item, 999))
        age = association_lookup.get(str(feature), {})
        theme = classify_feature_theme(str(feature))
        rows.append(
            {
                "feature": feature,
                "feature_kind": _feature_kind(str(feature)),
                "feature_label": _feature_label(str(feature)),
                "biology_theme": theme,
                "supporting_models": ",".join(models),
                "n_supporting_models": len(models),
                "top_model": str(best["model"]),
                "mean_selection_fraction": float(group["selection_fraction"].mean()),
                "max_abs_importance": float(group["abs_importance"].max()),
                "age_direction": age.get("direction", "not_tested"),
                "beta_per_10_years": age.get("beta_per_10_years", np.nan),
                "fdr": age.get("fdr", np.nan),
                "interpretation": _feature_interpretation(str(feature), theme, age.get("direction", "")),
                "caution": _feature_caution(str(feature), theme),
            }
        )
    result = pd.DataFrame(rows, columns=columns)
    result["_model_rank"] = result["top_model"].map(lambda item: model_rank.get(str(item), 999))
    result = result.sort_values(
        ["n_supporting_models", "mean_selection_fraction", "max_abs_importance", "_model_rank", "feature"],
        ascending=[False, False, False, True, True],
    )
    return result.drop(columns=["_model_rank"]).head(top_n).reset_index(drop=True)


def classify_feature_theme(feature: str) -> str:
    """Assign a broad manuscript theme from a feature name."""

    token = _feature_label(feature).lower()
    if any(term in token for term in ["ad_amyloid", "pd_synuclein"]):
        return "disease-linked module"
    if any(term in token for term in ["senescence", "sasp", "stressed", "dysfunctional"]):
        return "stress/senescence"
    if any(
        term in token
        for term in [
            "hbc",
            "basal",
            "suprabasal",
            "progenitor",
            "neuroblast",
            "cycling",
            "activation_injury",
        ]
    ):
        return "regenerative/progenitor epithelium"
    if any(term in token for term in ["neuron", "osn", "iosn", "mosn", "inp", "olfactory_transduction"]):
        return "neuronal lineage maturation"
    if any(
        term in token
        for term in [
            "antigen_presenting",
            "macrophage",
            "classical",
            "nonclassical",
            "cdc",
            "pdc",
            "dendritic",
            "cd56",
            "nk",
            "immune",
            "inflammation",
            "complement",
            "cytotoxic",
            "naive",
            "tcell",
            "cd4",
            "cd8",
            "plasma",
            "gcb",
        ]
    ):
        return "immune/inflammatory compartment"
    if any(
        term in token
        for term in [
            "bowman",
            "gland",
            "mucous",
            "serous",
            "goblet",
            "secretory",
            "sus",
            "sustentacular",
            "deuterosomal",
            "multiciliated",
            "ionocyte",
            "tuft",
            "club",
            "mv_type",
        ]
    ):
        return "supporting/secretory epithelium"
    if "ratio" in _feature_kind(feature):
        return "derived lineage ratio"
    return "unclassified feature"


def _association_lookup(associations: pd.DataFrame | None) -> dict[str, dict[str, Any]]:
    if associations is None or associations.empty or "feature" not in associations:
        return {}
    rows = {}
    for _, row in associations.iterrows():
        rows[str(row["feature"])] = {
            "direction": row.get("direction", "not_tested"),
            "beta_per_10_years": pd.to_numeric(row.get("beta_per_10_years"), errors="coerce"),
            "fdr": pd.to_numeric(row.get("fdr"), errors="coerce"),
        }
    return rows


def _feature_kind(feature: str) -> str:
    return feature.split("__", 1)[0] if "__" in feature else "feature"


def _feature_label(feature: str) -> str:
    return feature.split("__", 1)[1] if "__" in feature else feature


def _feature_interpretation(feature: str, theme: str, age_direction: str) -> str:
    label = _feature_label(feature).replace("_", " ")
    direction = str(age_direction)
    if direction in {"positive", "negative"}:
        trend = "increases with age" if direction == "positive" else "decreases with age"
        return f"{label} {trend} in donor-level association tests and maps to {theme}."
    if _feature_kind(feature) == "module_score":
        return f"{label} module contributes to ORA prediction and maps to {theme}."
    return f"{label} contributes to ORA prediction and maps to {theme}."


def _feature_caution(feature: str, theme: str) -> str:
    if theme == "disease-linked module":
        return "Exploratory module annotation; not evidence of disease prediction or causality."
    if theme == "immune/inflammatory compartment":
        return "Interpret with donor balance, cell yield, chemistry, and collection-method diagnostics."
    if _feature_kind(feature) == "ratio":
        return "Derived ratio; validate numerator and denominator cell-state stability separately."
    if _feature_kind(feature) == "module_score":
        return "Module score is supportive biology, not an independent external validation."
    return "Associational donor-level feature; do not interpret as measured lineage flux."


def _interpretation_columns() -> list[str]:
    return [
        "feature",
        "feature_kind",
        "feature_label",
        "biology_theme",
        "supporting_models",
        "n_supporting_models",
        "top_model",
        "mean_selection_fraction",
        "max_abs_importance",
        "age_direction",
        "beta_per_10_years",
        "fdr",
        "interpretation",
        "caution",
    ]
