"""MVP report generation for ORA analysis outputs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ora.utils import ensure_parent


FIGURE_NAMES = {
    "cohort": "mvp_cohort_overview.png",
    "performance": "mvp_model_performance.png",
    "associations": "mvp_top_age_associations.png",
    "predictions": "mvp_predicted_vs_age.png",
    "importance": "mvp_feature_importance.png",
}


def generate_mvp_report(
    *,
    manifest: pd.DataFrame,
    cohort_summary: pd.DataFrame,
    associations: pd.DataFrame,
    performance: pd.DataFrame,
    scores: pd.DataFrame,
    importance: pd.DataFrame,
    out_md: str | Path,
    figure_dir: str | Path,
    source: dict[str, Any] | None = None,
    paper_defaults: dict[str, Any] | None = None,
    schema: dict[str, Any] | None = None,
    top_n: int = 12,
) -> list[Path]:
    """Write a markdown report and supporting PNG figures."""

    figure_path = Path(figure_dir)
    figure_path.mkdir(parents=True, exist_ok=True)
    figure_paths = _write_figures(
        cohort_summary=cohort_summary,
        associations=associations,
        performance=performance,
        scores=scores,
        importance=importance,
        figure_dir=figure_path,
        top_n=top_n,
    )
    markdown = render_mvp_markdown(
        manifest=manifest,
        cohort_summary=cohort_summary,
        associations=associations,
        performance=performance,
        importance=importance,
        out_md=out_md,
        figure_paths=figure_paths,
        source=source or {},
        paper_defaults=paper_defaults or {},
        schema=schema or {},
        top_n=top_n,
    )
    report_path = ensure_parent(out_md)
    report_path.write_text(markdown, encoding="utf-8")
    return [report_path, *figure_paths.values()]


def load_schema(path: str | Path | None) -> dict[str, Any]:
    if not path:
        return {}
    schema_path = Path(path)
    if not schema_path.exists():
        return {}
    return json.loads(schema_path.read_text(encoding="utf-8"))


def rank_associations(associations: pd.DataFrame, top_n: int = 12) -> pd.DataFrame:
    """Return the top age associations by FDR, then p-value."""

    required = {"feature", "status", "p_value", "fdr", "beta_per_10_years"}
    if associations.empty or not required.issubset(associations.columns):
        return pd.DataFrame(columns=list(required))
    frame = associations[associations["status"].eq("ok")].copy()
    frame["p_value"] = pd.to_numeric(frame["p_value"], errors="coerce")
    frame["fdr"] = pd.to_numeric(frame["fdr"], errors="coerce")
    frame["beta_per_10_years"] = pd.to_numeric(frame["beta_per_10_years"], errors="coerce")
    frame = frame[np.isfinite(frame["p_value"]) & np.isfinite(frame["fdr"])]
    return frame.sort_values(["fdr", "p_value", "feature"]).head(top_n).reset_index(drop=True)


def best_predictive_model(performance: pd.DataFrame) -> pd.Series:
    """Return the best non-null model row by MAE, falling back to the overall best row."""

    if performance.empty:
        return pd.Series(dtype=object)
    frame = performance.copy()
    frame["mae"] = pd.to_numeric(frame["mae"], errors="coerce")
    candidates = frame[~frame["model"].eq("null_model") & frame["mae"].notna()]
    if candidates.empty:
        candidates = frame[frame["mae"].notna()]
    if candidates.empty:
        return pd.Series(dtype=object)
    return candidates.sort_values(["mae", "model"]).iloc[0]


def render_mvp_markdown(
    *,
    manifest: pd.DataFrame,
    cohort_summary: pd.DataFrame,
    associations: pd.DataFrame,
    performance: pd.DataFrame,
    importance: pd.DataFrame,
    out_md: str | Path,
    figure_paths: dict[str, Path],
    source: dict[str, Any],
    paper_defaults: dict[str, Any],
    schema: dict[str, Any],
    top_n: int = 12,
) -> str:
    """Render the markdown body for the MVP report."""

    report_path = Path(out_md)
    top_assoc = rank_associations(associations, top_n=top_n)
    best_model = best_predictive_model(performance)
    healthy_train = _usable_training_donors(manifest)
    lines = [
        "# Gateway ORA MVP Report",
        "",
        "## Source",
        "",
        f"- Dataset: {source.get('name', 'Gateway human olfactory epithelium atlas')}",
        f"- DOI: {source.get('doi', 'not recorded')}",
        f"- CELLxGENE collection: {source.get('collection_id', 'not recorded')}",
        f"- Matrix inspected: {_format_int(schema.get('n_obs'))} cells x {_format_int(schema.get('n_vars'))} genes",
        f"- Paper default cells: {_format_int(paper_defaults.get('cells_total'))}",
        f"- Embeddings available in export: {_format_list(schema.get('obsm_keys', []))}",
        "",
        "## Cohort",
        "",
        f"- Usable healthy training donors: {_format_int(healthy_train)}",
        f"- Donors with any age value: {_format_int(_donors_with_age(manifest))}",
        "",
        _markdown_table(
            cohort_summary,
            ["cohort", "donors", "samples", "cells", "median_age", "missing_age_samples", "lineage_cells", "mature_neurons"],
            max_rows=20,
        ),
        "",
        _figure_link(report_path, figure_paths.get("cohort"), "Cohort overview"),
        "",
        "## ORA Age Prediction",
        "",
        _model_summary_sentence(best_model),
        "",
        _markdown_table(performance, ["model", "n", "mae", "rmse", "r2", "spearman_r"], max_rows=20),
        "",
        _figure_link(report_path, figure_paths.get("performance"), "Model performance"),
        "",
        _figure_link(report_path, figure_paths.get("predictions"), "Predicted age versus chronological age"),
        "",
        "## Age Associations",
        "",
        f"- Association tests with status ok: {_format_int(int(associations['status'].eq('ok').sum())) if 'status' in associations else 'not recorded'}",
        f"- Top table is ranked by FDR, then p-value; beta is per 10 years.",
        "",
        _markdown_table(
            top_assoc,
            ["feature", "n", "beta_per_10_years", "p_value", "fdr", "direction"],
            max_rows=top_n,
        ),
        "",
        _figure_link(report_path, figure_paths.get("associations"), "Top age associations"),
        "",
        "## Feature Importance",
        "",
        _markdown_table(_top_importance(importance, top_n=10), ["model", "feature", "importance", "stability"], max_rows=20),
        "",
        _figure_link(report_path, figure_paths.get("importance"), "Top model features"),
        "",
        "## Interpretation Notes",
        "",
        "- This is a composition-only MVP; gene programs, pseudobulk DE, trajectory density, Milo, and cNMF are deferred commands.",
        "- ORA training is restricted to healthy donors with valid age.",
        "- Chemistry, collection method, site, and yield are treated as covariates or sensitivity variables rather than biological ORA features.",
        "- AD/PD donors are excluded from ORA training and reserved for later frozen-model projection.",
        "",
    ]
    return "\n".join(line for line in lines if line is not None)


def _write_figures(
    *,
    cohort_summary: pd.DataFrame,
    associations: pd.DataFrame,
    performance: pd.DataFrame,
    scores: pd.DataFrame,
    importance: pd.DataFrame,
    figure_dir: Path,
    top_n: int,
) -> dict[str, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = {key: figure_dir / name for key, name in FIGURE_NAMES.items()}
    _plot_cohort(cohort_summary, paths["cohort"], plt)
    _plot_performance(performance, paths["performance"], plt)
    _plot_associations(rank_associations(associations, top_n=top_n), paths["associations"], plt)
    _plot_predictions(scores, paths["predictions"], plt)
    _plot_importance(importance, paths["importance"], plt)
    return paths


def _plot_cohort(cohort_summary: pd.DataFrame, path: Path, plt: Any) -> None:
    frame = cohort_summary.copy()
    if frame.empty or "cohort" not in frame:
        _blank_figure(path, plt, "No cohort summary available")
        return
    cohorts = frame["cohort"].astype(str).tolist()
    donors = pd.to_numeric(frame.get("donors"), errors="coerce").fillna(0)
    cells = pd.to_numeric(frame.get("cells"), errors="coerce").fillna(0) / 1_000_000
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    axes[0].bar(cohorts, donors, color="#287c8e")
    axes[0].set_title("Donors")
    axes[0].set_ylabel("Count")
    axes[1].bar(cohorts, cells, color="#c2674f")
    axes[1].set_title("Cells")
    axes[1].set_ylabel("Millions")
    for ax in axes:
        ax.tick_params(axis="x", rotation=25)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Gateway Cohort Overview")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_performance(performance: pd.DataFrame, path: Path, plt: Any) -> None:
    if performance.empty or "model" not in performance:
        _blank_figure(path, plt, "No model performance available")
        return
    frame = performance.copy()
    frame["mae"] = pd.to_numeric(frame.get("mae"), errors="coerce")
    frame["spearman_r"] = pd.to_numeric(frame.get("spearman_r"), errors="coerce")
    labels = frame["model"].astype(str).str.replace("_", " ").tolist()
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    axes[0].bar(labels, frame["mae"], color="#4d7fb8")
    axes[0].set_title("MAE")
    axes[0].set_ylabel("Years")
    axes[1].bar(labels, frame["spearman_r"], color="#7a9f45")
    axes[1].set_title("Spearman r")
    axes[1].axhline(0, color="#555555", linewidth=0.8)
    for ax in axes:
        ax.tick_params(axis="x", rotation=25)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("ORA Model Performance")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_associations(top_assoc: pd.DataFrame, path: Path, plt: Any) -> None:
    if top_assoc.empty:
        _blank_figure(path, plt, "No age associations available")
        return
    frame = top_assoc.sort_values("beta_per_10_years")
    labels = [str(item).replace("__", ": ") for item in frame["feature"]]
    beta = pd.to_numeric(frame["beta_per_10_years"], errors="coerce")
    colors = np.where(beta >= 0, "#287c8e", "#c2674f")
    fig, ax = plt.subplots(figsize=(8, max(4, 0.34 * len(frame))), constrained_layout=True)
    ax.barh(labels, beta, color=colors)
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_xlabel("Beta per 10 years")
    ax.set_title("Top Age Associations")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_predictions(scores: pd.DataFrame, path: Path, plt: Any) -> None:
    if scores.empty or not {"model", "chronological_age", "ora"}.issubset(scores.columns):
        _blank_figure(path, plt, "No ORA predictions available")
        return
    frame = scores[scores["model"].isin(["elastic_net", "random_forest"])].copy()
    if frame.empty:
        frame = scores.copy()
    models = list(frame["model"].drop_duplicates())
    fig, axes = plt.subplots(1, len(models), figsize=(4.4 * len(models), 4), squeeze=False, constrained_layout=True)
    for ax, model in zip(axes.ravel(), models):
        sub = frame[frame["model"].eq(model)]
        x = pd.to_numeric(sub["chronological_age"], errors="coerce")
        y = pd.to_numeric(sub["ora"], errors="coerce")
        ax.scatter(x, y, s=20, alpha=0.72, color="#4d7fb8", edgecolor="none")
        low = float(np.nanmin([x.min(), y.min()]))
        high = float(np.nanmax([x.max(), y.max()]))
        ax.plot([low, high], [low, high], color="#555555", linewidth=0.9, linestyle="--")
        ax.set_title(str(model).replace("_", " "))
        ax.set_xlabel("Chronological age")
        ax.set_ylabel("Predicted age")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Donor-Level ORA Predictions")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_importance(importance: pd.DataFrame, path: Path, plt: Any) -> None:
    top = _top_importance(importance, top_n=10)
    if top.empty:
        _blank_figure(path, plt, "No feature importance available")
        return
    models = list(top["model"].drop_duplicates())
    fig, axes = plt.subplots(1, len(models), figsize=(4.8 * len(models), 4.2), squeeze=False, constrained_layout=True)
    for ax, model in zip(axes.ravel(), models):
        sub = top[top["model"].eq(model)].copy()
        sub["abs_importance"] = pd.to_numeric(sub["importance"], errors="coerce").abs()
        sub = sub.sort_values("abs_importance")
        labels = [str(item).replace("__", ": ") for item in sub["feature"]]
        ax.barh(labels, sub["abs_importance"], color="#7a9f45")
        ax.set_title(str(model).replace("_", " "))
        ax.set_xlabel("Absolute importance")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Top ORA Model Features")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _top_importance(importance: pd.DataFrame, top_n: int = 10) -> pd.DataFrame:
    if importance.empty or not {"model", "feature", "importance"}.issubset(importance.columns):
        return pd.DataFrame(columns=["model", "feature", "importance", "stability"])
    frame = importance[importance["model"].isin(["elastic_net", "random_forest"])].copy()
    if frame.empty:
        frame = importance.copy()
    frame["abs_importance"] = pd.to_numeric(frame["importance"], errors="coerce").abs()
    frame = frame[np.isfinite(frame["abs_importance"])]
    pieces = []
    for model, sub in frame.groupby("model", sort=False):
        pieces.append(sub.sort_values("abs_importance", ascending=False).head(top_n))
    if not pieces:
        return pd.DataFrame(columns=["model", "feature", "importance", "stability"])
    return pd.concat(pieces, ignore_index=True).drop(columns=["abs_importance"], errors="ignore")


def _blank_figure(path: Path, plt: Any, message: str) -> None:
    fig, ax = plt.subplots(figsize=(5, 2.5), constrained_layout=True)
    ax.text(0.5, 0.5, message, ha="center", va="center")
    ax.set_axis_off()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _markdown_table(frame: pd.DataFrame, columns: list[str], max_rows: int) -> str:
    available = [col for col in columns if col in frame.columns]
    if frame.empty or not available:
        return "_No rows available._"
    table = frame[available].head(max_rows).copy()
    for col in table.columns:
        if pd.api.types.is_numeric_dtype(table[col]):
            table[col] = table[col].map(_format_table_value)
        else:
            table[col] = table[col].fillna("").astype(str)
    return _render_markdown_table(table)


def _render_markdown_table(table: pd.DataFrame) -> str:
    headers = [str(col) for col in table.columns]
    rows = [[_escape_table_cell(value) for value in row] for row in table.to_numpy()]
    header = "| " + " | ".join(_escape_table_cell(col) for col in headers) + " |"
    divider = "| " + " | ".join("---" for _ in headers) + " |"
    body = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([header, divider, *body])


def _escape_table_cell(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _figure_link(report_path: Path, figure_path: Path | None, alt: str) -> str:
    if figure_path is None:
        return ""
    try:
        rel = figure_path.relative_to(report_path.parent)
    except ValueError:
        rel = Path("../figures") / figure_path.name
    return f"![{alt}]({rel.as_posix()})"


def _model_summary_sentence(best_model: pd.Series) -> str:
    if best_model.empty:
        return "_No model performance rows available._"
    return (
        f"Best non-null model by MAE: **{best_model.get('model')}** "
        f"(MAE {_format_float(best_model.get('mae'))} years, "
        f"Spearman r {_format_float(best_model.get('spearman_r'))})."
    )


def _usable_training_donors(manifest: pd.DataFrame) -> int | None:
    if "donor_id" not in manifest or "usable_for_ora_training" not in manifest:
        return None
    frame = manifest[manifest["usable_for_ora_training"].astype(bool)]
    if "age" in frame:
        frame = frame[frame["age"].notna()]
    return int(frame["donor_id"].nunique())


def _donors_with_age(manifest: pd.DataFrame) -> int | None:
    if "donor_id" not in manifest or "age" not in manifest:
        return None
    return int(manifest[manifest["age"].notna()]["donor_id"].nunique())


def _format_table_value(value: object) -> str:
    if pd.isna(value):
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number.is_integer():
        return f"{int(number):,}"
    if abs(number) >= 1000 and number.is_integer():
        return f"{int(number):,}"
    if abs(number) >= 10:
        return f"{number:.2f}"
    if abs(number) >= 0.01:
        return f"{number:.4f}"
    return f"{number:.3g}"


def _format_float(value: object) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "not recorded"
    if not np.isfinite(number):
        return "not recorded"
    return f"{number:.2f}"


def _format_int(value: object) -> str:
    if value is None:
        return "not recorded"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "not recorded"
    if not np.isfinite(number):
        return "not recorded"
    return f"{int(number):,}"


def _format_list(values: object) -> str:
    if not values:
        return "none recorded"
    if isinstance(values, (list, tuple)):
        return ", ".join(str(value) for value in values)
    return str(values)
