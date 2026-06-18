#!/usr/bin/env python3
"""Build polished manuscript-oriented figures from ORA result tables."""

from __future__ import annotations

import argparse
from pathlib import Path
from textwrap import wrap

import numpy as np
import pandas as pd


BG = "#fbfaf7"
INK = "#24262b"
MUTED = "#6f7782"
GRID = "#d8d2c6"
TEAL = "#167a7f"
BLUE = "#4267ac"
VERMILION = "#c6533f"
GOLD = "#d99b2b"
GREEN = "#5e8f4d"
PURPLE = "#7d6ba6"
GRAY = "#90969f"

THEME_COLORS = {
    "supporting/secretory epithelium": TEAL,
    "immune/inflammatory compartment": VERMILION,
    "stress/senescence response": GOLD,
    "neuronal-lineage state": BLUE,
    "regenerative/progenitor compartment": GREEN,
    "disease/external module": PURPLE,
    "other": GRAY,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables-dir", default="results/tables")
    parser.add_argument("--figures-dir", default="results/figures")
    args = parser.parse_args()

    tables = Path(args.tables_dir)
    figures = Path(args.figures_dir)
    figures.mkdir(parents=True, exist_ok=True)

    _setup_matplotlib()
    import matplotlib.pyplot as plt

    written = [
        _figure1_design(tables, figures, plt),
        _figure2_age_composition(tables, figures, plt),
        _figure3_modeling(tables, figures, plt),
        _figure4_feature_biology(tables, figures, plt),
        _figure5_external_ndd(tables, figures, plt),
        _figure6_de_latent(tables, figures, plt),
    ]
    print("Wrote manuscript figures:")
    for path in written:
        print(f"- {path}")


def _setup_matplotlib() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.facecolor": BG,
            "axes.facecolor": BG,
            "savefig.facecolor": BG,
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.7,
            "grid.alpha": 0.65,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _figure1_design(tables: Path, figures: Path, plt) -> Path:
    cohort = _read_table(tables / "cohort_summary.tsv")
    fig = plt.figure(figsize=(10.8, 6.3), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, width_ratios=[1.1, 1.2, 1.1], height_ratios=[1, 1])
    ax_donors = fig.add_subplot(grid[0, 0])
    ax_cells = fig.add_subplot(grid[1, 0])
    ax_flow = fig.add_subplot(grid[:, 1])
    ax_gates = fig.add_subplot(grid[:, 2])

    _panel_label(ax_donors, "A")
    _panel_label(ax_flow, "B")
    _panel_label(ax_gates, "C")

    focus = cohort[cohort["cohort"].isin(["healthy", "ad", "pd"])].copy()
    focus["cohort"] = pd.Categorical(focus["cohort"], ["healthy", "ad", "pd"], ordered=True)
    focus = focus.sort_values("cohort")
    labels = ["Healthy", "AD", "PD"]
    colors = [TEAL, VERMILION, GOLD]
    ax_donors.bar(labels, pd.to_numeric(focus["donors"], errors="coerce"), color=colors, width=0.65)
    ax_donors.set_title("Donor groups")
    ax_donors.set_ylabel("Donors")
    _annotate_bars(ax_donors, fmt="{:.0f}")

    ax_cells.bar(labels, pd.to_numeric(focus["cells"], errors="coerce") / 1_000_000, color=colors, width=0.65)
    ax_cells.set_title("Cells represented")
    ax_cells.set_ylabel("Million cells")
    _annotate_bars(ax_cells, fmt="{:.2f}")

    ax_flow.set_axis_off()
    flow = [
        ("Gateway atlas", "4.03M cells\n202 donors", TEAL),
        ("Healthy training", "187 age-usable donors\nDonor-level CV", BLUE),
        ("ORA axis", "Composition, lineage ratios,\ncurated modules", GREEN),
        ("Guarded extensions", "External mapping, NDD projection,\nDE audits, scaled scVI", GOLD),
    ]
    y = 0.78
    for idx, (title, body, color) in enumerate(flow):
        box_y = y - idx * 0.24
        _rounded_box(ax_flow, 0.08, box_y, 0.84, 0.15, title, body, color)
        if idx < len(flow) - 1:
            start = box_y - 0.01
            end = y - (idx + 1) * 0.24 + 0.16
            ax_flow.annotate(
                "",
                xy=(0.5, end),
                xytext=(0.5, start),
                arrowprops={"arrowstyle": "-|>", "color": MUTED, "lw": 1.2},
                xycoords=ax_flow.transAxes,
            )

    ax_gates.set_axis_off()
    gate_rows = [
        ("Primary", "Healthy ORA axis", "supported", TEAL),
        ("External", "GSE184117", "small-n mapped", GOLD),
        ("Disease", "AD/PD ORAA", "exploratory", VERMILION),
        ("DE", "edgeR + limma", "audited", BLUE),
        ("Latent", "250k + lineage scVI", "seed-gated", PURPLE),
    ]
    ax_gates.text(0, 0.97, "Claim gates", fontsize=12, fontweight="bold", transform=ax_gates.transAxes)
    for idx, (scope, item, status, color) in enumerate(gate_rows):
        yy = 0.83 - idx * 0.16
        ax_gates.scatter([0.04], [yy + 0.02], s=95, color=color, transform=ax_gates.transAxes, clip_on=False)
        ax_gates.text(0.12, yy + 0.05, scope, fontweight="bold", transform=ax_gates.transAxes)
        ax_gates.text(0.12, yy, item, transform=ax_gates.transAxes)
        ax_gates.text(
            0.12,
            yy - 0.055,
            status,
            color=color,
            fontweight="bold",
            fontsize=8,
            transform=ax_gates.transAxes,
        )

    fig.suptitle("ORA study design separates supported aging signal from exploratory extensions", fontsize=14, fontweight="bold")
    return _save(fig, figures / "manuscript_figure1_design")


def _figure2_age_composition(tables: Path, figures: Path, plt) -> Path:
    assoc = _read_table(tables / "age_cell_state_associations.tsv")
    interp = _read_table(tables / "ora_feature_interpretation.tsv")
    merged = assoc.merge(
        interp[["feature", "biology_theme"]],
        on="feature",
        how="left",
    )
    merged["biology_theme"] = merged["biology_theme"].fillna("other")
    merged["beta_per_10_years"] = pd.to_numeric(merged["beta_per_10_years"], errors="coerce")
    merged["fdr"] = pd.to_numeric(merged["fdr"], errors="coerce")
    composition = merged[
        merged["feature"].astype(str).str.startswith(("clr__", "prop__"))
        & merged["beta_per_10_years"].notna()
    ].copy()
    negative = composition.nsmallest(7, "beta_per_10_years")
    positive = composition.nlargest(7, "beta_per_10_years")
    top = pd.concat([negative, positive], ignore_index=True).drop_duplicates("feature").sort_values("beta_per_10_years")

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.4), constrained_layout=True, gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    _panel_label(ax, "A")
    labels = [_feature_label(item) for item in top["feature"]]
    colors = [_theme_color(theme) for theme in top["biology_theme"]]
    ax.barh(labels, top["beta_per_10_years"], color=colors)
    ax.axvline(0, color=INK, lw=0.9)
    ax.set_xlabel("Age association beta per 10 years")
    ax.set_title("Largest donor-level age associations")
    ax.grid(axis="x", visible=True)
    ax.grid(axis="y", visible=False)

    ax = axes[1]
    _panel_label(ax, "B")
    theme = interp.copy()
    theme["mean_selection_fraction"] = pd.to_numeric(theme["mean_selection_fraction"], errors="coerce")
    theme["max_abs_importance"] = pd.to_numeric(theme["max_abs_importance"], errors="coerce")
    summary = (
        theme.groupby("biology_theme", dropna=False)
        .agg(
            features=("feature", "count"),
            mean_selection=("mean_selection_fraction", "mean"),
            importance=("max_abs_importance", "mean"),
        )
        .reset_index()
        .sort_values("features")
    )
    y = np.arange(summary.shape[0])
    ax.barh(y, summary["features"], color=[_theme_color(t) for t in summary["biology_theme"]])
    ax.set_yticks(y, [_wrap_label(t, 23) for t in summary["biology_theme"]])
    ax.set_xlabel("Stable features")
    ax.set_title("Feature themes")
    for yi, row in enumerate(summary.itertuples()):
        ax.text(row.features + 0.15, yi, f"{row.features}", va="center", fontsize=8)

    fig.suptitle("Healthy olfactory aging is carried by interpretable epithelial, immune, and neuronal-state features", fontsize=14, fontweight="bold")
    return _save(fig, figures / "manuscript_figure2_age_composition")


def _figure3_modeling(tables: Path, figures: Path, plt) -> Path:
    repeated = _read_table(tables / "ora_repeated_cv_summary.tsv")
    augmented = _read_table(tables / "ora_augmented_candidate_repeated_cv_summary.tsv")
    permutation = _read_table(tables / "ora_permutation_empirical.tsv")
    calibration = _read_table(tables / "ora_calibration.tsv")

    fig = plt.figure(figsize=(11.6, 6.6), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)
    ax_perf = fig.add_subplot(grid[0, 0])
    ax_null = fig.add_subplot(grid[0, 1])
    ax_spear = fig.add_subplot(grid[1, 0])
    ax_cal = fig.add_subplot(grid[1, 1])

    _panel_label(ax_perf, "A")
    display = ["null_model", "ridge", "random_forest", "xgboost", "catboost"]
    perf = repeated[repeated["model"].isin(display)].copy()
    perf["model"] = pd.Categorical(perf["model"], display, ordered=True)
    perf = perf.sort_values("model")
    x = np.arange(perf.shape[0])
    ax_perf.errorbar(
        x,
        pd.to_numeric(perf["mae_mean"], errors="coerce"),
        yerr=[
            pd.to_numeric(perf["mae_mean"], errors="coerce") - pd.to_numeric(perf["mae_ci_low"], errors="coerce"),
            pd.to_numeric(perf["mae_ci_high"], errors="coerce") - pd.to_numeric(perf["mae_mean"], errors="coerce"),
        ],
        fmt="o",
        color=BLUE,
        ecolor=BLUE,
        capsize=3,
        ms=6,
    )
    ax_perf.set_xticks(x, [_model_label(m) for m in perf["model"]], rotation=25, ha="right")
    ax_perf.set_ylabel("MAE (years)")
    ax_perf.set_title("Repeated donor-level CV")

    _panel_label(ax_null, "B")
    perm = permutation[permutation["model"].isin(["catboost", "xgboost", "boosted_ensemble", "random_forest"])].copy()
    perm = perm.sort_values("observed_mae")
    y = np.arange(perm.shape[0])
    ax_null.barh(y, pd.to_numeric(perm["null_mae_mean"], errors="coerce"), color="#ded8cd", label="Shuffled-age null")
    ax_null.scatter(pd.to_numeric(perm["observed_mae"], errors="coerce"), y, color=VERMILION, s=48, label="Observed")
    ax_null.set_yticks(y, [_model_label(m) for m in perm["model"]])
    ax_null.set_xlabel("MAE (years)")
    ax_null.set_title("Observed models beat shuffled labels")
    ax_null.legend(loc="lower left")

    _panel_label(ax_spear, "C")
    aug = augmented.copy()
    aug = aug[aug["model"].isin(["hist_gradient_boosting", "xgboost", "catboost", "boosted_ensemble"])]
    aug = aug.sort_values("spearman_r_mean")
    y = np.arange(aug.shape[0])
    ax_spear.barh(y, pd.to_numeric(aug["spearman_r_mean"], errors="coerce"), color=GREEN)
    ax_spear.set_yticks(y, [_model_label(m) for m in aug["model"]])
    ax_spear.set_xlabel("Mean Spearman r")
    ax_spear.set_title("Candidate module-augmented models")

    _panel_label(ax_cal, "D")
    cal = calibration[calibration["model"].isin(["ridge", "random_forest", "xgboost", "catboost", "boosted_ensemble"])].copy()
    cal = cal.sort_values("calibration_slope_ora_on_age")
    y = np.arange(cal.shape[0])
    ax_cal.barh(y, pd.to_numeric(cal["calibration_slope_ora_on_age"], errors="coerce"), color=GOLD)
    ax_cal.axvline(1, color=INK, lw=0.9, ls="--", label="ideal")
    ax_cal.set_yticks(y, [_model_label(m) for m in cal["model"]])
    ax_cal.set_xlabel("ORA-on-age calibration slope")
    ax_cal.set_title("Predictions are under-dispersed")
    ax_cal.legend(loc="lower right")

    fig.suptitle("ORA is a reproducible relative aging axis, not an exact chronological-age clock", fontsize=14, fontweight="bold")
    return _save(fig, figures / "manuscript_figure3_modeling")


def _figure4_feature_biology(tables: Path, figures: Path, plt) -> Path:
    interp = _read_table(tables / "ora_feature_interpretation.tsv")
    deltas = _read_table(tables / "ora_feature_set_model_deltas.tsv")
    module_cov = _read_table(tables / "module_gene_coverage.tsv")

    fig = plt.figure(figsize=(11.4, 6.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.3, 1], height_ratios=[1, 1])
    ax_bubble = fig.add_subplot(grid[:, 0])
    ax_delta = fig.add_subplot(grid[0, 1])
    ax_cov = fig.add_subplot(grid[1, 1])

    _panel_label(ax_bubble, "A")
    work = interp.copy()
    work["max_abs_importance"] = pd.to_numeric(work["max_abs_importance"], errors="coerce")
    work["mean_selection_fraction"] = pd.to_numeric(work["mean_selection_fraction"], errors="coerce")
    work["beta_per_10_years"] = pd.to_numeric(work["beta_per_10_years"], errors="coerce")
    work = work.sort_values("max_abs_importance", ascending=False).head(18)
    y = np.arange(work.shape[0])
    size = 45 + 185 * work["mean_selection_fraction"].fillna(0)
    ax_bubble.scatter(
        work["beta_per_10_years"],
        y,
        s=size,
        color=[_theme_color(t) for t in work["biology_theme"]],
        alpha=0.9,
        edgecolor=BG,
        linewidth=0.8,
    )
    ax_bubble.axvline(0, color=INK, lw=0.9)
    ax_bubble.set_yticks(y, [_feature_label(f) for f in work["feature"]])
    ax_bubble.invert_yaxis()
    ax_bubble.set_xlabel("Age beta per 10 years")
    ax_bubble.set_title("Stable ORA features remain biologically legible")

    _panel_label(ax_delta, "B")
    if not deltas.empty and {"model", "delta_mae_mean"}.issubset(deltas.columns):
        d = deltas[~deltas["model"].eq("null_model")].copy()
        d["delta_mae_mean"] = pd.to_numeric(d["delta_mae_mean"], errors="coerce")
        d = d.sort_values("delta_mae_mean").head(8).sort_values("delta_mae_mean")
        colors = np.where(d["delta_mae_mean"] <= 0, TEAL, VERMILION)
        ax_delta.barh([_model_label(m) for m in d["model"]], d["delta_mae_mean"], color=colors)
        ax_delta.axvline(0, color=INK, lw=0.9)
        ax_delta.set_xlabel("Module MAE delta (years)")
    else:
        ax_delta.text(0.5, 0.5, "No module delta table", ha="center", va="center")
    ax_delta.set_title("Modules add modest predictive gain")

    _panel_label(ax_cov, "C")
    cov = module_cov.copy()
    if not cov.empty and {"module", "coverage_fraction"}.issubset(cov.columns):
        cov["coverage_fraction"] = pd.to_numeric(cov["coverage_fraction"], errors="coerce")
        cov = cov.sort_values("coverage_fraction").tail(10)
        ax_cov.barh([_wrap_label(str(m).replace("_", " "), 22) for m in cov["module"]], cov["coverage_fraction"], color=BLUE)
        ax_cov.set_xlim(0, 1.05)
        ax_cov.set_xlabel("Gene coverage")
    else:
        ax_cov.text(0.5, 0.5, "No module coverage table", ha="center", va="center")
    ax_cov.set_title("Curated module coverage")

    _theme_legend(fig, sorted(set(work["biology_theme"])))
    fig.suptitle("Biological interpretation comes from stable, audited features rather than black-box accuracy", fontsize=14, fontweight="bold")
    return _save(fig, figures / "manuscript_figure4_feature_biology")


def _figure5_external_ndd(tables: Path, figures: Path, plt) -> Path:
    concordance = _read_table(tables / "external_marker_age_concordance.tsv")
    ndd = _read_table(tables / "ndd_ora_projection_feature_comparison.tsv")
    perm = _read_table(tables / "ndd_label_permutation.tsv")
    evidence = _read_table(tables / "external_validation_evidence.tsv")

    fig = plt.figure(figsize=(11.6, 6.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)
    ax_conc = fig.add_subplot(grid[0, 0])
    ax_evid = fig.add_subplot(grid[1, 0])
    ax_ndd = fig.add_subplot(grid[0, 1])
    ax_perm = fig.add_subplot(grid[1, 1])

    _panel_label(ax_conc, "A")
    counts = concordance["concordance"].value_counts().reindex(["concordant", "discordant", "not_evaluable"]).fillna(0)
    ax_conc.bar(counts.index.str.replace("_", " "), counts.values, color=[TEAL, VERMILION, GRAY], width=0.62)
    ax_conc.set_ylabel("Mapped marker-feature rows")
    ax_conc.set_title("GSE184117 marker-age concordance")
    _annotate_bars(ax_conc, fmt="{:.0f}")

    _panel_label(ax_evid, "B")
    status_counts = evidence["validation_strength"].fillna(evidence["status"]).astype(str).value_counts().head(6)
    y = np.arange(status_counts.shape[0])
    ax_evid.barh(y, status_counts.values, color=GOLD)
    ax_evid.set_yticks(y, [_wrap_label(str(s).replace("_", " "), 30) for s in status_counts.index])
    ax_evid.set_xlabel("Evidence entries")
    ax_evid.set_title("Validation evidence stays claim-gated")

    _panel_label(ax_ndd, "C")
    models = ["xgboost", "catboost", "boosted_ensemble", "random_forest"]
    sub = ndd[ndd["model"].isin(models)].copy()
    sub["model"] = pd.Categorical(sub["model"], models, ordered=True)
    sub["augmented_mean_oraa"] = pd.to_numeric(sub["augmented_mean_oraa"], errors="coerce")
    pivot = sub.pivot_table(index="model", columns="disease_group", values="augmented_mean_oraa", aggfunc="mean", observed=False)
    x = np.arange(len(pivot.index))
    width = 0.34
    ax_ndd.bar(x - width / 2, pivot.get("ad", pd.Series(index=pivot.index, dtype=float)), width, color=VERMILION, label="AD")
    ax_ndd.bar(x + width / 2, pivot.get("pd", pd.Series(index=pivot.index, dtype=float)), width, color=GOLD, label="PD")
    ax_ndd.axhline(0, color=INK, lw=0.9)
    ax_ndd.set_xticks(x, [_model_label(m) for m in pivot.index], rotation=25, ha="right")
    ax_ndd.set_ylabel("Mean ORA acceleration")
    ax_ndd.set_title("Frozen NDD projections remain negative")
    ax_ndd.legend()

    _panel_label(ax_perm, "D")
    p = perm[perm["model"].isin(["xgboost", "catboost"])].copy()
    p["observed_difference_vs_reference"] = pd.to_numeric(p["observed_difference_vs_reference"], errors="coerce")
    p["null_ci_low"] = pd.to_numeric(p["null_ci_low"], errors="coerce")
    p["null_ci_high"] = pd.to_numeric(p["null_ci_high"], errors="coerce")
    p["label"] = p["model"].map(_model_label) + " " + p["disease_group"].str.upper()
    p = p.sort_values("observed_difference_vs_reference")
    y = np.arange(p.shape[0])
    ax_perm.hlines(y, p["null_ci_low"], p["null_ci_high"], color=GRAY, lw=5, alpha=0.7, label="permutation 95% interval")
    ax_perm.scatter(p["observed_difference_vs_reference"], y, color=BLUE, s=45, label="observed")
    ax_perm.axvline(0, color=INK, lw=0.9)
    ax_perm.set_yticks(y, p["label"])
    ax_perm.set_xlabel("Difference vs matched reference")
    ax_perm.set_title("Permutation frame for NDD projection")
    ax_perm.legend(loc="lower center", bbox_to_anchor=(0.5, -0.28), ncols=2)

    fig.suptitle("External and disease analyses support context, not primary disease claims", fontsize=14, fontweight="bold")
    return _save(fig, figures / "manuscript_figure5_external_ndd")


def _figure6_de_latent(tables: Path, figures: Path, plt) -> Path:
    edger = _read_table(tables / "pseudobulk_genomewide_de_audit.tsv")
    limma = _read_table(tables / "pseudobulk_genomewide_limma_voom_de_audit.tsv")
    matched_edger = _read_table(tables / "pseudobulk_genomewide_de_audit_matched_flex_v2_device.tsv")
    matched_limma = _read_table(tables / "pseudobulk_genomewide_limma_voom_de_audit_matched_flex_v2_device.tsv")
    scvi = _read_table(tables / "scvi_scaled_250k_validation.tsv")
    if scvi.empty:
        scvi = _read_table(tables / "scvi_pilot_validation.tsv")

    fig = plt.figure(figsize=(11.4, 6.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)
    ax_sig = fig.add_subplot(grid[0, 0])
    ax_flags = fig.add_subplot(grid[1, 0])
    ax_scvi = fig.add_subplot(grid[:, 1])

    _panel_label(ax_sig, "A")
    de_summary = _de_summary(edger, limma, matched_edger, matched_limma)
    pivot = de_summary.pivot_table(index="contrast", columns="engine", values="significant_rows", aggfunc="sum")
    x = np.arange(pivot.shape[0])
    engines = list(pivot.columns)
    width = 0.18
    palette = [VERMILION, BLUE, GOLD, TEAL]
    for idx, engine in enumerate(engines):
        ax_sig.bar(x + (idx - (len(engines) - 1) / 2) * width, pivot[engine], width, label=engine, color=palette[idx % len(palette)])
    ax_sig.set_xticks(x, [c.replace("_vs_healthy", "").upper() for c in pivot.index])
    ax_sig.set_ylabel("FDR-significant rows")
    ax_sig.set_title("DE significance is method- and match-sensitive")
    ax_sig.legend(ncols=2, fontsize=7)

    _panel_label(ax_flags, "B")
    flags = edger[["contrast", "is_sex_linked_initial_significant_rows", "is_hemoglobin_significant_rows", "is_immunoglobulin_significant_rows"]].copy()
    flags = flags.set_index("contrast")
    flag_labels = ["sex-linked", "hemoglobin", "immunoglobulin"]
    vals = flags.to_numpy(dtype=float)
    image = ax_flags.imshow(vals, aspect="auto", cmap="YlOrRd")
    ax_flags.set_xticks(np.arange(len(flag_labels)), flag_labels, rotation=20, ha="right")
    ax_flags.set_yticks(np.arange(flags.shape[0]), [c.replace("_vs_healthy", "").upper() for c in flags.index])
    ax_flags.set_title("Sentinel categories are explicitly audited")
    for (i, j), value in np.ndenumerate(vals):
        ax_flags.text(j, i, f"{int(value)}", ha="center", va="center", color=INK, fontsize=8)
    fig.colorbar(image, ax=ax_flags, shrink=0.75, label="Rows")

    _panel_label(ax_scvi, "C")
    ax_scvi.set_axis_off()
    ax_scvi.text(0, 0.98, "Latent-space readiness", fontsize=12, fontweight="bold", transform=ax_scvi.transAxes)
    rows = [
        ("CELLxGENE export", "X_umap only", GRAY),
        ("Scaled scVI", _scvi_detail(scvi, "pilot_h5ad"), TEAL),
        ("Embedding", _scvi_detail(scvi, "embedding_dimensions"), BLUE),
        ("Finite values", _scvi_detail(scvi, "embedding_finite"), GREEN),
        ("Donor/sample balance", _scvi_status(scvi, "metadata__sample_id"), GOLD),
        ("Trajectory/Milo/cNMF", _scvi_status(scvi, "claim_gate"), VERMILION),
    ]
    for idx, (label, value, color) in enumerate(rows):
        yy = 0.84 - idx * 0.13
        ax_scvi.scatter([0.03], [yy + 0.025], s=85, color=color, transform=ax_scvi.transAxes, clip_on=False)
        ax_scvi.text(0.09, yy + 0.045, label, fontweight="bold", transform=ax_scvi.transAxes)
        ax_scvi.text(0.09, yy - 0.005, _wrap_label(value, 42), transform=ax_scvi.transAxes, fontsize=8)
    ax_scvi.text(
        0,
        0.04,
        "Interpretation: scVI is now technically feasible, but biological latent claims remain gated.",
        color=MUTED,
        fontsize=8,
        transform=ax_scvi.transAxes,
    )

    fig.suptitle("Genome-wide DE and latent-space analyses are presented with audit gates", fontsize=14, fontweight="bold")
    return _save(fig, figures / "manuscript_figure6_de_latent")


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def _save(fig, stem: Path) -> Path:
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    import matplotlib.pyplot as plt

    plt.close(fig)
    return stem.with_suffix(".png")


def _panel_label(ax, label: str) -> None:
    ax.text(-0.08, 1.06, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top", ha="left")


def _rounded_box(ax, x: float, y: float, w: float, h: float, title: str, body: str, color: str) -> None:
    from matplotlib.patches import FancyBboxPatch

    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        transform=ax.transAxes,
        fc="#ffffff",
        ec=color,
        lw=1.8,
    )
    ax.add_patch(patch)
    ax.text(x + 0.035, y + h - 0.042, title, transform=ax.transAxes, fontweight="bold", color=color)
    ax.text(x + 0.035, y + 0.035, body, transform=ax.transAxes, fontsize=8, va="bottom")


def _annotate_bars(ax, fmt: str = "{:.1f}") -> None:
    for patch in ax.patches:
        height = patch.get_height()
        if not np.isfinite(height):
            continue
        ax.text(
            patch.get_x() + patch.get_width() / 2,
            height,
            fmt.format(height),
            ha="center",
            va="bottom",
            fontsize=8,
            color=INK,
        )


def _feature_label(feature: str) -> str:
    text = str(feature).replace("prop__", "prop: ").replace("clr__", "clr: ").replace("ratio__", "ratio: ")
    text = text.replace("_", " ")
    return _wrap_label(text, 30)


def _model_label(model: str) -> str:
    return str(model).replace("_", " ").title().replace("Xgboost", "XGBoost").replace("Ndd", "NDD")


def _wrap_label(value: str, width: int) -> str:
    return "\n".join(wrap(str(value), width=width, break_long_words=False))


def _theme_color(theme: str) -> str:
    return THEME_COLORS.get(str(theme), THEME_COLORS["other"])


def _theme_legend(fig, themes: list[str]) -> None:
    import matplotlib.lines as mlines

    handles = [
        mlines.Line2D([], [], color=_theme_color(theme), marker="o", linestyle="None", markersize=7, label=theme)
        for theme in themes
        if theme
    ]
    if handles:
        fig.legend(handles=handles, loc="lower center", ncols=3, fontsize=7, bbox_to_anchor=(0.5, -0.08))


def _de_summary(edger: pd.DataFrame, limma: pd.DataFrame, matched_edger: pd.DataFrame, matched_limma: pd.DataFrame) -> pd.DataFrame:
    pieces = []
    for label, frame in [
        ("edgeR all", edger),
        ("limma all", limma),
        ("edgeR matched", matched_edger),
        ("limma matched", matched_limma),
    ]:
        if frame.empty:
            continue
        sub = frame[["contrast", "significant_rows"]].copy()
        sub["engine"] = label
        pieces.append(sub)
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=["contrast", "significant_rows", "engine"])


def _scvi_detail(scvi: pd.DataFrame, check: str) -> str:
    if scvi.empty or "check" not in scvi:
        return "not available"
    sub = scvi[scvi["check"].eq(check)]
    if sub.empty:
        return "not available"
    return str(sub["detail"].iloc[0])


def _scvi_status(scvi: pd.DataFrame, check: str) -> str:
    if scvi.empty or "check" not in scvi:
        return "not available"
    sub = scvi[scvi["check"].eq(check)]
    if sub.empty:
        return "not available"
    return f"{sub['status'].iloc[0]}: {sub['detail'].iloc[0]}"


if __name__ == "__main__":
    main()
