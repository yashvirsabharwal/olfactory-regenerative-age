#!/usr/bin/env python3
"""Run ORA negative-control and technical-confound analyses."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.negative_controls import run_negative_controls
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", default="data/processed/ora_augmented_feature_matrix.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--scores", default="results/tables/augmented_donor_ora_scores.tsv")
    parser.add_argument("--ora-summary", default="results/tables/ora_augmented_candidate_repeated_cv_summary.tsv")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--n-shuffles", type=int, default=50)
    parser.add_argument("--performance-out", default=None)
    parser.add_argument("--comparison-out", default=None)
    parser.add_argument("--explainability-out", default=None)
    parser.add_argument("--figure-out", default=None)
    args = parser.parse_args()

    model_config = load_config(args.model_config)
    outputs = load_config(args.gateway_config).get("outputs", {})
    performance_out = args.performance_out or outputs.get(
        "ora_negative_control_performance_tsv",
        "results/tables/ora_negative_control_performance.tsv",
    )
    comparison_out = args.comparison_out or outputs.get(
        "ora_technical_baseline_comparison_tsv",
        "results/tables/ora_technical_baseline_comparison.tsv",
    )
    explainability_out = args.explainability_out or outputs.get(
        "ora_covariate_explainability_tsv",
        "results/tables/ora_covariate_explainability.tsv",
    )
    figure_out = args.figure_out or outputs.get(
        "extended_data_negative_controls_pdf",
        "results/figures/extended_data_negative_controls.pdf",
    )

    result = run_negative_controls(
        pd.read_csv(args.features, sep="\t"),
        pd.read_csv(args.manifest, sep="\t"),
        pd.read_csv(args.scores, sep="\t"),
        model_config,
        n_shuffles=args.n_shuffles,
    )
    comparison = _add_ora_reference(result.baseline_comparison, Path(args.ora_summary))
    result.performance.to_csv(ensure_parent(performance_out), sep="\t", index=False)
    comparison.to_csv(ensure_parent(comparison_out), sep="\t", index=False)
    result.covariate_explainability.to_csv(ensure_parent(explainability_out), sep="\t", index=False)
    _write_figure(result.performance, comparison, Path(figure_out))

    print(f"Wrote negative-control performance: {performance_out} ({result.performance.shape[0]} rows)")
    print(f"Wrote technical baseline comparison: {comparison_out} ({comparison.shape[0]} rows)")
    print(f"Wrote covariate explainability: {explainability_out} ({result.covariate_explainability.shape[0]} rows)")
    print(f"Wrote negative-control figure: {figure_out}")


def _write_figure(performance: pd.DataFrame, comparison: pd.DataFrame, figure_out: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ensure_parent(figure_out)
    ok = performance[performance["status"].eq("ok")].copy()
    fig, ax = plt.subplots(figsize=(8.8, 5.3), constrained_layout=True)
    if ok.empty:
        ax.text(0.5, 0.5, "No negative-control rows available", ha="center", va="center")
        ax.axis("off")
    else:
        observed_controls = [
            ("biological_ridge_cv", "Biological"),
            ("technical_only_ridge_cv", "Technical"),
            ("yield_only_ridge_cv", "Yield"),
            ("null_mean_cv", "Null"),
        ]
        labels = []
        values = []
        colors = []
        for control, label in observed_controls:
            frame = ok[ok["control"].eq(control)]
            if frame.empty:
                continue
            labels.append(label)
            values.append(float(frame["mae"].iloc[0]))
            colors.append("#006d77" if control == "biological_ridge_cv" else "#9aa0a6")
        shuffle = ok[ok["control"].eq("age_shuffle_within_technical_strata")]
        if not shuffle.empty:
            labels.append("Shuffled")
            values.append(float(pd.to_numeric(shuffle["mae"], errors="coerce").mean()))
            colors.append("#c77d29")
        x = np.arange(len(labels))
        ax.bar(x, values, color=colors, edgecolor="#263238", linewidth=0.5)
        if not shuffle.empty and labels[-1] == "Shuffled":
            low = float(pd.to_numeric(shuffle["mae"], errors="coerce").quantile(0.025))
            high = float(pd.to_numeric(shuffle["mae"], errors="coerce").quantile(0.975))
            ax.errorbar(
                x[-1],
                values[-1],
                yerr=[[values[-1] - low], [high - values[-1]]],
                fmt="none",
                ecolor="#263238",
                capsize=3,
                linewidth=0.9,
            )
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=0)
        ax.set_ylabel("Cross-validated MAE (years)")
        ax.set_title("ORA negative controls and technical baselines")
        ax.grid(axis="y", color="#d8dee4", linewidth=0.6, alpha=0.8)
        ax.spines[["top", "right"]].set_visible(False)
        if not comparison.empty and "shuffle_empirical_p_mae_le_observed" in comparison.columns:
            shuffle_row = comparison[comparison["comparison"].eq("age_shuffle_within_technical_strata")]
            if not shuffle_row.empty:
                p_value = shuffle_row["shuffle_empirical_p_mae_le_observed"].iloc[0]
                ax.text(
                    0.01,
                    0.02,
                    f"Age-shuffle empirical p={p_value:.3f} for MAE <= biological ridge",
                    transform=ax.transAxes,
                    fontsize=8,
                    color="#4b5563",
                )
    fig.savefig(figure_out)
    fig.savefig(figure_out.with_suffix(".png"), dpi=220)
    plt.close(fig)


def _add_ora_reference(comparison: pd.DataFrame, summary_path: Path) -> pd.DataFrame:
    output = comparison.copy()
    output["best_observed_ora_model"] = ""
    output["best_observed_ora_mae"] = np.nan
    output["delta_mae_vs_best_observed_ora"] = np.nan
    output["interpretation_vs_best_observed_ora"] = ""
    if not summary_path.exists() or output.empty:
        return output
    summary = pd.read_csv(summary_path, sep="\t")
    mae_col = "mae_mean" if "mae_mean" in summary.columns else "mae"
    if mae_col not in summary.columns or "model" not in summary.columns:
        return output
    best = summary.assign(_mae=pd.to_numeric(summary[mae_col], errors="coerce")).sort_values("_mae").head(1)
    if best.empty or not np.isfinite(float(best["_mae"].iloc[0])):
        return output
    best_model = str(best["model"].iloc[0])
    best_mae = float(best["_mae"].iloc[0])
    output["best_observed_ora_model"] = best_model
    output["best_observed_ora_mae"] = best_mae
    output["delta_mae_vs_best_observed_ora"] = pd.to_numeric(output["mae"], errors="coerce") - best_mae
    output["interpretation_vs_best_observed_ora"] = np.where(
        pd.to_numeric(output["mae"], errors="coerce") <= best_mae,
        "control_matches_or_beats_best_observed_ora",
        "control_worse_than_best_observed_ora",
    )
    return output


if __name__ == "__main__":
    main()
