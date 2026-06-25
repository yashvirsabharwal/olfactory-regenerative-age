#!/usr/bin/env python3
"""Run donor-level compositional age models for ORA."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.compositional import run_compositional_age_model
from ora.config import load_config
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", default="data/processed/donor_cell_state_counts.tsv")
    parser.add_argument("--manifest", default="data/processed/cohort_manifest.tsv")
    parser.add_argument("--age-associations", default="results/tables/age_cell_state_associations.tsv")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--sensitivity-out", default=None)
    parser.add_argument("--figure-out", default=None)
    parser.add_argument("--min-scenario-donors", type=int, default=30)
    parser.add_argument("--min-nonzero-donors", type=int, default=5)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    summary_out = args.summary_out or outputs.get(
        "compositional_age_model_summary_tsv",
        "results/tables/compositional_age_model_summary.tsv",
    )
    sensitivity_out = args.sensitivity_out or outputs.get(
        "compositional_age_model_sensitivity_tsv",
        "results/tables/compositional_age_model_sensitivity.tsv",
    )
    figure_out = args.figure_out or outputs.get(
        "extended_data_compositional_model_pdf",
        "results/figures/extended_data_compositional_model.pdf",
    )

    age_associations = None
    age_path = Path(args.age_associations)
    if age_path.exists():
        age_associations = pd.read_csv(age_path, sep="\t")

    result = run_compositional_age_model(
        pd.read_csv(args.counts, sep="\t"),
        pd.read_csv(args.manifest, sep="\t"),
        age_associations=age_associations,
        min_scenario_donors=args.min_scenario_donors,
        min_nonzero_donors=args.min_nonzero_donors,
    )
    result.summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    result.sensitivity.to_csv(ensure_parent(sensitivity_out), sep="\t", index=False)
    _write_figure(result.summary, Path(figure_out))

    ok = result.summary[result.summary["status"].eq("ok")] if "status" in result.summary.columns else result.summary
    supported = int(ok.get("supported_by_primary_and_sensitivity", pd.Series(dtype=bool)).fillna(False).sum())
    significant = int(ok.get("significant_q_0_10", pd.Series(dtype=bool)).fillna(False).sum())
    print(f"Wrote compositional summary: {summary_out} ({result.summary.shape[0]} cell states)")
    print(f"Wrote compositional sensitivity: {sensitivity_out} ({result.sensitivity.shape[0]} rows)")
    print(f"Wrote compositional figure: {figure_out}")
    print(f"Primary q<=0.10 cell states: {significant}; stable in sensitivity: {supported}")


def _write_figure(summary: pd.DataFrame, figure_out: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ensure_parent(figure_out)
    plot = summary[summary.get("status", "").eq("ok")].copy() if "status" in summary.columns else summary.copy()
    if plot.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No compositional model rows available", ha="center", va="center")
        ax.axis("off")
    else:
        plot["q_value"] = pd.to_numeric(plot["q_value"], errors="coerce")
        plot["abs_effect"] = pd.to_numeric(plot["age_beta_per_10_years"], errors="coerce").abs()
        plot = plot.sort_values(["q_value", "abs_effect"], ascending=[True, False]).head(18)
        plot = plot.iloc[::-1].copy()
        colors = np.where(
            pd.to_numeric(plot["q_value"], errors="coerce").le(0.10),
            "#006d77",
            "#9aa0a6",
        )
        fig, ax = plt.subplots(figsize=(8.8, 6.2), constrained_layout=True)
        y = np.arange(plot.shape[0])
        x = pd.to_numeric(plot["age_beta_per_10_years"], errors="coerce").to_numpy(dtype=float)
        xerr = pd.to_numeric(plot["standard_error_per_10_years"], errors="coerce").to_numpy(dtype=float)
        ax.barh(y, x, color=colors, edgecolor="#263238", linewidth=0.4)
        ax.errorbar(x, y, xerr=xerr, fmt="none", ecolor="#263238", elinewidth=0.8, capsize=2)
        ax.axvline(0, color="#263238", linewidth=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([_pretty_label(value) for value in plot["cell_state"]], fontsize=8)
        ax.set_xlabel("CLR age coefficient per 10 years")
        ax.set_title("Donor-level compositional age model")
        ax.grid(axis="x", color="#d8dee4", linewidth=0.6, alpha=0.8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(
            0.01,
            0.02,
            "Teal: primary q<=0.10; whiskers: 1 SE",
            transform=ax.transAxes,
            fontsize=8,
            color="#4b5563",
        )
    fig.savefig(figure_out)
    png_out = figure_out.with_suffix(".png")
    fig.savefig(png_out, dpi=220)
    plt.close(fig)


def _pretty_label(value: object) -> str:
    text = str(value).replace("_", " ")
    if len(text) <= 34:
        return text
    return text[:31] + "..."


if __name__ == "__main__":
    main()
