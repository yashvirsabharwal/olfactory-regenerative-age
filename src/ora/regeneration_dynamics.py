"""Feasibility and pilot analyses for regeneration-dynamics claims."""

from __future__ import annotations

from importlib import util
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

from .utils import ensure_parent


LINEAGE_FINE_ORDER = {
    "Quiescent_HBC": 0,
    "Cycling_HBC": 1,
    "Suprabasal": 2,
    "Early_INP": 3,
    "Late_INP": 4,
    "Early_iOSN": 5,
    "Late_iOSN": 6,
    "Early_mature_mOSN": 7,
    "Fully_mature_mOSN": 8,
    "Stressed_mOSN": 8,
}

LINEAGE_COARSE_ORDER = {
    "Resp_HBC": 0,
    "Olf_INPs": 3,
    "Olf_iOSNs": 5,
    "Olf_mOSNs": 8,
}


def audit_h5ad_dynamics_inputs(paths: list[str | Path]) -> pd.DataFrame:
    """Inspect local H5AD files for fate/velocity method prerequisites."""

    rows = []
    for path_like in paths:
        path = Path(path_like)
        if not path.exists():
            rows.append(
                {
                    "h5ad_path": str(path),
                    "exists": False,
                    "n_obs": np.nan,
                    "n_vars": np.nan,
                    "layers": "",
                    "obsm": "",
                    "has_raw": False,
                    "has_spliced": False,
                    "has_unspliced": False,
                    "has_velocity": False,
                    "has_x_scvi": False,
                    "has_x_umap": False,
                    "has_lineage_labels": False,
                    "has_batch_metadata": False,
                    "notes": "missing file",
                }
            )
            continue

        with h5py.File(path, "r") as handle:
            layers = sorted(str(key) for key in handle.get("layers", {}).keys())
            obsm = sorted(str(key) for key in handle.get("obsm", {}).keys())
            obs_keys = set(str(key) for key in handle.get("obs", {}).keys())
            n_obs, n_vars = _h5ad_shape(handle)
            has_raw = "raw" in handle
            has_spliced = "spliced" in layers
            has_unspliced = "unspliced" in layers
            has_velocity = any(key in layers for key in ("velocity", "Ms", "Mu", "ambiguous"))
            has_lineage_labels = bool({"fine_celltype", "coarse_celltype"} & obs_keys)
            has_batch_metadata = bool({"sample_id", "donor_id", "flex_version", "device_guided"} <= obs_keys)
            rows.append(
                {
                    "h5ad_path": str(path),
                    "exists": True,
                    "n_obs": n_obs,
                    "n_vars": n_vars,
                    "layers": ",".join(layers),
                    "obsm": ",".join(obsm),
                    "has_raw": has_raw,
                    "has_spliced": has_spliced,
                    "has_unspliced": has_unspliced,
                    "has_velocity": has_velocity,
                    "has_x_scvi": "X_scvi" in obsm,
                    "has_x_umap": "X_umap" in obsm,
                    "has_lineage_labels": has_lineage_labels,
                    "has_batch_metadata": has_batch_metadata,
                    "notes": _input_note(
                        has_spliced=has_spliced,
                        has_unspliced=has_unspliced,
                        has_lineage_labels=has_lineage_labels,
                        has_x_scvi="X_scvi" in obsm,
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_dynamics_feasibility(audit: pd.DataFrame) -> pd.DataFrame:
    """Classify each candidate dynamics method from observed input support."""

    has_velocity_input = bool((audit["has_spliced"] & audit["has_unspliced"]).any())
    has_lineage_scvi = bool((audit["has_lineage_labels"] & audit["has_x_scvi"]).any())
    has_lineage_expression = bool(audit["has_lineage_labels"].any())
    palantir_installed = util.find_spec("palantir") is not None
    cellrank_installed = util.find_spec("cellrank") is not None
    scvelo_installed = util.find_spec("scvelo") is not None

    rows = [
        {
            "method": "RNA velocity / scVelo",
            "status": "no_go",
            "installed": scvelo_installed,
            "input_support": has_velocity_input,
            "claim_status": "blocked",
            "rationale": "No inspected H5AD contains both spliced and unspliced layers.",
        },
        {
            "method": "CellRank velocity kernel",
            "status": "no_go",
            "installed": cellrank_installed,
            "input_support": has_velocity_input,
            "claim_status": "blocked",
            "rationale": "CellRank velocity kernels require velocity-compatible spliced/unspliced inputs.",
        },
        {
            "method": "Palantir",
            "status": "install_needed" if not palantir_installed else "feasible_exploratory",
            "installed": palantir_installed,
            "input_support": has_lineage_expression,
            "claim_status": "exploratory",
            "rationale": "Lineage labels and expression matrices exist, but Palantir is not a project dependency.",
        },
        {
            "method": "Scanpy diffusion pseudotime",
            "status": "feasible_exploratory" if has_lineage_scvi else "blocked",
            "installed": util.find_spec("scanpy") is not None,
            "input_support": has_lineage_scvi,
            "claim_status": "exploratory",
            "rationale": "Existing lineage scVI embedding can support a topology sanity check rooted in HBC states.",
        },
    ]
    return pd.DataFrame(rows)


def run_scanpy_dpt_pilot(
    *,
    h5ad_path: str | Path,
    max_cells: int = 40000,
    seed: int = 17,
    n_neighbors: int = 30,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run a guarded diffusion-pseudotime pilot on basal-to-neuronal lineage cells."""

    try:
        import scanpy as sc
    except ModuleNotFoundError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError("scanpy is required for the DPT pilot") from exc

    adata = sc.read_h5ad(h5ad_path)
    if "fine_celltype" not in adata.obs:
        raise KeyError("fine_celltype column is required for regeneration-dynamics pilot")
    if "X_scvi" not in adata.obsm:
        raise KeyError("X_scvi embedding is required for regeneration-dynamics pilot")

    fine = adata.obs["fine_celltype"].astype(str)
    keep = fine.isin(LINEAGE_FINE_ORDER)
    adata = adata[keep].copy()
    if adata.n_obs == 0:
        raise ValueError("No recognized basal-to-neuronal lineage cells found")
    if adata.n_obs > max_cells:
        adata = _stratified_obs_sample(adata, "fine_celltype", max_cells=max_cells, seed=seed)

    adata.obs["expected_lineage_rank"] = (
        adata.obs["fine_celltype"].astype(str).map(LINEAGE_FINE_ORDER).astype(float)
    )
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_scvi", random_state=seed)
    sc.tl.diffmap(adata)
    root_idx = _choose_root_index(adata)
    adata.uns["iroot"] = int(root_idx)
    sc.tl.dpt(adata, n_dcs=min(10, adata.obsm["X_diffmap"].shape[1] - 1))

    cell_scores = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "fine_celltype": adata.obs["fine_celltype"].astype(str).to_numpy(),
            "coarse_celltype": adata.obs.get("coarse_celltype", pd.Series(index=adata.obs_names, dtype=str))
            .astype(str)
            .to_numpy(),
            "donor_id": adata.obs.get("donor_id", pd.Series(index=adata.obs_names, dtype=str))
            .astype(str)
            .to_numpy(),
            "flex_version": adata.obs.get("flex_version", pd.Series(index=adata.obs_names, dtype=str))
            .astype(str)
            .to_numpy(),
            "device_guided": adata.obs.get("device_guided", pd.Series(index=adata.obs_names, dtype=str))
            .astype(str)
            .to_numpy(),
            "condition": adata.obs.get("condition", pd.Series(index=adata.obs_names, dtype=str))
            .astype(str)
            .to_numpy(),
            "expected_lineage_rank": adata.obs["expected_lineage_rank"].to_numpy(dtype=float),
            "dpt_pseudotime": adata.obs["dpt_pseudotime"].to_numpy(dtype=float),
        }
    )
    summary = summarize_dpt_pilot(cell_scores, root_cell_id=str(adata.obs_names[root_idx]))
    return summary, cell_scores


def summarize_dpt_pilot(cell_scores: pd.DataFrame, *, root_cell_id: str) -> pd.DataFrame:
    """Summarize DPT ordering overall, by state, and by technical strata."""

    rows = []
    rows.append(_association_row(cell_scores, "overall", "all", root_cell_id=root_cell_id))
    for column in ("flex_version", "device_guided", "condition"):
        if column in cell_scores:
            for value, group in cell_scores.groupby(column, dropna=False):
                rows.append(_association_row(group, column, str(value), root_cell_id=root_cell_id))

    state = (
        cell_scores.groupby("fine_celltype", observed=True)
        .agg(
            n_cells=("dpt_pseudotime", "size"),
            expected_lineage_rank=("expected_lineage_rank", "median"),
            median_dpt=("dpt_pseudotime", "median"),
            mean_dpt=("dpt_pseudotime", "mean"),
        )
        .reset_index()
        .sort_values(["expected_lineage_rank", "fine_celltype"])
    )
    for _, row in state.iterrows():
        rows.append(
            {
                "summary_type": "fine_celltype",
                "stratum": row["fine_celltype"],
                "n_cells": int(row["n_cells"]),
                "root_cell_id": root_cell_id,
                "spearman_r": np.nan,
                "p_value": np.nan,
                "expected_rank_min": float(row["expected_lineage_rank"]),
                "expected_rank_max": float(row["expected_lineage_rank"]),
                "median_dpt": float(row["median_dpt"]),
                "mean_dpt": float(row["mean_dpt"]),
                "claim_status": "state_summary",
            }
        )
    return pd.DataFrame(rows)


def render_regeneration_dynamics_report(
    *,
    audit: pd.DataFrame,
    feasibility: pd.DataFrame,
    dpt_summary: pd.DataFrame | None = None,
) -> str:
    """Render a concise Markdown report for the dynamics gate."""

    dpt_row = None
    if dpt_summary is not None and not dpt_summary.empty:
        overall = dpt_summary[
            (dpt_summary["summary_type"] == "overall") & (dpt_summary["stratum"] == "all")
        ]
        if not overall.empty:
            dpt_row = overall.iloc[0]

    lines = [
        "# Regeneration Dynamics Feasibility",
        "",
        "Status: exploratory topology pilot complete; velocity-based fate inference is not supported by current public H5AD layers.",
        "",
        "## Input Audit",
        "",
    ]
    for _, row in audit.iterrows():
        lines.append(
            "- "
            f"`{row['h5ad_path']}`: {int(row['n_obs']) if pd.notna(row['n_obs']) else 'NA'} cells, "
            f"{int(row['n_vars']) if pd.notna(row['n_vars']) else 'NA'} genes, "
            f"layers=`{row['layers'] or 'none'}`, obsm=`{row['obsm'] or 'none'}`."
        )
    lines.extend(
        [
            "",
            "## Method Decision",
            "",
            "- RNA velocity/scVelo and CellRank velocity kernels are no-go because no inspected object contains both `spliced` and `unspliced` layers.",
            "- CellRank/Palantir-style fate claims remain inappropriate without temporal, perturbational, velocity, or spatial direction evidence.",
            "- Scanpy diffusion pseudotime is allowed only as an exploratory topology sanity check rooted in HBC states.",
            "",
        ]
    )
    if dpt_row is not None:
        caveat = _state_ordering_caveat(dpt_summary)
        lines.extend(
            [
                "## DPT Pilot",
                "",
                f"- Overall Spearman correlation between expected lineage rank and DPT pseudotime: {float(dpt_row['spearman_r']):.3f} "
                f"(n={int(dpt_row['n_cells'])}).",
                f"- State-order caveat: {caveat}",
                "- Interpretation: use as guarded support that the lineage embedding contains a basal-to-neuronal ordering; do not call it lineage flux or fate transition.",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "## DPT Pilot",
                "",
                "- No pilot summary was produced.",
                "",
            ]
        )

    lines.extend(
        [
            "## Claim Status",
            "",
            "Allowed claim: basal, progenitor, immature OSN, and mature OSN states show an exploratory pseudotime-like ordering in the lineage embedding.",
            "",
            "Prohibited claim: cross-sectional data demonstrate regeneration trajectory speed, lineage flux, or cell fate transition rates.",
            "",
            "## Feasibility Table",
            "",
            _markdown_table(feasibility),
            "",
        ]
    )
    return "\n".join(lines)


def write_regeneration_dynamics_outputs(
    *,
    audit: pd.DataFrame,
    feasibility: pd.DataFrame,
    report: str,
    audit_out: str | Path,
    feasibility_out: str | Path,
    report_out: str | Path,
    dpt_summary: pd.DataFrame | None = None,
    dpt_cells: pd.DataFrame | None = None,
    summary_out: str | Path | None = None,
    cells_out: str | Path | None = None,
    figure_pdf: str | Path | None = None,
    figure_png: str | Path | None = None,
) -> None:
    """Write regeneration-dynamics tables, report, and optional figure."""

    ensure_parent(audit_out)
    audit.to_csv(audit_out, sep="\t", index=False)
    ensure_parent(feasibility_out)
    feasibility.to_csv(feasibility_out, sep="\t", index=False)
    ensure_parent(report_out)
    Path(report_out).write_text(report, encoding="utf-8")
    if dpt_summary is not None and summary_out is not None:
        ensure_parent(summary_out)
        dpt_summary.to_csv(summary_out, sep="\t", index=False)
    if dpt_cells is not None and cells_out is not None:
        ensure_parent(cells_out)
        dpt_cells.to_csv(cells_out, sep="\t", index=False)
    if dpt_cells is not None and (figure_pdf or figure_png):
        _plot_dpt_by_state(dpt_cells, figure_pdf=figure_pdf, figure_png=figure_png)


def _h5ad_shape(handle: h5py.File) -> tuple[int | float, int | float]:
    if "X" not in handle:
        return np.nan, np.nan
    x = handle["X"]
    if hasattr(x, "shape") and len(x.shape) == 2:
        return int(x.shape[0]), int(x.shape[1])
    shape = x.attrs.get("shape")
    if shape is not None and len(shape) == 2:
        return int(shape[0]), int(shape[1])
    if {"data", "indices", "indptr"} <= set(x.keys()):
        n_obs = len(x["indptr"]) - 1
        n_vars = int(x.attrs.get("shape", [np.nan, np.nan])[1])
        return n_obs, n_vars
    return np.nan, np.nan


def _markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "_No rows._"
    columns = list(frame.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in frame.iterrows():
        values = [str(row[col]).replace("|", "/") for col in columns]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def _state_ordering_caveat(dpt_summary: pd.DataFrame) -> str:
    states = dpt_summary[dpt_summary["summary_type"].eq("fine_celltype")].copy()
    if states.empty:
        return "state-level medians were not available."
    states = states.sort_values(["expected_rank_min", "stratum"])
    medians = pd.to_numeric(states["median_dpt"], errors="coerce").to_numpy(dtype=float)
    labels = states["stratum"].astype(str).tolist()
    dips = []
    for idx in range(1, len(medians)):
        if np.isfinite(medians[idx]) and np.isfinite(medians[idx - 1]) and medians[idx] < medians[idx - 1]:
            dips.append(f"{labels[idx]} median is below {labels[idx - 1]}")
    if not dips:
        return "state medians are monotonic across the expected ordering."
    return "; ".join(dips[:3]) + "."


def _input_note(
    *,
    has_spliced: bool,
    has_unspliced: bool,
    has_lineage_labels: bool,
    has_x_scvi: bool,
) -> str:
    notes = []
    if not (has_spliced and has_unspliced):
        notes.append("velocity layers absent")
    if has_lineage_labels and has_x_scvi:
        notes.append("lineage scVI topology pilot feasible")
    elif has_lineage_labels:
        notes.append("lineage labels present")
    else:
        notes.append("lineage labels absent")
    return "; ".join(notes)


def _stratified_obs_sample(adata, column: str, *, max_cells: int, seed: int):
    rng = np.random.default_rng(seed)
    obs = adata.obs.reset_index(names="_obs_name")
    per_group = max(1, max_cells // max(1, obs[column].nunique()))
    selected = []
    for _, group in obs.groupby(column, observed=True):
        take = min(len(group), per_group)
        selected.extend(rng.choice(group["_obs_name"].to_numpy(), size=take, replace=False).tolist())
    if len(selected) < max_cells:
        remaining = np.setdiff1d(obs["_obs_name"].to_numpy(), np.asarray(selected), assume_unique=False)
        take = min(len(remaining), max_cells - len(selected))
        if take:
            selected.extend(rng.choice(remaining, size=take, replace=False).tolist())
    return adata[selected].copy()


def _choose_root_index(adata) -> int:
    fine = adata.obs["fine_celltype"].astype(str)
    candidates = np.flatnonzero(fine.isin(["Quiescent_HBC", "Cycling_HBC"]).to_numpy())
    if len(candidates) == 0:
        candidates = np.flatnonzero(fine.map(LINEAGE_FINE_ORDER).fillna(99).to_numpy() == 0)
    if len(candidates) == 0:
        return 0
    diffmap = np.asarray(adata.obsm["X_diffmap"])
    root_space = diffmap[candidates, 1:] if diffmap.shape[1] > 1 else diffmap[candidates, :]
    center = np.nanmedian(root_space, axis=0)
    distance = np.linalg.norm(root_space - center, axis=1)
    return int(candidates[int(np.nanargmin(distance))])


def _association_row(group: pd.DataFrame, summary_type: str, stratum: str, *, root_cell_id: str) -> dict:
    valid = group[["expected_lineage_rank", "dpt_pseudotime"]].dropna()
    if len(valid) >= 3 and valid["expected_lineage_rank"].nunique() >= 2:
        corr = stats.spearmanr(valid["expected_lineage_rank"], valid["dpt_pseudotime"])
        r = float(corr.statistic)
        p = float(corr.pvalue)
    else:
        r = np.nan
        p = np.nan
    return {
        "summary_type": summary_type,
        "stratum": stratum,
        "n_cells": int(len(group)),
        "root_cell_id": root_cell_id,
        "spearman_r": r,
        "p_value": p,
        "expected_rank_min": float(valid["expected_lineage_rank"].min()) if len(valid) else np.nan,
        "expected_rank_max": float(valid["expected_lineage_rank"].max()) if len(valid) else np.nan,
        "median_dpt": float(valid["dpt_pseudotime"].median()) if len(valid) else np.nan,
        "mean_dpt": float(valid["dpt_pseudotime"].mean()) if len(valid) else np.nan,
        "claim_status": "exploratory_topology",
    }


def _plot_dpt_by_state(
    cell_scores: pd.DataFrame,
    *,
    figure_pdf: str | Path | None,
    figure_png: str | Path | None,
) -> None:
    import matplotlib.pyplot as plt

    ordered = [
        state
        for state, _ in sorted(LINEAGE_FINE_ORDER.items(), key=lambda item: (item[1], item[0]))
        if state in set(cell_scores["fine_celltype"])
    ]
    data = [
        cell_scores.loc[cell_scores["fine_celltype"].eq(state), "dpt_pseudotime"].dropna().to_numpy()
        for state in ordered
    ]
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    ax.boxplot(data, labels=ordered, showfliers=False, patch_artist=True)
    ax.set_ylabel("DPT pseudotime")
    ax.set_xlabel("Expected basal-to-neuronal state order")
    ax.set_title("Exploratory regeneration-lineage pseudotime pilot")
    ax.tick_params(axis="x", rotation=35, labelsize=8)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    for out in (figure_pdf, figure_png):
        if out:
            ensure_parent(out)
            fig.savefig(out, dpi=220)
    plt.close(fig)
