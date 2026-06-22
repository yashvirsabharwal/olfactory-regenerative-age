"""Count export and parity summaries for neighborhood differential abundance."""

from __future__ import annotations

import numpy as np
import pandas as pd


def export_neighborhood_count_inputs(
    memberships: pd.DataFrame,
    donor_metadata: pd.DataFrame,
    *,
    donor_column: str = "donor_id",
    age_column: str = "age",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Export neighborhood-by-donor counts and matching donor design metadata."""

    _require_columns(memberships, ["neighborhood_id", donor_column], "memberships")
    _require_columns(donor_metadata, [donor_column, age_column], "donor metadata")
    donors = donor_metadata.copy()
    donors[donor_column] = donors[donor_column].astype(str)
    donors[age_column] = pd.to_numeric(donors[age_column], errors="coerce")
    donors = donors[donors[age_column].notna()].drop_duplicates(donor_column).copy()

    members = memberships[["neighborhood_id", donor_column]].copy()
    members[donor_column] = members[donor_column].astype(str)
    members["neighborhood_id"] = pd.to_numeric(members["neighborhood_id"], errors="coerce").astype("Int64")
    members = members[members["neighborhood_id"].notna() & members[donor_column].isin(set(donors[donor_column]))].copy()
    members["neighborhood_id"] = members["neighborhood_id"].astype(int)
    if members.empty:
        raise ValueError("No membership rows remain after donor filtering.")

    counts = (
        members.groupby(["neighborhood_id", donor_column], observed=True)
        .size()
        .unstack(fill_value=0)
        .sort_index()
        .astype(int)
    )
    donor_order = [donor for donor in donors[donor_column].astype(str).tolist() if donor in counts.columns]
    counts = counts.reindex(columns=donor_order, fill_value=0)
    counts_out = counts.reset_index()
    design = donors.set_index(donor_column).loc[donor_order].reset_index()
    design["age_scaled"] = _zscore(design[age_column].to_numpy(dtype=float))
    if "total_cells" in design:
        design["log_total_cells"] = np.log1p(pd.to_numeric(design["total_cells"], errors="coerce").fillna(0).to_numpy(dtype=float))
    summary = pd.DataFrame(
        [
            {"metric": "neighborhoods", "value": int(counts.shape[0]), "detail": "rows in count matrix"},
            {"metric": "donors", "value": int(counts.shape[1]), "detail": "columns in count matrix"},
            {"metric": "membership_rows", "value": int(members.shape[0]), "detail": "filtered membership rows"},
            {"metric": "total_counts", "value": int(counts.to_numpy().sum()), "detail": "sum of neighborhood-by-donor counts"},
        ]
    )
    return counts_out, design, summary


def summarize_edger_parity(
    python_da: pd.DataFrame,
    edger_da: pd.DataFrame,
    *,
    run_name: str,
    fdr_threshold: float = 0.10,
    top_n: int = 100,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize concordance between Python donor-logit DA and edgeR count DA."""

    _require_columns(python_da, ["neighborhood_id", "age_coef", "age_fdr"], "Python DA")
    _require_columns(edger_da, ["neighborhood_id", "logFC", "FDR"], "edgeR DA")
    py = python_da[["neighborhood_id", "age_coef", "age_pvalue", "age_fdr", "top_fine_celltype", "top_coarse_celltype", "status"]].copy()
    edge = edger_da[["neighborhood_id", "logFC", "logCPM", "F", "PValue", "FDR"]].copy()
    merged = py.merge(edge, on="neighborhood_id", how="inner")
    merged["run"] = run_name
    merged["python_significant"] = merged["status"].astype(str).eq("tested") & pd.to_numeric(merged["age_fdr"], errors="coerce").lt(fdr_threshold)
    merged["edger_significant"] = pd.to_numeric(merged["FDR"], errors="coerce").lt(fdr_threshold)
    merged["python_direction"] = _direction(merged["age_coef"])
    merged["edger_direction"] = _direction(merged["logFC"])
    merged["direction_agrees"] = merged["python_direction"].eq(merged["edger_direction"])
    merged["python_rank"] = pd.to_numeric(merged["age_pvalue"], errors="coerce").rank(method="min", na_option="bottom")
    merged["edger_rank"] = pd.to_numeric(merged["PValue"], errors="coerce").rank(method="min", na_option="bottom")

    top_py = merged.nsmallest(top_n, "python_rank")
    top_edge = merged.nsmallest(top_n, "edger_rank")
    sig_py = merged[merged["python_significant"]]
    sig_edge = merged[merged["edger_significant"]]
    summary = pd.DataFrame(
        [
            {"run": run_name, "metric": "neighborhoods_compared", "value": int(merged.shape[0]), "detail": "shared neighborhoods"},
            {"run": run_name, "metric": "python_fdr_lt_0_10", "value": int(sig_py.shape[0]), "detail": "Python DA significant neighborhoods"},
            {"run": run_name, "metric": "edger_fdr_lt_0_10", "value": int(sig_edge.shape[0]), "detail": "edgeR QL significant neighborhoods"},
            {"run": run_name, "metric": "significant_overlap", "value": int((merged["python_significant"] & merged["edger_significant"]).sum()), "detail": "neighborhoods significant in both"},
            {"run": run_name, "metric": "python_sig_direction_agreement", "value": _mean_bool(sig_py["direction_agrees"]), "detail": "direction agreement among Python-significant neighborhoods"},
            {"run": run_name, "metric": "edger_sig_direction_agreement", "value": _mean_bool(sig_edge["direction_agrees"]), "detail": "direction agreement among edgeR-significant neighborhoods"},
            {"run": run_name, "metric": f"top_{top_n}_python_direction_agreement", "value": _mean_bool(top_py["direction_agrees"]), "detail": "direction agreement among top Python neighborhoods"},
            {"run": run_name, "metric": f"top_{top_n}_edger_direction_agreement", "value": _mean_bool(top_edge["direction_agrees"]), "detail": "direction agreement among top edgeR neighborhoods"},
            {"run": run_name, "metric": "signed_effect_spearman", "value": _corr(merged["age_coef"], merged["logFC"]), "detail": "Spearman correlation of signed age effects"},
            {"run": run_name, "metric": "rank_spearman", "value": _corr(merged["python_rank"], merged["edger_rank"]), "detail": "Spearman correlation of p-value ranks"},
        ]
    )
    return merged, summary


def _zscore(values: np.ndarray) -> np.ndarray:
    sd = float(np.nanstd(values))
    if sd == 0 or not np.isfinite(sd):
        return np.zeros_like(values, dtype=float)
    return (values - float(np.nanmean(values))) / sd


def _direction(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    return pd.Series(np.where(numeric.gt(0), "positive", np.where(numeric.lt(0), "negative", "zero")), index=values.index)


def _mean_bool(values: pd.Series) -> float:
    if values.empty:
        return float("nan")
    return float(values.mean())


def _corr(left: pd.Series, right: pd.Series) -> float:
    frame = pd.DataFrame({"left": pd.to_numeric(left, errors="coerce"), "right": pd.to_numeric(right, errors="coerce")}).dropna()
    if frame.shape[0] < 3 or frame["left"].nunique() < 2 or frame["right"].nunique() < 2:
        return float("nan")
    return float(frame["left"].corr(frame["right"], method="spearman"))


def _require_columns(frame: pd.DataFrame, columns: list[str], context: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{context} missing required columns: {', '.join(missing)}")
