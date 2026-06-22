"""Age-bin robustness summaries for Milo-style neighborhood memberships."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class AgeBinConfig:
    run_name: str = "neighborhood_run"
    age_column: str = "age"
    donor_column: str = "donor_id"
    fdr_threshold: float = 0.10
    bins: tuple[tuple[str, float, float], ...] = (
        ("lt45", 0.0, 45.0),
        ("45_59", 45.0, 60.0),
        ("60_74", 60.0, 75.0),
        ("75_plus", 75.0, np.inf),
    )


def summarize_neighborhood_age_bins(
    memberships: pd.DataFrame,
    donor_metadata: pd.DataFrame,
    *,
    da_table: pd.DataFrame | None = None,
    config: AgeBinConfig | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize neighborhood membership intensity across donor age bins.

    Membership rows can overlap across neighborhoods, so the fraction reported here
    is a donor-normalized membership-participation fraction. It is intended as a
    directionality robustness check, not as a replacement for the DA regression.
    """

    cfg = config or AgeBinConfig()
    _require_columns(memberships, ["neighborhood_id", cfg.donor_column], "memberships")
    _require_columns(donor_metadata, [cfg.donor_column, cfg.age_column], "donor metadata")

    donors = _prepare_donors(donor_metadata, cfg)
    members = memberships[["neighborhood_id", cfg.donor_column]].copy()
    members[cfg.donor_column] = members[cfg.donor_column].astype(str)
    members["neighborhood_id"] = pd.to_numeric(members["neighborhood_id"], errors="coerce").astype("Int64")
    members = members[members["neighborhood_id"].notna()].copy()
    members["neighborhood_id"] = members["neighborhood_id"].astype(int)
    members = members[members[cfg.donor_column].isin(donors.index)].copy()
    if members.empty:
        empty = _neighborhood_columns()
        return pd.DataFrame(columns=empty), _summary_table(cfg, donors, pd.DataFrame(columns=empty))

    donor_totals = members[cfg.donor_column].value_counts().rename("donor_membership_rows")
    observed_donors = donors.loc[donor_totals.index.intersection(donors.index)].copy()
    counts = members.groupby(["neighborhood_id", cfg.donor_column], observed=True).size().rename("membership_rows").reset_index()
    counts = counts.merge(donor_totals, left_on=cfg.donor_column, right_index=True, how="left")
    counts = counts.merge(
        donors[[cfg.age_column, "age_bin", "age_bin_order"]],
        left_on=cfg.donor_column,
        right_index=True,
        how="left",
    )
    counts["membership_fraction"] = (counts["membership_rows"].astype(float) + 0.5) / (
        counts["donor_membership_rows"].astype(float) + 1.0
    )
    counts["membership_logit_fraction"] = _logit(counts["membership_fraction"].to_numpy(dtype=float))
    bin_long = _summarize_bins(counts, cfg)
    neighborhoods = _summarize_neighborhoods(bin_long, cfg)
    if da_table is not None and not da_table.empty:
        neighborhoods = _merge_da(neighborhoods, da_table, cfg)
    summary = _summary_table(cfg, observed_donors, neighborhoods)
    return neighborhoods, summary


def _prepare_donors(donor_metadata: pd.DataFrame, cfg: AgeBinConfig) -> pd.DataFrame:
    donors = donor_metadata.copy()
    donors[cfg.donor_column] = donors[cfg.donor_column].astype(str)
    donors[cfg.age_column] = pd.to_numeric(donors[cfg.age_column], errors="coerce")
    donors = donors[donors[cfg.age_column].notna()].copy()
    donors["age_bin"] = pd.Categorical(
        _assign_bins(donors[cfg.age_column].to_numpy(dtype=float), cfg),
        categories=[label for label, _, _ in cfg.bins],
        ordered=True,
    )
    donors = donors[donors["age_bin"].notna()].copy()
    donors["age_bin_order"] = donors["age_bin"].cat.codes.astype(int)
    return donors.set_index(cfg.donor_column, drop=False)


def _assign_bins(ages: np.ndarray, cfg: AgeBinConfig) -> list[str | float]:
    labels: list[str | float] = []
    for age in ages:
        label: str | float = np.nan
        for candidate, lower, upper in cfg.bins:
            if age >= lower and age < upper:
                label = candidate
                break
        labels.append(label)
    return labels


def _summarize_bins(counts: pd.DataFrame, cfg: AgeBinConfig) -> pd.DataFrame:
    grouped = (
        counts.groupby(["neighborhood_id", "age_bin", "age_bin_order"], observed=True)
        .agg(
            n_donors=(cfg.donor_column, "nunique"),
            total_membership_rows=("membership_rows", "sum"),
            mean_membership_fraction=("membership_fraction", "mean"),
            median_membership_fraction=("membership_fraction", "median"),
            mean_logit_fraction=("membership_logit_fraction", "mean"),
            median_logit_fraction=("membership_logit_fraction", "median"),
        )
        .reset_index()
    )
    grouped["run"] = cfg.run_name
    return grouped[
        [
            "run",
            "neighborhood_id",
            "age_bin",
            "age_bin_order",
            "n_donors",
            "total_membership_rows",
            "mean_membership_fraction",
            "median_membership_fraction",
            "mean_logit_fraction",
            "median_logit_fraction",
        ]
    ]


def _summarize_neighborhoods(bin_long: pd.DataFrame, cfg: AgeBinConfig) -> pd.DataFrame:
    rows = []
    for neighborhood_id, group in bin_long.groupby("neighborhood_id", observed=True):
        ordered = group.sort_values("age_bin_order")
        by_label = ordered.set_index("age_bin", drop=False)
        values = ordered["median_logit_fraction"].to_numpy(dtype=float)
        orders = ordered["age_bin_order"].to_numpy(dtype=float)
        first_label = str(ordered["age_bin"].iloc[0])
        last_label = str(ordered["age_bin"].iloc[-1])
        young = _bin_value(by_label, first_label)
        old = _bin_value(by_label, last_label)
        old_minus_young = old - young if np.isfinite(young) and np.isfinite(old) else np.nan
        trend = _spearman(orders, values)
        rows.append(
            {
                "run": cfg.run_name,
                "neighborhood_id": int(neighborhood_id),
                "n_age_bins_observed": int(group["age_bin"].nunique()),
                "min_bin_donors": int(group["n_donors"].min()),
                "max_bin_donors": int(group["n_donors"].max()),
                "youngest_bin": first_label,
                "oldest_bin": last_label,
                "youngest_median_logit_fraction": young,
                "oldest_median_logit_fraction": old,
                "old_minus_young_median_logit_fraction": old_minus_young,
                "age_bin_spearman": trend,
                "age_bin_direction": _direction(old_minus_young),
                "age_bin_profile": _profile(group),
            }
        )
    result = pd.DataFrame(rows, columns=_neighborhood_columns())
    return result.sort_values("neighborhood_id").reset_index(drop=True)


def _merge_da(neighborhoods: pd.DataFrame, da_table: pd.DataFrame, cfg: AgeBinConfig) -> pd.DataFrame:
    keep_cols = [
        "neighborhood_id",
        "top_fine_celltype",
        "top_coarse_celltype",
        "age_coef",
        "age_fdr",
        "status",
    ]
    da = da_table[[col for col in keep_cols if col in da_table.columns]].copy()
    merged = neighborhoods.merge(da, on="neighborhood_id", how="left")
    merged["age_coef"] = pd.to_numeric(merged.get("age_coef"), errors="coerce")
    merged["age_fdr"] = pd.to_numeric(merged.get("age_fdr"), errors="coerce")
    merged["is_age_associated_fdr_0_10"] = merged["status"].astype(str).eq("tested") & merged["age_fdr"].lt(cfg.fdr_threshold)
    merged["regression_direction"] = np.where(merged["age_coef"].gt(0), "positive", np.where(merged["age_coef"].lt(0), "negative", "zero"))
    merged["bin_agrees_with_regression"] = (
        merged["is_age_associated_fdr_0_10"]
        & merged["regression_direction"].isin(["negative", "positive"])
        & merged["regression_direction"].eq(merged["age_bin_direction"])
    )
    return merged


def _summary_table(cfg: AgeBinConfig, donors: pd.DataFrame, neighborhoods: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {"run": cfg.run_name, "metric": "donors_with_age_bins", "value": int(donors.shape[0]), "detail": "donors assigned to configured age bins"},
    ]
    for label in [label for label, _, _ in cfg.bins]:
        count = int(donors["age_bin"].astype(str).eq(label).sum()) if not donors.empty else 0
        rows.append({"run": cfg.run_name, "metric": f"donors_{label}", "value": count, "detail": "donors in age bin"})
    if not neighborhoods.empty:
        sig = neighborhoods.get("is_age_associated_fdr_0_10", pd.Series(False, index=neighborhoods.index)).fillna(False).astype(bool)
        neg_sig = sig & neighborhoods.get("regression_direction", "").astype(str).eq("negative")
        pos_sig = sig & neighborhoods.get("regression_direction", "").astype(str).eq("positive")
        agrees = neighborhoods.get("bin_agrees_with_regression", pd.Series(False, index=neighborhoods.index)).fillna(False).astype(bool)
        rows.extend(
            [
                {"run": cfg.run_name, "metric": "neighborhoods_total", "value": int(neighborhoods.shape[0]), "detail": "neighborhoods with membership-bin summaries"},
                {"run": cfg.run_name, "metric": "age_fdr_lt_0_10", "value": int(sig.sum()), "detail": "DA-significant neighborhoods in merged table"},
                {"run": cfg.run_name, "metric": "negative_sig_bin_agreement", "value": int((neg_sig & agrees).sum()), "detail": "negative DA neighborhoods also lower in oldest than youngest age bin"},
                {"run": cfg.run_name, "metric": "positive_sig_bin_agreement", "value": int((pos_sig & agrees).sum()), "detail": "positive DA neighborhoods also higher in oldest than youngest age bin"},
            ]
        )
    return pd.DataFrame(rows, columns=["run", "metric", "value", "detail"])


def _bin_value(by_label: pd.DataFrame, label: str) -> float:
    if label not in by_label.index:
        return float("nan")
    return float(by_label.loc[label, "median_logit_fraction"])


def _profile(group: pd.DataFrame) -> str:
    ordered = group.sort_values("age_bin_order")
    return ";".join(f"{row.age_bin}:{row.median_logit_fraction:.4g}:n={int(row.n_donors)}" for row in ordered.itertuples())


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.nanstd(y[ok]) == 0:
        return float("nan")
    return float(pd.Series(x[ok]).corr(pd.Series(y[ok]), method="spearman"))


def _direction(value: float) -> str:
    if not np.isfinite(value):
        return "not_available"
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "zero"


def _logit(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(values, 1e-12, 1.0 - 1e-12)
    return np.log(clipped / (1.0 - clipped))


def _require_columns(frame: pd.DataFrame, columns: list[str], context: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{context} missing required columns: {', '.join(missing)}")


def _neighborhood_columns() -> list[str]:
    return [
        "run",
        "neighborhood_id",
        "n_age_bins_observed",
        "min_bin_donors",
        "max_bin_donors",
        "youngest_bin",
        "oldest_bin",
        "youngest_median_logit_fraction",
        "oldest_median_logit_fraction",
        "old_minus_young_median_logit_fraction",
        "age_bin_spearman",
        "age_bin_direction",
        "age_bin_profile",
    ]
