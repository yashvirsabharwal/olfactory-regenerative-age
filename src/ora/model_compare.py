"""Model comparison helpers for ORA feature-set benchmarks."""

from __future__ import annotations

import pandas as pd


REQUIRED_SUMMARY_COLUMNS = {
    "model",
    "mae_mean",
    "rmse_mean",
    "r2_mean",
    "spearman_r_mean",
}


def rank_feature_set_summaries(
    summaries: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Stack repeated-CV summaries and rank all feature-set/model combinations."""

    frames = []
    for feature_set, summary in summaries.items():
        _validate_summary(summary, label=feature_set)
        frame = summary.copy()
        frame.insert(0, "feature_set", feature_set)
        frame["mae_rank_within_feature_set"] = frame["mae_mean"].rank(method="min", ascending=True).astype(int)
        frames.append(frame)

    ranked = pd.concat(frames, ignore_index=True)
    ranked["mae_rank_overall"] = ranked["mae_mean"].rank(method="min", ascending=True).astype(int)
    ranked["is_best_within_feature_set"] = ranked["mae_rank_within_feature_set"].eq(1)
    ranked["is_best_overall"] = ranked["mae_rank_overall"].eq(1)
    return ranked.sort_values(["mae_mean", "rmse_mean", "model", "feature_set"]).reset_index(drop=True)


def compare_feature_set_deltas(
    base_summary: pd.DataFrame,
    augmented_summary: pd.DataFrame,
    *,
    base_label: str = "composition",
    augmented_label: str = "composition_plus_modules",
) -> pd.DataFrame:
    """Compare repeated-CV metrics for matched models across two feature sets."""

    _validate_summary(base_summary, label=base_label)
    _validate_summary(augmented_summary, label=augmented_label)
    base = base_summary.set_index("model")
    augmented = augmented_summary.set_index("model")
    shared_models = sorted(set(base.index).intersection(set(augmented.index)))
    rows = []
    for model in shared_models:
        base_row = base.loc[model]
        augmented_row = augmented.loc[model]
        rows.append(
            {
                "model": model,
                "base_feature_set": base_label,
                "augmented_feature_set": augmented_label,
                "base_mae_mean": float(base_row["mae_mean"]),
                "augmented_mae_mean": float(augmented_row["mae_mean"]),
                "delta_mae_mean": float(augmented_row["mae_mean"] - base_row["mae_mean"]),
                "base_rmse_mean": float(base_row["rmse_mean"]),
                "augmented_rmse_mean": float(augmented_row["rmse_mean"]),
                "delta_rmse_mean": float(augmented_row["rmse_mean"] - base_row["rmse_mean"]),
                "base_r2_mean": float(base_row["r2_mean"]),
                "augmented_r2_mean": float(augmented_row["r2_mean"]),
                "delta_r2_mean": float(augmented_row["r2_mean"] - base_row["r2_mean"]),
                "base_spearman_r_mean": float(base_row["spearman_r_mean"]),
                "augmented_spearman_r_mean": float(augmented_row["spearman_r_mean"]),
                "delta_spearman_r_mean": float(augmented_row["spearman_r_mean"] - base_row["spearman_r_mean"]),
                "mae_improved": bool(augmented_row["mae_mean"] < base_row["mae_mean"]),
            }
        )
    return pd.DataFrame(rows).sort_values(["delta_mae_mean", "model"]).reset_index(drop=True)


def _validate_summary(summary: pd.DataFrame, *, label: str) -> None:
    missing = sorted(REQUIRED_SUMMARY_COLUMNS.difference(summary.columns))
    if missing:
        raise ValueError(f"{label} summary is missing required columns: {', '.join(missing)}")
