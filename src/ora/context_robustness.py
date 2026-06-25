"""Leave-one-context-out robustness checks for donor-level ORA models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .age_model import (
    _performance_row,
    biological_feature_columns,
    fit_model_predictions,
    fit_preprocessor,
    model_names_from_config,
    summarize_repeated_performance,
    transform_preprocessor,
)


@dataclass
class LeaveContextOutResult:
    feasibility: pd.DataFrame
    performance: pd.DataFrame
    summary: pd.DataFrame
    scores: pd.DataFrame
    feature_stability: pd.DataFrame


def run_leave_context_out(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    model_config: dict[str, Any],
    *,
    contexts: list[str] | None = None,
    repeats: int | None = None,
    min_train_donors: int = 40,
    min_test_donors: int = 5,
) -> LeaveContextOutResult:
    """Train on all but one context level and evaluate the held-out level."""

    model_config = dict(model_config or {})
    repeats = max(1, int(repeats or model_config.get("outer_cv_repeats", 1)))
    donor_meta = _donor_manifest(manifest)
    donor_meta = _add_yield_contexts(donor_meta)
    context_specs = build_context_splits(
        donor_meta,
        contexts=contexts,
        min_train_donors=min_train_donors,
        min_test_donors=min_test_donors,
    )
    data = features.set_index("donor_id").join(donor_meta.set_index("donor_id"), how="inner", rsuffix="_meta").reset_index()
    performance_rows = []
    score_rows = []
    importance_rows = []
    feasibility_rows = []
    base_seed = int(model_config.get("random_seed", 42))

    for spec in context_specs:
        feasibility_rows.append(spec.as_row())
        if spec.status != "ok":
            continue
        train_mask = data["donor_id"].astype(str).isin(spec.train_donors)
        test_mask = data["donor_id"].astype(str).isin(spec.test_donors)
        train = data[train_mask].copy().reset_index(drop=True)
        test = data[test_mask].copy().reset_index(drop=True)
        for repeat in range(repeats):
            repeat_config = dict(model_config)
            repeat_config["random_seed"] = base_seed + repeat
            try:
                repeat_perf, repeat_scores, repeat_importance = _fit_context_models(
                    train,
                    test,
                    repeat_config,
                    context=spec.context,
                    level=spec.level,
                    repeat=repeat,
                )
            except Exception as exc:  # pragma: no cover - real-run defensive capture
                row = spec.as_row()
                row["status"] = "failed"
                row["error"] = str(exc)
                feasibility_rows.append(row)
                break
            performance_rows.append(repeat_perf)
            score_rows.append(repeat_scores)
            if not repeat_importance.empty:
                importance_rows.append(repeat_importance)

    performance = pd.concat(performance_rows, ignore_index=True) if performance_rows else pd.DataFrame()
    scores = pd.concat(score_rows, ignore_index=True) if score_rows else pd.DataFrame()
    feature_stability = _summarize_context_feature_stability(
        pd.concat(importance_rows, ignore_index=True) if importance_rows else pd.DataFrame()
    )
    summary = summarize_leave_context_performance(performance)
    return LeaveContextOutResult(
        feasibility=pd.DataFrame(feasibility_rows),
        performance=performance,
        summary=summary,
        scores=scores,
        feature_stability=feature_stability,
    )


def summarize_leave_context_performance(performance: pd.DataFrame) -> pd.DataFrame:
    """Summarize leave-context-out performance across repeats."""

    if performance.empty:
        return pd.DataFrame()
    rows = []
    for (context, level), frame in performance.groupby(["context", "level"], observed=True, sort=False):
        summary = summarize_repeated_performance(frame)
        summary.insert(0, "level", level)
        summary.insert(0, "context", context)
        summary["test_donors"] = int(frame["test_donors"].median())
        summary["train_donors"] = int(frame["train_donors"].median())
        rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def build_context_splits(
    donor_meta: pd.DataFrame,
    *,
    contexts: list[str] | None = None,
    min_train_donors: int = 40,
    min_test_donors: int = 5,
) -> list["_ContextSpec"]:
    """Create leave-one-context split specifications."""

    contexts = contexts or ["site", "chemistry", "collection_method", "sex", "race_ethnicity", "yield_bin"]
    trainable = donor_meta[_boolean_series(donor_meta.get("usable_for_ora_training")) & donor_meta["age"].notna()].copy()
    specs: list[_ContextSpec] = []
    for context in contexts:
        if context not in trainable:
            specs.append(_ContextSpec(context, "missing", "skipped_missing_context", [], [], "Context column is unavailable."))
            continue
        values = _context_values(trainable[context])
        if values.nunique(dropna=False) < 2:
            level = str(values.iloc[0]) if not values.empty else "missing"
            specs.append(_ContextSpec(context, level, "skipped_single_level", [], [], "Context has fewer than two levels."))
            continue
        for level in sorted(values.unique()):
            test = trainable.loc[values.eq(level), "donor_id"].astype(str).tolist()
            train = trainable.loc[~values.eq(level), "donor_id"].astype(str).tolist()
            status = "ok"
            note = ""
            if len(test) < min_test_donors:
                status = "too_few_test_donors"
                note = f"Need at least {min_test_donors} held-out donors."
            elif len(train) < min_train_donors:
                status = "too_few_train_donors"
                note = f"Need at least {min_train_donors} training donors."
            specs.append(_ContextSpec(context, str(level), status, train, test, note))
    return specs


def _fit_context_models(
    train: pd.DataFrame,
    test: pd.DataFrame,
    model_config: dict[str, Any],
    *,
    context: str,
    level: str,
    repeat: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feature_cols = biological_feature_columns(train, model_config)
    if not feature_cols:
        raise ValueError("No biological numeric feature columns available for context model.")
    max_missing = float(model_config.get("missingness_max_fraction", 0.30))
    missing_fraction = train[feature_cols].isna().mean()
    feature_cols = [col for col in feature_cols if missing_fraction[col] <= max_missing]
    if not feature_cols:
        raise ValueError("All feature columns exceeded the missingness threshold.")

    prep = fit_preprocessor(train[feature_cols])
    x_train = transform_preprocessor(train[feature_cols], prep)
    x_test = transform_preprocessor(test[feature_cols], prep)
    y_train = train["age"].astype(float).to_numpy()
    y_test = test["age"].astype(float).to_numpy()
    perf_rows = []
    score_rows = []
    importance_rows = []
    for model_name in model_names_from_config(model_config):
        pred, importance, backend = fit_model_predictions(model_name, x_train, y_train, x_test, model_config)
        row = _performance_row(model_name, y_test, pred, backend)
        row.update(
            {
                "repeat": repeat,
                "context": context,
                "level": level,
                "train_donors": int(train["donor_id"].nunique()),
                "test_donors": int(test["donor_id"].nunique()),
                "n_features": len(feature_cols),
            }
        )
        perf_rows.append(row)
        score_rows.append(
            pd.DataFrame(
                {
                    "repeat": repeat,
                    "context": context,
                    "level": level,
                    "donor_id": test["donor_id"].astype(str).to_numpy(),
                    "model": model_name,
                    "chronological_age": y_test,
                    "ora": pred,
                    "error": pred - y_test,
                    "abs_error": np.abs(pred - y_test),
                }
            )
        )
        if importance is not None:
            importance_rows.append(
                pd.DataFrame(
                    {
                        "repeat": repeat,
                        "context": context,
                        "level": level,
                        "model": model_name,
                        "feature": feature_cols,
                        "importance": importance,
                    }
                )
            )
    return (
        pd.DataFrame(perf_rows),
        pd.concat(score_rows, ignore_index=True),
        pd.concat(importance_rows, ignore_index=True) if importance_rows else pd.DataFrame(),
    )


def _summarize_context_feature_stability(feature_importance: pd.DataFrame) -> pd.DataFrame:
    if feature_importance.empty:
        return pd.DataFrame()
    frame = feature_importance.copy()
    frame["importance"] = pd.to_numeric(frame["importance"], errors="coerce").fillna(0.0)
    frame["selected"] = frame["importance"].abs().gt(1e-12)
    summary = (
        frame.groupby(["context", "level", "model", "feature"], observed=True)
        .agg(
            mean_importance=("importance", "mean"),
            sd_importance=("importance", "std"),
            selection_fraction=("selected", "mean"),
            repeats=("repeat", "nunique"),
        )
        .reset_index()
    )
    summary["abs_mean_importance"] = summary["mean_importance"].abs()
    return summary.sort_values(
        ["context", "level", "model", "selection_fraction", "abs_mean_importance"],
        ascending=[True, True, True, False, False],
    ).reset_index(drop=True)


def _donor_manifest(manifest: pd.DataFrame) -> pd.DataFrame:
    donor_meta = manifest.sort_values(["donor_id", "sample_id"] if "sample_id" in manifest else ["donor_id"])
    return donor_meta.drop_duplicates("donor_id").copy()


def _add_yield_contexts(donor_meta: pd.DataFrame) -> pd.DataFrame:
    output = donor_meta.copy()
    if "total_cells" not in output:
        output["yield_bin"] = "missing"
        return output
    total_cells = pd.to_numeric(output["total_cells"], errors="coerce")
    median = total_cells.median()
    output["yield_bin"] = np.where(total_cells >= median, "high_yield", "low_yield")
    output.loc[total_cells.isna(), "yield_bin"] = "missing"
    return output


def _context_values(values: pd.Series) -> pd.Series:
    output = values.astype("string").str.strip().fillna("missing")
    output = output.mask(output.str.lower().isin(["", "nan", "none", "na", "<na>", "unknown"]), "missing")
    return output.astype(str)


def _boolean_series(values: pd.Series | None) -> pd.Series:
    if values is None:
        return pd.Series(dtype=bool)
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y", "t"})


@dataclass
class _ContextSpec:
    context: str
    level: str
    status: str
    train_donors: list[str]
    test_donors: list[str]
    note: str = ""

    def as_row(self) -> dict[str, object]:
        return {
            "context": self.context,
            "level": self.level,
            "status": self.status,
            "train_donors": len(self.train_donors),
            "test_donors": len(self.test_donors),
            "note": self.note,
        }
