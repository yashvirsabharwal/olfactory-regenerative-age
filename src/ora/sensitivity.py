"""Sensitivity analysis helpers for ORA modeling."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd

from .age_model import train_ora_models


@dataclass
class SensitivityResult:
    scenarios: pd.DataFrame
    performance: pd.DataFrame
    scores: pd.DataFrame


def default_ora_scenarios(manifest: pd.DataFrame, min_cell_thresholds: list[int] | None = None) -> list[dict[str, object]]:
    """Build conservative ORA sensitivity scenarios from available manifest columns."""

    scenarios: list[dict[str, object]] = [{"scenario": "baseline", "filter_type": "all", "filter_value": "all"}]
    for col in ["chemistry", "collection_method"]:
        if col not in manifest:
            continue
        for value in sorted(manifest[col].dropna().astype(str).unique()):
            if value and value != "unknown":
                scenarios.append({"scenario": f"{col}__{value}", "filter_type": col, "filter_value": value})
    for threshold in min_cell_thresholds or [500, 1000, 5000, 10000]:
        scenarios.append(
            {
                "scenario": f"min_total_cells__{threshold}",
                "filter_type": "min_total_cells",
                "filter_value": int(threshold),
            }
        )
    if {"chemistry", "collection_method"}.issubset(manifest.columns):
        scenarios.append(
            {
                "scenario": "matched_flex_v2_device",
                "filter_type": "compound",
                "filter_value": "chemistry=flex_v2;collection_method=device",
            }
        )
    if "total_cells" in manifest:
        scenarios.extend(
            [
                {
                    "scenario": "exclude_yield_extremes__10pct",
                    "filter_type": "exclude_yield_extremes",
                    "filter_value": 0.10,
                },
                {
                    "scenario": "exclude_yield_extremes__20pct",
                    "filter_type": "exclude_yield_extremes",
                    "filter_value": 0.20,
                },
            ]
        )
    scenarios.append({"scenario": "healthy_only", "filter_type": "disease_group", "filter_value": "healthy"})
    return scenarios


def run_ora_sensitivity(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    model_config: dict[str, Any],
    *,
    min_cell_thresholds: list[int] | None = None,
    min_train_donors: int = 20,
) -> SensitivityResult:
    """Run ORA age-model sensitivity scenarios with donor/sample filters."""

    scenarios = default_ora_scenarios(manifest, min_cell_thresholds=min_cell_thresholds)
    scenario_rows = []
    performance_rows = []
    score_rows = []
    for spec in scenarios:
        scenario = str(spec["scenario"])
        subset_manifest = filter_manifest_for_scenario(manifest, spec)
        donor_ids = set(subset_manifest["donor_id"].astype(str))
        subset_features = features[features["donor_id"].astype(str).isin(donor_ids)].copy()
        healthy_train = subset_manifest[
            subset_manifest.get("usable_for_ora_training", False).astype(bool) & subset_manifest["age"].notna()
        ]["donor_id"].nunique()
        row = {
            "scenario": scenario,
            "filter_type": spec["filter_type"],
            "filter_value": spec["filter_value"],
            "donors": subset_manifest["donor_id"].nunique(),
            "samples": subset_manifest["sample_id"].nunique() if "sample_id" in subset_manifest else subset_manifest.shape[0],
            "healthy_train_donors": healthy_train,
            "ad_donors": _count_group(subset_manifest, "ad"),
            "pd_donors": _count_group(subset_manifest, "pd"),
            "total_cells": int(pd.to_numeric(subset_manifest.get("total_cells"), errors="coerce").fillna(0).sum()),
            "status": "ok" if healthy_train >= min_train_donors else "too_few_training_donors",
        }
        scenario_rows.append(row)
        if row["status"] != "ok":
            continue
        try:
            result = train_ora_models(subset_features, subset_manifest, model_config)
        except Exception as exc:  # pragma: no cover - defensive status capture for real runs
            scenario_rows[-1]["status"] = "failed"
            scenario_rows[-1]["error"] = str(exc)
            continue
        perf = result.performance.copy()
        perf.insert(0, "scenario", scenario)
        perf["filter_type"] = spec["filter_type"]
        perf["filter_value"] = spec["filter_value"]
        perf["healthy_train_donors"] = healthy_train
        performance_rows.append(perf)
        scores = result.predictions.copy()
        scores.insert(0, "scenario", scenario)
        score_rows.append(scores)

    return SensitivityResult(
        scenarios=pd.DataFrame(scenario_rows),
        performance=pd.concat(performance_rows, ignore_index=True) if performance_rows else pd.DataFrame(),
        scores=pd.concat(score_rows, ignore_index=True) if score_rows else pd.DataFrame(),
    )


def filter_manifest_for_scenario(manifest: pd.DataFrame, scenario: dict[str, object]) -> pd.DataFrame:
    """Apply one sensitivity scenario to a manifest."""

    filter_type = str(scenario["filter_type"])
    value = scenario["filter_value"]
    if filter_type == "all":
        return manifest.copy()
    if filter_type == "min_total_cells":
        total_cells = pd.to_numeric(manifest.get("total_cells"), errors="coerce").fillna(0)
        return manifest[total_cells >= int(value)].copy()
    if filter_type == "exclude_yield_extremes":
        total_cells = pd.to_numeric(manifest.get("total_cells"), errors="coerce")
        if total_cells.notna().sum() < 4:
            return manifest.copy()
        fraction = float(value)
        low = total_cells.quantile(fraction)
        high = total_cells.quantile(1.0 - fraction)
        return manifest[total_cells.between(low, high, inclusive="both")].copy()
    if filter_type == "compound":
        subset = manifest.copy()
        for token in str(value).split(";"):
            if "=" not in token:
                continue
            col, expected = token.split("=", 1)
            if col not in subset:
                return subset.iloc[0:0].copy()
            subset = subset[subset[col].astype(str).eq(expected)].copy()
        return subset
    if filter_type not in manifest:
        return manifest.iloc[0:0].copy()
    return manifest[manifest[filter_type].astype(str).eq(str(value))].copy()


def _count_group(manifest: pd.DataFrame, group: str) -> int:
    if "disease_group" not in manifest:
        return 0
    return int(manifest.loc[manifest["disease_group"].astype(str).eq(group), "donor_id"].nunique())
