"""Donor-level compositional age modeling."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np
import pandas as pd

from .aggregate import DEFAULT_PSEUDOCOUNT
from .stats import bh_fdr
from .utils import normalize_token


DEFAULT_COVARIATES = ("sex", "chemistry", "collection_method", "site", "log10_total_cells")
DEFAULT_MIN_SCENARIO_DONORS = 30
DEFAULT_MIN_NONZERO_DONORS = 5


@dataclass
class CompositionalModelResult:
    """Container for primary and sensitivity compositional model outputs."""

    summary: pd.DataFrame
    sensitivity: pd.DataFrame


def run_compositional_age_model(
    counts: pd.DataFrame,
    manifest: pd.DataFrame,
    *,
    age_associations: pd.DataFrame | None = None,
    pseudocount: float = DEFAULT_PSEUDOCOUNT,
    min_scenario_donors: int = DEFAULT_MIN_SCENARIO_DONORS,
    min_nonzero_donors: int = DEFAULT_MIN_NONZERO_DONORS,
    covariates: tuple[str, ...] = DEFAULT_COVARIATES,
) -> CompositionalModelResult:
    """Fit CLR composition-by-age models and scenario sensitivities."""

    count_matrix = _donor_count_matrix(counts)
    donor_meta = _donor_metadata(manifest)
    joined = donor_meta.join(count_matrix.sum(axis=1).rename("observed_total_cells"), how="inner")
    joined["age"] = pd.to_numeric(joined["age"], errors="coerce")
    joined["log10_total_cells"] = np.log10(pd.to_numeric(joined["observed_total_cells"], errors="coerce") + 1.0)

    scenarios = _build_scenarios(joined, min_scenario_donors=min_scenario_donors)
    rows: list[dict[str, Any]] = []
    for scenario in scenarios:
        scenario_meta = joined.loc[scenario["mask"]].copy()
        scenario_counts = count_matrix.loc[scenario_meta.index]
        rows.extend(
            _fit_scenario(
                scenario=scenario,
                counts=scenario_counts,
                meta=scenario_meta,
                pseudocount=pseudocount,
                min_nonzero_donors=min_nonzero_donors,
                covariates=covariates,
            )
        )

    sensitivity = pd.DataFrame(rows)
    if sensitivity.empty:
        sensitivity = pd.DataFrame(
            columns=[
                "scenario",
                "scenario_label",
                "cell_state",
                "ora_clr_feature",
                "n_donors",
                "age_beta_per_10_years",
                "standard_error_per_10_years",
                "p_value",
                "q_value",
                "age_direction",
                "status",
            ]
        )
    else:
        sensitivity["q_value"] = np.nan
        for scenario_name, index in sensitivity.groupby("scenario").groups.items():
            p_values = sensitivity.loc[index, "p_value"].to_numpy(dtype=float)
            sensitivity.loc[index, "q_value"] = bh_fdr(p_values)
        sensitivity = _merge_age_associations(sensitivity, age_associations)
        sensitivity = sensitivity.sort_values(
            ["scenario_rank", "q_value", "p_value", "cell_state"], na_position="last"
        ).reset_index(drop=True)

    summary = _primary_summary(sensitivity)
    return CompositionalModelResult(summary=summary, sensitivity=sensitivity)


def _donor_count_matrix(counts: pd.DataFrame) -> pd.DataFrame:
    required = {"donor_id", "fine_cell_type", "cell_count"}
    missing = sorted(required.difference(counts.columns))
    if missing:
        raise KeyError(f"Cell count table missing columns: {missing}")
    matrix = (
        counts.assign(cell_count=pd.to_numeric(counts["cell_count"], errors="coerce").fillna(0.0))
        .groupby(["donor_id", "fine_cell_type"], observed=True)["cell_count"]
        .sum()
        .unstack(fill_value=0.0)
        .sort_index(axis=1)
    )
    matrix.index = matrix.index.astype(str)
    return matrix


def _donor_metadata(manifest: pd.DataFrame) -> pd.DataFrame:
    if "donor_id" not in manifest.columns:
        raise KeyError("Manifest missing required column: donor_id")
    if "age" not in manifest.columns:
        raise KeyError("Manifest missing required column: age")
    sort_cols = [col for col in ["donor_id", "sample_id"] if col in manifest.columns]
    meta = manifest.sort_values(sort_cols).drop_duplicates("donor_id").copy()
    meta["donor_id"] = meta["donor_id"].astype(str)
    return meta.set_index("donor_id", drop=True)


def _build_scenarios(meta: pd.DataFrame, *, min_scenario_donors: int) -> list[dict[str, Any]]:
    primary_mask = _analysis_mask(meta)
    scenarios: list[dict[str, Any]] = [
        {
            "scenario": "primary_all_healthy",
            "scenario_label": "Primary healthy donors",
            "scenario_rank": 0,
            "mask": primary_mask,
            "scenario_status": "ok" if int(primary_mask.sum()) >= min_scenario_donors else "too_few_samples",
        }
    ]

    if "passes_strict_ora_training_rule" in meta.columns:
        strict_mask = primary_mask & _as_bool(meta["passes_strict_ora_training_rule"])
        scenarios.append(
            {
                "scenario": "strict_threshold",
                "scenario_label": "Strict cell-yield threshold",
                "scenario_rank": 1,
                "mask": strict_mask,
                "scenario_status": "ok" if int(strict_mask.sum()) >= min_scenario_donors else "too_few_samples",
            }
        )

    if {"chemistry", "collection_method"}.issubset(meta.columns):
        candidates = meta.loc[primary_mask, ["chemistry", "collection_method"]].copy()
        candidates["chemistry"] = _clean_category(candidates["chemistry"])
        candidates["collection_method"] = _clean_category(candidates["collection_method"])
        combo_counts = candidates.groupby(["chemistry", "collection_method"], dropna=False).size()
        if not combo_counts.empty:
            chemistry, collection_method = combo_counts.sort_values(ascending=False).index[0]
            combo_mask = (
                primary_mask
                & _clean_category(meta["chemistry"]).eq(chemistry)
                & _clean_category(meta["collection_method"]).eq(collection_method)
            )
            label = f"Single chemistry/collection: {chemistry}/{collection_method}"
            scenarios.append(
                {
                    "scenario": f"single_{_scenario_token(chemistry)}_{_scenario_token(collection_method)}",
                    "scenario_label": label,
                    "scenario_rank": 2,
                    "mask": combo_mask,
                    "scenario_status": "ok" if int(combo_mask.sum()) >= min_scenario_donors else "too_few_samples",
                }
            )
    return scenarios


def _analysis_mask(meta: pd.DataFrame) -> pd.Series:
    mask = meta["age"].notna()
    if "usable_for_ora_training" in meta.columns:
        mask = mask & _as_bool(meta["usable_for_ora_training"])
    elif "passes_primary_ora_training_rule" in meta.columns:
        mask = mask & _as_bool(meta["passes_primary_ora_training_rule"])
    elif "is_healthy" in meta.columns:
        mask = mask & _as_bool(meta["is_healthy"])
    return mask.fillna(False)


def _fit_scenario(
    *,
    scenario: dict[str, Any],
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    pseudocount: float,
    min_nonzero_donors: int,
    covariates: tuple[str, ...],
) -> list[dict[str, Any]]:
    n_donors = int(meta.shape[0])
    states = [
        col
        for col in counts.columns
        if int((pd.to_numeric(counts[col], errors="coerce").fillna(0.0) > 0).sum()) >= min_nonzero_donors
    ]
    if not states:
        return [
            _scenario_empty_row(
                scenario,
                n_donors=n_donors,
                cell_state="",
                status="no_eligible_cell_states",
                covariates_used="",
            )
        ]

    design = _design_matrix(meta, covariates=covariates)
    scenario_status = scenario.get("scenario_status", "ok")
    if scenario_status != "ok":
        return [
            _scenario_empty_row(
                scenario,
                n_donors=n_donors,
                cell_state=state,
                status=str(scenario_status),
                covariates_used=";".join(design["covariates_used"]),
            )
            for state in states
        ]

    x = design["matrix"]
    if x.shape[0] <= x.shape[1] + 1:
        return [
            _scenario_empty_row(
                scenario,
                n_donors=n_donors,
                cell_state=state,
                status="underdetermined",
                covariates_used=";".join(design["covariates_used"]),
            )
            for state in states
        ]

    clr = _clr_transform(counts[states], pseudocount=pseudocount)
    totals = counts.sum(axis=1).replace(0, np.nan)
    rows: list[dict[str, Any]] = []
    for state in states:
        y = pd.to_numeric(clr[state], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(y) & np.isfinite(x).all(axis=1)
        state_counts = pd.to_numeric(counts[state], errors="coerce").fillna(0.0)
        if int(valid.sum()) <= x.shape[1] + 1:
            status = "underdetermined"
            rows.append(
                _scenario_empty_row(
                    scenario,
                    n_donors=int(valid.sum()),
                    cell_state=state,
                    status=status,
                    covariates_used=";".join(design["covariates_used"]),
                )
            )
            continue
        result = _ols(y[valid], x[valid])
        age_idx = design["names"].index("age")
        age_coef = float(result["coef"][age_idx] * 10.0)
        age_se = float(result["se"][age_idx] * 10.0)
        z_value = float(result["coef"][age_idx] / result["se"][age_idx]) if result["se"][age_idx] > 0 else np.nan
        rows.append(
            {
                "scenario": scenario["scenario"],
                "scenario_label": scenario["scenario_label"],
                "scenario_rank": scenario["scenario_rank"],
                "cell_state": state,
                "ora_clr_feature": _clr_feature_name(state),
                "n_donors": int(valid.sum()),
                "nonzero_donors": int((state_counts > 0).sum()),
                "prevalence": float((state_counts > 0).mean()),
                "total_cell_count": float(state_counts.sum()),
                "mean_fraction": float((state_counts / totals).replace([np.inf, -np.inf], np.nan).mean()),
                "age_beta_per_10_years": age_coef,
                "standard_error_per_10_years": age_se,
                "z_value": z_value,
                "p_value": _two_sided_normal_p(z_value),
                "q_value": np.nan,
                "age_direction": _direction(age_coef),
                "covariates_used": ";".join(design["covariates_used"]),
                "method": "pseudocount_clr_ols",
                "compositional_baseline": "donor_geometric_mean_clr",
                "status": "ok",
            }
        )
    return rows


def _clr_transform(counts: pd.DataFrame, *, pseudocount: float) -> pd.DataFrame:
    adjusted = counts.astype(float) + pseudocount
    proportions = adjusted.div(adjusted.sum(axis=1), axis=0)
    log_values = np.log(proportions)
    return log_values.sub(log_values.mean(axis=1), axis=0)


def _design_matrix(meta: pd.DataFrame, *, covariates: tuple[str, ...]) -> dict[str, Any]:
    frame = meta.copy()
    age = pd.to_numeric(frame["age"], errors="coerce").to_numpy(dtype=float)
    x_parts = [np.ones(frame.shape[0]), age]
    names = ["intercept", "age"]
    covariates_used: list[str] = []

    for covariate in covariates:
        if covariate not in frame.columns:
            continue
        series = frame[covariate]
        if pd.api.types.is_numeric_dtype(series):
            numeric = pd.to_numeric(series, errors="coerce")
            if numeric.notna().sum() == 0:
                continue
            numeric = numeric.fillna(numeric.median())
            if numeric.nunique(dropna=True) <= 1:
                continue
            x_parts.append(numeric.to_numpy(dtype=float))
            names.append(covariate)
            covariates_used.append(covariate)
            continue

        categorical = _clean_category(series)
        if categorical.nunique(dropna=True) <= 1:
            continue
        dummies = pd.get_dummies(categorical, prefix=covariate, drop_first=True, dtype=float)
        added = False
        for dummy_col in dummies.columns:
            if dummies[dummy_col].nunique(dropna=True) <= 1:
                continue
            x_parts.append(dummies[dummy_col].to_numpy(dtype=float))
            names.append(str(dummy_col))
            added = True
        if added:
            covariates_used.append(covariate)

    x = np.vstack(x_parts).T
    return {"matrix": x, "names": names, "covariates_used": covariates_used}


def _merge_age_associations(sensitivity: pd.DataFrame, age_associations: pd.DataFrame | None) -> pd.DataFrame:
    output = sensitivity.copy()
    output["ora_age_association_beta_per_10_years"] = np.nan
    output["ora_age_association_fdr"] = np.nan
    output["ora_age_association_direction"] = ""
    output["direction_concordant_with_ora_age_association"] = pd.Series(
        [pd.NA] * output.shape[0],
        dtype="boolean",
        index=output.index,
    )
    if age_associations is None or age_associations.empty or "feature" not in age_associations.columns:
        return output

    assoc = age_associations.drop_duplicates("feature").set_index("feature")
    for idx, row in output.iterrows():
        feature = row.get("ora_clr_feature")
        if feature not in assoc.index:
            continue
        assoc_row = assoc.loc[feature]
        assoc_beta = pd.to_numeric(pd.Series([assoc_row.get("beta_per_10_years")]), errors="coerce").iloc[0]
        assoc_fdr = pd.to_numeric(pd.Series([assoc_row.get("fdr")]), errors="coerce").iloc[0]
        assoc_direction = str(assoc_row.get("direction", ""))
        output.at[idx, "ora_age_association_beta_per_10_years"] = assoc_beta
        output.at[idx, "ora_age_association_fdr"] = assoc_fdr
        output.at[idx, "ora_age_association_direction"] = assoc_direction
        if np.isfinite(float(row.get("age_beta_per_10_years", np.nan))) and np.isfinite(float(assoc_beta)):
            output.at[idx, "direction_concordant_with_ora_age_association"] = bool(
                np.sign(float(row["age_beta_per_10_years"])) == np.sign(float(assoc_beta))
            )
    return output


def _primary_summary(sensitivity: pd.DataFrame) -> pd.DataFrame:
    if sensitivity.empty or "scenario" not in sensitivity.columns:
        return sensitivity.copy()
    primary = sensitivity[sensitivity["scenario"].eq("primary_all_healthy")].copy()
    if primary.empty:
        return primary

    stability = []
    for cell_state, group in sensitivity[sensitivity["status"].eq("ok")].groupby("cell_state", dropna=False):
        sensitivity_group = group[~group["scenario"].eq("primary_all_healthy")]
        primary_rows = group[group["scenario"].eq("primary_all_healthy")]
        primary_direction = primary_rows["age_direction"].iloc[0] if not primary_rows.empty else ""
        directions = [
            f"{row.scenario}:{row.age_direction}"
            for row in sensitivity_group.itertuples()
            if getattr(row, "age_direction", "") in {"positive", "negative"}
        ]
        stable = bool(directions) and all(item.rsplit(":", 1)[-1] == primary_direction for item in directions)
        stability.append(
            {
                "cell_state": cell_state,
                "sensitivity_ok_scenarios": int(sensitivity_group.shape[0]),
                "sensitivity_directions": "|".join(directions),
                "directionally_stable_in_sensitivity": stable,
            }
        )
    stability_frame = pd.DataFrame(stability)
    if not stability_frame.empty:
        primary = primary.merge(stability_frame, on="cell_state", how="left")
    primary["significant_q_0_10"] = pd.to_numeric(primary["q_value"], errors="coerce") <= 0.10
    primary["supported_by_primary_and_sensitivity"] = (
        primary["significant_q_0_10"].fillna(False) & primary["directionally_stable_in_sensitivity"].fillna(False)
    )
    return primary.sort_values(["q_value", "p_value", "cell_state"], na_position="last").reset_index(drop=True)


def _scenario_empty_row(
    scenario: dict[str, Any],
    *,
    n_donors: int,
    cell_state: str,
    status: str,
    covariates_used: str,
) -> dict[str, Any]:
    return {
        "scenario": scenario["scenario"],
        "scenario_label": scenario["scenario_label"],
        "scenario_rank": scenario["scenario_rank"],
        "cell_state": cell_state,
        "ora_clr_feature": _clr_feature_name(cell_state) if cell_state else "",
        "n_donors": n_donors,
        "nonzero_donors": np.nan,
        "prevalence": np.nan,
        "total_cell_count": np.nan,
        "mean_fraction": np.nan,
        "age_beta_per_10_years": np.nan,
        "standard_error_per_10_years": np.nan,
        "z_value": np.nan,
        "p_value": np.nan,
        "q_value": np.nan,
        "age_direction": "",
        "covariates_used": covariates_used,
        "method": "pseudocount_clr_ols",
        "compositional_baseline": "donor_geometric_mean_clr",
        "status": status,
    }


def _ols(y: np.ndarray, x: np.ndarray) -> dict[str, np.ndarray]:
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ coef
    df = max(x.shape[0] - x.shape[1], 1)
    sigma2 = float((resid @ resid) / df)
    xtx_inv = np.linalg.pinv(x.T @ x)
    se = np.sqrt(np.clip(np.diag(xtx_inv) * sigma2, 0, np.inf))
    return {"coef": coef, "se": se}


def _two_sided_normal_p(value: float) -> float:
    if not np.isfinite(value):
        return np.nan
    return math.erfc(abs(float(value)) / math.sqrt(2.0))


def _direction(value: float) -> str:
    if not np.isfinite(value):
        return ""
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "zero"


def _clr_feature_name(cell_state: object) -> str:
    text = normalize_token(cell_state)
    text = text.replace("/", " ").replace("+", " plus ").replace(" ", "_")
    return f"clr__{text}"


def _scenario_token(value: object) -> str:
    token = normalize_token(value).replace(" ", "_").replace("/", "_")
    return token or "unknown"


def _clean_category(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("unknown").str.strip().replace("", "unknown")


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    normalized = series.astype("string").str.strip().str.lower()
    return normalized.isin({"true", "1", "yes", "y"})
