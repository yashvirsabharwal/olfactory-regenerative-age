"""Regeneration-module metadata, age associations, and ORA correlations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from .regeneration_axis import compute_feature_age_associations
from .utils import normalize_token


@dataclass(frozen=True)
class RegenerationModuleMetadata:
    module: str
    module_feature: str
    description: str
    theme: str
    expected_age_direction: str
    source: str
    citation: str


DEFAULT_ANALYSIS_SETS = (
    ("primary", "usable_for_ora_training"),
    ("strict", "passes_strict_ora_training_rule"),
)

DEFAULT_ADJUSTMENT_COVARIATES = (
    "sex",
    "race_ethnicity",
    "chemistry",
    "collection_method",
    "total_cells",
    "lineage_cells",
)


def parse_regeneration_module_metadata(config: dict) -> pd.DataFrame:
    """Parse module-level metadata from the regeneration gene-set YAML."""

    rows = []
    for module, spec in config.get("gene_sets", {}).items():
        if not isinstance(spec, dict):
            spec = {"genes": spec}
        rows.append(
            {
                "module": str(module),
                "module_feature": f"module_score__{_slugify(module)}",
                "description": str(spec.get("description", "")),
                "theme": str(spec.get("theme", "")),
                "expected_age_direction": str(spec.get("expected_age_direction", "unknown")),
                "source": str(spec.get("source", "")),
                "citation": str(spec.get("citation", "")),
            }
        )
    return pd.DataFrame(rows)


def build_regeneration_module_age_associations(
    *,
    donor_module_features: pd.DataFrame,
    manifest: pd.DataFrame,
    module_metadata: pd.DataFrame,
    coverage: pd.DataFrame | None = None,
    analysis_sets: tuple[tuple[str, str], ...] = DEFAULT_ANALYSIS_SETS,
    adjustment_covariates: tuple[str, ...] = DEFAULT_ADJUSTMENT_COVARIATES,
) -> pd.DataFrame:
    """Test donor-level module score association with age across cohort masks."""

    metadata = _metadata_lookup(module_metadata)
    coverage_lookup = _coverage_lookup(coverage)
    rows = []
    for analysis_set, mask_col in analysis_sets:
        if mask_col not in manifest:
            continue
        masked_manifest = manifest[["donor_id", "age", mask_col]].copy()
        masked_manifest["usable_for_ora_training"] = masked_manifest[mask_col].fillna(False).astype(bool)
        unadjusted = compute_feature_age_associations(donor_module_features, masked_manifest)
        adjusted_manifest = manifest.copy()
        adjusted_manifest["usable_for_ora_training"] = manifest[mask_col].fillna(False).astype(bool)
        adjusted = _covariate_adjusted_feature_age_associations(
            donor_module_features,
            adjusted_manifest,
            adjustment_covariates=adjustment_covariates,
        )
        for model_type, age in [("unadjusted", unadjusted), ("covariate_adjusted", adjusted)]:
            for feature, assoc in age.iterrows():
                rows.append(
                    _module_age_row(
                        analysis_set=analysis_set,
                        model_type=model_type,
                        feature=feature,
                        assoc=assoc,
                        metadata=metadata,
                        coverage_lookup=coverage_lookup,
                    )
                )
    result = pd.DataFrame(rows, columns=_age_columns())
    result["fdr"] = np.nan
    for _, idx in result.groupby(["analysis_set", "model_type"], observed=True).groups.items():
        p_values = pd.to_numeric(result.loc[idx, "p_value"], errors="coerce")
        ok = p_values.notna()
        if ok.any():
            result.loc[p_values.index[ok], "fdr"] = _benjamini_hochberg(
                p_values.loc[ok].to_numpy(dtype=float)
            )
    result["observed_vs_expected"] = result.apply(
        lambda row: _observed_vs_expected(
            row.get("expected_age_direction", "unknown"),
            row.get("direction", "not_tested"),
            row.get("fdr", np.nan),
        ),
        axis=1,
    )
    return result


def build_regeneration_module_ora_correlations(
    *,
    donor_module_features: pd.DataFrame,
    ora_feature_matrix: pd.DataFrame,
    module_metadata: pd.DataFrame,
    manifest: pd.DataFrame | None = None,
    feature_map: pd.DataFrame | None = None,
    analysis_sets: tuple[tuple[str, str], ...] = DEFAULT_ANALYSIS_SETS,
) -> pd.DataFrame:
    """Correlate regeneration modules with ORA composition/module features."""

    metadata = _metadata_lookup(module_metadata)
    feature_lookup = _feature_lookup(feature_map)
    module_cols = [
        col
        for col in donor_module_features.columns
        if col != "donor_id" and pd.api.types.is_numeric_dtype(donor_module_features[col])
    ]
    ora_cols = [
        col
        for col in ora_feature_matrix.columns
        if col != "donor_id" and pd.api.types.is_numeric_dtype(ora_feature_matrix[col])
    ]
    overlapping = set(module_cols).intersection(ora_cols)
    module_join_cols = {col: f"{col}__regen" for col in overlapping}
    ora_join_cols = {col: f"{col}__ora" for col in overlapping}
    joined = donor_module_features.rename(columns=module_join_cols).merge(
        ora_feature_matrix.rename(columns=ora_join_cols),
        on="donor_id",
        how="inner",
    )
    if manifest is not None and "donor_id" in manifest:
        keep_cols = ["donor_id", *[col for _, col in analysis_sets if col in manifest]]
        joined = joined.merge(manifest[keep_cols], on="donor_id", how="left")
    rows = []
    for analysis_set, mask_col in analysis_sets:
        frame = joined.copy()
        if mask_col in frame:
            frame = frame[frame[mask_col].fillna(False).astype(bool)].copy()
        for module_feature in module_cols:
            module = _module_from_feature(module_feature)
            meta = metadata.get(module, {})
            module_frame_col = module_join_cols.get(module_feature, module_feature)
            module_values = pd.to_numeric(frame[module_frame_col], errors="coerce")
            for feature in ora_cols:
                feature_frame_col = ora_join_cols.get(feature, feature)
                feature_values = pd.to_numeric(frame[feature_frame_col], errors="coerce")
                mask = module_values.notna() & feature_values.notna()
                n = int(mask.sum())
                if n < 8 or module_values.loc[mask].nunique() < 2 or feature_values.loc[mask].nunique() < 2:
                    r_value = np.nan
                    p_value = np.nan
                    direction = "not_tested"
                else:
                    r_value, p_value = stats.pearsonr(module_values.loc[mask], feature_values.loc[mask])
                    direction = "positive" if r_value > 0 else "negative" if r_value < 0 else "flat"
                feature_meta = feature_lookup.get(feature, {})
                rows.append(
                    {
                        "analysis_set": analysis_set,
                        "module": module,
                        "module_feature": module_feature,
                        "module_theme": meta.get("theme", ""),
                        "feature": feature,
                        "feature_theme": feature_meta.get("primary_theme", ""),
                        "specificity_class": feature_meta.get("specificity_class", ""),
                        "n": n,
                        "pearson_r": r_value,
                        "abs_pearson_r": abs(r_value) if pd.notna(r_value) else np.nan,
                        "p_value": p_value,
                        "direction": direction,
                    }
                )
    result = pd.DataFrame(rows, columns=_correlation_columns())
    result["fdr"] = np.nan
    for _, idx in result.groupby("analysis_set", observed=True).groups.items():
        p_values = pd.to_numeric(result.loc[idx, "p_value"], errors="coerce")
        ok = p_values.notna()
        if ok.any():
            result.loc[p_values.index[ok], "fdr"] = _benjamini_hochberg(
                p_values.loc[ok].to_numpy(dtype=float)
            )
    result["correlation_rank_within_module"] = (
        result.groupby(["analysis_set", "module"], observed=True)["abs_pearson_r"]
        .rank(method="min", ascending=False)
        .astype("Int64")
    )
    return result.sort_values(
        ["analysis_set", "module", "correlation_rank_within_module", "feature"],
        ascending=[True, True, True, True],
    ).reset_index(drop=True)


def _module_age_row(
    *,
    analysis_set: str,
    model_type: str,
    feature: str,
    assoc: pd.Series,
    metadata: dict[str, dict],
    coverage_lookup: dict[str, dict],
) -> dict[str, object]:
    module = _module_from_feature(feature)
    meta = metadata.get(module, {})
    coverage_row = coverage_lookup.get(module, {})
    return {
        "analysis_set": analysis_set,
        "model_type": model_type,
        "module": module,
        "module_feature": feature,
        "theme": meta.get("theme", ""),
        "description": meta.get("description", ""),
        "expected_age_direction": meta.get("expected_age_direction", "unknown"),
        "source": meta.get("source", ""),
        "citation": meta.get("citation", ""),
        "n_requested": coverage_row.get("n_requested", np.nan),
        "n_present": coverage_row.get("n_present", np.nan),
        "coverage_fraction": coverage_row.get("coverage_fraction", np.nan),
        "n": assoc.get("n", np.nan),
        "beta_per_10_years": assoc.get("beta_per_10_years", np.nan),
        "standard_error": assoc.get("standard_error", np.nan),
        "p_value": assoc.get("p_value", np.nan),
        "fdr": assoc.get("fdr", np.nan),
        "direction": assoc.get("direction", "not_tested"),
        "status": assoc.get("status", "not_tested"),
        "adjustment_covariates": assoc.get("adjustment_covariates", ""),
        "observed_vs_expected": _observed_vs_expected(
            meta.get("expected_age_direction", "unknown"),
            assoc.get("direction", "not_tested"),
            assoc.get("fdr", np.nan),
        ),
    }


def _covariate_adjusted_feature_age_associations(
    feature_matrix: pd.DataFrame,
    manifest: pd.DataFrame,
    *,
    adjustment_covariates: tuple[str, ...],
) -> pd.DataFrame:
    try:
        import statsmodels.api as sm
    except ModuleNotFoundError:
        return _not_tested_adjusted_rows(feature_matrix, "statsmodels_missing").set_index("feature", drop=False)

    manifest_cols = ["donor_id", "age", "usable_for_ora_training"]
    covariates = [col for col in adjustment_covariates if col in manifest.columns]
    manifest_cols.extend(covariates)
    frame = feature_matrix.merge(manifest[manifest_cols], on="donor_id", how="left")
    frame = frame[frame["usable_for_ora_training"].fillna(False).astype(bool)].copy()
    frame["age"] = pd.to_numeric(frame["age"], errors="coerce")
    rows = []
    for feature in _numeric_feature_columns(feature_matrix):
        values = pd.to_numeric(frame[feature], errors="coerce")
        work = frame.loc[values.notna() & frame["age"].notna(), ["age", *covariates]].copy()
        work["feature_value"] = values.loc[work.index].to_numpy(dtype=float)
        covariate_frame, used_covariates = _design_covariates(work, covariates)
        model_frame = pd.concat([work[["age", "feature_value"]], covariate_frame], axis=1).dropna()
        n = int(model_frame.shape[0])
        if n < 8 or model_frame["feature_value"].nunique(dropna=True) < 2:
            rows.append(_adjusted_age_row(feature, n, "not_enough_data", used_covariates))
            continue
        x = model_frame.drop(columns=["feature_value"])
        x = sm.add_constant(x.astype(float), has_constant="add")
        if n <= x.shape[1] + 2:
            rows.append(_adjusted_age_row(feature, n, "too_few_degrees_of_freedom", used_covariates))
            continue
        try:
            fit = sm.OLS(model_frame["feature_value"].astype(float), x).fit()
        except Exception:
            rows.append(_adjusted_age_row(feature, n, "fit_failed", used_covariates))
            continue
        if "age" not in fit.params:
            rows.append(_adjusted_age_row(feature, n, "age_term_missing", used_covariates))
            continue
        beta = float(fit.params["age"] * 10.0)
        se = float(fit.bse["age"] * 10.0)
        p_value = float(fit.pvalues["age"])
        direction = "positive" if beta > 0 else "negative" if beta < 0 else "flat"
        rows.append(
            {
                "feature": feature,
                "n": n,
                "beta_per_10_years": beta,
                "standard_error": se,
                "p_value": p_value,
                "direction": direction,
                "status": "ok",
                "adjustment_covariates": ";".join(used_covariates),
            }
        )
    return pd.DataFrame(rows).set_index("feature", drop=False)


def _design_covariates(frame: pd.DataFrame, covariates: list[str]) -> tuple[pd.DataFrame, list[str]]:
    pieces = []
    used = []
    for covariate in covariates:
        series = frame[covariate]
        numeric = pd.to_numeric(series, errors="coerce")
        if numeric.notna().sum() == series.notna().sum() and numeric.nunique(dropna=True) > 1:
            transformed = np.log1p(numeric.clip(lower=0)) if covariate.endswith("cells") else numeric
            pieces.append(transformed.rename(covariate))
            used.append(covariate)
            continue
        categorical = series.astype("string").fillna("unknown")
        if categorical.nunique(dropna=True) <= 1:
            continue
        dummies = pd.get_dummies(categorical, prefix=covariate, drop_first=True, dtype=float)
        if not dummies.empty:
            pieces.append(dummies)
            used.append(covariate)
    if not pieces:
        return pd.DataFrame(index=frame.index), []
    return pd.concat(pieces, axis=1), used


def _adjusted_age_row(
    feature: str,
    n: int,
    status: str,
    used_covariates: list[str],
) -> dict[str, object]:
    return {
        "feature": feature,
        "n": n,
        "beta_per_10_years": np.nan,
        "standard_error": np.nan,
        "p_value": np.nan,
        "direction": "not_tested",
        "status": status,
        "adjustment_covariates": ";".join(used_covariates),
    }


def _not_tested_adjusted_rows(feature_matrix: pd.DataFrame, status: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            _adjusted_age_row(feature, 0, status, [])
            for feature in _numeric_feature_columns(feature_matrix)
        ]
    )


def _metadata_lookup(module_metadata: pd.DataFrame) -> dict[str, dict]:
    if module_metadata is None or module_metadata.empty:
        return {}
    lookup = {}
    for _, row in module_metadata.iterrows():
        data = row.to_dict()
        module = str(row["module"])
        lookup[module] = data
        lookup[_slugify(module)] = data
    return lookup


def _coverage_lookup(coverage: pd.DataFrame | None) -> dict[str, dict]:
    if coverage is None or coverage.empty or "module" not in coverage:
        return {}
    lookup = {}
    for _, row in coverage.iterrows():
        data = row.to_dict()
        module = str(row["module"])
        lookup[module] = data
        lookup[_slugify(module)] = data
    return lookup


def _feature_lookup(feature_map: pd.DataFrame | None) -> dict[str, dict]:
    if feature_map is None or feature_map.empty or "feature" not in feature_map:
        return {}
    return {str(row["feature"]): row.to_dict() for _, row in feature_map.iterrows()}


def _module_from_feature(feature: str) -> str:
    return str(feature).removeprefix("module_score__")


def _numeric_feature_columns(feature_matrix: pd.DataFrame) -> list[str]:
    return [
        col
        for col in feature_matrix.columns
        if col != "donor_id" and pd.api.types.is_numeric_dtype(feature_matrix[col])
    ]


def _slugify(value: object) -> str:
    return normalize_token(value).replace("/", " ").replace("+", " plus ").replace(" ", "_")


def _observed_vs_expected(expected: object, observed: object, fdr: object) -> str:
    expected = str(expected)
    observed = str(observed)
    if observed not in {"positive", "negative"}:
        return "not_tested"
    fdr_value = pd.to_numeric(fdr, errors="coerce")
    if pd.isna(fdr_value) or fdr_value >= 0.05:
        return "observed_not_fdr_significant"
    if expected not in {"positive", "negative"}:
        return "no_directional_prior"
    return "aligned" if expected == observed else "opposite"


def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    n = len(ranked)
    adjusted = np.empty(n, dtype=float)
    cumulative = 1.0
    for idx in range(n - 1, -1, -1):
        rank = idx + 1
        cumulative = min(cumulative, ranked[idx] * n / rank)
        adjusted[idx] = cumulative
    result = np.empty(n, dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def _age_columns() -> list[str]:
    return [
        "analysis_set",
        "model_type",
        "module",
        "module_feature",
        "theme",
        "description",
        "expected_age_direction",
        "source",
        "citation",
        "n_requested",
        "n_present",
        "coverage_fraction",
        "n",
        "beta_per_10_years",
        "standard_error",
        "p_value",
        "fdr",
        "direction",
        "status",
        "adjustment_covariates",
        "observed_vs_expected",
    ]


def _correlation_columns() -> list[str]:
    return [
        "analysis_set",
        "module",
        "module_feature",
        "module_theme",
        "feature",
        "feature_theme",
        "specificity_class",
        "n",
        "pearson_r",
        "abs_pearson_r",
        "p_value",
        "direction",
    ]
