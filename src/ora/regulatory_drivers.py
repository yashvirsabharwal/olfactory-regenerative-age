"""First-pass regulatory-driver scoring from genomewide pseudobulk counts."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .regeneration_modules import DEFAULT_ADJUSTMENT_COVARIATES
from .utils import ensure_parent, normalize_token


DEFAULT_ANALYSIS_SETS = (
    ("primary", "usable_for_ora_training"),
    ("strict", "passes_strict_ora_training_rule"),
)


def parse_driver_metadata(config: dict) -> pd.DataFrame:
    """Parse curated regulatory-driver target-set metadata."""

    rows = []
    for driver_set, spec in config.get("gene_sets", {}).items():
        rows.append(
            {
                "driver_set": str(driver_set),
                "driver": str(spec.get("driver", driver_set)),
                "driver_class": str(spec.get("driver_class", "")),
                "target_theme": str(spec.get("target_theme", "")),
                "expected_age_direction": str(spec.get("expected_age_direction", "unknown")),
                "source": str(spec.get("source", "")),
                "citation": str(spec.get("citation", "")),
                "description": str(spec.get("description", "")),
                "genes": tuple(str(gene) for gene in spec.get("genes", []) if str(gene).strip()),
            }
        )
    return pd.DataFrame(rows)


def score_regulatory_driver_activity(
    *,
    counts_path: str | Path,
    metadata: pd.DataFrame,
    driver_metadata: pd.DataFrame,
    chunksize: int = 1000,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score curated driver target programs as mean log1p CPM in pseudobulk groups."""

    pseudobulk_ids = metadata["pseudobulk_id"].astype(str).tolist()
    target_rows = _read_target_gene_rows(counts_path, driver_metadata, chunksize=chunksize)
    coverage = _driver_gene_coverage(driver_metadata, target_rows)
    if target_rows.empty:
        return _empty_activity(metadata, driver_metadata), _empty_donor_activity(), coverage

    library = metadata.set_index("pseudobulk_id").loc[pseudobulk_ids, "sum_n_counts"].to_numpy(dtype=float)
    library = np.where(library > 0, library, np.nan)
    values = target_rows[pseudobulk_ids].to_numpy(dtype=float)
    log_cpm = np.log1p((values / library.reshape(1, -1)) * 1_000_000.0)
    token_to_positions = _token_positions(target_rows["gene_symbol"])

    activity_rows = []
    for _, driver in driver_metadata.iterrows():
        positions = []
        for gene in driver["genes"]:
            positions.extend(token_to_positions.get(normalize_token(gene), []))
        positions = sorted(set(positions))
        score = np.nanmean(log_cpm[positions, :], axis=0) if positions else np.full(len(pseudobulk_ids), np.nan)
        frame = metadata.copy()
        frame["driver_set"] = driver["driver_set"]
        frame["driver"] = driver["driver"]
        frame["driver_class"] = driver["driver_class"]
        frame["target_theme"] = driver["target_theme"]
        frame["activity_score"] = score
        frame["lineage_group"] = frame.apply(_lineage_group, axis=1)
        activity_rows.append(frame)
    activity = pd.concat(activity_rows, ignore_index=True)
    donor_activity = aggregate_driver_activity_by_donor_lineage(activity)
    return activity, donor_activity, coverage


def aggregate_driver_activity_by_donor_lineage(activity: pd.DataFrame) -> pd.DataFrame:
    """Aggregate pseudobulk driver activity to donor x lineage groups."""

    rows = []
    group_cols = ["driver_set", "driver", "driver_class", "target_theme", "donor_id", "lineage_group"]
    for keys, group in activity.groupby(group_cols, observed=True, dropna=False):
        rows.append(_weighted_activity_row(keys, group, group_cols))
    all_cols = ["driver_set", "driver", "driver_class", "target_theme", "donor_id"]
    for keys, group in activity.groupby(all_cols, observed=True, dropna=False):
        row = _weighted_activity_row((*keys, "all_lineages"), group, [*all_cols, "lineage_group"])
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["driver_set", "donor_id", "lineage_group"]).reset_index(drop=True)


def build_driver_age_associations(
    *,
    donor_activity: pd.DataFrame,
    manifest: pd.DataFrame,
    driver_metadata: pd.DataFrame,
    analysis_sets: tuple[tuple[str, str], ...] = DEFAULT_ANALYSIS_SETS,
    covariates: tuple[str, ...] = DEFAULT_ADJUSTMENT_COVARIATES,
) -> pd.DataFrame:
    """Run covariate-adjusted driver_activity ~ age models by lineage group."""

    rows = []
    metadata = _driver_lookup(driver_metadata)
    for analysis_set, mask_col in analysis_sets:
        if mask_col not in manifest:
            continue
        for (driver_set, lineage_group), group in donor_activity.groupby(
            ["driver_set", "lineage_group"],
            observed=True,
        ):
            merged = group.merge(manifest, on="donor_id", how="left")
            merged = merged[merged[mask_col].fillna(False).astype(bool)].copy()
            assoc = _ols_age_association(merged, covariates=covariates)
            meta = metadata.get(str(driver_set), {})
            rows.append(
                {
                    "analysis_set": analysis_set,
                    "driver_set": driver_set,
                    "driver": meta.get("driver", driver_set),
                    "driver_class": meta.get("driver_class", ""),
                    "target_theme": meta.get("target_theme", ""),
                    "lineage_group": lineage_group,
                    "expected_age_direction": meta.get("expected_age_direction", "unknown"),
                    **assoc,
                    "observed_vs_expected": _observed_vs_expected(
                        meta.get("expected_age_direction", "unknown"),
                        assoc.get("direction", "not_tested"),
                        assoc.get("fdr", np.nan),
                    ),
                }
            )
    result = pd.DataFrame(rows)
    result["fdr"] = np.nan
    for _, idx in result.groupby(["analysis_set", "lineage_group"], observed=True).groups.items():
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
    return result.sort_values(
        ["analysis_set", "lineage_group", "fdr", "driver_set"],
        na_position="last",
    ).reset_index(drop=True)


def build_driver_ora_correlations(
    *,
    donor_activity: pd.DataFrame,
    ora_scores: pd.DataFrame,
    driver_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Correlate donor-level driver activity with averaged ORA/ORAA scores."""

    scores = _summarize_ora_scores(ora_scores)
    metadata = _driver_lookup(driver_metadata)
    rows = []
    for (driver_set, lineage_group), group in donor_activity.groupby(
        ["driver_set", "lineage_group"],
        observed=True,
    ):
        merged = group.merge(scores, on="donor_id", how="inner")
        meta = metadata.get(str(driver_set), {})
        for metric in ["ora", "oraa"]:
            if metric not in merged:
                continue
            assoc = _pearson_association(merged["activity_score"], merged[metric])
            rows.append(
                {
                    "driver_set": driver_set,
                    "driver": meta.get("driver", driver_set),
                    "driver_class": meta.get("driver_class", ""),
                    "target_theme": meta.get("target_theme", ""),
                    "lineage_group": lineage_group,
                    "score_metric": metric,
                    **assoc,
                }
            )
    result = pd.DataFrame(rows)
    result["fdr"] = np.nan
    for _, idx in result.groupby(["score_metric", "lineage_group"], observed=True).groups.items():
        p_values = pd.to_numeric(result.loc[idx, "p_value"], errors="coerce")
        ok = p_values.notna()
        if ok.any():
            result.loc[p_values.index[ok], "fdr"] = _benjamini_hochberg(
                p_values.loc[ok].to_numpy(dtype=float)
            )
    result["abs_pearson_r"] = pd.to_numeric(result["pearson_r"], errors="coerce").abs()
    return result.sort_values(
        ["score_metric", "lineage_group", "fdr", "abs_pearson_r"],
        ascending=[True, True, True, False],
        na_position="last",
    ).reset_index(drop=True)


def build_regulatory_driver_map(
    *,
    driver_metadata: pd.DataFrame,
    coverage: pd.DataFrame,
    age_associations: pd.DataFrame,
    ora_correlations: pd.DataFrame,
) -> pd.DataFrame:
    """Build one ranked summary row per curated regulatory driver."""

    rows = []
    coverage_lookup = {str(row["driver_set"]): row.to_dict() for _, row in coverage.iterrows()}
    for _, driver in driver_metadata.iterrows():
        driver_set = str(driver["driver_set"])
        age = age_associations[
            age_associations["driver_set"].eq(driver_set)
            & age_associations["analysis_set"].eq("primary")
        ].copy()
        corr = ora_correlations[
            ora_correlations["driver_set"].eq(driver_set)
            & ora_correlations["score_metric"].eq("oraa")
        ].copy()
        top_age = age.sort_values(["fdr", "abs_beta_per_10_years"], ascending=[True, False]).head(1)
        top_corr = corr.sort_values(["fdr", "abs_pearson_r"], ascending=[True, False]).head(1)
        cov = coverage_lookup.get(driver_set, {})
        rows.append(
            {
                "driver_set": driver_set,
                "driver": driver["driver"],
                "driver_class": driver["driver_class"],
                "target_theme": driver["target_theme"],
                "expected_age_direction": driver["expected_age_direction"],
                "description": driver["description"],
                "source": driver["source"],
                "citation": driver["citation"],
                "n_requested": cov.get("n_requested", np.nan),
                "n_present": cov.get("n_present", np.nan),
                "coverage_fraction": cov.get("coverage_fraction", np.nan),
                "top_age_lineage_group": _first_value(top_age, "lineage_group"),
                "top_age_beta_per_10_years": _first_value(top_age, "beta_per_10_years"),
                "top_age_fdr": _first_value(top_age, "fdr"),
                "top_age_direction": _first_value(top_age, "direction"),
                "top_oraa_lineage_group": _first_value(top_corr, "lineage_group"),
                "top_oraa_pearson_r": _first_value(top_corr, "pearson_r"),
                "top_oraa_fdr": _first_value(top_corr, "fdr"),
                "driver_priority_score": _priority_score(top_age, top_corr),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["driver_priority_score", "coverage_fraction", "driver_set"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def write_regulatory_driver_figure(
    driver_map: pd.DataFrame,
    *,
    pdf_out: str | Path,
    png_out: str | Path | None = None,
) -> None:
    """Plot the highest-priority first-pass regulatory drivers."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402

    plot = driver_map.sort_values("driver_priority_score", ascending=False).head(12).iloc[::-1]
    labels = plot["driver"].astype(str)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.barh(labels, plot["driver_priority_score"], color="#4c78a8")
    ax.set_xlabel("First-pass priority score")
    ax.set_ylabel("Curated driver")
    ax.set_title("Regulatory-driver hypotheses from pseudobulk target programs")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(ensure_parent(pdf_out))
    if png_out is not None:
        fig.savefig(ensure_parent(png_out), dpi=200)
    plt.close(fig)


def render_regulatory_driver_feasibility(
    *,
    driver_map: pd.DataFrame,
    coverage: pd.DataFrame,
) -> str:
    """Render the M5.3 feasibility note."""

    top = driver_map.head(8)
    lines = [
        "# Regulatory Driver Feasibility",
        "",
        "Updated: 2026-06-25",
        "",
        "## Decision",
        "",
        "Primary local method: curated TF/pathway target programs scored from genomewide pseudobulk logCPM by donor and cell state. This satisfies the first-pass driver-hypothesis layer without requiring new package installs or per-cell regulon inference.",
        "",
        "Deferred heavy methods: decoupler with DoRothEA/CollecTRI or PROGENy, pySCENIC/SCENIC+, and chromatin/motif workflows. These remain useful upgrades, but they require additional databases/packages and stronger compute/runtime validation.",
        "",
        "## Local Tool Status",
        "",
        "- `scanpy`: installed.",
        "- `decoupler`, `pyscenic`, `gseapy`, `omnipath`: not installed in the current `.venv` at feasibility time.",
        "",
        "## Lineage Groups",
        "",
        "- HBC",
        "- GBC/INP",
        "- Immature OSN",
        "- Mature OSN",
        "- Sustentacular",
        "- Immune",
        "- Respiratory/secretory",
        "- All lineages",
        "",
        "## Gene Coverage",
        "",
        f"- Driver target sets: {coverage.shape[0]}",
        f"- Minimum coverage fraction: {coverage['coverage_fraction'].min():.3f}",
        f"- Mean coverage fraction: {coverage['coverage_fraction'].mean():.3f}",
        "",
        "## Top Ranked First-Pass Drivers",
        "",
        "| Driver | Theme | Coverage | Top age lineage | Top age FDR | Top ORAA lineage | Top ORAA FDR |",
        "| --- | --- | ---: | --- | ---: | --- | ---: |",
    ]
    for _, row in top.iterrows():
        lines.append(
            "| "
            f"{row['driver']} | {row['target_theme']} | {row['coverage_fraction']:.3f} | "
            f"{row['top_age_lineage_group']} | {_format_float(row['top_age_fdr'])} | "
            f"{row['top_oraa_lineage_group']} | {_format_float(row['top_oraa_fdr'])} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Guardrail",
            "",
            "These rows are ranked driver hypotheses from curated target programs. They are not causal perturbation evidence and should not be described as measured regulons until an external prior method such as decoupler/CollecTRI or SCENIC-style inference is run and sensitivity-checked.",
            "",
        ]
    )
    return "\n".join(lines)


def write_regulatory_driver_outputs(
    *,
    activity: pd.DataFrame,
    donor_activity: pd.DataFrame,
    coverage: pd.DataFrame,
    driver_map: pd.DataFrame,
    age_associations: pd.DataFrame,
    ora_correlations: pd.DataFrame,
    feasibility_note: str,
    activity_out: str | Path,
    donor_activity_out: str | Path,
    coverage_out: str | Path,
    driver_map_out: str | Path,
    age_out: str | Path,
    ora_out: str | Path,
    feasibility_out: str | Path,
    figure_pdf: str | Path,
    figure_png: str | Path | None = None,
) -> None:
    activity.to_csv(ensure_parent(activity_out), sep="\t", index=False)
    donor_activity.to_csv(ensure_parent(donor_activity_out), sep="\t", index=False)
    coverage.to_csv(ensure_parent(coverage_out), sep="\t", index=False)
    driver_map.to_csv(ensure_parent(driver_map_out), sep="\t", index=False)
    age_associations.to_csv(ensure_parent(age_out), sep="\t", index=False)
    ora_correlations.to_csv(ensure_parent(ora_out), sep="\t", index=False)
    ensure_parent(feasibility_out).write_text(feasibility_note + "\n", encoding="utf-8")
    write_regulatory_driver_figure(driver_map, pdf_out=figure_pdf, png_out=figure_png)


def _read_target_gene_rows(
    counts_path: str | Path,
    driver_metadata: pd.DataFrame,
    *,
    chunksize: int,
) -> pd.DataFrame:
    target_tokens = {
        normalize_token(gene)
        for genes in driver_metadata["genes"]
        for gene in genes
    }
    rows = []
    for chunk in pd.read_csv(counts_path, sep="\t", chunksize=chunksize):
        keep = chunk["gene_symbol"].map(lambda gene: normalize_token(gene) in target_tokens)
        if keep.any():
            rows.append(chunk.loc[keep].copy())
    if not rows:
        return pd.DataFrame()
    frame = pd.concat(rows, ignore_index=True).copy()
    frame["_token"] = frame["gene_symbol"].map(normalize_token)
    return frame.drop_duplicates("_token").drop(columns=["_token"]).reset_index(drop=True)


def _driver_gene_coverage(driver_metadata: pd.DataFrame, target_rows: pd.DataFrame) -> pd.DataFrame:
    present_tokens = set(target_rows["gene_symbol"].map(normalize_token)) if not target_rows.empty else set()
    rows = []
    for _, driver in driver_metadata.iterrows():
        genes = list(driver["genes"])
        present = [gene for gene in genes if normalize_token(gene) in present_tokens]
        missing = [gene for gene in genes if normalize_token(gene) not in present_tokens]
        rows.append(
            {
                "driver_set": driver["driver_set"],
                "driver": driver["driver"],
                "target_theme": driver["target_theme"],
                "n_requested": len(genes),
                "n_present": len(present),
                "coverage_fraction": len(present) / len(genes) if genes else np.nan,
                "present_genes": ",".join(present),
                "missing_genes": ",".join(missing),
            }
        )
    return pd.DataFrame(rows)


def _token_positions(symbols: pd.Series) -> dict[str, list[int]]:
    positions: dict[str, list[int]] = {}
    for idx, symbol in enumerate(symbols):
        positions.setdefault(normalize_token(symbol), []).append(idx)
    return positions


def _lineage_group(row: pd.Series) -> str:
    fine = normalize_token(row.get("fine_cell_type", ""))
    coarse = normalize_token(row.get("coarse_cell_type", ""))
    text = f"{coarse} {fine}"
    if "hbc" in text:
        return "hbc"
    if "inp" in text or "suprabasal" in text or "progenitor" in text:
        return "gbc_inp"
    if "iosn" in text or "immature" in text:
        return "immature_osn"
    if "mosn" in text or "mature" in text:
        return "mature_osn"
    if "sus" in text or "sustentacular" in text:
        return "sustentacular"
    if any(term in text for term in ["dendritic", "tcell", "bcell", "nk", "macrophage", "plasma", "cd4", "cd8"]):
        return "immune"
    if any(term in text for term in ["bowman", "gland", "goblet", "ciliated", "resp", "ionocyte", "tuft", "club", "mv"]):
        return "respiratory_secretory"
    return "other"


def _weighted_activity_row(keys: tuple, group: pd.DataFrame, group_cols: list[str]) -> dict[str, object]:
    weights = pd.to_numeric(group["n_cells"], errors="coerce").fillna(0.0).to_numpy(dtype=float)
    values = pd.to_numeric(group["activity_score"], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(values) & (weights > 0)
    score = np.average(values[valid], weights=weights[valid]) if valid.any() else np.nan
    row = dict(zip(group_cols, keys, strict=True))
    row["n_pseudobulk_groups"] = int(group.shape[0])
    row["n_cells"] = int(weights.sum())
    row["activity_score"] = score
    return row


def _ols_age_association(frame: pd.DataFrame, *, covariates: tuple[str, ...]) -> dict[str, object]:
    try:
        import statsmodels.api as sm
    except ModuleNotFoundError:
        return _not_tested_assoc(frame, "statsmodels_missing", [])

    work = frame.copy()
    work["age"] = pd.to_numeric(work["age"], errors="coerce")
    work["activity_score"] = pd.to_numeric(work["activity_score"], errors="coerce")
    covariate_cols = [col for col in covariates if col in work.columns]
    cov_frame, used = _design_covariates(work, covariate_cols)
    model_frame = pd.concat([work[["age", "activity_score"]], cov_frame], axis=1).dropna()
    n = int(model_frame.shape[0])
    if n < 8 or model_frame["activity_score"].nunique(dropna=True) < 2:
        return _not_tested_assoc(model_frame, "not_enough_data", used)
    x = model_frame.drop(columns=["activity_score"])
    x = sm.add_constant(x.astype(float), has_constant="add")
    if n <= x.shape[1] + 2:
        return _not_tested_assoc(model_frame, "too_few_degrees_of_freedom", used)
    try:
        fit = sm.OLS(model_frame["activity_score"].astype(float), x).fit()
    except Exception:
        return _not_tested_assoc(model_frame, "fit_failed", used)
    beta = float(fit.params["age"] * 10.0)
    se = float(fit.bse["age"] * 10.0)
    direction = "positive" if beta > 0 else "negative" if beta < 0 else "flat"
    return {
        "n": n,
        "beta_per_10_years": beta,
        "abs_beta_per_10_years": abs(beta),
        "standard_error": se,
        "p_value": float(fit.pvalues["age"]),
        "direction": direction,
        "status": "ok",
        "adjustment_covariates": ";".join(used),
    }


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
        pieces.append(pd.get_dummies(categorical, prefix=covariate, drop_first=True, dtype=float))
        used.append(covariate)
    if not pieces:
        return pd.DataFrame(index=frame.index), []
    return pd.concat(pieces, axis=1), used


def _not_tested_assoc(frame: pd.DataFrame, status: str, used: list[str]) -> dict[str, object]:
    return {
        "n": int(frame.shape[0]),
        "beta_per_10_years": np.nan,
        "abs_beta_per_10_years": np.nan,
        "standard_error": np.nan,
        "p_value": np.nan,
        "direction": "not_tested",
        "status": status,
        "adjustment_covariates": ";".join(used),
    }


def _pearson_association(x: pd.Series, y: pd.Series) -> dict[str, object]:
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    mask = x.notna() & y.notna()
    n = int(mask.sum())
    if n < 8 or x.loc[mask].nunique() < 2 or y.loc[mask].nunique() < 2:
        return {"n": n, "pearson_r": np.nan, "p_value": np.nan, "direction": "not_tested"}
    r_value, p_value = stats.pearsonr(x.loc[mask], y.loc[mask])
    direction = "positive" if r_value > 0 else "negative" if r_value < 0 else "flat"
    return {"n": n, "pearson_r": float(r_value), "p_value": float(p_value), "direction": direction}


def _summarize_ora_scores(ora_scores: pd.DataFrame) -> pd.DataFrame:
    scores = ora_scores.copy()
    if "model" in scores:
        scores = scores[~scores["model"].astype(str).eq("null_model")].copy()
    metrics = [col for col in ["ora", "oraa"] if col in scores.columns]
    return scores.groupby("donor_id", observed=True)[metrics].mean().reset_index()


def _driver_lookup(driver_metadata: pd.DataFrame) -> dict[str, dict[str, object]]:
    return {str(row["driver_set"]): row.to_dict() for _, row in driver_metadata.iterrows()}


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
        cumulative = min(cumulative, ranked[idx] * n / (idx + 1))
        adjusted[idx] = cumulative
    result = np.empty(n, dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def _first_value(frame: pd.DataFrame, column: str) -> object:
    if frame.empty or column not in frame:
        return np.nan
    return frame.iloc[0][column]


def _priority_score(age: pd.DataFrame, corr: pd.DataFrame) -> float:
    score = 0.0
    if not age.empty and pd.notna(age.iloc[0].get("fdr")):
        score += -np.log10(max(float(age.iloc[0]["fdr"]), 1e-300))
    if not corr.empty and pd.notna(corr.iloc[0].get("fdr")):
        score += -np.log10(max(float(corr.iloc[0]["fdr"]), 1e-300))
    return float(score)


def _format_float(value: object) -> str:
    number = pd.to_numeric(value, errors="coerce")
    if pd.isna(number):
        return ""
    if number < 0.001:
        return f"{number:.2e}"
    return f"{number:.3f}"


def _empty_activity(metadata: pd.DataFrame, driver_metadata: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, driver in driver_metadata.iterrows():
        frame = metadata.copy()
        frame["driver_set"] = driver["driver_set"]
        frame["driver"] = driver["driver"]
        frame["driver_class"] = driver["driver_class"]
        frame["target_theme"] = driver["target_theme"]
        frame["activity_score"] = np.nan
        frame["lineage_group"] = frame.apply(_lineage_group, axis=1)
        rows.append(frame)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _empty_donor_activity() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "driver_set",
            "driver",
            "driver_class",
            "target_theme",
            "donor_id",
            "lineage_group",
            "n_pseudobulk_groups",
            "n_cells",
            "activity_score",
        ]
    )
