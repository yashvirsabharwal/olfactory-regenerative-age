"""Curated ligand-receptor niche signaling from pseudobulk expression."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .regeneration_modules import DEFAULT_ADJUSTMENT_COVARIATES, DEFAULT_ANALYSIS_SETS
from .regulatory_drivers import (
    _benjamini_hochberg,
    _first_value,
    _format_float,
    _ols_age_association,
    _pearson_association,
    _summarize_ora_scores,
)
from .utils import ensure_parent, normalize_token


def parse_niche_interactions(config: dict) -> pd.DataFrame:
    """Parse curated ligand-receptor interaction metadata."""

    rows = []
    for interaction_id, spec in config.get("interactions", {}).items():
        ligand_genes = tuple(str(gene) for gene in spec.get("ligand_genes", []) if str(gene).strip())
        receptor_genes = tuple(str(gene) for gene in spec.get("receptor_genes", []) if str(gene).strip())
        rows.append(
            {
                "interaction_id": str(interaction_id),
                "family": str(spec.get("family", interaction_id)),
                "ligand_genes": ligand_genes,
                "receptor_genes": receptor_genes,
                "sender_groups": tuple(str(group) for group in spec.get("sender_groups", [])),
                "receiver_groups": tuple(str(group) for group in spec.get("receiver_groups", [])),
                "expected_age_direction": str(spec.get("expected_age_direction", "unknown")),
                "source": str(spec.get("source", "")),
                "citation": str(spec.get("citation", "")),
                "rationale": str(spec.get("rationale", "")),
            }
        )
    return pd.DataFrame(rows)


def score_niche_interactions(
    *,
    counts_path: str | Path,
    metadata: pd.DataFrame,
    interactions: pd.DataFrame,
    chunksize: int = 1000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score donor-level sender ligand x receiver receptor interaction activity."""

    metadata = metadata.copy()
    pseudobulk_ids = metadata["pseudobulk_id"].astype(str).tolist()
    metadata["niche_group"] = metadata.apply(_niche_group, axis=1)
    target_rows = _read_target_gene_rows(counts_path, interactions, chunksize=chunksize)
    coverage = _interaction_gene_coverage(interactions, target_rows)
    if target_rows.empty:
        return pd.DataFrame(), coverage

    library = metadata.set_index("pseudobulk_id").loc[pseudobulk_ids, "sum_n_counts"].to_numpy(dtype=float)
    library = np.where(library > 0, library, np.nan)
    values = target_rows[pseudobulk_ids].to_numpy(dtype=float)
    log_cpm = np.log1p((values / library.reshape(1, -1)) * 1_000_000.0)
    token_to_positions = _token_positions(target_rows["gene_symbol"])

    rows = []
    base = metadata[
        [
            "pseudobulk_id",
            "donor_id",
            "sample_id",
            "disease_group",
            "coarse_cell_type",
            "fine_cell_type",
            "niche_group",
            "n_cells",
        ]
    ].copy()
    for _, interaction in interactions.iterrows():
        ligand_score = _gene_set_score(log_cpm, token_to_positions, interaction["ligand_genes"])
        receptor_score = _gene_set_score(log_cpm, token_to_positions, interaction["receptor_genes"])
        frame = base.copy()
        frame["ligand_score"] = ligand_score
        frame["receptor_score"] = receptor_score
        grouped = _aggregate_scores_by_donor_group(frame)
        for donor_id, donor_group in grouped.groupby("donor_id", observed=True):
            by_group = {str(row["niche_group"]): row for _, row in donor_group.iterrows()}
            for sender_group in interaction["sender_groups"]:
                sender = by_group.get(str(sender_group))
                if sender is None:
                    continue
                for receiver_group in interaction["receiver_groups"]:
                    receiver = by_group.get(str(receiver_group))
                    if receiver is None:
                        continue
                    sender_ligand = float(sender["ligand_score"])
                    receiver_receptor = float(receiver["receptor_score"])
                    rows.append(
                        {
                            "interaction_id": interaction["interaction_id"],
                            "family": interaction["family"],
                            "ligand_genes": ",".join(interaction["ligand_genes"]),
                            "receptor_genes": ",".join(interaction["receptor_genes"]),
                            "sender_group": str(sender_group),
                            "receiver_group": str(receiver_group),
                            "donor_id": donor_id,
                            "sender_ligand_score": sender_ligand,
                            "receiver_receptor_score": receiver_receptor,
                            "interaction_score": float(np.sqrt(max(sender_ligand, 0.0) * max(receiver_receptor, 0.0))),
                            "interaction_min_score": float(min(sender_ligand, receiver_receptor)),
                            "n_sender_cells": int(sender["n_cells"]),
                            "n_receiver_cells": int(receiver["n_cells"]),
                            "sender_pseudobulk_groups": int(sender["n_pseudobulk_groups"]),
                            "receiver_pseudobulk_groups": int(receiver["n_pseudobulk_groups"]),
                            "expected_age_direction": interaction["expected_age_direction"],
                            "source": interaction["source"],
                            "citation": interaction["citation"],
                            "rationale": interaction["rationale"],
                        }
                    )
    return pd.DataFrame(rows), coverage


def build_niche_age_associations(
    *,
    donor_scores: pd.DataFrame,
    manifest: pd.DataFrame,
    interactions: pd.DataFrame,
    analysis_sets: tuple[tuple[str, str], ...] = DEFAULT_ANALYSIS_SETS,
    covariates: tuple[str, ...] = DEFAULT_ADJUSTMENT_COVARIATES,
) -> pd.DataFrame:
    """Test donor-level interaction scores against age."""

    metadata = _interaction_lookup(interactions)
    rows = []
    group_cols = ["interaction_id", "sender_group", "receiver_group"]
    for analysis_set, mask_col in analysis_sets:
        if mask_col not in manifest:
            continue
        for keys, group in donor_scores.groupby(group_cols, observed=True):
            merged = group.merge(manifest, on="donor_id", how="left")
            merged = merged[merged[mask_col].fillna(False).astype(bool)].copy()
            merged["activity_score"] = pd.to_numeric(merged["interaction_score"], errors="coerce")
            assoc = _ols_age_association(merged, covariates=covariates)
            meta = metadata.get(str(keys[0]), {})
            rows.append(
                {
                    "analysis_set": analysis_set,
                    "interaction_id": keys[0],
                    "family": meta.get("family", keys[0]),
                    "sender_group": keys[1],
                    "receiver_group": keys[2],
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
    if result.empty:
        return result
    result["fdr"] = np.nan
    for _, idx in result.groupby(["analysis_set", "receiver_group"], observed=True).groups.items():
        p_values = pd.to_numeric(result.loc[idx, "p_value"], errors="coerce")
        ok = p_values.notna()
        if ok.any():
            result.loc[p_values.index[ok], "fdr"] = _benjamini_hochberg(p_values.loc[ok].to_numpy(dtype=float))
    result["observed_vs_expected"] = result.apply(
        lambda row: _observed_vs_expected(
            row.get("expected_age_direction", "unknown"),
            row.get("direction", "not_tested"),
            row.get("fdr", np.nan),
        ),
        axis=1,
    )
    return result.sort_values(
        ["analysis_set", "receiver_group", "fdr", "abs_beta_per_10_years"],
        ascending=[True, True, True, False],
        na_position="last",
    ).reset_index(drop=True)


def build_niche_ora_associations(
    *,
    donor_scores: pd.DataFrame,
    ora_scores: pd.DataFrame,
) -> pd.DataFrame:
    """Correlate donor-level niche scores with averaged ORA/ORAA scores."""

    scores = _summarize_ora_scores(ora_scores)
    rows = []
    group_cols = ["interaction_id", "family", "sender_group", "receiver_group"]
    for keys, group in donor_scores.groupby(group_cols, observed=True):
        merged = group.merge(scores, on="donor_id", how="inner")
        for metric in ["ora", "oraa"]:
            if metric not in merged:
                continue
            assoc = _pearson_association(merged["interaction_score"], merged[metric])
            rows.append(
                {
                    "interaction_id": keys[0],
                    "family": keys[1],
                    "sender_group": keys[2],
                    "receiver_group": keys[3],
                    "score_metric": metric,
                    **assoc,
                }
            )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr"] = np.nan
    for _, idx in result.groupby(["score_metric", "receiver_group"], observed=True).groups.items():
        p_values = pd.to_numeric(result.loc[idx, "p_value"], errors="coerce")
        ok = p_values.notna()
        if ok.any():
            result.loc[p_values.index[ok], "fdr"] = _benjamini_hochberg(p_values.loc[ok].to_numpy(dtype=float))
    result["abs_pearson_r"] = pd.to_numeric(result["pearson_r"], errors="coerce").abs()
    return result.sort_values(
        ["score_metric", "receiver_group", "fdr", "abs_pearson_r"],
        ascending=[True, True, True, False],
        na_position="last",
    ).reset_index(drop=True)


def build_niche_priority_table(
    *,
    interactions: pd.DataFrame,
    coverage: pd.DataFrame,
    age_associations: pd.DataFrame,
    ora_associations: pd.DataFrame,
) -> pd.DataFrame:
    """Build ranked niche-driver hypotheses per sender/receiver pair."""

    metadata = _interaction_lookup(interactions)
    coverage_lookup = {str(row["interaction_id"]): row.to_dict() for _, row in coverage.iterrows()}
    pairs = (
        age_associations[["interaction_id", "family", "sender_group", "receiver_group"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    rows = []
    for _, pair in pairs.iterrows():
        interaction_id = str(pair["interaction_id"])
        sender = str(pair["sender_group"])
        receiver = str(pair["receiver_group"])
        age = age_associations[
            age_associations["interaction_id"].eq(interaction_id)
            & age_associations["sender_group"].eq(sender)
            & age_associations["receiver_group"].eq(receiver)
            & age_associations["analysis_set"].eq("primary")
        ].copy()
        corr = ora_associations[
            ora_associations["interaction_id"].eq(interaction_id)
            & ora_associations["sender_group"].eq(sender)
            & ora_associations["receiver_group"].eq(receiver)
            & ora_associations["score_metric"].eq("oraa")
        ].copy()
        top_age = age.sort_values(["fdr", "abs_beta_per_10_years"], ascending=[True, False]).head(1)
        top_corr = corr.sort_values(["fdr", "abs_pearson_r"], ascending=[True, False]).head(1)
        cov = coverage_lookup.get(interaction_id, {})
        meta = metadata.get(interaction_id, {})
        rows.append(
            {
                "interaction_id": interaction_id,
                "family": meta.get("family", pair["family"]),
                "ligand_genes": ",".join(meta.get("ligand_genes", ())),
                "receptor_genes": ",".join(meta.get("receptor_genes", ())),
                "sender_group": sender,
                "receiver_group": receiver,
                "expected_age_direction": meta.get("expected_age_direction", "unknown"),
                "coverage_fraction": cov.get("coverage_fraction", np.nan),
                "top_age_beta_per_10_years": _first_value(top_age, "beta_per_10_years"),
                "top_age_fdr": _first_value(top_age, "fdr"),
                "top_age_direction": _first_value(top_age, "direction"),
                "top_oraa_pearson_r": _first_value(top_corr, "pearson_r"),
                "top_oraa_fdr": _first_value(top_corr, "fdr"),
                "lineage_neighborhood_support": _lineage_neighborhood_support(receiver),
                "priority_score": _niche_priority_score(top_age, top_corr, cov.get("coverage_fraction", np.nan), receiver),
                "rationale": meta.get("rationale", ""),
                "claim_status": "hypothesis_only",
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["priority_score", "coverage_fraction", "interaction_id", "sender_group", "receiver_group"],
        ascending=[False, False, True, True, True],
    ).reset_index(drop=True)


def render_niche_signaling_report(
    *,
    interactions: pd.DataFrame,
    coverage: pd.DataFrame,
    priority: pd.DataFrame,
) -> str:
    """Render a concise M5.4 feasibility/interpretation note."""

    top = priority.head(10)
    lines = [
        "# Niche Signaling Feasibility",
        "",
        "Updated: 2026-06-25",
        "",
        "## Decision",
        "",
        "Primary local method: curated ligand-receptor pseudobulk scoring by donor, sender group, and receiver group. This is a transparent hypothesis layer, not a full LIANA/NicheNet/CellChat/CellPhoneDB inference run.",
        "",
        "## Scope",
        "",
        f"- Curated interaction families: {interactions['family'].nunique()}",
        f"- Interaction definitions: {interactions.shape[0]}",
        f"- Minimum gene coverage: {coverage['coverage_fraction'].min():.3f}",
        f"- Mean gene coverage: {coverage['coverage_fraction'].mean():.3f}",
        "",
        "## Top Ranked Niche Hypotheses",
        "",
        "| Family | Sender | Receiver | Coverage | Age beta/10y | Age FDR | ORAA r | ORAA FDR |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for _, row in top.iterrows():
        lines.append(
            "| "
            f"{row['family']} | {row['sender_group']} | {row['receiver_group']} | "
            f"{_format_float(row['coverage_fraction'])} | {_format_float(row['top_age_beta_per_10_years'])} | "
            f"{_format_float(row['top_age_fdr'])} | {_format_float(row['top_oraa_pearson_r'])} | "
            f"{_format_float(row['top_oraa_fdr'])} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Guardrail",
            "",
            "These results nominate sender-receiver signaling hypotheses. They do not prove physical proximity, receptor activation, or causal niche control. Spatial transcriptomics, perturbation, or a full LR method with permutation/background testing would be required to promote the claims.",
            "",
        ]
    )
    return "\n".join(lines)


def write_niche_signaling_outputs(
    *,
    donor_scores: pd.DataFrame,
    coverage: pd.DataFrame,
    age_associations: pd.DataFrame,
    ora_associations: pd.DataFrame,
    priority: pd.DataFrame,
    report: str,
    donor_scores_out: str | Path,
    coverage_out: str | Path,
    age_out: str | Path,
    ora_out: str | Path,
    priority_out: str | Path,
    report_out: str | Path,
    figure_pdf: str | Path,
    figure_png: str | Path | None = None,
) -> None:
    donor_scores.to_csv(ensure_parent(donor_scores_out), sep="\t", index=False)
    coverage.to_csv(ensure_parent(coverage_out), sep="\t", index=False)
    age_associations.to_csv(ensure_parent(age_out), sep="\t", index=False)
    ora_associations.to_csv(ensure_parent(ora_out), sep="\t", index=False)
    priority.to_csv(ensure_parent(priority_out), sep="\t", index=False)
    ensure_parent(report_out).write_text(report + "\n", encoding="utf-8")
    write_niche_signaling_figure(priority, pdf_out=figure_pdf, png_out=figure_png)


def write_niche_signaling_figure(
    priority: pd.DataFrame,
    *,
    pdf_out: str | Path,
    png_out: str | Path | None = None,
) -> None:
    """Plot top niche hypotheses."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402

    plot = priority.head(12).iloc[::-1].copy()
    labels = plot.apply(
        lambda row: f"{row['family']}: {row['sender_group']}->{row['receiver_group']}",
        axis=1,
    )
    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    ax.barh(labels, plot["priority_score"], color="#527a3a")
    ax.set_xlabel("Niche hypothesis priority score")
    ax.set_ylabel("Curated ligand-receptor pair")
    ax.set_title("First-pass niche signaling hypotheses")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(ensure_parent(pdf_out))
    if png_out is not None:
        fig.savefig(ensure_parent(png_out), dpi=220)
    plt.close(fig)


def _read_target_gene_rows(
    counts_path: str | Path,
    interactions: pd.DataFrame,
    *,
    chunksize: int,
) -> pd.DataFrame:
    target_tokens = {
        normalize_token(gene)
        for _, row in interactions.iterrows()
        for genes in (row["ligand_genes"], row["receptor_genes"])
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


def _interaction_gene_coverage(interactions: pd.DataFrame, target_rows: pd.DataFrame) -> pd.DataFrame:
    present_tokens = set(target_rows["gene_symbol"].map(normalize_token)) if not target_rows.empty else set()
    rows = []
    for _, interaction in interactions.iterrows():
        ligands = list(interaction["ligand_genes"])
        receptors = list(interaction["receptor_genes"])
        ligand_present = [gene for gene in ligands if normalize_token(gene) in present_tokens]
        receptor_present = [gene for gene in receptors if normalize_token(gene) in present_tokens]
        n_requested = len(ligands) + len(receptors)
        n_present = len(ligand_present) + len(receptor_present)
        rows.append(
            {
                "interaction_id": interaction["interaction_id"],
                "family": interaction["family"],
                "n_ligands_requested": len(ligands),
                "n_ligands_present": len(ligand_present),
                "n_receptors_requested": len(receptors),
                "n_receptors_present": len(receptor_present),
                "n_requested": n_requested,
                "n_present": n_present,
                "coverage_fraction": n_present / n_requested if n_requested else np.nan,
                "present_ligands": ",".join(ligand_present),
                "missing_ligands": ",".join(gene for gene in ligands if gene not in ligand_present),
                "present_receptors": ",".join(receptor_present),
                "missing_receptors": ",".join(gene for gene in receptors if gene not in receptor_present),
            }
        )
    return pd.DataFrame(rows)


def _token_positions(symbols: pd.Series) -> dict[str, list[int]]:
    positions: dict[str, list[int]] = {}
    for idx, symbol in enumerate(symbols):
        positions.setdefault(normalize_token(symbol), []).append(idx)
    return positions


def _gene_set_score(log_cpm: np.ndarray, token_to_positions: dict[str, list[int]], genes: tuple[str, ...]) -> np.ndarray:
    positions = []
    for gene in genes:
        positions.extend(token_to_positions.get(normalize_token(gene), []))
    positions = sorted(set(positions))
    if not positions:
        return np.full(log_cpm.shape[1], np.nan)
    return np.nanmean(log_cpm[positions, :], axis=0)


def _aggregate_scores_by_donor_group(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (donor_id, niche_group), group in frame.groupby(["donor_id", "niche_group"], observed=True):
        weights = pd.to_numeric(group["n_cells"], errors="coerce").fillna(0.0).to_numpy(dtype=float)
        valid_weights = weights > 0
        row = {
            "donor_id": donor_id,
            "niche_group": niche_group,
            "n_pseudobulk_groups": int(group.shape[0]),
            "n_cells": int(weights.sum()),
        }
        for column in ("ligand_score", "receptor_score"):
            values = pd.to_numeric(group[column], errors="coerce").to_numpy(dtype=float)
            valid = np.isfinite(values) & valid_weights
            row[column] = float(np.average(values[valid], weights=weights[valid])) if valid.any() else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _niche_group(row: pd.Series) -> str:
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
    if any(term in text for term in ["dendritic", "tcell", "bcell", "nk", "macrophage", "plasma", "cd4", "cd8", "mono"]):
        return "immune"
    if any(term in text for term in ["bowman", "gland", "goblet", "ciliated", "resp", "ionocyte", "tuft", "club", "mv", "secretory"]):
        return "respiratory_secretory"
    return "other"


def _interaction_lookup(interactions: pd.DataFrame) -> dict[str, dict[str, object]]:
    return {str(row["interaction_id"]): row.to_dict() for _, row in interactions.iterrows()}


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


def _lineage_neighborhood_support(receiver_group: str) -> str:
    if receiver_group == "immature_osn":
        return "strongest_existing_milo_age_neighborhood_is_early_iosn"
    if receiver_group in {"hbc", "gbc_inp", "mature_osn"}:
        return "within_lineage_milo_scope"
    if receiver_group == "sustentacular":
        return "support_epithelium_not_lineage_neighborhood"
    return "not_in_lineage_scope"


def _niche_priority_score(age: pd.DataFrame, corr: pd.DataFrame, coverage: object, receiver_group: str) -> float:
    score = 0.0
    if not age.empty and pd.notna(age.iloc[0].get("fdr")):
        score += -np.log10(max(float(age.iloc[0]["fdr"]), 1e-300))
    if not corr.empty and pd.notna(corr.iloc[0].get("fdr")):
        score += -np.log10(max(float(corr.iloc[0]["fdr"]), 1e-300))
    cov = pd.to_numeric(coverage, errors="coerce")
    if pd.notna(cov):
        score += float(cov)
    if receiver_group in {"hbc", "gbc_inp", "immature_osn", "mature_osn"}:
        score += 0.5
    return float(score)
