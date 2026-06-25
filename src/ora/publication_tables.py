"""Build compact manuscript-facing tables from ORA result artifacts."""

from __future__ import annotations

from datetime import date
from pathlib import Path

import pandas as pd


def build_publication_tables(tables_dir: str | Path) -> dict[str, pd.DataFrame]:
    """Return publication-ready summary tables from generated ORA outputs."""

    root = Path(tables_dir)
    return {
        "manuscript_table_cohort": _cohort_table(_read(root / "cohort_summary.tsv")),
        "manuscript_table_model_card": _model_card_table(_read(root / "ora_model_card.tsv")),
        "manuscript_table_external_validation_strength": _external_table(_read(root / "external_validation_evidence.tsv")),
        "manuscript_table_latent_neighborhood_gates": _latent_neighborhood_table(root),
        "manuscript_table_de_audit_summary": _de_audit_table(root),
        "manuscript_table_ndd_guardrails": _ndd_table(_read(root / "ndd_ora_projection_summary.tsv")),
    }


def render_publication_table_index(tables: dict[str, pd.DataFrame]) -> str:
    """Render a Markdown index describing the publication table bundle."""

    lines = [
        "# ORA Publication Table Bundle",
        "",
        f"Updated: {date.today().isoformat()}",
        "",
        "These compact tables are generated from the real ORA result artifacts and are intended for manuscript, extended-data, and supplement assembly.",
        "",
        "| Table | Rows | Purpose |",
        "| --- | ---: | --- |",
    ]
    purposes = {
        "manuscript_table_cohort": "Cohort and cell-count summary.",
        "manuscript_table_model_card": "Compact ORA model-card metrics and limitations.",
        "manuscript_table_external_validation_strength": "External validation readiness and claim strength.",
        "manuscript_table_latent_neighborhood_gates": "Full 4M scVI, Milo-style, edgeR, and MiloR claim gates.",
        "manuscript_table_de_audit_summary": "Genome-wide DE audit and sentinel summary.",
        "manuscript_table_ndd_guardrails": "NDD projection guardrails by model and disease group.",
    }
    for name, table in tables.items():
        lines.append(f"| `{name}.tsv` | {table.shape[0]} | {purposes.get(name, 'Publication summary table.')} |")
    lines.extend(
        [
            "",
            "Claim rule: tables summarize evidence strength; they do not upgrade guarded analyses to primary claims.",
            "",
        ]
    )
    return "\n".join(lines)


def _cohort_table(cohort: pd.DataFrame) -> pd.DataFrame:
    if cohort.empty:
        return pd.DataFrame(columns=["cohort", "donors", "cells", "samples", "publication_role"])
    keep = cohort.copy()
    for col in ["donors", "cells", "samples"]:
        if col in keep:
            keep[col] = pd.to_numeric(keep[col], errors="coerce")
    keep["publication_role"] = keep["cohort"].map(
        {
            "healthy": "ORA training and healthy-aging primary claim",
            "ad": "held-out exploratory NDD projection",
            "pd": "held-out exploratory NDD projection",
        }
    ).fillna("context")
    columns = [col for col in ["cohort", "donors", "cells", "samples", "publication_role"] if col in keep]
    return keep[columns].sort_values("cohort").reset_index(drop=True)


def _model_card_table(model_card: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "model",
        "feature_set",
        "role",
        "n",
        "repeats",
        "backend",
        "backend_package",
        "backend_version",
        "fallback_used",
        "mae_mean",
        "mae_ci_low",
        "mae_ci_high",
        "spearman_r_mean",
        "calibration_slope",
        "permutation_p_mae",
        "limitations",
    ]
    if model_card.empty:
        return pd.DataFrame(columns=columns)
    table = model_card[[col for col in columns if col in model_card]].copy()
    role_order = {"preferred_benchmark": 0, "secondary_benchmark": 1, "interpretable_baseline": 2, "negative_control": 3}
    table["_role_order"] = table["role"].map(role_order).fillna(9)
    table["mae_mean"] = pd.to_numeric(table["mae_mean"], errors="coerce")
    table = table.sort_values(["_role_order", "feature_set", "mae_mean"]).drop(columns="_role_order")
    return table.head(12).reset_index(drop=True)


def _external_table(evidence: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "dataset_id",
        "accession",
        "evidence_type",
        "feature_level",
        "readiness_class",
        "validation_strength",
        "n_samples",
        "n_donors",
        "supports_primary_claim",
        "limitation",
        "next_action",
    ]
    if evidence.empty:
        return pd.DataFrame(columns=columns)
    table = evidence[[col for col in columns if col in evidence]].copy()
    table["supports_primary_claim"] = table["supports_primary_claim"].fillna("no")
    return table.sort_values(["dataset_id", "evidence_type"]).reset_index(drop=True)


def _latent_neighborhood_table(root: Path) -> pd.DataFrame:
    rows = []
    scvi = _read(root / "scvi_embedding_claim_gates.tsv")
    for _, row in scvi.iterrows():
        rows.append(
            {
                "analysis": row.get("model"),
                "category": "scVI embedding",
                "scale": row.get("role"),
                "primary_metric": f"{_fmt(row.get('cells'))} cells; fine purity {_fmt(row.get('fine_label_purity'), 3)}",
                "claim_gate": row.get("claim_gate"),
                "interpretation": "Full 4M model is primary; smaller models are sensitivity anchors.",
            }
        )
    summary_specs = [
        ("python_lineage_full", "Milo-style DA", root / "milo_full_4m_lineage_summary.tsv", "age_fdr_lt_0_10", "neighborhoods_tested", "secondary"),
        ("python_lineage_matched", "Milo-style DA", root / "milo_full_4m_lineage_matched_summary.tsv", "age_fdr_lt_0_10", "neighborhoods_tested", "guarded_exact_neighborhood"),
        ("edger_lineage_full", "edgeR parity", root / "milo_full_4m_lineage_edger_parity_summary.tsv", "edger_fdr_lt_0_10", "neighborhoods_compared", "directionality_sensitivity"),
        ("edger_lineage_matched", "edgeR parity", root / "milo_full_4m_lineage_matched_edger_parity_summary.tsv", "edger_fdr_lt_0_10", "neighborhoods_compared", "directionality_sensitivity"),
        ("milor_lineage_subset", "official MiloR subset", root / "milor_lineage_subset_summary.tsv", "fdr_lt_0_10", "neighborhoods", "implementation_sensitivity"),
        ("milor_lineage_matched_subset", "official MiloR subset", root / "milor_lineage_matched_subset_summary.tsv", "fdr_lt_0_10", "neighborhoods", "implementation_sensitivity"),
    ]
    for name, category, path, sig_metric, tested_metric, gate in summary_specs:
        table = _read(path)
        sig = _metric(table, sig_metric)
        tested = _metric(table, tested_metric)
        rows.append(
            {
                "analysis": name,
                "category": category,
                "scale": _scale_label(name),
                "primary_metric": f"{_fmt(sig)} / {_fmt(tested)} FDR<0.10 neighborhoods",
                "claim_gate": gate,
                "interpretation": _latent_interpretation(name),
            }
        )
    return pd.DataFrame(rows)


def _de_audit_table(root: Path) -> pd.DataFrame:
    rows = []
    for label, path in [
        ("edgeR_all", root / "pseudobulk_genomewide_de_audit.tsv"),
        ("edgeR_matched", root / "pseudobulk_genomewide_de_audit_matched_flex_v2_device.tsv"),
        ("limma_all", root / "pseudobulk_genomewide_limma_voom_de_audit.tsv"),
        ("limma_matched", root / "pseudobulk_genomewide_limma_voom_de_audit_matched_flex_v2_device.tsv"),
    ]:
        table = _read(path)
        for _, row in table.iterrows():
            rows.append(
                {
                    "engine_context": label,
                    "contrast": row.get("contrast"),
                    "tested_rows": row.get("tested_rows"),
                    "significant_rows": row.get("significant_rows"),
                    "sex_linked_significant_rows": row.get("is_sex_linked_initial_significant_rows"),
                    "hemoglobin_significant_rows": row.get("is_hemoglobin_significant_rows"),
                    "immunoglobulin_significant_rows": row.get("is_immunoglobulin_significant_rows"),
                    "claim_gate": "hypothesis_generation_with_audit",
                }
            )
    return pd.DataFrame(rows)


def _ndd_table(ndd: pd.DataFrame) -> pd.DataFrame:
    columns = ["model", "disease_group", "donors", "mean_age", "mean_ora", "mean_oraa", "sd_oraa", "claim_gate"]
    if ndd.empty:
        return pd.DataFrame(columns=columns)
    table = ndd[ndd["disease_group"].isin(["ad", "pd"])].copy()
    table["claim_gate"] = "exploratory_small_n_flex_v2_device"
    return table[[col for col in columns if col in table]].sort_values(["model", "disease_group"]).reset_index(drop=True)


def _read(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def _metric(summary: pd.DataFrame, metric: str) -> float:
    if summary.empty or "metric" not in summary or "value" not in summary:
        return 0.0
    rows = summary[summary["metric"].astype(str).eq(metric)]
    if rows.empty:
        return 0.0
    return float(pd.to_numeric(rows["value"], errors="coerce").fillna(0).iloc[0])


def _fmt(value: object, digits: int = 0) -> str:
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(number):
        return "NA"
    if digits:
        return f"{float(number):.{digits}f}"
    return f"{int(round(float(number))):,}"


def _scale_label(name: str) -> str:
    if "matched" in name:
        return "matched FLEX v2/device"
    if "subset" in name:
        return "stratified subset"
    return "all healthy donors"


def _latent_interpretation(name: str) -> str:
    if name == "python_lineage_matched":
        return "Single matched Early iOSN exact-neighborhood result; keep narrow."
    if name.startswith("milor"):
        return "Official MiloR supports broad lineage-neighborhood structure but not dominant Early iOSN replication."
    if name.startswith("edger"):
        return "Same-neighborhood count-model parity supports signed direction."
    return "Full-scale Milo-style map supports broad secondary lineage-neighborhood remodeling."
