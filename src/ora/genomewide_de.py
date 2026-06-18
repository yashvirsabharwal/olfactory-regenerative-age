"""Summaries for genome-wide pseudobulk DE results."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


SEX_LINKED_SYMBOLS = {
    "DDX3Y",
    "EIF1AY",
    "KDM5D",
    "NLGN4Y",
    "PRKY",
    "RPS4Y1",
    "TBL1Y",
    "TMSB4Y",
    "TTTY15",
    "TXLNGY",
    "USP9Y",
    "UTY",
    "XIST",
    "ZFY",
}

AUDIT_FLAG_COLUMNS = [
    "is_sex_linked_initial",
    "is_mitochondrial",
    "is_ribosomal",
    "is_hemoglobin",
    "is_immunoglobulin",
]


def summarize_genomewide_de(
    de_path: str | Path,
    summary_path: str | Path,
    *,
    fdr_threshold: float = 0.05,
    top_n: int = 100,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return compact summary and top-hit tables for genome-wide DE output."""

    de = pd.read_csv(de_path, sep="\t")
    run_summary = pd.read_csv(summary_path, sep="\t")
    de["fdr"] = pd.to_numeric(de["fdr"], errors="coerce")
    de["p_value"] = pd.to_numeric(de["p_value"], errors="coerce")
    de["log2fc"] = pd.to_numeric(de["log2fc"], errors="coerce")
    de = add_gene_audit_flags(de)
    ok_summary = run_summary[run_summary["status"].eq("ok")].copy() if "status" in run_summary else run_summary.copy()
    sig = de[de["fdr"].lt(fdr_threshold)].copy()
    rows = []
    for contrast, frame in de.groupby("contrast", observed=True):
        sig_frame = sig[sig["contrast"].eq(contrast)]
        rows.append(
            {
                "contrast": contrast,
                "tested_rows": int(frame.shape[0]),
                "tested_genes": int(frame["gene_symbol"].nunique()),
                "tested_cell_states": int(frame["fine_cell_type"].nunique()),
                "ok_cell_state_models": int(ok_summary[ok_summary["contrast"].eq(contrast)]["fine_cell_type"].nunique())
                if "contrast" in ok_summary
                else pd.NA,
                "fdr_threshold": fdr_threshold,
                "significant_rows": int(sig_frame.shape[0]),
                "significant_genes": int(sig_frame["gene_symbol"].nunique()),
                "significant_cell_states": int(sig_frame["fine_cell_type"].nunique()),
                "sex_linked_significant_rows": int(sig_frame["is_sex_linked_initial"].sum()),
            }
        )
    top_hits = de.sort_values(["fdr", "p_value", "contrast", "fine_cell_type", "gene_symbol"]).head(top_n).reset_index(drop=True)
    return pd.DataFrame(rows), top_hits


def add_gene_audit_flags(de: pd.DataFrame) -> pd.DataFrame:
    """Add sentinel gene-category flags used to triage genome-wide DE hits."""

    frame = de.copy()
    symbols = frame.get("gene_symbol", pd.Series("", index=frame.index)).fillna("").astype(str).str.upper()
    frame["is_sex_linked_initial"] = symbols.isin(SEX_LINKED_SYMBOLS)
    frame["is_mitochondrial"] = symbols.str.startswith(("MT-", "MT."))
    frame["is_ribosomal"] = symbols.str.match(r"^(RPL|RPS|MRPL|MRPS)[0-9A-Z-]*$", na=False)
    frame["is_hemoglobin"] = symbols.str.match(r"^HB[ABDEGMQZ][0-9A-Z-]*$", na=False)
    frame["is_immunoglobulin"] = symbols.str.match(r"^(IGH|IGK|IGL)[A-Z0-9-]*$", na=False)
    return frame


def audit_genomewide_de(
    de_path: str | Path,
    run_summary_path: str | Path,
    metadata_path: str | Path | None = None,
    manifest_path: str | Path | None = None,
    *,
    fdr_threshold: float = 0.05,
    min_case_donors: int = 3,
    min_control_donors: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Summarize donor balance, sentinel hit classes, and matched-DE feasibility."""

    de = add_gene_audit_flags(pd.read_csv(de_path, sep="\t"))
    run_summary = pd.read_csv(run_summary_path, sep="\t")
    de["fdr"] = pd.to_numeric(de.get("fdr"), errors="coerce")
    significant = de[de["fdr"].lt(fdr_threshold)].copy()
    audit_rows = []
    for contrast, frame in de.groupby("contrast", observed=True):
        sig = significant[significant["contrast"].eq(contrast)]
        row = {
            "contrast": contrast,
            "tested_rows": int(frame.shape[0]),
            "significant_rows": int(sig.shape[0]),
            "fdr_threshold": fdr_threshold,
        }
        for flag in AUDIT_FLAG_COLUMNS:
            row[f"{flag}_significant_rows"] = int(sig.get(flag, pd.Series(dtype=bool)).astype(bool).sum())
            row[f"{flag}_significant_genes"] = (
                int(sig.loc[sig.get(flag, pd.Series(False, index=sig.index)).astype(bool), "gene_symbol"].nunique())
                if "gene_symbol" in sig
                else 0
            )
        audit_rows.append(row)

    donor_balance = summarize_donor_balance(
        run_summary,
        min_case_donors=min_case_donors,
        min_control_donors=min_control_donors,
    )
    matched = (
        matched_de_feasibility(
            metadata_path,
            manifest_path=manifest_path,
            min_case_donors=min_case_donors,
            min_control_donors=min_case_donors,
        )
        if metadata_path
        else pd.DataFrame()
    )
    return pd.DataFrame(audit_rows), donor_balance, matched


def summarize_donor_balance(
    run_summary: pd.DataFrame,
    *,
    min_case_donors: int = 3,
    min_control_donors: int = 10,
) -> pd.DataFrame:
    """Classify per-contrast/cell-state donor balance before interpreting DE."""

    if run_summary is None or run_summary.empty:
        return pd.DataFrame()
    frame = run_summary.copy()
    frame["n_case"] = pd.to_numeric(frame.get("n_case"), errors="coerce").fillna(0).astype(int)
    frame["n_control"] = pd.to_numeric(frame.get("n_control"), errors="coerce").fillna(0).astype(int)
    frame["balance_status"] = "ok"
    frame.loc[frame["n_case"].lt(min_case_donors), "balance_status"] = "low_case_donors"
    frame.loc[frame["n_control"].lt(min_control_donors), "balance_status"] = "low_control_donors"
    frame.loc[
        frame["n_case"].lt(min_case_donors) & frame["n_control"].lt(min_control_donors),
        "balance_status",
    ] = "low_case_and_control_donors"
    columns = ["contrast", "fine_cell_type", "n_case", "n_control", "n_genes_tested", "status", "balance_status"]
    for col in columns:
        if col not in frame:
            frame[col] = pd.NA
    return frame[columns].sort_values(["contrast", "balance_status", "fine_cell_type"]).reset_index(drop=True)


def matched_de_feasibility(
    metadata_path: str | Path,
    *,
    manifest_path: str | Path | None = None,
    chemistry: str = "flex_v2",
    collection_method: str = "device",
    min_case_donors: int = 3,
    min_control_donors: int = 3,
) -> pd.DataFrame:
    """Report cell states with enough matched FLEX/device donors for DE reruns."""

    metadata = pd.read_csv(metadata_path, sep="\t")
    if manifest_path and Path(manifest_path).exists():
        manifest_cols = ["donor_id", "chemistry", "collection_method"]
        manifest = pd.read_csv(manifest_path, sep="\t")
        available = [col for col in manifest_cols if col in manifest.columns]
        if set(available) == set(manifest_cols):
            donor_manifest = manifest[manifest_cols].drop_duplicates("donor_id")
            metadata = metadata.drop(columns=[col for col in ["chemistry", "collection_method"] if col in metadata])
            metadata = metadata.merge(donor_manifest, on="donor_id", how="left")
    required = {"donor_id", "disease_group", "fine_cell_type", "chemistry", "collection_method"}
    if metadata.empty or not required.issubset(metadata.columns):
        return pd.DataFrame()
    frame = metadata[
        metadata["chemistry"].astype(str).eq(chemistry)
        & metadata["collection_method"].astype(str).eq(collection_method)
        & metadata["disease_group"].astype(str).isin(["healthy", "ad", "pd"])
    ].copy()
    rows = []
    for (state, disease), group in frame[frame["disease_group"].isin(["ad", "pd"])].groupby(
        ["fine_cell_type", "disease_group"],
        observed=True,
    ):
        healthy = frame[frame["fine_cell_type"].eq(state) & frame["disease_group"].eq("healthy")]
        n_case = int(group["donor_id"].nunique())
        n_control = int(healthy["donor_id"].nunique())
        rows.append(
            {
                "contrast": f"{disease}_vs_healthy",
                "fine_cell_type": state,
                "chemistry": chemistry,
                "collection_method": collection_method,
                "n_case": n_case,
                "n_matched_healthy": n_control,
                "ready_for_matched_de": bool(n_case >= min_case_donors and n_control >= min_control_donors),
            }
        )
    return pd.DataFrame(rows).sort_values(["contrast", "ready_for_matched_de", "fine_cell_type"], ascending=[True, False, True])
