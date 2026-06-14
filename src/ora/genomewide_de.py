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
    de["is_sex_linked_initial"] = de["gene_symbol"].astype(str).isin(SEX_LINKED_SYMBOLS)
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
