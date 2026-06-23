#!/usr/bin/env python3
"""Reporting, package-check, and provenance command group."""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.manuscript import build_model_card
from ora.provenance import command_manifest_table, output_provenance_table
from ora.publication_tables import build_publication_tables, render_publication_table_index
from ora.utils import ensure_parent


EXPECTED_EXTENDED_FIGURES = [
    "extended_data_figure1_model_card.pdf",
    "extended_data_figure2_external_evidence.pdf",
    "extended_data_figure3_scvi_validation.pdf",
    "extended_data_figure4_de_audit.pdf",
    "extended_data_figure5_latent_robustness.pdf",
    "extended_data_figure6_ndd_guardrails.pdf",
]

EXPECTED_MANUSCRIPT_TABLES = [
    "manuscript_table_cohort.tsv",
    "manuscript_table_model_card.tsv",
    "manuscript_table_external_validation_strength.tsv",
    "manuscript_table_latent_neighborhood_gates.tsv",
    "manuscript_table_de_audit_summary.tsv",
    "manuscript_table_ndd_guardrails.tsv",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_model_card(subparsers)
    _add_publication_tables(subparsers)
    _add_manuscript_check(subparsers)
    _add_output_provenance(subparsers)
    args = parser.parse_args()
    args.func(args)


def _add_model_card(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("model-card")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--feature-set-comparison", default=None)
    parser.add_argument("--calibration", default=None)
    parser.add_argument("--permutation", default=None)
    parser.add_argument("--nested-tuning", default=None)
    parser.add_argument("--stacking", default=None)
    parser.add_argument("--out", default=None)
    parser.set_defaults(func=_model_card)


def _model_card(args: argparse.Namespace) -> None:
    config = load_config(args.gateway_config)
    outputs = config.get("outputs", {})
    paths = {
        "feature_set_comparison": args.feature_set_comparison
        or outputs.get("ora_feature_set_model_comparison_tsv", "results/tables/ora_feature_set_model_comparison.tsv"),
        "calibration": args.calibration or outputs.get("ora_calibration_tsv", "results/tables/ora_calibration.tsv"),
        "permutation": args.permutation or outputs.get("ora_permutation_empirical_tsv", "results/tables/ora_permutation_empirical.tsv"),
        "nested_tuning": args.nested_tuning or outputs.get("ora_nested_tuning_summary_tsv", "results/tables/ora_nested_tuning_summary.tsv"),
        "stacking": args.stacking or outputs.get("ora_stacking_summary_tsv", "results/tables/ora_stacking_summary.tsv"),
    }
    card = build_model_card(
        feature_set_comparison=_read_optional(paths["feature_set_comparison"]),
        calibration=_read_optional(paths["calibration"]),
        permutation=_read_optional(paths["permutation"]),
        nested_tuning=_read_optional(paths["nested_tuning"]),
        stacking=_read_optional(paths["stacking"]),
    )
    out_path = args.out or outputs.get("ora_model_card_tsv", "results/tables/ora_model_card.tsv")
    card.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(f"Wrote ORA model card: {out_path} ({card.shape[0]} rows)")


def _add_publication_tables(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("publication-tables")
    parser.add_argument("--tables-dir", default="results/tables")
    parser.add_argument("--out-dir", default="results/tables")
    parser.add_argument("--index-out", default="results/reports/publication_tables.md")
    parser.set_defaults(func=_publication_tables)


def _publication_tables(args: argparse.Namespace) -> None:
    tables = build_publication_tables(args.tables_dir)
    out_dir = Path(args.out_dir)
    for name, table in tables.items():
        table.to_csv(ensure_parent(out_dir / f"{name}.tsv"), sep="\t", index=False)
    Path(args.index_out).write_text(render_publication_table_index(tables), encoding="utf-8")
    print(f"Wrote {len(tables)} publication tables to {out_dir}")
    print(f"Wrote publication table index: {args.index_out}")


def _add_manuscript_check(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("manuscript-check")
    parser.add_argument("--manuscript", default="manuscript/main.tex")
    parser.add_argument("--bib", default="manuscript/references.bib")
    parser.add_argument("--figures-dir", default="results/figures")
    parser.add_argument("--tables-dir", default="results/tables")
    parser.add_argument("--out", default="results/reports/manuscript_package_check.tsv")
    parser.add_argument("--markdown-out", default="results/reports/manuscript_package_check.md")
    parser.set_defaults(func=_manuscript_check)


def _manuscript_check(args: argparse.Namespace) -> None:
    rows = check_package(
        manuscript=Path(args.manuscript),
        bib=Path(args.bib),
        figures_dir=Path(args.figures_dir),
        tables_dir=Path(args.tables_dir),
    )
    table = pd.DataFrame(rows, columns=["check", "status", "detail"])
    _write_tsv(table, Path(args.out))
    _write_text(render_markdown(table), Path(args.markdown_out))
    failures = table["status"].eq("fail").sum()
    blocked = table["status"].eq("blocked").sum()
    print(f"Wrote manuscript package check: {args.out} ({len(table)} checks; {failures} fail; {blocked} blocked)")
    print(f"Wrote manuscript package audit: {args.markdown_out}")
    if failures:
        raise SystemExit(1)


def _add_output_provenance(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("output-provenance")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--command-manifest", default="configs/command_manifest.yaml")
    parser.add_argument("--checksum-max-bytes", type=int, default=100 * 1024 * 1024)
    parser.add_argument("--command-out", default=None)
    parser.add_argument("--provenance-out", default=None)
    parser.set_defaults(func=_output_provenance)


def _output_provenance(args: argparse.Namespace) -> None:
    gateway_config = load_config(args.gateway_config)
    command_config = load_config(args.command_manifest)
    outputs = gateway_config.get("outputs", {})
    command_out = args.command_out or outputs.get("command_manifest_tsv", "results/reports/command_manifest.tsv")
    provenance_out = args.provenance_out or outputs.get("output_provenance_tsv", "results/reports/output_provenance.tsv")
    commands = command_manifest_table(command_config)
    provenance = output_provenance_table(command_config, checksum_max_bytes=args.checksum_max_bytes)
    commands.to_csv(ensure_parent(command_out), sep="\t", index=False)
    provenance.to_csv(ensure_parent(provenance_out), sep="\t", index=False)
    missing = int((~provenance["exists"].astype(bool)).sum()) if "exists" in provenance else 0
    print(f"Wrote command manifest: {command_out} ({commands.shape[0]} stages)")
    print(f"Wrote output provenance: {provenance_out} ({provenance.shape[0]} outputs; {missing} missing)")


def check_package(
    *,
    manuscript: Path,
    bib: Path,
    figures_dir: Path,
    tables_dir: Path,
) -> list[dict[str, str]]:
    tex = manuscript.read_text(encoding="utf-8") if manuscript.exists() else ""
    bib_text = bib.read_text(encoding="utf-8") if bib.exists() else ""
    rows: list[dict[str, str]] = []
    rows.append(_existence_row("manuscript_source", manuscript))
    rows.append(_existence_row("bibliography_source", bib))
    rows.append(_citation_row(tex, bib_text))
    rows.extend(_main_figure_rows(tex, figures_dir))
    rows.extend(_expected_file_rows("extended_data_figure", figures_dir, EXPECTED_EXTENDED_FIGURES))
    rows.extend(_table_rows(tables_dir))
    rows.append(_tex_engine_row())
    return rows


def render_markdown(table: pd.DataFrame) -> str:
    lines = [
        "# Manuscript Package Check",
        "",
        "Updated: 2026-06-23",
        "",
        "This audit is generated by `make manuscript-check`. A `blocked` TeX row means the local machine lacks a TeX engine; it is not a manuscript asset failure.",
        "",
        "| Check | Status | Detail |",
        "| --- | --- | --- |",
    ]
    for _, row in table.iterrows():
        lines.append(f"| `{row['check']}` | {row['status']} | {row['detail']} |")
    lines.append("")
    return "\n".join(lines)


def _citation_row(tex: str, bib_text: str) -> dict[str, str]:
    cite_keys = _citation_keys(tex)
    bib_keys = set(re.findall(r"@\w+\s*\{\s*([^,]+)", bib_text))
    missing = sorted(cite_keys - bib_keys)
    unused = sorted(bib_keys - cite_keys)
    if missing:
        return {
            "check": "citation_keys",
            "status": "fail",
            "detail": "missing from bibliography: " + ", ".join(missing),
        }
    return {
        "check": "citation_keys",
        "status": "pass",
        "detail": f"{len(cite_keys)} cited keys resolve; {len(unused)} unused bibliography keys.",
    }


def _citation_keys(tex: str) -> set[str]:
    keys: set[str] = set()
    pattern = re.compile(r"\\cite\w*\*?(?:\[[^\]]*\]){0,2}\{([^}]+)\}")
    for match in pattern.finditer(tex):
        keys.update(key.strip() for key in match.group(1).split(",") if key.strip())
    return keys


def _main_figure_rows(tex: str, figures_dir: Path) -> list[dict[str, str]]:
    refs = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex)
    rows = [_existence_row(f"main_figure:{ref}", figures_dir / ref) for ref in refs]
    if not refs:
        rows.append({"check": "main_figures", "status": "fail", "detail": "No includegraphics references found."})
    return rows


def _expected_file_rows(prefix: str, base_dir: Path, filenames: list[str]) -> list[dict[str, str]]:
    return [_existence_row(f"{prefix}:{name}", base_dir / name) for name in filenames]


def _table_rows(tables_dir: Path) -> list[dict[str, str]]:
    rows = []
    for name in EXPECTED_MANUSCRIPT_TABLES:
        path = tables_dir / name
        row = _existence_row(f"manuscript_table:{name}", path)
        if path.exists():
            n_lines = sum(1 for _ in path.open(encoding="utf-8"))
            row["detail"] += f"; {max(n_lines - 1, 0)} data rows"
        rows.append(row)
    return rows


def _tex_engine_row() -> dict[str, str]:
    latexmk = shutil.which("latexmk")
    pdflatex = shutil.which("pdflatex")
    bibtex = shutil.which("bibtex")
    if latexmk:
        return {"check": "tex_engine", "status": "pass", "detail": f"latexmk available at {latexmk}"}
    if pdflatex and bibtex:
        return {"check": "tex_engine", "status": "pass", "detail": "pdflatex and bibtex available"}
    return {
        "check": "tex_engine",
        "status": "blocked",
        "detail": "No local latexmk or pdflatex+bibtex; run `make manuscript` in a TeX-enabled environment.",
    }


def _existence_row(check: str, path: Path) -> dict[str, str]:
    if not path.exists():
        return {"check": check, "status": "fail", "detail": f"missing: {path}"}
    size = path.stat().st_size
    if size <= 0:
        return {"check": check, "status": "fail", "detail": f"empty file: {path}"}
    return {"check": check, "status": "pass", "detail": f"{path} ({size} bytes)"}


def _read_optional(path: str | None) -> pd.DataFrame | None:
    if not path:
        return None
    candidate = Path(path)
    if not candidate.exists():
        return None
    return pd.read_csv(candidate, sep="\t")


def _write_tsv(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)


def _write_text(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
