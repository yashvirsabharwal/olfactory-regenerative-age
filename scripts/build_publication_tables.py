#!/usr/bin/env python3
"""Build compact ORA publication tables."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.publication_tables import build_publication_tables, render_publication_table_index
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables-dir", default="results/tables")
    parser.add_argument("--out-dir", default="results/tables")
    parser.add_argument("--index-out", default="results/reports/publication_tables.md")
    args = parser.parse_args()

    tables = build_publication_tables(args.tables_dir)
    out_dir = Path(args.out_dir)
    for name, table in tables.items():
        table.to_csv(ensure_parent(out_dir / f"{name}.tsv"), sep="\t", index=False)
    Path(args.index_out).write_text(render_publication_table_index(tables), encoding="utf-8")
    print(f"Wrote {len(tables)} publication tables to {out_dir}")
    print(f"Wrote publication table index: {args.index_out}")


if __name__ == "__main__":
    main()
