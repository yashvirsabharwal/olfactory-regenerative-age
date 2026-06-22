#!/usr/bin/env python3
"""Compare scVI validation runs into publication claim gates."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.scvi_comparison import compare_scvi_validation_runs
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--validation",
        action="append",
        nargs=2,
        metavar=("MODEL", "TSV"),
        required=True,
        help="Model label and validation TSV path. Can be repeated.",
    )
    parser.add_argument("--summary-out", default="results/tables/scvi_embedding_claim_gates.tsv")
    parser.add_argument("--markers-out", default="results/tables/scvi_embedding_marker_concordance.tsv")
    parser.add_argument("--note-out", default="docs/scvi_embedding_comparison.md")
    args = parser.parse_args()

    validation_paths = {model: path for model, path in args.validation}
    summary, markers, note = compare_scvi_validation_runs(validation_paths)
    summary.to_csv(ensure_parent(args.summary_out), sep="\t", index=False)
    markers.to_csv(ensure_parent(args.markers_out), sep="\t", index=False)
    Path(args.note_out).write_text(note, encoding="utf-8")
    print(f"Wrote scVI embedding claim gates: {args.summary_out} ({summary.shape[0]} rows)")
    print(f"Wrote scVI marker concordance: {args.markers_out} ({markers.shape[0]} rows)")
    print(f"Wrote scVI comparison note: {args.note_out}")


if __name__ == "__main__":
    main()
