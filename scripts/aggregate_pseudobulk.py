#!/usr/bin/env python3
"""Deferred heavy-stage command for donor x cell-type pseudobulk aggregation."""

from __future__ import annotations

import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.parse_args()
    raise SystemExit(
        "Pseudobulk aggregation is deferred from the composition-MVP scaffold. "
        "Use backed/chunked count aggregation before running DESeq2 or limma."
    )


if __name__ == "__main__":
    main()

