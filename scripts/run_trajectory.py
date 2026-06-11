#!/usr/bin/env python3
"""Deferred heavy-stage command for lineage pseudotime and bottleneck analysis."""

from __future__ import annotations

import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.parse_args()
    raise SystemExit(
        "Trajectory analysis is deferred from the MVP. "
        "Use the Gateway scANVI latent representation when implementing this stage."
    )


if __name__ == "__main__":
    main()

