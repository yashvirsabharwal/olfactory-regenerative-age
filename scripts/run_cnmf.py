#!/usr/bin/env python3
"""Deferred external-compute command for cNMF program discovery."""

from __future__ import annotations

import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.parse_args()
    raise SystemExit("cNMF is an external-compute stretch goal and is not part of the MVP scaffold.")


if __name__ == "__main__":
    main()

