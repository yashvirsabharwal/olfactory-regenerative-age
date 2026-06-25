#!/usr/bin/env python3
"""Build the publication release artifact manifest."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.release import build_release_manifest, render_release_manifest_markdown
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--command-manifest", default="configs/command_manifest.yaml")
    parser.add_argument("--checksum-max-bytes", type=int, default=100 * 1024 * 1024)
    parser.add_argument("--out", default=None)
    parser.add_argument("--markdown-out", default=None)
    args = parser.parse_args()

    gateway_config = load_config(args.config)
    command_config = load_config(args.command_manifest)
    outputs = gateway_config.get("outputs", {})
    out_path = args.out or outputs.get(
        "release_artifact_manifest_tsv",
        "results/reports/release_artifact_manifest.tsv",
    )
    markdown_path = args.markdown_out or outputs.get(
        "release_artifact_manifest_md",
        "results/reports/release_artifact_manifest.md",
    )
    manifest = build_release_manifest(
        command_config,
        checksum_max_bytes=args.checksum_max_bytes,
    )
    manifest.to_csv(ensure_parent(out_path), sep="\t", index=False)
    ensure_parent(markdown_path).write_text(render_release_manifest_markdown(manifest))
    missing_required = int(
        (
            manifest["required_for_review"].astype(bool)
            & ~manifest["artifact_status"].astype(str).isin(["present", "archived"])
        ).sum()
    )
    print(f"Wrote release artifact manifest: {out_path} ({manifest.shape[0]} rows)")
    print(f"Wrote release artifact summary: {markdown_path} ({missing_required} required missing/deferred)")


if __name__ == "__main__":
    main()
