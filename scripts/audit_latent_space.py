#!/usr/bin/env python3
"""Audit local and portal latent embeddings for Gateway trajectory readiness."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config, project_path
from ora.latent import (
    fetch_cellxgene_collection_assets,
    inspect_h5ad_obsm,
    latent_readiness_summary,
    render_latent_space_plan,
)
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--skip-portal", action="store_true")
    parser.add_argument("--local-audit-out", default=None)
    parser.add_argument("--portal-assets-out", default=None)
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--plan-out", default=None)
    parser.add_argument("--portal-json-out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    source = config.get("source", {})
    h5ad_path = project_path(args.h5ad or source.get("h5ad_path", "data/raw/gateway.h5ad"))
    local_audit_out = args.local_audit_out or outputs.get(
        "latent_space_local_audit_tsv",
        "results/tables/latent_space_local_audit.tsv",
    )
    portal_assets_out = args.portal_assets_out or outputs.get(
        "latent_space_portal_assets_tsv",
        "results/tables/latent_space_portal_assets.tsv",
    )
    summary_out = args.summary_out or outputs.get(
        "latent_space_readiness_tsv",
        "results/tables/latent_space_readiness.tsv",
    )
    plan_out = args.plan_out or outputs.get("latent_space_plan_md", "results/reports/latent_space_recovery_plan.md")
    portal_json_out = args.portal_json_out or outputs.get(
        "latent_space_portal_json",
        "results/reports/latent_space_portal_collection.json",
    )

    local_audit = inspect_h5ad_obsm(h5ad_path, config.get("embedding_priority", []))
    collection_id = str(source.get("collection_id", ""))
    if args.skip_portal or not collection_id:
        portal_assets, portal_payload = _skipped_portal_assets(collection_id), {}
    else:
        portal_assets, portal_payload = fetch_cellxgene_collection_assets(collection_id)
    summary = latent_readiness_summary(local_audit, portal_assets)
    plan = render_latent_space_plan(summary, local_audit, portal_assets)

    local_audit.to_csv(ensure_parent(local_audit_out), sep="\t", index=False)
    portal_assets.to_csv(ensure_parent(portal_assets_out), sep="\t", index=False)
    summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    ensure_parent(plan_out).write_text(plan, encoding="utf-8")
    if portal_payload:
        ensure_parent(portal_json_out).write_text(json.dumps(portal_payload, indent=2, sort_keys=True), encoding="utf-8")
    print(f"Wrote local latent audit: {local_audit_out} ({local_audit.shape[0]} rows)")
    print(f"Wrote portal latent asset audit: {portal_assets_out} ({portal_assets.shape[0]} rows)")
    print(f"Wrote latent readiness summary: {summary_out}")
    print(f"Wrote latent-space recovery plan: {plan_out}")


def _skipped_portal_assets(collection_id: str):
    import pandas as pd

    return pd.DataFrame(
        [
            {
                "source": "cellxgene_api",
                "collection_id": collection_id,
                "dataset_id": "",
                "dataset_version_id": "",
                "dataset_title": "",
                "asset_filetype": "",
                "asset_filesize_bytes": 0,
                "asset_url": "",
                "portal_embeddings": "",
                "default_embedding": "",
                "status": "skipped",
                "notes": "Portal fetch skipped by CLI flag.",
            }
        ]
    )


if __name__ == "__main__":
    main()
