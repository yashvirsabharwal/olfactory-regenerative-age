"""Latent-space readiness audits for trajectory and neighborhood analyses."""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import h5py
import pandas as pd


LATENT_READY_KEYS = {"X_scANVI", "X_scanvi", "X_scvi", "X_pca", "X_harmony", "X_totalVI"}
VISUALIZATION_ONLY_KEYS = {"X_umap", "X_tsne", "X_draw_graph_fa"}


def inspect_h5ad_obsm(path: str | Path, embedding_priority: list[str] | tuple[str, ...] = ()) -> pd.DataFrame:
    """Inspect `.obsm` keys and shapes without loading the full H5AD matrix."""

    h5ad_path = Path(path)
    if not h5ad_path.exists():
        return pd.DataFrame(
            [
                {
                    "source": "local_h5ad",
                    "embedding_key": "",
                    "n_cells": 0,
                    "n_dimensions": 0,
                    "priority_rank": pd.NA,
                    "readiness": "missing_h5ad",
                    "recommended_use": "none",
                    "notes": f"Missing H5AD: {h5ad_path}",
                }
            ]
        )
    priority = {str(key): idx + 1 for idx, key in enumerate(embedding_priority)}
    rows = []
    with h5py.File(h5ad_path, "r") as handle:
        obsm = handle.get("obsm")
        if obsm is None or len(obsm.keys()) == 0:
            rows.append(
                {
                    "source": "local_h5ad",
                    "embedding_key": "",
                    "n_cells": 0,
                    "n_dimensions": 0,
                    "priority_rank": pd.NA,
                    "readiness": "no_embeddings",
                    "recommended_use": "none",
                    "notes": "No obsm group or no obsm keys found.",
                }
            )
        else:
            for key in sorted(obsm.keys()):
                shape = _hdf5_shape(obsm[key])
                n_cells = int(shape[0]) if len(shape) >= 1 else 0
                n_dimensions = int(shape[1]) if len(shape) >= 2 else 0
                readiness, recommended_use = classify_embedding(key, n_dimensions)
                rows.append(
                    {
                        "source": "local_h5ad",
                        "embedding_key": str(key),
                        "n_cells": n_cells,
                        "n_dimensions": n_dimensions,
                        "priority_rank": priority.get(str(key), pd.NA),
                        "readiness": readiness,
                        "recommended_use": recommended_use,
                        "notes": _embedding_note(key, n_dimensions),
                    }
                )
    return pd.DataFrame(rows).sort_values(["priority_rank", "embedding_key"], na_position="last").reset_index(drop=True)


def classify_embedding(key: str, n_dimensions: int) -> tuple[str, str]:
    """Classify whether an embedding is acceptable for latent-space analyses."""

    if key in LATENT_READY_KEYS and n_dimensions >= 5:
        return "usable_latent", "trajectory_milo_neighbors_after_validation"
    if key in VISUALIZATION_ONLY_KEYS:
        return "visualization_only", "display_only"
    if n_dimensions >= 10:
        return "candidate_latent_unknown", "requires_validation"
    return "insufficient_dimensions", "display_or_diagnostics_only"


def fetch_cellxgene_collection_assets(collection_id: str, timeout: int = 30) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Fetch CELLxGENE collection asset and embedding metadata."""

    url = f"https://api.cellxgene.cziscience.com/curation/v1/collections/{collection_id}"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        return (
            pd.DataFrame(
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
                        "status": "fetch_failed",
                        "notes": str(exc),
                    }
                ]
            ),
            {},
        )
    rows = []
    for dataset in payload.get("datasets", []):
        embeddings = ",".join(str(item) for item in dataset.get("embeddings", []))
        default_embedding = str(dataset.get("default_embedding", ""))
        assets = dataset.get("assets", []) or [{}]
        for asset in assets:
            rows.append(
                {
                    "source": "cellxgene_api",
                    "collection_id": collection_id,
                    "dataset_id": dataset.get("dataset_id", ""),
                    "dataset_version_id": dataset.get("dataset_version_id", ""),
                    "dataset_title": dataset.get("title", ""),
                    "asset_filetype": asset.get("filetype", ""),
                    "asset_filesize_bytes": int(asset.get("filesize") or 0),
                    "asset_url": asset.get("url", ""),
                    "portal_embeddings": embeddings,
                    "default_embedding": default_embedding,
                    "status": "ok",
                    "notes": "",
                }
            )
    if not rows:
        rows.append(
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
                "status": "no_datasets",
                "notes": "Collection fetched but no datasets were listed.",
            }
        )
    return pd.DataFrame(rows), payload


def latent_readiness_summary(local_audit: pd.DataFrame, portal_assets: pd.DataFrame | None = None) -> pd.DataFrame:
    """Create a one-row publication-readiness summary from latent audit tables."""

    local = local_audit.copy() if local_audit is not None else pd.DataFrame()
    portal = portal_assets.copy() if portal_assets is not None else pd.DataFrame()
    usable = (
        local[local["readiness"].isin(["usable_latent", "candidate_latent_unknown"])]
        if "readiness" in local
        else pd.DataFrame()
    )
    visual = (
        local[local["readiness"].eq("visualization_only")]
        if "readiness" in local
        else pd.DataFrame()
    )
    portal_embeddings = sorted(
        {
            item
            for value in portal.get("portal_embeddings", pd.Series(dtype=object)).fillna("").astype(str)
            for item in value.split(",")
            if item
        }
    )
    status = "ready_for_latent_validation" if not usable.empty else "latent_recompute_required"
    recommendation = (
        "Validate existing non-UMAP embedding before trajectory, Milo, or cNMF."
        if not usable.empty
        else "Search author assets or recompute a non-UMAP latent representation before trajectory, Milo, or cNMF."
    )
    return pd.DataFrame(
        [
            {
                "status": status,
                "local_embeddings": ",".join(
                    str(item) for item in local.get("embedding_key", pd.Series(dtype=object)).dropna() if str(item)
                )
                or "none",
                "portal_embeddings": ",".join(portal_embeddings) or "none",
                "usable_local_embeddings": ",".join(usable.get("embedding_key", pd.Series(dtype=object)).astype(str))
                or "none",
                "visualization_only_embeddings": ",".join(visual.get("embedding_key", pd.Series(dtype=object)).astype(str))
                or "none",
                "recommendation": recommendation,
            }
        ]
    )


def render_latent_space_plan(
    summary: pd.DataFrame,
    local_audit: pd.DataFrame,
    portal_assets: pd.DataFrame | None = None,
) -> str:
    """Render a Markdown latent-space recovery plan from audit outputs."""

    row = summary.iloc[0] if not summary.empty else pd.Series(dtype=object)
    status = str(row.get("status", "latent_recompute_required"))
    local_embeddings = str(row.get("local_embeddings", "")) or "none recorded"
    portal_embeddings = str(row.get("portal_embeddings", "")) or "none recorded"
    usable = str(row.get("usable_local_embeddings", "")) or "none"
    recommendation = str(row.get("recommendation", "Recompute a non-UMAP latent representation."))
    portal_asset_lines = _portal_asset_lines(portal_assets)
    local_lines = _local_embedding_lines(local_audit)
    return f"""# Latent-Space Recovery Plan

Status: `{status}`

## Current Export

- Local H5AD embeddings: {local_embeddings}
- CELLxGENE portal-reported embeddings: {portal_embeddings}
- Usable non-UMAP latent embeddings found locally: {usable}
- UMAP alone is not acceptable for trajectory, Milo, or cNMF claims.

## CELLxGENE Asset Audit

{portal_asset_lines}

## Local Embedding Audit

{local_lines}

## Recommendation

{recommendation}

## Next Steps

1. Ask the Gateway authors or source repository for original `X_scANVI` or `X_scvi` coordinates if available outside CELLxGENE.
2. If unavailable, recompute a latent representation from the H5AD using a documented Scanpy/scvi-tools workflow.
3. Use HVG selection within olfactory epithelial/neuronal subsets, preserve donor/sample metadata, and include chemistry, collection method, and donor/sample batch covariates where supported.
4. Validate the latent space with marker continuity, donor/chemistry mixing diagnostics, and negative-control checks before running pseudotime, Milo, or cNMF.
5. Add report sections only after the latent-space validation table exists.

## Claim Guardrail

Trajectory and neighborhood findings must be described as deferred until a non-UMAP latent space is available and validated.
"""


def _hdf5_shape(node: h5py.Dataset | h5py.Group) -> tuple[int, ...]:
    if isinstance(node, h5py.Dataset):
        return tuple(int(dim) for dim in node.shape)
    if isinstance(node, h5py.Group):
        for name in ["data", "X", "values"]:
            child = node.get(name)
            if isinstance(child, h5py.Dataset):
                return tuple(int(dim) for dim in child.shape)
    return ()


def _embedding_note(key: str, n_dimensions: int) -> str:
    if key in VISUALIZATION_ONLY_KEYS:
        return "Visualization embedding; do not use as the basis for trajectory or neighborhood claims."
    if key in LATENT_READY_KEYS and n_dimensions >= 5:
        return "Candidate latent embedding; validate before downstream analyses."
    if n_dimensions < 5:
        return "Too few dimensions for robust neighborhood/trajectory analyses."
    return "Unknown embedding; validate provenance and batch behavior before use."


def _portal_asset_lines(portal_assets: pd.DataFrame | None) -> str:
    if portal_assets is None or portal_assets.empty:
        return "- No portal asset metadata recorded."
    lines = []
    for _, row in portal_assets.iterrows():
        if str(row.get("status", "")) != "ok":
            lines.append(f"- CELLxGENE API status: `{row.get('status')}` ({row.get('notes')}).")
            continue
        lines.append(
            "- "
            f"{row.get('dataset_title', 'dataset')}: asset `{row.get('asset_filetype', '')}` "
            f"with portal embeddings `{row.get('portal_embeddings', '') or 'none'}`."
        )
    return "\n".join(lines)


def _local_embedding_lines(local_audit: pd.DataFrame) -> str:
    if local_audit is None or local_audit.empty:
        return "- No local embedding audit rows recorded."
    lines = []
    for _, row in local_audit.iterrows():
        key = str(row.get("embedding_key", "")) or "(none)"
        lines.append(
            "- "
            f"`{key}`: {row.get('n_cells')} cells x {row.get('n_dimensions')} dimensions; "
            f"readiness `{row.get('readiness')}`."
        )
    return "\n".join(lines)
