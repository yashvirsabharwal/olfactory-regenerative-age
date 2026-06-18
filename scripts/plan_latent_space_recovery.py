#!/usr/bin/env python3
"""Write the latent-space recovery plan for deferred trajectory/neighborhood work."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.reporting import load_schema
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--schema", default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    schema_path = args.schema or outputs.get("schema_json", "results/reports/h5ad_schema.json")
    out_path = args.out or outputs.get("latent_space_plan_md", "docs/latent_space_recovery_plan.md")
    schema = load_schema(schema_path)
    embeddings = [str(item) for item in schema.get("obsm_keys", [])]
    preferred = [key for key in config.get("embedding_priority", []) if key in embeddings]
    has_latent = any(key in embeddings for key in ["X_scANVI", "X_scvi", "X_pca"])
    text = _render_plan(embeddings, preferred, has_latent)
    ensure_parent(out_path).write_text(text, encoding="utf-8")
    print(f"Wrote latent-space recovery plan: {out_path}")


def _render_plan(embeddings: list[str], preferred: list[str], has_latent: bool) -> str:
    status = "ready_for_latent_work" if has_latent else "latent_recompute_required"
    return f"""# Latent-Space Recovery Plan

Status: `{status}`

## Current Export

- Available embeddings: {", ".join(embeddings) if embeddings else "none recorded"}
- Preferred usable embeddings found: {", ".join(preferred) if preferred else "none"}
- UMAP alone is not acceptable for trajectory, Milo, or cNMF claims.

## Next Steps

1. Search the Gateway portal, paper supplement, and author-provided assets for original `X_scANVI` or `X_scvi` coordinates.
2. If unavailable, recompute a latent representation from the H5AD using a documented Scanpy/scvi-tools workflow.
3. Use HVG selection within olfactory epithelial/neuronal subsets, preserve donor/sample metadata, and include chemistry, collection method, and donor/sample batch covariates where supported.
4. Validate the latent space with marker continuity, donor/chemistry mixing diagnostics, and negative-control checks before running pseudotime, Milo, or cNMF.
5. Add report sections only after the latent-space validation table exists.

## Claim Guardrail

Trajectory and neighborhood findings must be described as deferred until a non-UMAP latent space is available and validated.
"""


if __name__ == "__main__":
    main()
