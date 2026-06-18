# Latent-Space Recovery Plan

Status: `latent_recompute_required`

## Current Export

- Local H5AD embeddings: X_umap
- CELLxGENE portal-reported embeddings: X_umap
- Usable non-UMAP latent embeddings found locally: none
- UMAP alone is not acceptable for trajectory, Milo, or cNMF claims.

## CELLxGENE Asset Audit

- Crownlands Gateway — Olfactory + Respiratory Single-Cell Atlas (v1): asset `H5AD` with portal embeddings `X_umap`.

## Local Embedding Audit

- `X_umap`: 4028275 cells x 2 dimensions; readiness `visualization_only`.

## Recommendation

Search author assets or recompute a non-UMAP latent representation before trajectory, Milo, or cNMF.

## Next Steps

1. Ask the Gateway authors or source repository for original `X_scANVI` or `X_scvi` coordinates if available outside CELLxGENE.
2. If unavailable, recompute a latent representation from the H5AD using a documented Scanpy/scvi-tools workflow.
3. Use HVG selection within olfactory epithelial/neuronal subsets, preserve donor/sample metadata, and include chemistry, collection method, and donor/sample batch covariates where supported.
4. Validate the latent space with marker continuity, donor/chemistry mixing diagnostics, and negative-control checks before running pseudotime, Milo, or cNMF.
5. Add report sections only after the latent-space validation table exists.

## Claim Guardrail

Trajectory and neighborhood findings must be described as deferred until a non-UMAP latent space is available and validated.
