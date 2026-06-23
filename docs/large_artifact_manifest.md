# Large Artifact Manifest

Updated: 2026-06-23

Purpose: make large local and deferred artifacts reviewable without committing them to Git. Generated analysis data remain ignored by `.gitignore`; this manifest records what exists locally, what is deferred, and what must be restored or regenerated before final submission packaging.

## Policy

| Rule | Decision |
| --- | --- |
| Git storage | Do not commit H5AD files, model directories, raw Gateway data, GEO archives, or large generated result tables. |
| Reviewer package | Provide checksums and storage/access notes for every heavyweight artifact needed to reproduce manuscript-facing tables/figures. |
| Remote storage | Before submission, move final heavyweight artifacts to controlled remote storage or an institutional repository and record stable URIs here. |
| Checksum standard | Use SHA-256 for local and remote artifact identity. |
| Regeneration record | Every heavyweight artifact must map to a Make target or command-manifest entry. |

## Present Local Heavy Artifacts

| Artifact | Size | SHA-256 | Role | Regeneration command |
| --- | ---: | --- | --- | --- |
| `data/processed/gateway_scvi_stratified_250k.h5ad` | 641M | `47c707caa6b923596df7d4cbc71afab9f1b1b023631c433bbe421304dd77560c` | 250k scVI atlas used for scaled validation and external scANVI reference input. | `make scvi-scaled-250k` |
| `data/processed/gateway_scvi_lineage_basal_neural_100k.h5ad` | 251M | `6ec288b8318fea4422e8fbaea77313e5e60d48f5471388d022723fb5bba7f7ed` | Lineage-focused latent validation and pilot neighborhood support. | `make scvi-lineage-basal-neural` |
| `data/processed/gse184117_scanvi_mapped.h5ad` | 240M | `91f9aa512832d440d55ae8abf1021fc700c5dc9778cd6094b922fb34437c9eea` | GSE184117 scANVI/scArches mapped external validation object. | `make external-scanvi-reference-map` |
| `data/processed/gateway_scvi_pilot_25k.h5ad` | 30M | `29d6a37d4ed6fd591ae23780dd133f059abd1512c6877f3aac834741fb7ed` | Early pilot latent sanity artifact. | `make scvi-pilot` |
| `data/processed/gse184117_marker_mapped.h5ad` | 2.2G | `06c46a102e6877f3aac834741fb3747b7e800d045002cf09264b5bc3d276b4ad` | GSE184117 marker-reference mapped external validation object. | `make external-gse184117-mapped` |
| `results/models/gateway_scanvi_reference/model.pt` | 11M | `88a2c97f6f63affe30840a3ea43dcc433d0c17acd57b35840fb7ea932b7fd017` | Gateway scANVI reference model used for GSE184117 query transfer. | `make external-scanvi-reference-map` |

## Deferred Or Missing Heavy Artifacts

These are intentionally not present in the current local checkout. Their manuscript-facing result summaries exist, but the upstream H5AD/model artifacts need restoration, remote storage, or regeneration before the final reproducibility claim is closed.

| Artifact | Current state | Why it matters | Required action before submission |
| --- | --- | --- | --- |
| `data/processed/gateway_scvi_stratified_250k_seed23.h5ad` | Missing/deferred | Seed-stability input for scaled scVI comparison. | Restore or rerun `make scvi-scaled-250k-seed23`; record checksum and storage URI. |
| `results/models/gateway_scvi_stratified_250k_seed23` | Missing/deferred | Model directory paired with the seed-23 250k atlas. | Restore or rerun with the H5AD artifact. |
| `data/processed/gateway_scvi_stratified_500k.h5ad` | Missing/deferred | Large fallback scaled atlas. | Restore only if this fallback is cited; otherwise mark as deferred in the final package. |
| `results/models/gateway_scvi_stratified_500k` | Missing/deferred | Model directory for 500k fallback. | Restore only if the 500k run is cited. |
| `data/processed/gateway_scvi_stratified_1m.h5ad` | Missing/deferred | Large fallback scaled atlas. | Restore only if this fallback is cited; otherwise keep deferred. |
| `results/models/gateway_scvi_stratified_1m` | Missing/deferred | Model directory for 1M fallback. | Restore only if the 1M run is cited. |
| `data/processed/gateway_hvg3003_4m.h5ad` | Missing/deferred | Chunked 3003-gene all-cell substrate for the reduced 4M scVI model. | Restore or rerun `make scvi-reduced-4m`; record checksum and storage URI. |
| `data/processed/gateway_scvi_full_4m.h5ad` | Missing/deferred | Direct full 4M scVI attempt. | Restore only if cited; otherwise keep as deferred alternative. |
| `results/models/gateway_scvi_full_4m` | Missing/deferred | Model directory for direct full 4M attempt. | Restore only if cited. |
| `data/processed/gateway_scvi_full_4m_safe.h5ad` | Missing/deferred | Reduced-gene memory-safe full 4M alternative. | Restore only if cited. |
| `results/models/gateway_scvi_full_4m_safe` | Missing/deferred | Model directory for memory-safe full 4M alternative. | Restore only if cited. |
| `data/processed/gateway_scvi_full_4m_reduced.h5ad` | Missing/deferred | Primary all-cell reduced scVI substrate for publication-scale Milo-style outputs. | Restore or rerun `make scvi-full-4m-reduced`; record checksum and storage URI. |
| `results/models/gateway_scvi_full_4m_reduced` | Missing/deferred | Model directory for the primary all-cell reduced scVI substrate. | Restore or rerun with the H5AD artifact. |

## Present Result Tables From Deferred Runs

| Result family | Current local state | Review stance |
| --- | --- | --- |
| `results/tables/scvi_full_4m_reduced_validation.tsv` and scVI claim gates | Present | Can support current manuscript tables, but final package needs the upstream 4M H5AD/model checksum. |
| `results/tables/milo_full_4m_*` summaries, neighborhoods, age-bin, edgeR parity, and program tables | Present | Supports latent/neighborhood claim, with the 4M substrate listed as deferred until restored. |
| `results/tables/milor_lineage*_subset_*` summaries and DA tables | Present | Supports official MiloR subset sensitivity; final package needs environment and upstream substrate notes. |

## Final Packaging Checklist

| Check | Required result |
| --- | --- |
| Local checksum refresh | Recompute `shasum -a 256` after final reruns or restores. |
| Remote location | Add stable storage URI for every H5AD/model artifact used in final figures or tables. |
| Provenance refresh | Run `PYTHON=.venv/bin/python make output-provenance`; require zero missing non-deferred outputs. |
| Manuscript consistency | Ensure `docs/run_hierarchy.md`, `docs/journal_acceptance_tracker.md`, and data/code availability text all point to this manifest. |

