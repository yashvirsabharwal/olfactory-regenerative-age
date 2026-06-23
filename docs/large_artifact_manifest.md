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

## Present Remote Compute Artifacts

Remote location: `sabharwaly2@mia.ninds.nih.gov:/home/sabharwaly2/olfactory`. These paths are current compute storage locations, not stable publication archive URIs. Before submission, copy or register the final artifacts in controlled long-term storage and preserve the checksums below.

| Artifact | Remote path | Size | SHA-256 | Role | Regeneration command |
| --- | --- | ---: | --- | --- | --- |
| `data/processed/gateway_hvg3003_4m.h5ad` | `/home/sabharwaly2/olfactory/data/processed/gateway_hvg3003_4m.h5ad` | 11,296,884,912 bytes | `91cde68121d9609d993761fa56639a7dfb9cf90863d5726bf48317d9e892d40b` | Chunked 3,003-gene all-cell substrate for primary full 4M reduced scVI. | `make scvi-reduced-4m` |
| `data/processed/gateway_scvi_full_4m_reduced.h5ad` | `/home/sabharwaly2/olfactory/data/processed/gateway_scvi_full_4m_reduced.h5ad` | 11,722,478,874 bytes | `fd195a86703909ba2331bdc32e817d500321d641dd8d06809ac4c26ea259a3cd` | Primary all-cell reduced scVI substrate for publication-scale latent/neighborhood analyses. | `make scvi-full-4m-reduced` |
| `results/models/gateway_scvi_full_4m_reduced/model.pt` | `/home/sabharwaly2/olfactory/results/models/gateway_scvi_full_4m_reduced/model.pt` | 44,237,347 bytes | `38b88b26bc0f82dc89dcf8da22d5c2a4c8f130d0809fed75653b26fc9698425f` | Model weights paired with the primary full 4M reduced scVI H5AD. | `make scvi-full-4m-reduced` |
| `data/processed/gateway_scvi_stratified_250k_seed23.h5ad` | `/home/sabharwaly2/olfactory/data/processed/gateway_scvi_stratified_250k_seed23.h5ad` | 675,514,610 bytes | `42eb3ce87140d60ba14c358dc8ef225c39bd03e064e65681dc3dadb5830cdbe9` | Seed-stability atlas for scaled scVI comparison. | `make scvi-scaled-250k-seed23` |
| `results/models/gateway_scvi_stratified_250k_seed23/model.pt` | `/home/sabharwaly2/olfactory/results/models/gateway_scvi_stratified_250k_seed23/model.pt` | 10,249,827 bytes | `2d42db0b4ce1707eab0695cf9a2b3adddf07a9965697932059ec64118d8f6505` | Model weights paired with the seed-23 250k atlas. | `make scvi-scaled-250k-seed23` |

## Deferred Or Missing Heavy Artifacts

These are intentionally not present in the current local checkout and were not found in the current `mia` project artifact scan. They are alternate or fallback runs unless explicitly cited.

| Artifact | Current state | Why it matters | Required action before submission |
| --- | --- | --- | --- |
| `data/processed/gateway_scvi_stratified_500k.h5ad` | Missing/deferred | Large fallback scaled atlas. | Restore only if this fallback is cited; otherwise mark as deferred in the final package. |
| `results/models/gateway_scvi_stratified_500k` | Missing/deferred | Model directory for 500k fallback. | Restore only if the 500k run is cited. |
| `data/processed/gateway_scvi_stratified_1m.h5ad` | Missing/deferred | Large fallback scaled atlas. | Restore only if this fallback is cited; otherwise keep deferred. |
| `results/models/gateway_scvi_stratified_1m` | Missing/deferred | Model directory for 1M fallback. | Restore only if the 1M run is cited. |
| `data/processed/gateway_scvi_full_4m.h5ad` | Missing/deferred | Direct full 4M scVI attempt. | Restore only if cited; otherwise keep as deferred alternative. |
| `results/models/gateway_scvi_full_4m` | Missing/deferred | Model directory for direct full 4M attempt. | Restore only if cited. |
| `data/processed/gateway_scvi_full_4m_safe.h5ad` | Missing/deferred | Reduced-gene memory-safe full 4M alternative. | Restore only if cited. |
| `results/models/gateway_scvi_full_4m_safe` | Missing/deferred | Model directory for memory-safe full 4M alternative. | Restore only if cited. |

## Present Result Tables From Deferred Runs

| Result family | Current local state | Review stance |
| --- | --- | --- |
| `results/tables/scvi_full_4m_reduced_validation.tsv` and scVI claim gates | Present locally and on `mia` | Supports current manuscript tables; upstream 4M H5AD/model checksums are recorded above. |
| `results/tables/milo_full_4m_*` summaries, neighborhoods, age-bin, edgeR parity, and program tables | Present locally and on `mia` | Supports latent/neighborhood claim; primary 4M substrate is identified on `mia` until archived. |
| `results/tables/milor_lineage*_subset_*` summaries and DA tables | Present locally and on `mia` | Supports official MiloR subset sensitivity; environment and rerun notes are in `docs/manuscript_rerun_profile.md` and `docs/run_hierarchy.md`. |

## Final Packaging Checklist

| Check | Required result |
| --- | --- |
| Local checksum refresh | Recompute `shasum -a 256` after final reruns or restores. |
| Remote location | Replace current `mia` compute paths with stable storage URIs for every H5AD/model artifact used in final figures or tables. |
| Provenance refresh | Run `PYTHON=.venv/bin/python make output-provenance`; require zero missing non-deferred outputs. |
| Manuscript consistency | Ensure `docs/run_hierarchy.md`, `docs/journal_acceptance_tracker.md`, and data/code availability text all point to this manifest. |
