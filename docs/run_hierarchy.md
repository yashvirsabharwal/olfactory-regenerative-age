# ORA Run Hierarchy And Publication Roles

Updated: 2026-06-22

This document standardizes how to describe the different scVI, scANVI, and neighborhood runs. The goal is to keep the manuscript and tracker consistent: large all-cell analyses carry the main latent/neighborhood evidence, while smaller runs are labeled as engineering checks, reference-mapping anchors, or sensitivity analyses.

## Core Principle

Do not mix cell counts as if they are interchangeable. Each run has a role. A 25k or 100k run can validate code paths and lineage plausibility, but it is not the same claim as a full 4,028,275-cell run. Conversely, an official MiloR subset can be a useful implementation check without needing to replace the full-scale Python Milo-style map.

## Run Roles

| Run or artifact | Cell count | Data source | Role | Publication use | Status |
| --- | ---: | --- | --- | --- | --- |
| scVI pilot | 25,000 | Gateway subset | Engineering smoke test for dependency, memory, and output shape | Not manuscript evidence | Complete |
| Lineage-focused scVI | 100,000 | Gateway lineage-enriched subset | Sensitivity anchor for basal/progenitor/neuronal structure | Supplementary validation only | Complete |
| Stratified scVI seed 13 | 250,000 | Gateway stratified subset | Scalable latent reference and scANVI/scArches anchor | Reference-mapping and sensitivity support | Complete |
| Stratified scVI seed 23 | 250,000 | Gateway stratified subset | Seed sensitivity for the scaled reference | Supplementary validation support | Complete |
| 500k / 1M scVI targets | 500,000 to 1,000,000 | Gateway stratified subset | Optional stress-test targets | No claim unless run and validated | Optional/deferred |
| Full reduced Gateway scVI | 4,028,275 | Gateway all cells, 3,003 HVG/marker genes | Primary latent substrate | Main latent-atlas method and input to full-scale neighborhoods | Complete |
| Full 4M Python Milo-style neighborhoods | 4,028,275 latent substrate; 20,000 neighborhoods | Full reduced Gateway scVI | Primary publication-scale neighborhood map | Secondary mechanistic layer, labeled "Milo-style" | Complete |
| Full 4M exact-neighborhood edgeR parity | Same 20,000 Python neighborhoods | Full reduced Gateway scVI memberships | Statistical-core sensitivity on the same neighborhoods | Extended Data or supplement supporting directionality | Complete |
| Official MiloR lineage subset | 100,000 | Stratified subset from full 4M latent space | Independent neighborhood-construction sensitivity | Extended Data or supplement, not the main map | Complete |
| Official MiloR matched lineage subset | 75,000 | Matched FLEX v2/device subset from full 4M latent space | Independent matched sensitivity | Extended Data or supplement; narrows Early iOSN language | Complete |
| GSE184117 scANVI/scArches mapping | 59,656 query cells | External presbyosmia/normosmia samples | External reference-mapping sanity check | Small-n external support only | Complete |

## Publication Rules

- The primary latent analysis is the all-cell reduced Gateway scVI model trained on 4,028,275 cells.
- The primary neighborhood map is the full-scale Python Milo-style workflow on that all-cell latent substrate.
- The manuscript must call the Python workflow "Milo-style" rather than official MiloR.
- Official MiloR is useful and should be kept, but it is a subset sensitivity. It does not need to be rerun on all 4M cells for publication unless a reviewer specifically requires official-MiloR-only discovery.
- The exact-neighborhood edgeR parity is the cleanest statistical-core sensitivity because it tests the same full-scale Python neighborhoods with a negative-binomial count model.
- The official MiloR subset changes both neighborhood construction and the statistical model. It answers a different question: whether age-associated structure is visible under the canonical MiloR framework on a tractable, stratified sample of the same latent atlas.
- The 25k, 100k, 250k, 500k, and 1M labels should never be described as progressively better versions of the same final result. They are pilot, sensitivity, reference, or optional stress-test runs.

## Current Interpretation

The consistent manuscript framing is:

- Full 4M scVI establishes a validated latent substrate for the Gateway atlas.
- Full 4M Milo-style neighborhoods show broad age-associated latent-neighborhood structure.
- Exact-neighborhood edgeR parity supports the signed direction of the full-scale Python neighborhood results, including the matched Early iOSN hit.
- Official MiloR subset sensitivity confirms broad age-associated HBC/suprabasal/sustentacular lineage-neighborhood structure, but it does not independently make matched Early iOSN the dominant official-MiloR finding.
- Therefore, the safest publication language is a broad secondary claim about age-associated lineage-neighborhood remodeling, with Early iOSN described as a narrow exact-neighborhood result supported by edgeR, age-bin directionality, and curated program enrichment.

## Answer To The 4M MiloR Question

We do not need official MiloR on all 4M cells for publication. We already have the publication-scale 4M latent substrate and the full-scale 20,000-neighborhood Python Milo-style map. Official MiloR on 75k/100k cells is valuable as an implementation sensitivity because it checks the canonical method on a stratified subset from the same full latent space. It should strengthen the paper as a transparent sensitivity analysis, not redefine the primary result.
