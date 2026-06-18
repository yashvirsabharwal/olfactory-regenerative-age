# Full 4M Milo-Style Neighborhood Results

Updated: 2026-06-18

## Scope

These are publication-scale Milo-style donor-level neighborhood differential-abundance runs on the all-cell reduced Gateway scVI atlas:

- Input latent atlas: `data/processed/gateway_scvi_full_4m_reduced.h5ad`
- Cells represented: 4,028,275
- Latent dimensions: 10-dimensional `X_scvi`
- Seed neighborhoods per run: 20,000
- Nearest cells per neighborhood: 100
- Minimum donors per tested neighborhood: 30
- Donor-level model: logit neighborhood fraction ~ age + sex + chemistry + collection method + log donor cell yield
- Donor set: healthy age-annotated donors

This is a lightweight Python Milo-style workflow, not Bioconductor MiloR. Neighborhoods overlap, so counts of significant neighborhoods should be treated as a map of recurring local signal rather than independent discoveries.

## Results

| Run | Neighborhoods tested | FDR < 0.10 | Top neighborhood | Direction |
| --- | ---: | ---: | --- | --- |
| Full all-cell 4M | 20,000 | 3,673 | Mucous gland / respiratory secretory, FDR 1.62e-05 | Mostly negative age coefficients |
| Full 4M lineage-focused | 20,000 | 5,613 | Early iOSN / olfactory immature OSN, FDR 5.05e-04 | Strongly negative age coefficients |
| Full 4M secretory-focused | 20,000 | 285 | Mucous gland / respiratory secretory, FDR 8.78e-02 | Mixed direction within secretory-only denominator |

## Matched FLEX V2 / Device Sensitivity

The matched sensitivity restricts healthy donors to the same major chemistry and collection-method stratum as the AD/PD Gateway disease samples: FLEX v2 and device-guided collection. This is a much smaller healthy subset, with 27 donors, so the analysis is lower-powered but more technically controlled.

| Run | Neighborhoods tested | FDR < 0.10 | Top neighborhood | Direction |
| --- | ---: | ---: | --- | --- |
| Matched all-cell 4M | 20,000 | 10 | Naive CD8 / T cell, FDR 0.00524 | Mostly negative age coefficients |
| Matched 4M lineage-focused | 20,000 | 1 | Early iOSN / olfactory immature OSN, FDR 0.0427 | Negative age coefficient |
| Matched 4M secretory-focused | 20,000 | 0 | Olfactory sustentacular nominal top, FDR 0.437 | Not significant |

The matched all-cell significant neighborhoods were mostly Naive CD8/T-cell neighborhoods with negative age coefficients. One late iOSN neighborhood also survived in the matched all-cell analysis, and one Early iOSN neighborhood survived in the matched lineage-focused analysis. The secretory-only neighborhood signal did not survive matched correction.

## Biological Reading

The all-cell and lineage-focused full-scale runs now support a real neighborhood-level aging signal that was not visible in the 250k/100k pilots. The strongest recurring all-donor signal is reduced age-associated representation of regenerative neuronal-lineage neighborhoods, especially early iOSN, late iOSN, INP, and related HBC/suprabasal neighborhoods. The broad all-cell run also shows age-associated shifts in mucous gland/secretory, multiciliated, dendritic/T-cell, Bowman gland, sustentacular, and mature neuronal neighborhoods.

The matched FLEX v2/device sensitivity makes the interpretation more conservative. Under strict technical matching, the broad number of significant neighborhoods shrinks sharply, but a negative Early iOSN/iOSN signal remains detectable. This is the most manuscript-ready mechanistic result from the Milo-style analyses. The broad all-donor neighborhood map should be used to generate biological hypotheses and guide annotation, while strict matched results should govern the claim language.

The secretory-only run is weaker and directionally mixed because the denominator is restricted to secretory/sustentacular/glandular cells rather than all cells. It should be used as a compartment-specific sensitivity view, not as the primary secretory claim.

## Claim Status

Supported as an exploratory mechanistic layer:

- Full-scale scVI neighborhoods show age-associated abundance shifts after donor-level adjustment.
- The most coherent and technically matched lineage signal points toward reduced immature olfactory neuronal/regenerative neighborhoods with age.

Still gated before main-text promotion:

- marker/program annotation;
- age-bin robustness;
- comparison with official MiloR or a clearly documented reason to keep the Python implementation;
- replication against the 250k seed and lineage-focused sensitivity runs.

Use the matched secretory result as a guardrail: do not promote a secretory-specific differential-abundance claim from Milo-style neighborhoods unless follow-up marker/program or pseudotime analyses independently support it.
