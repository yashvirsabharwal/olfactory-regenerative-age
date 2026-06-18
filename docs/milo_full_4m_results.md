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

## Biological Reading

The all-cell and lineage-focused full-scale runs now support a real neighborhood-level aging signal that was not visible in the 250k/100k pilots. The strongest recurring signal is reduced age-associated representation of regenerative neuronal-lineage neighborhoods, especially early iOSN, late iOSN, INP, and related HBC/suprabasal neighborhoods. The broad all-cell run also shows age-associated shifts in mucous gland/secretory, multiciliated, dendritic/T-cell, Bowman gland, sustentacular, and mature neuronal neighborhoods.

The secretory-only run is weaker and directionally mixed because the denominator is restricted to secretory/sustentacular/glandular cells rather than all cells. It should be used as a compartment-specific sensitivity view, not as the primary secretory claim.

## Claim Status

Supported as an exploratory mechanistic layer:

- Full-scale scVI neighborhoods show age-associated abundance shifts after donor-level adjustment.
- The most coherent lineage signal points toward reduced immature olfactory neuronal/regenerative neighborhoods with age.

Still gated before main-text promotion:

- matched FLEX/device sensitivity;
- donor-yield and age-bin robustness;
- comparison with official MiloR or a clearly documented reason to keep the Python implementation;
- neighborhood marker/program annotation;
- replication against the 250k seed and lineage-focused sensitivity runs.
