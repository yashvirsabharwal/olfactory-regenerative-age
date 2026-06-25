# Niche Signaling Feasibility

Updated: 2026-06-25

## Decision

Primary local method: curated ligand-receptor pseudobulk scoring by donor, sender group, and receiver group. This is a transparent hypothesis layer, not a full LIANA/NicheNet/CellChat/CellPhoneDB inference run.

## Scope

- Curated interaction families: 12
- Interaction definitions: 12
- Minimum gene coverage: 0.889
- Mean gene coverage: 0.984

## Top Ranked Niche Hypotheses

| Family | Sender | Receiver | Coverage | Age beta/10y | Age FDR | ORAA r | ORAA FDR |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| IFN | immune | mature_osn | 0.923 | -8.81e-02 | 5.04e-05 | -4.19e-01 | 1.23e-07 |
| IFN | immune | hbc | 0.923 | -8.07e-02 | 6.67e-04 | -4.06e-01 | 5.13e-07 |
| CCL | immune | hbc | 1.000 | -2.93e-02 | 0.033 | -3.63e-01 | 8.95e-06 |
| IFN | immune | immature_osn | 0.923 | -9.02e-02 | 0.002 | -3.73e-01 | 4.01e-04 |
| TNF | immune | immature_osn | 1.000 | -9.98e-02 | 0.016 | -4.21e-01 | 9.34e-05 |
| IFN | immune | gbc_inp | 0.923 | -1.25e-01 | 3.30e-05 | -3.44e-01 | 0.050 |
| IFN | immune | sustentacular | 0.923 | -9.27e-02 | 0.001 | -3.22e-01 | 0.003 |
| MIF | immune | immature_osn | 1.000 | -1.06e-01 | 0.034 | -3.70e-01 | 4.01e-04 |
| MIF | sustentacular | immature_osn | 1.000 | -1.05e-01 | 0.006 | -3.18e-01 | 0.003 |
| CCL | respiratory_secretory | hbc | 1.000 | -1.41e-02 | 0.033 | -2.89e-01 | 7.26e-04 |

## Interpretation Guardrail

These results nominate sender-receiver signaling hypotheses. They do not prove physical proximity, receptor activation, or causal niche control. Spatial transcriptomics, perturbation, or a full LR method with permutation/background testing would be required to promote the claims.

