# Regulatory Driver Feasibility

Updated: 2026-06-25

## Decision

Primary local method: curated TF/pathway target programs scored from genomewide pseudobulk logCPM by donor and cell state. This satisfies the first-pass driver-hypothesis layer without requiring new package installs or per-cell regulon inference.

Deferred heavy methods: decoupler with DoRothEA/CollecTRI or PROGENy, pySCENIC/SCENIC+, and chromatin/motif workflows. These remain useful upgrades, but they require additional databases/packages and stronger compute/runtime validation.

## Local Tool Status

- `scanpy`: installed.
- `decoupler`, `pyscenic`, `gseapy`, `omnipath`: not installed in the current `.venv` at feasibility time.

## Lineage Groups

- HBC
- GBC/INP
- Immature OSN
- Mature OSN
- Sustentacular
- Immune
- Respiratory/secretory
- All lineages

## Gene Coverage

- Driver target sets: 15
- Minimum coverage fraction: 0.833
- Mean coverage fraction: 0.984

## Top Ranked First-Pass Drivers

| Driver | Theme | Coverage | Top age lineage | Top age FDR | Top ORAA lineage | Top ORAA FDR |
| --- | --- | ---: | --- | ---: | --- | ---: |
| WNT/CTNNB1 | neural_progenitor | 1.000 | hbc | 0.004 | all_lineages | 2.39e-04 |
| ASCL1 | neural_progenitor | 1.000 | hbc | 0.004 | hbc | 5.84e-04 |
| STAT/IRF | immune_inflammatory | 1.000 | gbc_inp | 0.003 | mature_osn | 0.002 |
| TP63 | basal_quiescence | 1.000 | gbc_inp | 0.001 | all_lineages | 0.005 |
| OMP/ADCY3 | mature_osn | 1.000 | gbc_inp | 1.41e-04 | immune | 0.101 |
| NOTCH/HES | basal_activation | 1.000 | gbc_inp | 0.001 | immune | 0.082 |
| YAP/TEAD | basal_activation | 1.000 | hbc | 0.023 | all_lineages | 0.020 |
| NF-kB/IL17/TNF | immune_inflammatory | 1.000 | immature_osn | 0.112 | respiratory_secretory | 0.005 |

## Interpretation Guardrail

These rows are ranked driver hypotheses from curated target programs. They are not causal perturbation evidence and should not be described as measured regulons until an external prior method such as decoupler/CollecTRI or SCENIC-style inference is run and sensitivity-checked.

