# GSE309325 Organoid Perturbation Adapter Status

Updated: 2026-06-25

## Scope

GSE309325 is a human nasal organoid scRNA-seq perturbation dataset containing mock samples and SARS-CoV-2 timepoints. It is olfactory-relevant because the organoid model includes olfactory and nasal respiratory epithelial compartments, but it is infection context, not aging or adult donor validation.

## Adapter Outputs

- Samples scored: 7
- Modules scored per sample: 24
- Perturbation contrasts: 96

## Sample QC

| Sample | Condition | Day | Cells | Selected genes present |
| --- | --- | ---: | ---: | ---: |
| GSM9265687 | mock | 0 | 5227 | 153 |
| GSM9265688 | mock | 0 | 8506 | 153 |
| GSM9265689 | mock | 0 | 8964 | 153 |
| GSM9265690 | sars_cov_2 | 1 | 8236 | 153 |
| GSM9265691 | sars_cov_2 | 2 | 6594 | 153 |
| GSM9265692 | sars_cov_2 | 7 | 8371 | 153 |
| GSM9265693 | sars_cov_2 | 11 | 6655 | 153 |

## Largest Module Shifts Versus Mock

| Sample | Day | Module | Direction | Delta vs mock | z vs mock |
| --- | ---: | --- | --- | ---: | ---: |
| GSM9265691 | 2 | interferon_response | increased | 1.036 | 351.1 |
| GSM9265690 | 1 | interferon_response | increased | 0.3265 | 110.6 |
| GSM9265691 | 2 | progenitor_neuroblast | increased | 0.2857 | 6.789 |
| GSM9265690 | 1 | immature_neuron | decreased | -0.2708 | -1.778 |
| GSM9265693 | 11 | immature_neuron | increased | 0.2653 | 1.742 |
| GSM9265690 | 1 | gap43_dcx_immature_osn | decreased | -0.2561 | -1.71 |
| GSM9265693 | 11 | gap43_dcx_immature_osn | increased | 0.2438 | 1.627 |
| GSM9265691 | 2 | ascl1_neurog_neurod_neurogenesis | increased | 0.2436 | 11 |
| GSM9265691 | 2 | ad_amyloid_tau | increased | 0.2253 | 6.075 |
| GSM9265691 | 2 | notch_fate_control | increased | 0.2065 | 5.159 |
| GSM9265691 | 2 | oxidative_stress_detox | increased | 0.1977 | 10.77 |
| GSM9265691 | 2 | pd_synuclein_mito | increased | 0.1878 | 6.908 |

## Claim Boundary

Use this as perturbation-context evidence only. Infected timepoints are single samples without matched biological replicates, so effect sizes are descriptive module shifts rather than inferential p-values.
