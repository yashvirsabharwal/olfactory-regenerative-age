# Cross-Tissue Age-Effect Estimates

Updated: 2026-06-25

## Scope

This analysis estimates donor-level age effects in selected public CELLxGENE nasal, bronchial, and lung comparator datasets. The primary claim gate uses adult donors only; fetal or pediatric lung stages are retained only as context.

## Comparator Assets

| Dataset | Tissue | Cells | Donors | Adult donors | Age range | Status | Notes |
| --- | --- | ---: | ---: | ---: | --- | --- | --- |
| Xu nasal scRNA-Seq Dataset | nasal | 34833 | 7 | 7 | 30.0-66.0 | ok | adult donor-level age metadata available |
| Xu bronchial scRNA-Seq Dataset | bronchial | 2075 | 17 | 17 | 51.0-87.0 | ok | adult donor-level age metadata available |
| LungMAP human lungs across age groups | lung | 46500 | 9 | 3 | 0.6-31.0 | ok | non-adult stages present; adult-only estimates are primary; adult donor count below primary age-effect threshold |

## Measured Adult Effects

| Dataset | Tissue | Tested adult effects |
| --- | --- | ---: |
| d2bc2703-a9b5-4949-b85d-61da7b9961b1 | bronchial | 36 |
| e3cc5b85-62ec-4add-bee3-73c8fd11e59d | nasal | 35 |

Nominal adult comparator effects with p<0.05: 3.
Adult comparator effects with within-dataset/scope FDR<0.05: 0.

## ORA Feature-Level Status

| Status | ORA features |
| --- | ---: |
| measured_adult_comparator_age_effect | 74 |
| expected_no_cross_tissue_cell_state_mapping | 27 |
| no_cross_tissue_feature_mapping | 9 |
| mapped_but_no_adult_comparator_effect | 6 |
| measured_airway_context_for_olfactory_feature | 4 |

## Interpretation Guardrails

- Adult nasal and bronchial estimates are donor-level comparator evidence for shared airway epithelial and immune programs.
- The LungMAP asset contains prenatal/child/adult stages, so its all-stage estimates are context-only unless an adult-only subset reaches the donor threshold.
- These estimates test cross-tissue sharing of feature families; they do not turn ORA into an absolute age clock or prove causal regeneration mechanisms.
