# Donor-Level Compositional Age Modeling Plan

## Purpose

ORA uses cell-state proportions, centered-log-ratio (CLR) features, and lineage
ratios. Because these inputs are compositional and sum-constrained, M1.2 adds a
formal donor-level compositional analysis that tests whether fine-cell-state
composition changes with age after available technical covariates are included.

## Primary response

The primary response is the donor-level fine-cell-state composition matrix from
`data/processed/donor_cell_state_counts.tsv`. Counts are aggregated by
`donor_id` and `fine_cell_type`, then transformed with a 0.5 pseudocount CLR:

`clr_ij = log(p_ij) - mean_j(log(p_ij))`

where `p_ij` is the pseudocount-adjusted fraction of cell state `j` in donor
`i`.

## Primary predictor

The primary predictor is continuous donor age in years. Reported coefficients
are scaled to a 10-year age interval.

## Covariates

The primary model includes available donor-level covariates when they have
usable variation in the scenario:

- `sex`
- `chemistry`
- `collection_method`
- `site`
- `log10_total_cells`

`site` is retained in the model specification but drops out automatically for
the current Gateway cohort because site is missing or invariant in the ORA
training subset.

## Compositional baseline

The primary baseline is the donor-specific CLR geometric mean across included
fine-cell states. This avoids choosing a single reference cell state and keeps
the model aligned with the ORA feature construction. A Bayesian multinomial
method such as scCODA remains a later optional upgrade, but the current fallback
is deliberately simple, reproducible, and implemented with the existing Python
stack.

## Scenarios

The workflow writes both primary and sensitivity rows:

- `primary_all_healthy`: donors usable for primary healthy ORA training.
- `strict_threshold`: primary donors that also pass strict total-cell,
  lineage-cell, and mature-neuron thresholds.
- `single_<chemistry>_<collection>`: the largest single chemistry and
  collection-method stratum, if it has at least the minimum scenario size.

Scenarios with too few donors are retained in the output table with a
`too_few_samples` status rather than being silently dropped.

## Directional comparison

Each modeled cell state is mapped back to its ORA CLR feature name. When
`results/tables/age_cell_state_associations.tsv` is available, the workflow
adds the ordinary age-association beta, FDR, and direction-concordance flag. The
goal is to identify cell states whose composition-aware age direction supports
the ORA feature direction and to expose discordant states.

## Artifacts

- `scripts/run_compositional_age_model.py`
- `results/tables/compositional_age_model_summary.tsv`
- `results/tables/compositional_age_model_sensitivity.tsv`
- `results/figures/extended_data_compositional_model.pdf`
- `results/figures/extended_data_compositional_model.png`

## Interpretation boundary

This model supports or qualifies donor-level compositional age associations. It
does not measure cell lineage flux, prove a universal biological clock, or
establish clinical utility.
