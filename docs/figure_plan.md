# ORA Figure Plan

Updated: 2026-06-16

## Main Figures

1. Cohort and claim-gated workflow: donor groups, healthy training split, AD/PD frozen projection, external sanity checks, DE audits, and latent-space gates.
2. Age-associated composition: top donor-level cell-state associations, lineage ratios, and theme labels that explain the regenerative-state signal.
3. ORA modeling: repeated-CV performance, shuffled-age null comparison, calibration caveat, and best-model residual behavior.
4. Stable biology and modules: top ORA feature themes, module coverage, module-augmented model deltas, and feature interpretation cautions.
5. External validation and NDD guardrails: GSE184117 marker-age concordance, evidence-ledger strength, NDD projection sensitivity, and label-permutation framing.
6. Genome-wide DE and latent readiness: edgeR/limma parity, sentinel-gene audit flags, matched-subset sensitivity, and scVI pilot validation.

## Figure Style Direction

- Use a clean off-white background, dark graphite text, restrained grid lines, and a mixed palette with teal, blue, vermilion, gold, and gray rather than a one-hue theme.
- Prefer compact multi-panel figures with direct labels and short subtitles over dense legends.
- Use consistent typography, high DPI, fixed panel dimensions, and figure captions that state the claim gate.
- Show uncertainty or validation status wherever possible: intervals, null bands, status pills, or audited/not-audited labels.
- Avoid decorative backgrounds. The figures should feel like a high-end computational biology paper, not a dashboard screenshot.

## Extended Data

- Full model-card table.
- Top ORA feature biological interpretation table with theme and caution fields.
- ORA residual diagnostics by sex, chemistry, collection method, site, and yield.
- NDD projection appendix with donor-level ORAA and matched FLEX v2/device reference.
- External validation readiness, evidence-ledger, and harmonization reports.
- scVI pilot validation and latent-space recovery plan.

## Deferred Figures

- External validation replication figures after `GSE184117` adapter outputs real donor features.
- Trajectory, Milo, or cNMF figures after a scaled or lineage-focused latent-space validation passes.
