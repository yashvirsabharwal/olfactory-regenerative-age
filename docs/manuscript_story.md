# ORA Manuscript Story and Novelty

Updated: 2026-06-18

## One-Sentence Claims

Primary claim: healthy human olfactory epithelial cell-state composition and curated module features encode a modest but reproducible regenerative aging axis.

Secondary claim: AD/PD olfactory donor projections are directionally stable but exploratory, hypothesis-generating observations rather than disease biomarkers.

## Best Manuscript Frame

The strongest current framing is an atlas reanalysis and reproducible methods/resource paper with an aging-biology centerpiece. The paper should not be framed as a clinical biological-age clock, a diagnostic NDD classifier, or a trajectory paper until true reference-mapped external validation and seed-stable latent mechanism analyses are complete.

## Novelty

- The analysis turns a large human olfactory epithelial atlas into an interpretable donor-level aging readout focused on regenerative tissue state, rather than another generic methylation, blood, or black-box aging clock.
- The feature space is biologically legible: cell-state proportions, centered-log-ratio composition, lineage ratios, and curated modules map to supporting/secretory epithelium, immune/inflammatory state, stress/senescence, neuronal maturation, and basal/progenitor biology.
- External validation is handled as a claim-gated evidence ledger. GSE184117 now has marker-reference and Gateway scANVI/scArches-mapped donor features; scANVI mapping confidence is high, but the cohort remains small-n and direction concordance is mixed. GSE151973 remains bulk marker context.
- Disease projection and genome-wide DE are explicitly demoted unless matched/sensitivity/audit gates pass, which makes the manuscript safer and more credible.

## Literature Context

- `GSE184117` / Oliva et al. is the most actionable external olfactory aging source. GEO describes single-cell 3' RNA-seq from 6 adult olfactory epithelium biopsies: 3 normosmic and 3 presbyosmic older subjects, plus one culture sample; its summary reports inflammation-associated changes in olfactory epithelial stem cells of presbyosmic patients. Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184117; PubMed/PMC: https://pubmed.ncbi.nlm.nih.gov/34990409/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC8843745/.
- `GSE151973` / Fodoulian et al. is useful as olfactory-versus-respiratory epithelial marker context, not donor-level ORA validation. GEO describes bulk RNA-seq from adult olfactory and respiratory epithelium biopsies; the paper used bulk and single-cell context to localize SARS-CoV-2 entry genes to sustentacular/respiratory ciliated states. Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151973; PubMed/PMC: https://pubmed.ncbi.nlm.nih.gov/33251489/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC7685946/.
- The aging-clock field is moving toward cell-type-aware and single-cell transcriptomic models. PubMed currently lists `scAgeClock: a single-cell transcriptome-based human aging clock model using gated multi-head attention neural networks` in NPJ Aging, 2026, DOI `10.1038/s41514-026-00379-5`. ORA should be positioned differently: less as a universal high-accuracy age predictor, more as a tissue-specific, interpretable regenerative-state axis.
- Human immune-aging clocks are also advancing; PubMed lists a 2026 Immunity paper, `Human immune aging clock identifies RUNX1 as a decelerator of T cell senescence`, DOI `10.1016/j.immuni.2026.02.007`. This supports the idea that interpretable cell-type aging models are timely, while also reinforcing that ORA's immune signals need confounding audits.

## Publication Readiness Verdict

On track for a conservative manuscript if framed as a reproducible atlas reanalysis/resource with a healthy-aging primary claim. The project is not yet ready for a strong disease-biomarker, external-validation, or trajectory claim.

Strengths:

- Large real-data base: 202 donors and more than 4 million cells.
- Healthy-only donor-level modeling is benchmarked, repeated, tuned, null-tested, and calibrated.
- Output provenance, command manifest, model cards, claim ledger, limitations, external evidence ledger, and feature interpretation table are now generated or documented.
- Genome-wide DE has edgeR/limma parity and sentinel-gene/donor-balance audits.

Main blockers:

- No large independent donor-feature-ready validation dataset yet.
- GSE184117 is useful and now scANVI/scArches mapped, but only n=3 versus n=3 with mixed feature-direction concordance.
- AD/PD projection has only 5 donors per disease group and is fully FLEX v2/device.
- A chunked reduced all-cell scVI atlas now trains on all 4,028,275 Gateway cells using a 3,000-gene HVG/marker feature set and validates finite 10-dimensional `X_scvi`, strong label purity, and acceptable technical mixing. Trajectory, Milo, and cNMF remain gated until full-model, 250k-seed, and lineage-focused marker diagnostics are reconciled.

## Draft Title

An Interpretable Olfactory Regenerative Aging Axis from Human Single-Cell Epithelial Composition

## Draft Abstract

Human olfactory epithelium is a continually renewing neuroepithelium, but whether its donor-level cellular composition captures aging-related regenerative state remains unclear. We reanalyzed a large human olfactory epithelial single-cell atlas to construct donor-level composition, lineage-ratio, and curated module features, then trained healthy-donor olfactory regenerative age models under donor-level cross-validation. Composition and module features encoded a modest but reproducible aging signal above shuffled-age null models, with stable features mapping to supporting/secretory epithelium, immune/inflammatory state, stress/senescence, neuronal maturation, and basal/progenitor biology. A chunked reduced scVI workflow trained an all-cell latent atlas across 4,028,275 Gateway cells, providing a validated substrate for guarded neighborhood and lineage follow-up. GSE184117 was converted into marker-reference and Gateway scANVI/scArches-mapped donor features, but small-n concordance remains guarded. Frozen healthy-trained models projected AD and PD olfactory donors toward negative residual age across model families, but small disease cohorts and matched-context diagnostics keep this observation exploratory. Genome-wide pseudobulk analyses, external dataset checks, and scaled latent analyses are reported with explicit donor-balance, sentinel-gene, cross-method, and validation-strength gates. Together, the analysis provides a reproducible, interpretable framework for studying olfactory epithelial regenerative aging while separating supported aging biology from deferred disease and trajectory claims.

## One-Page Significance Draft

Most biological-age models are optimized as broad predictors, often using blood, methylation, or high-dimensional molecular signatures that are difficult to translate into tissue biology. The olfactory epithelium offers a different opportunity: it is a human neuroepithelium with ongoing cell turnover, sensory-neuron replacement, immune exposure, and age-related functional decline. ORA asks whether donor-level cellular composition in this tissue carries an interpretable aging signal. The answer is yes, but modestly and with useful caveats. The strongest paper is therefore not "we built a perfect clock"; it is "we exposed a measurable regenerative-state axis in human olfactory epithelium and built the reproducible guardrails needed to interpret it."

The manuscript should lead with healthy aging, interpretable features, and reproducibility. NDD projection can be included as an appendix-style anchoring analysis, and genome-wide DE can be framed as audited hypothesis generation. External validation remains the key next upgrade.
