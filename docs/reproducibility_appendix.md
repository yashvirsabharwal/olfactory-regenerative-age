# ORA Reproducibility Appendix

Updated: 2026-06-22

## Environment

Use the local virtual environment for Python commands:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[dev]"
PYTHON=.venv/bin/python make test
```

The local R/Bioconductor environment for genome-wide edgeR, limma-voom, and optional DESeq2 work is configured through `R_ENV=.mamba/ora-r` and `RSCRIPT=$(HOME)/.local/bin/micromamba run -p .mamba/ora-r Rscript`.

The official MiloR subset sensitivity was run on `mia` with a user-level micromamba environment because source installation under system R 4.6 lacked `fontconfig`/`freetype` development headers:

```bash
mkdir -p ~/bin ~/tmp/micromamba
cd ~/tmp/micromamba
curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest -o micromamba.tar.bz2
tar -xjf micromamba.tar.bz2 bin/micromamba
cp bin/micromamba ~/bin/micromamba
~/bin/micromamba create -y -p ~/micromamba-envs/ora-milor -c conda-forge -c bioconda bioconductor-milor r-data.table r-optparse
```

## Core Real-Data Command Order

```bash
make inspect
make cohort
make aggregate
make features
make age-associations
make model-ora
make modules
make features-augmented
make model-ora-augmented
make project-ndd
make pseudobulk
make pseudobulk-covariate-de
make pseudobulk-genomewide
make pseudobulk-genomewide-qc
make pseudobulk-genomewide-edger
make pseudobulk-genomewide-de-summary
make pseudobulk-genomewide-de-audit
make pseudobulk-genomewide-edger-matched
make pseudobulk-genomewide-de-summary-matched
make pseudobulk-genomewide-de-audit-matched
make pseudobulk-genomewide-limma
make pseudobulk-genomewide-limma-de-summary
make pseudobulk-genomewide-limma-de-audit
make pseudobulk-genomewide-limma-matched
make pseudobulk-genomewide-limma-de-summary-matched
make pseudobulk-genomewide-limma-de-audit-matched
make external-validation
make external-gse184117-modules
make external-gse184117-markers
make external-marker-age-concordance
make external-evidence
make feature-interpretation
make latent-space-audit
make latent-space-recompute-plan
make scvi-pilot
make scvi-pilot-validation
make scvi-embedding-claim-gates
make model-card
make output-provenance
make report
```

The completed bounded scVI pilot used:

```bash
PYTHON=.venv/bin/python make scvi-pilot
PYTHON=.venv/bin/python make scvi-pilot-validation
```

The pilot output H5AD is local/ignored. Its validation table is tracked through the provenance manifest as `results/tables/scvi_pilot_validation.tsv`.

The completed scVI embedding comparison used:

```bash
PYTHON=.venv/bin/python make scvi-embedding-claim-gates
```

This writes the publication-facing latent claim gates to `results/tables/scvi_embedding_claim_gates.tsv`, marker-concordance gates to `results/tables/scvi_embedding_marker_concordance.tsv`, and the summary note `docs/scvi_embedding_comparison.md`.

The completed full 4M neighborhood parity checks used the remote full reduced scVI H5AD and membership tables:

```bash
PYTHON=.venv/bin/python make milo-full-4m-lineage-age-bins
PYTHON=.venv/bin/python make milo-full-4m-lineage-matched-age-bins
PYTHON=.venv/bin/python make milo-full-4m-lineage-edger-parity
PYTHON=.venv/bin/python make milo-full-4m-lineage-matched-edger-parity
PYTHON=.venv/bin/python make milor-lineage-subset-parity
PYTHON=.venv/bin/python make milor-lineage-matched-subset-parity
```

## Generated Artifact Policy

Large inputs, processed data, result tables, figures, and reports are ignored by Git. Use `make output-provenance` to write `results/reports/command_manifest.tsv` and `results/reports/output_provenance.tsv` with existence, size, modification time, and checksums for small artifacts.
