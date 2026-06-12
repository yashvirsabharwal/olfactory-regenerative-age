# Genome-Wide Pseudobulk DE Hooks

This workflow starts after `make pseudobulk-genomewide` exports:

- `data/processed/pseudobulk_genomewide_counts.tsv.gz`
- `data/processed/pseudobulk_genomewide_metadata.tsv`
- `data/processed/pseudobulk_genomewide_genes.tsv`
- `results/tables/pseudobulk_genomewide_summary.tsv`

The count matrix is gene-by-pseudobulk-group. Metadata rows match matrix columns by `pseudobulk_id`.

## Recommended Design

Run models within each `fine_cell_type` and compare disease groups only where donor counts are sufficient.

Base design:

```r
~ disease_group + age + sex + chemistry + collection_method
```

Add `site` only when available and not aliased with disease or chemistry:

```r
~ disease_group + age + sex + chemistry + collection_method + site
```

Primary contrasts:

```r
ad_vs_healthy
pd_vs_healthy
```

## edgeR Quasi-Likelihood Template

```r
library(edgeR)

counts <- read.delim("data/processed/pseudobulk_genomewide_counts.tsv.gz", check.names = FALSE)
meta <- read.delim("data/processed/pseudobulk_genomewide_metadata.tsv", check.names = FALSE)

rownames(counts) <- counts$gene_id
count_mat <- as.matrix(counts[, meta$pseudobulk_id, drop = FALSE])
meta$disease_group <- relevel(factor(meta$disease_group), ref = "healthy")

run_state <- function(state) {
  keep_samples <- meta$fine_cell_type == state
  y <- DGEList(counts = count_mat[, keep_samples, drop = FALSE], samples = meta[keep_samples, ])
  design <- model.matrix(~ disease_group + age + sex + chemistry + collection_method, data = y$samples)
  keep_genes <- filterByExpr(y, design = design)
  y <- calcNormFactors(y[keep_genes, , keep.lib.sizes = FALSE])
  fit <- glmQLFit(estimateDisp(y, design), design)
  qlf_ad <- glmQLFTest(fit, coef = "disease_groupad")
  qlf_pd <- glmQLFTest(fit, coef = "disease_grouppd")
  list(ad_vs_healthy = topTags(qlf_ad, n = Inf)$table, pd_vs_healthy = topTags(qlf_pd, n = Inf)$table)
}
```

## limma-voom Template

```r
library(edgeR)
library(limma)

counts <- read.delim("data/processed/pseudobulk_genomewide_counts.tsv.gz", check.names = FALSE)
meta <- read.delim("data/processed/pseudobulk_genomewide_metadata.tsv", check.names = FALSE)

rownames(counts) <- counts$gene_id
count_mat <- as.matrix(counts[, meta$pseudobulk_id, drop = FALSE])
meta$disease_group <- relevel(factor(meta$disease_group), ref = "healthy")

run_state <- function(state) {
  keep_samples <- meta$fine_cell_type == state
  design <- model.matrix(~ disease_group + age + sex + chemistry + collection_method, data = meta[keep_samples, ])
  y <- DGEList(counts = count_mat[, keep_samples, drop = FALSE])
  keep_genes <- filterByExpr(y, design = design)
  y <- calcNormFactors(y[keep_genes, , keep.lib.sizes = FALSE])
  v <- voom(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  list(
    ad_vs_healthy = topTable(fit, coef = "disease_groupad", number = Inf),
    pd_vs_healthy = topTable(fit, coef = "disease_grouppd", number = Inf)
  )
}
```

## DESeq2 Template

```r
library(DESeq2)

counts <- read.delim("data/processed/pseudobulk_genomewide_counts.tsv.gz", check.names = FALSE)
meta <- read.delim("data/processed/pseudobulk_genomewide_metadata.tsv", check.names = FALSE)

rownames(counts) <- counts$gene_id
count_mat <- round(as.matrix(counts[, meta$pseudobulk_id, drop = FALSE]))
meta$disease_group <- relevel(factor(meta$disease_group), ref = "healthy")

run_state <- function(state) {
  keep_samples <- meta$fine_cell_type == state
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat[, keep_samples, drop = FALSE],
    colData = meta[keep_samples, ],
    design = ~ age + sex + chemistry + collection_method + disease_group
  )
  keep_genes <- rowSums(counts(dds) >= 10) >= 3
  dds <- DESeq(dds[keep_genes, ])
  list(
    ad_vs_healthy = as.data.frame(results(dds, contrast = c("disease_group", "ad", "healthy"))),
    pd_vs_healthy = as.data.frame(results(dds, contrast = c("disease_group", "pd", "healthy")))
  )
}
```

## Sensitivity Hooks

Re-run the export and R models after stratifying or filtering metadata by:

- `chemistry`
- `collection_method`
- larger `--min-cells-per-group`
- larger `--min-donors-per-cell-state`
- healthy-only age contrasts
