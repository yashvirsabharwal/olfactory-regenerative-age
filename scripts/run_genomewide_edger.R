#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(optparse)
})

option_list <- list(
  make_option("--counts", default = "data/processed/pseudobulk_genomewide_counts.tsv.gz"),
  make_option("--metadata", default = "data/processed/pseudobulk_genomewide_metadata.tsv"),
  make_option("--manifest", default = "data/processed/cohort_manifest.tsv"),
  make_option("--out", default = "results/tables/pseudobulk_genomewide_edger.tsv.gz"),
  make_option("--summary-out", dest = "summary_out", default = "results/tables/pseudobulk_genomewide_edger_summary.tsv"),
  make_option("--contrasts", default = "ad:healthy,pd:healthy"),
  make_option("--min-donors", dest = "min_donors", type = "integer", default = 3),
  make_option("--min-cells", dest = "min_cells", type = "integer", default = 10)
)
opt <- parse_args(OptionParser(option_list = option_list))

ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

parse_contrasts <- function(value) {
  parts <- strsplit(value, ",", fixed = TRUE)[[1]]
  lapply(parts[nzchar(parts)], function(item) strsplit(item, ":", fixed = TRUE)[[1]])
}

first_non_missing <- function(x) {
  y <- x[!is.na(x) & x != ""]
  if (length(y) == 0) NA else y[[1]]
}

donor_state_counts <- function(count_mat, meta_state) {
  donor_ids <- meta_state$donor_id
  summed <- rowsum(t(count_mat[, meta_state$pseudobulk_id, drop = FALSE]), group = donor_ids, reorder = FALSE)
  t(summed)
}

donor_state_metadata <- function(meta_state, manifest) {
  base <- meta_state[, .(
    disease_group = first_non_missing(disease_group),
    fine_cell_type = first_non_missing(fine_cell_type),
    n_cells = sum(as.numeric(n_cells), na.rm = TRUE),
    sum_n_counts = sum(as.numeric(sum_n_counts), na.rm = TRUE)
  ), by = donor_id]
  donor_manifest <- unique(manifest[, .(
    donor_id,
    age,
    sex,
    race_ethnicity,
    chemistry,
    collection_method,
    site,
    total_cells
  )], by = "donor_id")
  merge(base, donor_manifest, by = "donor_id", all.x = TRUE)
}

build_design <- function(meta, case_group, covariates = c("age", "sex", "chemistry", "collection_method", "site")) {
  meta <- copy(meta)
  meta[, disease_group := factor(disease_group, levels = c("healthy", case_group))]
  usable_covariates <- c()
  for (covariate in covariates) {
    if (!covariate %in% names(meta)) next
    values <- meta[[covariate]]
    if (is.numeric(values) || is.integer(values)) {
      if (sum(!is.na(values)) >= 2 && length(unique(values[!is.na(values)])) > 1) {
        fill <- median(values, na.rm = TRUE)
        meta[[covariate]] <- ifelse(is.na(values), fill, values)
        usable_covariates <- c(usable_covariates, covariate)
      }
    } else {
      values <- as.character(values)
      if (length(unique(values[!is.na(values) & values != ""])) > 1) {
        meta[[covariate]] <- factor(ifelse(is.na(values) | values == "", "unknown", values))
        usable_covariates <- c(usable_covariates, covariate)
      }
    }
  }
  while (TRUE) {
    formula_text <- paste("~ disease_group", paste(usable_covariates, collapse = " + "), sep = ifelse(length(usable_covariates), " + ", ""))
    design <- model.matrix(as.formula(formula_text), data = meta)
    if (qr(design)$rank == ncol(design)) {
      return(list(design = design, covariates = usable_covariates, coefficient = paste0("disease_group", case_group)))
    }
    if (length(usable_covariates) == 0) {
      return(list(design = design[, qr(design)$pivot[seq_len(qr(design)$rank)], drop = FALSE], covariates = usable_covariates, coefficient = NA_character_))
    }
    usable_covariates <- usable_covariates[-length(usable_covariates)]
  }
}

run_one <- function(counts, gene_info, metadata, manifest, state, case_group, control_group, min_donors, min_cells) {
  meta_state <- metadata[fine_cell_type == state & disease_group %in% c(case_group, control_group) & n_cells >= min_cells]
  donor_meta <- donor_state_metadata(meta_state, manifest)
  donor_meta <- donor_meta[disease_group %in% c(case_group, control_group)]
  disease_values <- as.character(donor_meta$disease_group)
  n_case <- uniqueN(donor_meta[!is.na(disease_values) & disease_values == case_group]$donor_id)
  n_control <- uniqueN(donor_meta[!is.na(disease_values) & disease_values == control_group]$donor_id)
  if (n_case < min_donors || n_control < min_donors) {
    return(list(result = NULL, summary = data.table(
      contrast = paste0(case_group, "_vs_", control_group),
      fine_cell_type = state,
      n_case = n_case,
      n_control = n_control,
      n_genes_tested = 0L,
      covariates = "",
      status = "too_few_donors"
    )))
  }
  donor_counts <- donor_state_counts(counts, meta_state)
  common_donors <- intersect(colnames(donor_counts), donor_meta$donor_id)
  donor_counts <- donor_counts[, common_donors, drop = FALSE]
  donor_meta <- donor_meta[match(common_donors, donor_id)]
  design_info <- build_design(donor_meta, case_group)
  coefficient <- design_info$coefficient
  if (is.na(coefficient) || !coefficient %in% colnames(design_info$design)) {
    return(list(result = NULL, summary = data.table(
      contrast = paste0(case_group, "_vs_", control_group),
      fine_cell_type = state,
      n_case = n_case,
      n_control = n_control,
      n_genes_tested = 0L,
      covariates = paste(design_info$covariates, collapse = ","),
      status = "invalid_design"
    )))
  }
  y <- DGEList(counts = donor_counts, samples = donor_meta)
  keep <- filterByExpr(y, design = design_info$design)
  if (!any(keep)) {
    return(list(result = NULL, summary = data.table(
      contrast = paste0(case_group, "_vs_", control_group),
      fine_cell_type = state,
      n_case = n_case,
      n_control = n_control,
      n_genes_tested = 0L,
      covariates = paste(design_info$covariates, collapse = ","),
      status = "no_genes_after_filter"
    )))
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  fit <- glmQLFit(estimateDisp(y, design_info$design), design_info$design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = coefficient)
  table <- as.data.table(topTags(qlf, n = Inf)$table, keep.rownames = "gene_id")
  table <- merge(gene_info, table, by = "gene_id", all.y = TRUE)
  setnames(table, old = c("logFC", "logCPM", "PValue", "FDR"), new = c("log2fc", "logcpm", "p_value", "fdr"), skip_absent = TRUE)
  table[, `:=`(
    contrast = paste0(case_group, "_vs_", control_group),
    case_group = case_group,
    control_group = control_group,
    fine_cell_type = state,
    n_case = n_case,
    n_control = n_control,
    covariates = paste(design_info$covariates, collapse = ","),
    status = "ok"
  )]
  list(result = table[, .(
    contrast, case_group, control_group, fine_cell_type, gene_id, gene_symbol, gene_index,
    n_case, n_control, log2fc, logcpm, F, p_value, fdr, covariates, status
  )], summary = data.table(
    contrast = paste0(case_group, "_vs_", control_group),
    fine_cell_type = state,
    n_case = n_case,
    n_control = n_control,
    n_genes_tested = nrow(table),
    covariates = paste(design_info$covariates, collapse = ","),
    status = "ok"
  ))
}

message("Reading genome-wide counts: ", opt$counts)
counts_dt <- fread(opt$counts)
gene_info <- counts_dt[, .(gene_id, gene_symbol, gene_index)]
sample_cols <- setdiff(names(counts_dt), c("gene_id", "gene_symbol", "gene_index"))
counts <- as.matrix(counts_dt[, ..sample_cols])
rownames(counts) <- counts_dt$gene_id
storage.mode(counts) <- "integer"
rm(counts_dt)

metadata <- fread(opt$metadata)
manifest <- fread(opt$manifest)
metadata <- metadata[pseudobulk_id %in% sample_cols]
setkey(metadata, pseudobulk_id)
metadata <- metadata[sample_cols]
metadata[, disease_group := as.character(disease_group)]
metadata[, fine_cell_type := as.character(fine_cell_type)]

contrasts <- parse_contrasts(opt$contrasts)
states <- sort(unique(metadata$fine_cell_type))
results <- list()
summaries <- list()
counter <- 1L
for (contrast in contrasts) {
  case_group <- contrast[[1]]
  control_group <- contrast[[2]]
  for (state in states) {
    message("Running ", case_group, "_vs_", control_group, " / ", state)
    fit <- run_one(counts, gene_info, metadata, manifest, state, case_group, control_group, opt$min_donors, opt$min_cells)
    if (!is.null(fit$result)) {
      results[[length(results) + 1L]] <- fit$result
    }
    summaries[[length(summaries) + 1L]] <- fit$summary
    counter <- counter + 1L
  }
}

ensure_parent(opt$out)
ensure_parent(opt$summary_out)
if (length(results)) {
  fwrite(rbindlist(results, use.names = TRUE), opt$out, sep = "\t")
} else {
  fwrite(data.table(), opt$out, sep = "\t")
}
fwrite(rbindlist(summaries, use.names = TRUE), opt$summary_out, sep = "\t")
message("Wrote edgeR results: ", opt$out)
message("Wrote edgeR summary: ", opt$summary_out)
