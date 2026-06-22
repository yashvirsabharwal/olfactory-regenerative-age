#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
})

parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  if (length(args) == 0) return(out)
  idx <- 1
  while (idx <= length(args)) {
    key <- sub("^--", "", args[[idx]])
    key <- gsub("-", "_", key)
    if (!key %in% names(out) || idx == length(args)) {
      stop("Invalid or incomplete argument: ", args[[idx]])
    }
    out[[key]] <- args[[idx + 1]]
    idx <- idx + 2
  }
  out
}

opt <- parse_args(list(
  counts = "data/processed/milo_full_4m_lineage_neighborhood_counts.tsv",
  design = "data/processed/milo_full_4m_lineage_neighborhood_design.tsv",
  out = "results/tables/milo_full_4m_lineage_edger_parity.tsv",
  summary_out = "results/tables/milo_full_4m_lineage_edger_parity_run_summary.tsv"
))

ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

add_covariate_if_usable <- function(design_frame, covariates, column, numeric = FALSE) {
  if (!column %in% names(design_frame)) return(list(design_frame = design_frame, covariates = covariates))
  values <- design_frame[[column]]
  if (numeric) {
    values <- suppressWarnings(as.numeric(values))
    if (sum(!is.na(values)) >= 3 && length(unique(values[!is.na(values)])) > 1) {
      fill <- median(values, na.rm = TRUE)
      design_frame[[column]] <- ifelse(is.na(values), fill, values)
      return(list(design_frame = design_frame, covariates = c(covariates, column)))
    }
  } else {
    values <- as.character(values)
    usable <- values[!is.na(values) & values != ""]
    if (length(unique(usable)) > 1) {
      design_frame[[column]] <- factor(ifelse(is.na(values) | values == "", "unknown", values))
      return(list(design_frame = design_frame, covariates = c(covariates, column)))
    }
  }
  list(design_frame = design_frame, covariates = covariates)
}

build_design <- function(sample_data) {
  sample_data <- as.data.frame(sample_data)
  sample_data$age_scaled <- as.numeric(sample_data$age_scaled)
  covariates <- c()
  update <- add_covariate_if_usable(sample_data, covariates, "sex")
  sample_data <- update$design_frame
  covariates <- update$covariates
  update <- add_covariate_if_usable(sample_data, covariates, "chemistry")
  sample_data <- update$design_frame
  covariates <- update$covariates
  update <- add_covariate_if_usable(sample_data, covariates, "collection_method")
  sample_data <- update$design_frame
  covariates <- update$covariates
  update <- add_covariate_if_usable(sample_data, covariates, "log_total_cells", numeric = TRUE)
  sample_data <- update$design_frame
  covariates <- update$covariates
  while (TRUE) {
    formula_text <- paste("~ age_scaled", paste(covariates, collapse = " + "), sep = ifelse(length(covariates), " + ", ""))
    design <- model.matrix(as.formula(formula_text), data = sample_data)
    if (qr(design)$rank == ncol(design)) {
      return(list(design = design, covariates = covariates))
    }
    if (length(covariates) == 0) {
      keep <- qr(design)$pivot[seq_len(qr(design)$rank)]
      return(list(design = design[, keep, drop = FALSE], covariates = covariates))
    }
    covariates <- covariates[-length(covariates)]
  }
}

message("Reading neighborhood counts: ", opt$counts)
counts_dt <- read.delim(opt$counts, check.names = FALSE)
neighborhood_ids <- counts_dt$neighborhood_id
donor_cols <- setdiff(names(counts_dt), "neighborhood_id")
counts <- as.matrix(counts_dt[, donor_cols, drop = FALSE])
rownames(counts) <- as.character(neighborhood_ids)
storage.mode(counts) <- "integer"
rm(counts_dt)

design_dt <- read.delim(opt$design, check.names = FALSE)
design_dt$donor_id <- as.character(design_dt$donor_id)
common_donors <- intersect(donor_cols, design_dt$donor_id)
counts <- counts[, common_donors, drop = FALSE]
design_dt <- design_dt[match(common_donors, design_dt$donor_id), , drop = FALSE]
design_info <- build_design(design_dt)

message("Fitting edgeR QL model: ", nrow(counts), " neighborhoods x ", ncol(counts), " donors")
y <- DGEList(counts = counts, samples = design_dt)
y <- calcNormFactors(y)
fit <- glmQLFit(estimateDisp(y, design_info$design), design_info$design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "age_scaled")
result <- topTags(qlf, n = Inf, sort.by = "none")$table
result <- data.frame(neighborhood_id = as.integer(rownames(result)), result, check.names = FALSE)
result$covariates <- paste(design_info$covariates, collapse = ",")
result$n_donors <- ncol(counts)

summary <- data.frame(
  metric = c("neighborhoods", "donors", "fdr_lt_0_10", "covariates"),
  value = c(nrow(result), ncol(counts), sum(result$FDR < 0.10, na.rm = TRUE), paste(design_info$covariates, collapse = ",")),
  detail = c("edgeR rows tested", "donors in count matrix", "edgeR QL neighborhoods at BH FDR < 0.10", "covariates retained after rank checks"),
  check.names = FALSE
)

ensure_parent(opt$out)
write.table(result, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
ensure_parent(opt$summary_out)
write.table(summary, opt$summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote edgeR parity DA: ", opt$out, " (", nrow(result), " rows)")
message("Wrote edgeR parity summary: ", opt$summary_out)
