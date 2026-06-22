#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(miloR)
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
  out$k <- as.integer(out$k)
  out$prop <- as.numeric(out$prop)
  out
}

opt <- parse_args(list(
  metadata = "data/processed/milor_lineage_subset_metadata.tsv",
  embedding = "data/processed/milor_lineage_subset_embedding.tsv",
  run_name = "lineage_full_subset",
  k = "100",
  prop = "0.05",
  out = "results/tables/milor_lineage_subset_da.tsv",
  summary_out = "results/tables/milor_lineage_subset_summary.tsv"
))

ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

top_label <- function(values) {
  values <- as.character(values)
  values[is.na(values) | values == ""] <- "unknown"
  tab <- sort(table(values), decreasing = TRUE)
  if (length(tab) == 0) return(list(label = "unknown", fraction = NA_real_))
  list(label = names(tab)[[1]], fraction = as.numeric(tab[[1]]) / sum(tab))
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

build_design <- function(metadata) {
  design <- unique(metadata[, c("donor_id", "age", "sex", "chemistry", "collection_method", "total_cells"), drop = FALSE])
  design$age <- as.numeric(design$age)
  design$age_scaled <- as.numeric(scale(design$age))
  if ("total_cells" %in% names(design)) {
    design$log_total_cells <- log1p(as.numeric(design$total_cells))
  }
  covariates <- c()
  for (column in c("sex", "chemistry", "collection_method")) {
    update <- add_covariate_if_usable(design, covariates, column)
    design <- update$design_frame
    covariates <- update$covariates
  }
  update <- add_covariate_if_usable(design, covariates, "log_total_cells", numeric = TRUE)
  design <- update$design_frame
  covariates <- update$covariates
  while (TRUE) {
    formula_text <- paste("~ age_scaled", paste(covariates, collapse = " + "), sep = ifelse(length(covariates), " + ", ""))
    design_matrix <- model.matrix(as.formula(formula_text), data = design)
    if (qr(design_matrix)$rank == ncol(design_matrix)) {
      rownames(design) <- design$donor_id
      return(list(design = design, formula = as.formula(formula_text), covariates = covariates))
    }
    if (length(covariates) == 0) {
      rownames(design) <- design$donor_id
      return(list(design = design, formula = as.formula("~ age_scaled"), covariates = covariates))
    }
    covariates <- covariates[-length(covariates)]
  }
}

message("Reading MiloR metadata: ", opt$metadata)
metadata <- read.delim(opt$metadata, check.names = FALSE)
message("Reading MiloR embedding: ", opt$embedding)
embedding_dt <- read.delim(opt$embedding, check.names = FALSE)
common <- intersect(metadata$cell_id, embedding_dt$cell_id)
metadata <- metadata[match(common, metadata$cell_id), , drop = FALSE]
embedding_dt <- embedding_dt[match(common, embedding_dt$cell_id), , drop = FALSE]
embedding_cols <- setdiff(names(embedding_dt), "cell_id")
embedding <- as.matrix(embedding_dt[, embedding_cols, drop = FALSE])
storage.mode(embedding) <- "double"
rownames(embedding) <- embedding_dt$cell_id

counts <- Matrix(0, nrow = 1, ncol = nrow(metadata), sparse = TRUE)
rownames(counts) <- "dummy"
colnames(counts) <- metadata$cell_id
sce <- SingleCellExperiment(assays = list(counts = counts), colData = as.data.frame(metadata))
reducedDim(sce, "X_scvi") <- embedding

design_info <- build_design(metadata)
message("Building Milo graph: k=", opt$k, "; prop=", opt$prop)
milo <- Milo(sce)
milo <- buildGraph(milo, k = opt$k, d = ncol(embedding), reduced.dim = "X_scvi")
milo <- makeNhoods(milo, prop = opt$prop, k = opt$k, d = ncol(embedding), refined = TRUE, reduced_dims = "X_scvi")
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample = "donor_id")
milo <- calcNhoodDistance(milo, d = ncol(embedding), reduced.dim = "X_scvi")
message("Testing Milo neighborhoods with design: ", deparse(design_info$formula))
da <- as.data.frame(testNhoods(milo, design = design_info$formula, design.df = design_info$design, reduced.dim = "X_scvi"))
da$milor_neighborhood_id <- seq_len(nrow(da))

nh <- nhoods(milo)
annotation_rows <- vector("list", ncol(nh))
for (idx in seq_len(ncol(nh))) {
  cell_idx <- which(nh[, idx] > 0)
  fine <- top_label(metadata$fine_celltype[cell_idx])
  coarse <- top_label(metadata$coarse_celltype[cell_idx])
  annotation_rows[[idx]] <- data.frame(
    milor_neighborhood_id = idx,
    n_cells = length(cell_idx),
    top_fine_celltype = fine$label,
    top_fine_fraction = fine$fraction,
    top_coarse_celltype = coarse$label,
    top_coarse_fraction = coarse$fraction,
    check.names = FALSE
  )
}
annotation <- do.call(rbind, annotation_rows)
da <- merge(da, annotation, by = "milor_neighborhood_id", all.x = TRUE)
da$run <- opt$run_name
da$n_cells_exported <- nrow(metadata)
da$n_donors <- length(unique(metadata$donor_id))
da$covariates <- paste(design_info$covariates, collapse = ",")

fdr_col <- if ("SpatialFDR" %in% names(da)) "SpatialFDR" else if ("FDR" %in% names(da)) "FDR" else NA_character_
summary <- data.frame(
  run = opt$run_name,
  metric = c("cells", "donors", "neighborhoods", "fdr_lt_0_10", "covariates"),
  value = c(nrow(metadata), length(unique(metadata$donor_id)), nrow(da), if (!is.na(fdr_col)) sum(da[[fdr_col]] < 0.10, na.rm = TRUE) else NA, paste(design_info$covariates, collapse = ",")),
  detail = c("cells used in official MiloR subset", "donors used in official MiloR subset", "MiloR neighborhoods tested", "MiloR neighborhoods at FDR < 0.10", "covariates retained after rank checks"),
  check.names = FALSE
)

ensure_parent(opt$out)
write.table(da, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
ensure_parent(opt$summary_out)
write.table(summary, opt$summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote MiloR parity DA: ", opt$out, " (", nrow(da), " rows)")
message("Wrote MiloR parity summary: ", opt$summary_out)
