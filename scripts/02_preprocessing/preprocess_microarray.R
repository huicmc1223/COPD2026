###############################################################################
# preprocess_microarray.R — Multi-platform microarray preprocessing
#
# Datasets & platforms:
#   Affymetrix HG-U133_Plus_2 (GPL570):
#     GSE24206, GSE21369, GSE5058, GSE8545, GSE20257, GSE11906
#   Affymetrix HuGene-1_0-st (GPL6244):
#     GSE110147, GSE32537
#   Agilent 4x44K / 8x60K (GPL6480 / GPL14550):
#     GSE47460, GSE38974
#   Illumina HT-12 V4 (GPL10558):
#     GSE76925
#
# Pipeline per platform:
#   Affy → oligo::rma()  |  Agilent → limma quantile  |  Illumina → lumi
# Then: probe→gene mapping → cross-platform ComBat → PCA QC
#
# Run: Rscript scripts/02_preprocessing/preprocess_microarray.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(limma)
  library(sva)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(AnnotationDbi)
})

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Affymetrix datasets (RMA normalisation)
# ══════════════════════════════════════════════════════════════════════════════

preprocess_affy <- function(geo_id, cel_dir, annotation_db = "hgu133plus2.db") {
  message("\n── Affymetrix: ", geo_id, " (", annotation_db, ") ──")
  
  if (!requireNamespace("oligo", quietly = TRUE)) stop("oligo not installed")
  
  # Read CEL files
  cel_files <- list.files(cel_dir, pattern = "\\.CEL(\\.gz)?$",
                          full.names = TRUE, ignore.case = TRUE)
  if (length(cel_files) == 0) {
    message("  ⚠ No CEL files found in ", cel_dir)
    return(NULL)
  }
  message("  Found ", length(cel_files), " CEL files")
  
  raw_data <- oligo::read.celfiles(cel_files)
  eset <- oligo::rma(raw_data)
  
  # Probe → Gene symbol mapping
  if (requireNamespace(annotation_db, quietly = TRUE)) {
    db <- get(annotation_db)
    probe_ids <- rownames(exprs(eset))
    gene_map <- AnnotationDbi::select(db, keys = probe_ids,
                                      columns = "SYMBOL", keytype = "PROBEID")
    gene_map <- gene_map[!is.na(gene_map$SYMBOL), ]
    gene_map <- gene_map[!duplicated(gene_map$PROBEID), ]
    
    # Collapse probes to genes (keep highest IQR)
    expr_mat <- exprs(eset)
    expr_mat <- expr_mat[gene_map$PROBEID, ]
    rownames(expr_mat) <- gene_map$SYMBOL
    
    # For duplicated genes, keep the probe with highest IQR
    if (any(duplicated(rownames(expr_mat)))) {
      iqrs <- apply(expr_mat, 1, IQR)
      expr_mat <- expr_mat[order(iqrs, decreasing = TRUE), ]
      expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
    }
    
    message("  Gene-level: ", nrow(expr_mat), " genes")
  } else {
    expr_mat <- exprs(eset)
    message("  ⚠ Annotation DB not available; using probe IDs")
  }
  
  # Save
  out_dir <- file.path(proc_array, geo_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(expr_mat, file.path(out_dir, "expr_rma_gene.csv"))
  saveRDS(eset, file.path(out_dir, "eset.rds"))
  message("  ✔ Saved: expr_rma_gene.csv")
  
  return(expr_mat)
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Agilent datasets (limma quantile normalisation)
# ══════════════════════════════════════════════════════════════════════════════

preprocess_agilent <- function(geo_id) {
  message("\n── Agilent: ", geo_id, " ──")
  
  # For Agilent, getGEO gives processed data in ExpressionSet
  eset_path <- file.path(raw_array, geo_id)
  rds_files <- list.files(eset_path, pattern = "\\.rds$", full.names = TRUE)
  
  # Try to load from getGEO result
  gse_files <- list.files(eset_path, pattern = "series_matrix",
                          full.names = TRUE)
  if (length(gse_files) > 0) {
    # Already have series matrix — use it
    gse <- GEOquery::getGEO(filename = gse_files[1])
    expr_mat <- exprs(gse)
  } else {
    # Re-download
    gse <- GEOquery::getGEO(geo_id, destdir = eset_path)
    if (is.list(gse)) gse <- gse[[1]]
    expr_mat <- exprs(gse)
  }
  
  # Quantile normalise (if not already)
  if (max(expr_mat, na.rm = TRUE) > 100) {
    # Likely raw; log2 + quantile
    expr_mat[expr_mat <= 0] <- NA
    expr_mat <- log2(expr_mat)
  }
  expr_mat <- limma::normalizeBetweenArrays(expr_mat, method = "quantile")
  
  # Map probes to gene symbols using fData
  fd <- fData(gse)
  symbol_col <- grep("symbol|gene.symbol", colnames(fd),
                     ignore.case = TRUE, value = TRUE)
  if (length(symbol_col) > 0) {
    gene_symbols <- fd[[symbol_col[1]]]
    valid <- !is.na(gene_symbols) & gene_symbols != ""
    expr_mat <- expr_mat[valid, ]
    rownames(expr_mat) <- gene_symbols[valid]
    
    # Remove duplicate genes (keep highest IQR)
    if (any(duplicated(rownames(expr_mat)))) {
      iqrs <- apply(expr_mat, 1, IQR, na.rm = TRUE)
      expr_mat <- expr_mat[order(iqrs, decreasing = TRUE), ]
      expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
    }
  }
  
  message("  Gene-level: ", nrow(expr_mat), " genes × ", ncol(expr_mat),
          " samples")
  
  # Save
  out_dir <- file.path(proc_array, geo_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(expr_mat, file.path(out_dir, "expr_quantile_gene.csv"))
  message("  ✔ Saved: expr_quantile_gene.csv")
  
  return(expr_mat)
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Illumina dataset (lumi / neqc)
# ══════════════════════════════════════════════════════════════════════════════

preprocess_illumina <- function(geo_id) {
  message("\n── Illumina: ", geo_id, " ──")
  
  eset_path <- file.path(raw_array, geo_id)
  gse <- GEOquery::getGEO(geo_id, destdir = eset_path, GSEMatrix = TRUE)
  if (is.list(gse)) gse <- gse[[1]]
  
  expr_mat <- exprs(gse)
  
  # Background correction + quantile normalisation
  if (max(expr_mat, na.rm = TRUE) > 100) {
    expr_mat[expr_mat <= 0] <- 1
    expr_mat <- log2(expr_mat)
  }
  expr_mat <- limma::normalizeBetweenArrays(expr_mat, method = "quantile")
  
  # Map to gene symbols
  fd <- fData(gse)
  symbol_col <- grep("symbol|gene.symbol", colnames(fd),
                     ignore.case = TRUE, value = TRUE)
  if (length(symbol_col) > 0) {
    gene_symbols <- fd[[symbol_col[1]]]
    valid <- !is.na(gene_symbols) & gene_symbols != ""
    expr_mat <- expr_mat[valid, ]
    rownames(expr_mat) <- gene_symbols[valid]
    
    if (any(duplicated(rownames(expr_mat)))) {
      iqrs <- apply(expr_mat, 1, IQR, na.rm = TRUE)
      expr_mat <- expr_mat[order(iqrs, decreasing = TRUE), ]
      expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
    }
  }
  
  message("  Gene-level: ", nrow(expr_mat), " genes × ", ncol(expr_mat),
          " samples")
  
  out_dir <- file.path(proc_array, geo_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(expr_mat, file.path(out_dir, "expr_quantile_gene.csv"))
  message("  ✔ Saved: expr_quantile_gene.csv")
  
  return(expr_mat)
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — Process all available datasets
# ══════════════════════════════════════════════════════════════════════════════

results_array <- list()

# Affymetrix (check if CEL files exist)
for (geo_id in c("GSE24206", "GSE21369", "GSE32537", "GSE110147")) {
  cel_dir <- file.path(raw_array, geo_id)
  annot <- ifelse(geo_id %in% c("GSE32537", "GSE110147"),
                  "hugene10stv1.db", "hgu133plus2.db")
  if (dir.exists(cel_dir)) {
    results_array[[geo_id]] <- preprocess_affy(geo_id, cel_dir, annot)
  } else {
    message("⚠ ", geo_id, " directory not found — skipping")
  }
}

# Agilent
for (geo_id in c("GSE47460", "GSE38974")) {
  results_array[[geo_id]] <- tryCatch(
    preprocess_agilent(geo_id),
    error = function(e) { message("⚠ ", geo_id, ": ", e$message); NULL }
  )
}

# Illumina
results_array[["GSE76925"]] <- tryCatch(
  preprocess_illumina("GSE76925"),
  error = function(e) { message("⚠ GSE76925: ", e$message); NULL }
)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Cross-platform integration (common genes + ComBat)
# ══════════════════════════════════════════════════════════════════════════════

available <- results_array[!sapply(results_array, is.null)]

if (length(available) >= 2) {
  message("\n── Cross-platform integration ──")
  
  # Find common genes across all platforms
  gene_lists <- lapply(available, rownames)
  common_genes <- Reduce(intersect, gene_lists)
  message("  Common genes across ", length(available), " datasets: ",
          length(common_genes))
  
  # Merge
  merged <- do.call(cbind, lapply(available, function(x) x[common_genes, ]))
  batch_labels <- unlist(lapply(names(available), function(nm) {
    rep(nm, ncol(available[[nm]]))
  }))
  
  # ComBat batch correction
  corrected <- combat_correct(merged, factor(batch_labels))
  
  # Save
  write.csv(corrected, file.path(proc_array, "cross_platform_combat.csv"))
  message("  ✔ Cross-platform batch-corrected matrix saved (",
          nrow(corrected), " genes × ", ncol(corrected), " samples)")
}

# ══════════════════════════════════════════════════════════════════════════════
message("\n✔ Microarray preprocessing complete")
save_session_info(report_dir)
