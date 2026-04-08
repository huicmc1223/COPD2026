###############################################################################
# preprocess_bulk_rnaseq.R — Bulk RNA-seq preprocessing pipeline
#
# Datasets processed:
#   GSE57148  — COPD (Korean, GOLD 1-2), raw counts
#   GSE134692 — IPF (Caucasian), raw counts
#   GSE124685 — IPF (Caucasian), raw counts
#   GSE213001 — IPF (Australian), raw counts [Priority 3]
#
# Pipeline: raw counts → filter low-expression → DESeq2 → VST → PCA QC
# Cross-cohort: ComBat-seq batch correction (within disease)
#
# Run: Rscript scripts/02_preprocessing/preprocess_bulk_rnaseq.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
  library(ggplot2)
  library(patchwork)
})

# ══════════════════════════════════════════════════════════════════════════════
# Helper: DESeq2 preprocessing for a single dataset
# ══════════════════════════════════════════════════════════════════════════════

preprocess_one_bulk <- function(geo_id, counts_file, meta_file,
                                condition_col = "condition",
                                min_count = 10, min_samples = 3) {
  message("\n── Processing ", geo_id, " ──")
  
  # Read counts
  counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  meta   <- read.csv(meta_file, row.names = 1, stringsAsFactors = FALSE)
  
  # Align samples
  common <- intersect(colnames(counts), rownames(meta))
  counts <- counts[, common]
  meta   <- meta[common, , drop = FALSE]
  
  message("  Samples: ", ncol(counts), "  Genes: ", nrow(counts))
  
  # Filter low-expression genes
  keep <- rowSums(counts >= min_count) >= min_samples
  counts <- counts[keep, ]
  message("  After filtering: ", nrow(counts), " genes retained")
  
  # DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData   = meta,
    design    = as.formula(paste0("~ ", condition_col))
  )
  dds <- DESeq(dds)
  
  # VST transformation (for visualisation & batch correction)
  vst_data <- vst(dds, blind = FALSE)
  
  # PCA QC plot
  pca_res <- plotPCA(vst_data, intgroup = condition_col, returnData = TRUE)
  pct_var <- round(100 * attr(pca_res, "percentVar"), 1)
  
  p_pca <- ggplot(pca_res, aes(PC1, PC2, colour = .data[[condition_col]])) +
    geom_point(size = 2.5, alpha = 0.8) +
    xlab(paste0("PC1 (", pct_var[1], "%)")) +
    ylab(paste0("PC2 (", pct_var[2], "%)")) +
    labs(title = paste0(geo_id, " — PCA (VST)")) +
    scale_color_manual(values = col_condition) +
    theme_nature()
  
  # Save outputs
  out_dir <- file.path(proc_bulk, geo_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  saveRDS(dds, file.path(out_dir, "dds.rds"))
  saveRDS(vst_data, file.path(out_dir, "vst.rds"))
  write.csv(assay(vst_data), file.path(out_dir, "vst_matrix.csv"))
  write.csv(counts(dds, normalized = TRUE),
            file.path(out_dir, "normalized_counts.csv"))
  save_nature_fig(p_pca, paste0(geo_id, "_pca_qc"), width = 4, height = 3.5,
                  dir = file.path(report_dir))
  
  message("  ✔ Saved: dds.rds, vst.rds, vst_matrix.csv, PCA plot")
  
  return(list(dds = dds, vst = vst_data, pca_plot = p_pca))
}

# ══════════════════════════════════════════════════════════════════════════════
# Process each dataset
# NOTE: Adjust file paths after actual download. The paths below assume
# GEO supplementary files have been extracted to the expected directories.
# ══════════════════════════════════════════════════════════════════════════════

# --- GSE57148 (COPD) ---
# After download, expect: data/raw/bulk_rnaseq/GSE57148/counts.csv + meta.csv
# The actual file names depend on GEO supplementary format — adapt as needed.
gse57148_counts <- file.path(raw_bulk, "GSE57148", "counts.csv")
gse57148_meta   <- file.path(raw_bulk, "GSE57148", "phenodata_1.csv")

if (file.exists(gse57148_counts)) {
  res_57148 <- preprocess_one_bulk("GSE57148", gse57148_counts, gse57148_meta)
} else {
  message("⚠ GSE57148 counts not found — run download_priority1.R first")
}

# --- GSE134692 (IPF) ---
gse134692_counts <- file.path(raw_bulk, "GSE134692", "counts.csv")
gse134692_meta   <- file.path(raw_bulk, "GSE134692", "phenodata_1.csv")

if (file.exists(gse134692_counts)) {
  res_134692 <- preprocess_one_bulk("GSE134692", gse134692_counts, gse134692_meta)
} else {
  message("⚠ GSE134692 counts not found — run download_priority1.R first")
}

# --- GSE124685 (IPF) ---
gse124685_counts <- file.path(raw_bulk, "GSE124685", "counts.csv")
gse124685_meta   <- file.path(raw_bulk, "GSE124685", "phenodata_1.csv")

if (file.exists(gse124685_counts)) {
  res_124685 <- preprocess_one_bulk("GSE124685", gse124685_counts, gse124685_meta)
} else {
  message("⚠ GSE124685 counts not found — run download_priority1.R first")
}

# ══════════════════════════════════════════════════════════════════════════════
# Cross-cohort batch correction (IPF bulk: GSE134692 + GSE124685)
# ══════════════════════════════════════════════════════════════════════════════

if (exists("res_134692") && exists("res_124685")) {
  message("\n── Cross-cohort batch correction: IPF bulk RNA-seq ──")
  
  # Merge VST matrices on common genes
  vst1 <- assay(res_134692$vst)
  vst2 <- assay(res_124685$vst)
  common_genes <- intersect(rownames(vst1), rownames(vst2))
  
  merged_vst <- cbind(vst1[common_genes, ], vst2[common_genes, ])
  batch <- factor(c(rep("GSE134692", ncol(vst1)),
                     rep("GSE124685", ncol(vst2))))
  
  # Condition labels
  cond1 <- colData(res_134692$dds)$condition
  cond2 <- colData(res_124685$dds)$condition
  condition <- factor(c(cond1, cond2))
  
  # ComBat
  corrected <- combat_correct(merged_vst, batch,
                              mod = model.matrix(~ condition))
  
  # PCA diagnostic
  p_batch <- plot_batch_pca(merged_vst, corrected, batch, condition,
                            "IPF Bulk RNA-seq")
  save_nature_fig(p_batch, "ipf_bulk_batch_correction_pca",
                  width = 7, height = 6, dir = report_dir)
  
  # Save
  write.csv(corrected, file.path(proc_bulk, "ipf_bulk_combat_vst.csv"))
  message("  ✔ IPF bulk batch-corrected VST saved")
}

# ══════════════════════════════════════════════════════════════════════════════
message("\n✔ Bulk RNA-seq preprocessing complete")
save_session_info(report_dir)
