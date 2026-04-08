###############################################################################
# copd_deg_bulk.R — TMEM164 differential expression across COPD bulk/array
#
# Step 10: Per-cohort DEG → Meta-analysis (random effects) → Forest plot
#
# Input:  processed bulk/array matrices from 02_preprocessing/
# Output: results/tables/copd_tmem164_deg_meta.csv
#         results/figures/copd_tmem164_forest.pdf
#
# Run: Rscript scripts/03_copd_analysis/copd_deg_bulk.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(metafor)
  library(ggplot2)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# Helper: extract TMEM164 stats from a DESeq2 result
# ══════════════════════════════════════════════════════════════════════════════

extract_tmem164_deseq2 <- function(dds_path, geo_id) {
  dds <- readRDS(dds_path)
  res <- results(dds, contrast = c("condition", "COPD", "Control"))
  
  if (TARGET_GENE %in% rownames(res)) {
    row <- res[TARGET_GENE, ]
    data.frame(
      geo_id   = geo_id,
      type     = "bulk_rnaseq",
      log2FC   = row$log2FoldChange,
      se       = row$lfcSE,
      pvalue   = row$pvalue,
      padj     = row$padj,
      baseMean = row$baseMean,
      stringsAsFactors = FALSE
    )
  } else {
    message("  ⚠ ", TARGET_GENE, " not found in ", geo_id)
    NULL
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# Helper: extract TMEM164 stats from limma array result
# ══════════════════════════════════════════════════════════════════════════════

extract_tmem164_limma <- function(expr_path, meta_path, geo_id,
                                   condition_col = "condition") {
  expr <- as.matrix(read.csv(expr_path, row.names = 1, check.names = FALSE))
  meta <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)
  
  common <- intersect(colnames(expr), rownames(meta))
  expr <- expr[, common]
  meta <- meta[common, , drop = FALSE]
  
  if (!(TARGET_GENE %in% rownames(expr))) {
    message("  ⚠ ", TARGET_GENE, " not found in ", geo_id)
    return(NULL)
  }
  
  # limma
  design <- model.matrix(~ factor(meta[[condition_col]]))
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  
  if (TARGET_GENE %in% rownames(tt)) {
    row <- tt[TARGET_GENE, ]
    data.frame(
      geo_id = geo_id,
      type   = "microarray",
      log2FC = row$logFC,
      se     = sqrt(row$s2.post) / sqrt(ncol(expr)),  # approximate SE
      pvalue = row$P.Value,
      padj   = row$adj.P.Val,
      baseMean = mean(expr[TARGET_GENE, ]),
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# Collect results from all COPD cohorts
# ══════════════════════════════════════════════════════════════════════════════

deg_results <- list()

# --- Bulk RNA-seq: GSE57148 ---
dds_path <- file.path(proc_bulk, "GSE57148", "dds.rds")
if (file.exists(dds_path)) {
  deg_results[["GSE57148"]] <- extract_tmem164_deseq2(dds_path, "GSE57148")
}

# --- Microarray: GSE76925 ---
expr_path <- file.path(proc_array, "GSE76925", "expr_quantile_gene.csv")
meta_path <- file.path(raw_array, "GSE76925", "phenodata_1.csv")
if (file.exists(expr_path) && file.exists(meta_path)) {
  deg_results[["GSE76925"]] <- extract_tmem164_limma(expr_path, meta_path,
                                                      "GSE76925")
}

# --- Microarray: GSE47460 (COPD subset) ---
expr_path <- file.path(proc_array, "GSE47460", "expr_quantile_gene.csv")
meta_path <- file.path(raw_array, "GSE47460", "phenodata_1.csv")
if (file.exists(expr_path) && file.exists(meta_path)) {
  deg_results[["GSE47460"]] <- extract_tmem164_limma(expr_path, meta_path,
                                                      "GSE47460")
}

# --- Microarray: GSE38974 (Priority 2) ---
expr_path <- file.path(proc_array, "GSE38974", "expr_quantile_gene.csv")
meta_path <- file.path(raw_array, "GSE38974", "phenodata_1.csv")
if (file.exists(expr_path) && file.exists(meta_path)) {
  deg_results[["GSE38974"]] <- extract_tmem164_limma(expr_path, meta_path,
                                                      "GSE38974")
}

# Combine
deg_df <- do.call(rbind, deg_results[!sapply(deg_results, is.null)])

if (is.null(deg_df) || nrow(deg_df) == 0) {
  message("⚠ No DEG results available — ensure preprocessing is complete")
  quit(save = "no")
}

# Save
write.csv(deg_df, file.path(table_dir, "copd_tmem164_deg_percohort.csv"),
          row.names = FALSE)
message("Per-cohort results: ", nrow(deg_df), " datasets")
print(deg_df)

# ══════════════════════════════════════════════════════════════════════════════
# Meta-analysis (random effects model)
# ══════════════════════════════════════════════════════════════════════════════

if (nrow(deg_df) >= 2) {
  message("\n── Meta-analysis ──")
  
  ma <- rma(yi = deg_df$log2FC,
            sei = deg_df$se,
            slab = deg_df$geo_id,
            method = "REML")
  
  # Summary
  ma_summary <- data.frame(
    pooled_log2FC = coef(ma),
    se            = ma$se,
    ci_lower      = ma$ci.lb,
    ci_upper      = ma$ci.ub,
    pvalue        = ma$pval,
    I2            = ma$I2,
    tau2          = ma$tau2
  )
  write.csv(ma_summary, file.path(table_dir, "copd_tmem164_meta_summary.csv"),
            row.names = FALSE)
  
  message("  Pooled log2FC: ", round(ma_summary$pooled_log2FC, 3),
          " [", round(ma_summary$ci_lower, 3), ", ",
          round(ma_summary$ci_upper, 3), "]")
  message("  P-value: ", signif(ma_summary$pvalue, 3),
          "  I²: ", round(ma_summary$I2, 1), "%")
  
  # ── Forest plot ──
  pdf(file.path(fig_dir, "copd_tmem164_forest.pdf"),
      width = 5, height = 3 + nrow(deg_df) * 0.3)
  forest(ma, xlab = paste0(TARGET_GENE, " log2 Fold Change (COPD vs Control)"),
         header = TRUE, cex = 0.75, col = "#E64B35", border = "#E64B35")
  dev.off()
  
  png(file.path(fig_dir, "copd_tmem164_forest.png"),
      width = 5, height = 3 + nrow(deg_df) * 0.3, units = "in", res = 300)
  forest(ma, xlab = paste0(TARGET_GENE, " log2 Fold Change (COPD vs Control)"),
         header = TRUE, cex = 0.75, col = "#E64B35", border = "#E64B35")
  dev.off()
  
  message("  ✔ Forest plot saved")
}

message("\n✔ COPD bulk/array DEG analysis complete")
save_session_info(report_dir)
