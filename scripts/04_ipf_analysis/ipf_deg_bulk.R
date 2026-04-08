###############################################################################
# ipf_deg_bulk.R — TMEM164 differential expression across IPF bulk/array
#
# Step 15: Per-cohort DEG → Meta-analysis → Forest plot (IPF-specific)
#
# Datasets: GSE134692, GSE124685, GSE47460, GSE32537, GSE24206, GSE110147
#
# Run: Rscript scripts/04_ipf_analysis/ipf_deg_bulk.R
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
# Collect TMEM164 stats from each IPF cohort
# ══════════════════════════════════════════════════════════════════════════════

extract_tmem164 <- function(path, geo_id, type = "array") {
  if (!file.exists(path)) {
    message("  ⚠ ", geo_id, " not found")
    return(NULL)
  }
  
  if (type == "deseq2") {
    dds <- readRDS(path)
    res <- results(dds, contrast = c("condition", "IPF", "Control"))
    if (TARGET_GENE %in% rownames(res)) {
      row <- res[TARGET_GENE, ]
      return(data.frame(geo_id = geo_id, type = "bulk_rnaseq",
                        log2FC = row$log2FoldChange, se = row$lfcSE,
                        pvalue = row$pvalue, padj = row$padj))
    }
  } else {
    expr <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
    meta_path <- sub("expr_.*\\.csv", "../../raw/microarray/", path)
    # For simplicity, use the expression matrix directly with limma
    if (TARGET_GENE %in% rownames(expr)) {
      # Simple t-test approach (meta column extraction depends on curation)
      vals <- expr[TARGET_GENE, ]
      # Placeholder: will be refined per-dataset after metadata curation
      return(data.frame(geo_id = geo_id, type = "microarray",
                        log2FC = NA, se = NA,
                        pvalue = NA, padj = NA))
    }
  }
  NULL
}

deg_results <- list()

# Bulk RNA-seq
for (geo_id in c("GSE134692", "GSE124685")) {
  dds_path <- file.path(proc_bulk, geo_id, "dds.rds")
  if (file.exists(dds_path)) {
    dds <- readRDS(dds_path)
    res <- results(dds, contrast = c("condition", "IPF", "Control"))
    if (TARGET_GENE %in% rownames(res)) {
      row <- res[TARGET_GENE, ]
      deg_results[[geo_id]] <- data.frame(
        geo_id = geo_id, type = "bulk_rnaseq",
        log2FC = row$log2FoldChange, se = row$lfcSE,
        pvalue = row$pvalue, padj = row$padj
      )
    }
  }
}

# Microarrays — limma analysis
for (geo_id in c("GSE47460", "GSE32537", "GSE24206", "GSE110147")) {
  expr_path <- file.path(proc_array, geo_id)
  expr_files <- list.files(expr_path, pattern = "expr_.*gene\\.csv",
                           full.names = TRUE)
  meta_path <- file.path(raw_array, geo_id, "phenodata_1.csv")
  
  if (length(expr_files) > 0 && file.exists(meta_path)) {
    expr <- as.matrix(read.csv(expr_files[1], row.names = 1,
                               check.names = FALSE))
    meta <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)
    
    common <- intersect(colnames(expr), rownames(meta))
    if (length(common) > 0 && TARGET_GENE %in% rownames(expr)) {
      expr <- expr[, common]
      meta <- meta[common, , drop = FALSE]
      
      # Detect condition column
      cond_col <- grep("condition|disease|group", colnames(meta),
                       ignore.case = TRUE, value = TRUE)[1]
      if (!is.na(cond_col)) {
        design <- model.matrix(~ factor(meta[[cond_col]]))
        fit <- lmFit(expr, design)
        fit <- eBayes(fit)
        tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
        
        if (TARGET_GENE %in% rownames(tt)) {
          row <- tt[TARGET_GENE, ]
          deg_results[[geo_id]] <- data.frame(
            geo_id = geo_id, type = "microarray",
            log2FC = row$logFC,
            se = row$logFC / row$t,  # SE ≈ logFC / t-statistic
            pvalue = row$P.Value, padj = row$adj.P.Val
          )
        }
      }
    }
  }
}

# Combine
deg_df <- rbindlist(deg_results, fill = TRUE)
deg_df <- deg_df[!is.na(log2FC)]

if (nrow(deg_df) == 0) {
  message("⚠ No IPF DEG results — ensure preprocessing is complete")
  quit(save = "no")
}

write.csv(deg_df, file.path(table_dir, "ipf_tmem164_deg_percohort.csv"),
          row.names = FALSE)
message("IPF per-cohort results: ", nrow(deg_df), " datasets")

# ══════════════════════════════════════════════════════════════════════════════
# Meta-analysis
# ══════════════════════════════════════════════════════════════════════════════

if (nrow(deg_df) >= 2) {
  ma <- rma(yi = deg_df$log2FC, sei = deg_df$se,
            slab = deg_df$geo_id, method = "REML")
  
  ma_summary <- data.frame(
    pooled_log2FC = coef(ma), se = ma$se,
    ci_lower = ma$ci.lb, ci_upper = ma$ci.ub,
    pvalue = ma$pval, I2 = ma$I2, tau2 = ma$tau2
  )
  write.csv(ma_summary, file.path(table_dir, "ipf_tmem164_meta_summary.csv"),
            row.names = FALSE)
  
  message("  Pooled log2FC: ", round(ma_summary$pooled_log2FC, 3))
  
  # Forest plot
  pdf(file.path(fig_dir, "ipf_tmem164_forest.pdf"),
      width = 5, height = 3 + nrow(deg_df) * 0.3)
  forest(ma, xlab = paste0(TARGET_GENE, " log2FC (IPF vs Control)"),
         header = TRUE, cex = 0.75, col = "#F39B7F", border = "#F39B7F")
  dev.off()
  
  png(file.path(fig_dir, "ipf_tmem164_forest.png"),
      width = 5, height = 3 + nrow(deg_df) * 0.3, units = "in", res = 300)
  forest(ma, xlab = paste0(TARGET_GENE, " log2FC (IPF vs Control)"),
         header = TRUE, cex = 0.75, col = "#F39B7F", border = "#F39B7F")
  dev.off()
}

message("\n✔ IPF bulk/array DEG analysis complete")
save_session_info(report_dir)
