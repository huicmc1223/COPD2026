###############################################################################
# fig1_overview.R — Figure 1: Study overview & TMEM164 pan-cohort expression
#
# Panels:
#   A. Study design schematic (placeholder annotation)
#   B. TMEM164 expression across bulk cohorts (COPD + IPF + Control)
#   C. Forest plot — meta-analysis summary
#   D. Volcano plot — top DEG cohort
#
# Run: Rscript scripts/06_figures/fig1_overview.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# Panel B: TMEM164 expression across bulk cohorts
# ══════════════════════════════════════════════════════════════════════════════

collect_tmem164_expr <- function() {
  result <- list()
  
  for (geo_id in c("GSE57148", "GSE76925", "GSE47460",
                    "GSE134692", "GSE124685", "GSE32537")) {
    paths <- c(file.path(proc_bulk, geo_id, "vst_matrix.csv"),
               file.path(proc_array, geo_id, "expr_quantile_gene.csv"),
               file.path(proc_array, geo_id, "expr_rma_gene.csv"))
    expr_path <- paths[file.exists(paths)][1]
    meta_path <- file.path(raw_array, geo_id, "phenodata_1.csv")
    
    if (is.na(expr_path) || !file.exists(meta_path)) next
    
    expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
    meta <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)
    common <- intersect(colnames(expr), rownames(meta))
    
    if (!(TARGET_GENE %in% rownames(expr)) || length(common) < 5) next
    
    cond_col <- grep("condition|disease|group", colnames(meta),
                     ignore.case = TRUE, value = TRUE)[1]
    if (is.na(cond_col)) next
    
    meta <- meta[common, , drop = FALSE]
    result[[geo_id]] <- data.frame(
      geo_id = geo_id,
      sample = common,
      expr = as.numeric(expr[TARGET_GENE, common]),
      condition = meta[[cond_col]]
    )
  }
  
  if (length(result) > 0) rbindlist(result) else NULL
}

expr_df <- collect_tmem164_expr()

if (!is.null(expr_df) && nrow(expr_df) > 0) {
  # Standardise condition labels
  expr_df[, condition := fcase(
    grepl("COPD|copd", condition), "COPD",
    grepl("IPF|ipf|fibrosis", condition, ignore.case = TRUE), "IPF",
    default = "Control"
  )]
  
  p_B <- ggplot(expr_df, aes(x = geo_id, y = expr, fill = condition)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
    scale_fill_manual(values = c(Control = "#00A087", COPD = "#4DBBD5",
                                  IPF = "#F39B7F")) +
    labs(x = NULL, y = paste0(TARGET_GENE, " expression"),
         title = paste0(TARGET_GENE, " across cohorts")) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  save_nature_fig(p_B, "fig1_B_tmem164_cohorts", width = 5, height = 3.5)
}

# ══════════════════════════════════════════════════════════════════════════════
# Panel C: Forest plot (already generated — copy reference)
# ══════════════════════════════════════════════════════════════════════════════

message("Panel C: Forest plots → see copd_tmem164_forest.pdf / ipf_tmem164_forest.pdf")

# ══════════════════════════════════════════════════════════════════════════════
# Panel D: Volcano plot — largest COPD cohort
# ══════════════════════════════════════════════════════════════════════════════

# Re-use per-cohort DEG results if available
for (geo_id in c("GSE76925", "GSE57148")) {
  dds_path <- file.path(proc_bulk, geo_id, "dds.rds")
  if (!file.exists(dds_path)) {
    # Try limma results
    res_path <- file.path(table_dir, paste0(geo_id, "_limma_results.csv"))
    if (file.exists(res_path)) {
      res <- fread(res_path)
      if ("logFC" %in% colnames(res)) {
        colnames(res)[colnames(res) == "logFC"] <- "log2FoldChange"
        colnames(res)[colnames(res) == "adj.P.Val"] <- "padj"
        colnames(res)[colnames(res) == "P.Value"] <- "pvalue"
      }
    } else {
      next
    }
  } else {
    library(DESeq2)
    dds <- readRDS(dds_path)
    res <- as.data.frame(results(dds, contrast = c("condition", "COPD", "Control")))
    res$gene <- rownames(res)
  }
  
  if (!("log2FoldChange" %in% colnames(res))) next
  
  res <- as.data.table(res)
  res[, sig := fcase(
    padj < PADJ_CUTOFF & log2FoldChange > LFC_CUTOFF, "Up",
    padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF, "Down",
    default = "NS"
  )]
  
  # Label TMEM164
  res[, label := ifelse(gene == TARGET_GENE, TARGET_GENE, "")]
  
  p_D <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue),
                          colour = sig)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_manual(values = c(Up = "#E64B35", Down = "#4DBBD5",
                                   NS = "grey80")) +
    geom_text(data = res[label != ""], aes(label = label),
              colour = "black", size = 2, nudge_y = 0.5) +
    geom_point(data = res[label != ""], colour = "black", size = 1.5) +
    labs(x = "log2FC", y = "-log10(p-value)",
         title = paste0(geo_id, ": COPD vs Control")) +
    theme_nature() + NoLegend()
  
  save_nature_fig(p_D, "fig1_D_volcano", width = 3.5, height = 3.5)
  break
}

# ══════════════════════════════════════════════════════════════════════════════
# Composite (placeholder: panels will be assembled in Illustrator/Inkscape)
# ══════════════════════════════════════════════════════════════════════════════

if (exists("p_B") && exists("p_D")) {
  p_composite <- p_B / p_D + plot_annotation(tag_levels = list(c("B", "D")))
  save_nature_fig(p_composite, "fig1_composite", width = 7.2, height = 7)
}

message("\n✔ Figure 1 panels generated")
save_session_info(report_dir)
