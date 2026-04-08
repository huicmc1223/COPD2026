###############################################################################
# cross_disease_comparison.R — COPD vs IPF cross-disease TMEM164 analysis
#
# Step 20: Shared & divergent patterns between COPD and IPF
#
# Run: Rscript scripts/05_integrative/cross_disease_comparison.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# 1. Meta-analysis comparison: forest plot overlay
# ══════════════════════════════════════════════════════════════════════════════

copd_deg <- fread(file.path(table_dir, "copd_tmem164_deg_percohort.csv"))
ipf_deg  <- fread(file.path(table_dir, "ipf_tmem164_deg_percohort.csv"))

if (nrow(copd_deg) > 0 && nrow(ipf_deg) > 0) {
  copd_deg$disease <- "COPD"
  ipf_deg$disease  <- "IPF"
  all_deg <- rbind(copd_deg, ipf_deg, fill = TRUE)
  all_deg <- all_deg[!is.na(log2FC) & !is.na(se)]
  
  # Combined forest-style dot plot
  p_forest <- ggplot(all_deg, aes(x = log2FC, y = reorder(geo_id, log2FC),
                                   colour = disease)) +
    geom_errorbarh(aes(xmin = log2FC - 1.96 * se, xmax = log2FC + 1.96 * se),
                   height = 0.2) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c(COPD = "#4DBBD5", IPF = "#F39B7F")) +
    labs(x = paste0(TARGET_GENE, " log2FC (Disease vs Control)"),
         y = "Dataset", title = paste0(TARGET_GENE, " across diseases")) +
    theme_nature()
  
  save_nature_fig(p_forest, "cross_disease_deg_comparison",
                  width = 5, height = 4)
  
  write.csv(all_deg,
            file.path(table_dir, "cross_disease_deg_combined.csv"),
            row.names = FALSE)
  message("✔ Cross-disease DEG comparison saved")
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. Cell-type–level comparison (scRNA-seq)
# ══════════════════════════════════════════════════════════════════════════════

copd_sc <- tryCatch(fread(file.path(table_dir, "copd_tmem164_sc_deg.csv")),
                    error = function(e) NULL)
ipf_sc  <- tryCatch(fread(file.path(table_dir, "ipf_tmem164_sc_deg.csv")),
                    error = function(e) NULL)

if (!is.null(copd_sc) && !is.null(ipf_sc)) {
  copd_sc$disease <- "COPD"
  ipf_sc$disease  <- "IPF"
  all_sc <- rbind(copd_sc, ipf_sc, fill = TRUE)
  
  # Heatmap: cell types × diseases, value = avg_log2FC
  mat <- dcast(all_sc, celltype ~ disease, value.var = "avg_log2FC")
  mat_m <- as.matrix(mat[, -1, drop = FALSE])
  rownames(mat_m) <- mat$celltype
  
  pdf(file.path(fig_dir, "cross_disease_celltype_heatmap.pdf"),
      width = 3.5, height = 5)
  pheatmap(mat_m,
           color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
           breaks = seq(-2, 2, length.out = 101),
           display_numbers = TRUE, number_format = "%.2f",
           fontsize = 7, cluster_cols = FALSE,
           main = paste0(TARGET_GENE, " log2FC by cell type"))
  dev.off()
  
  # Divergence plot
  common_ct <- mat_m[rowSums(!is.na(mat_m)) == 2, , drop = FALSE]
  if (nrow(common_ct) > 0) {
    df_plot <- data.frame(
      celltype = rownames(common_ct),
      COPD = common_ct[, "COPD"],
      IPF  = common_ct[, "IPF"]
    )
    df_long <- melt(as.data.table(df_plot), id.vars = "celltype",
                    variable.name = "disease", value.name = "log2FC")
    
    p_bar <- ggplot(df_long, aes(x = reorder(celltype, -log2FC),
                                  y = log2FC, fill = disease)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c(COPD = "#4DBBD5", IPF = "#F39B7F")) +
      geom_hline(yintercept = 0) +
      labs(x = "Cell type", y = paste0(TARGET_GENE, " log2FC"),
           title = "Shared vs divergent cell-type regulation") +
      theme_nature() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    save_nature_fig(p_bar, "cross_disease_celltype_barplot",
                    width = 5, height = 3.5)
  }
  
  message("✔ Cross-disease cell-type comparison saved")
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Pathway overlap: Venn/UpSet
# ══════════════════════════════════════════════════════════════════════════════

# Compare GSEA leading-edge genes from COPD vs IPF
copd_gsea <- tryCatch(fread(file.path(table_dir, "copd_gsea_kegg.csv")),
                      error = function(e) NULL)
ipf_gsea <- tryCatch(fread(file.path(table_dir, "ipf_tmem164_pathway_correlations.csv")),
                     error = function(e) NULL)

if (!is.null(copd_gsea) && !is.null(ipf_gsea)) {
  # Enriched pathways in COPD
  copd_pathways <- unique(copd_gsea[padj < 0.05, Description])
  
  # Will be expanded when full GSEA for IPF is run
  message("COPD enriched pathways: ", length(copd_pathways))
}

message("\n✔ Cross-disease integrative analysis complete")
save_session_info(report_dir)
