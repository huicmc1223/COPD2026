###############################################################################
# ipf_marker_interaction.R — TMEM164 vs IPF core markers
#
# Step 19: SPP1, MMP7, MUC5B, KRT17, ACTA2 interaction analysis
# (Mirrors copd_marker_interaction.R with IPF-specific focus)
#
# Run: Rscript scripts/04_ipf_analysis/ipf_marker_interaction.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(data.table)
})

ipf_interaction_markers <- c(TARGET_GENE, "SPP1", "MMP7", "MUC5B",
                              "KRT17", "ACTA2")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Bulk-level Spearman correlations
# ══════════════════════════════════════════════════════════════════════════════

cor_results <- list()

for (geo_id in c("GSE134692", "GSE124685", "GSE47460", "GSE32537")) {
  # Try multiple processed paths
  paths <- c(file.path(proc_bulk, geo_id, "vst_matrix.csv"),
             file.path(proc_array, geo_id, "expr_quantile_gene.csv"),
             file.path(proc_array, geo_id, "expr_rma_gene.csv"))
  
  expr_path <- paths[file.exists(paths)][1]
  if (is.na(expr_path)) next
  
  expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
  genes_avail <- intersect(ipf_interaction_markers, rownames(expr))
  
  if (!(TARGET_GENE %in% genes_avail)) next
  
  sub_expr <- t(expr[genes_avail, ])
  
  for (g in setdiff(genes_avail, TARGET_GENE)) {
    ct <- cor.test(sub_expr[, TARGET_GENE], sub_expr[, g],
                   method = "spearman")
    cor_results[[paste0(geo_id, "_", g)]] <- data.frame(
      dataset = geo_id, gene = g,
      rho = ct$estimate, pvalue = ct$p.value
    )
  }
}

if (length(cor_results) > 0) {
  cor_df <- rbindlist(cor_results)
  cor_df$padj <- p.adjust(cor_df$pvalue, method = "BH")
  write.csv(cor_df,
            file.path(table_dir, "ipf_tmem164_marker_correlations.csv"),
            row.names = FALSE)
  
  # Heatmap
  cor_wide <- dcast(cor_df, gene ~ dataset, value.var = "rho")
  cor_mat <- as.matrix(cor_wide[, -1, drop = FALSE])
  rownames(cor_mat) <- cor_wide$gene
  
  pdf(file.path(fig_dir, "ipf_tmem164_marker_correlation_heatmap.pdf"),
      width = 5, height = 3)
  pheatmap(cor_mat,
           color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
           breaks = seq(-1, 1, length.out = 101),
           display_numbers = TRUE, number_format = "%.2f",
           fontsize = 7,
           main = paste0(TARGET_GENE, " vs IPF markers (Spearman ρ)"))
  dev.off()
  
  message("✔ IPF bulk correlations saved")
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. scRNA-seq co-expression
# ══════════════════════════════════════════════════════════════════════════════

seu_path <- file.path(proc_sc, "IPF_sc_annotated.rds")
if (file.exists(seu_path)) {
  library(Seurat)
  seu <- readRDS(seu_path)
  
  genes_sc <- intersect(ipf_interaction_markers, rownames(seu))
  
  if (TARGET_GENE %in% genes_sc && length(genes_sc) >= 2) {
    expr_data <- FetchData(seu, vars = c(genes_sc, "celltype_predicted",
                                          "condition"))
    
    # Scatter plots
    scatter_plots <- list()
    for (g in setdiff(genes_sc, TARGET_GENE)) {
      p <- ggplot(expr_data, aes_string(x = TARGET_GENE, y = g)) +
        geom_point(aes(colour = condition), size = 0.1, alpha = 0.2) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.5,
                    aes(colour = condition)) +
        scale_color_manual(values = col_condition) +
        labs(title = paste0(TARGET_GENE, " vs ", g)) +
        theme_nature(base_size = 6) +
        theme(legend.position = "none")
      scatter_plots[[g]] <- p
    }
    
    p_combined <- wrap_plots(scatter_plots, ncol = 3)
    save_nature_fig(p_combined, "ipf_tmem164_marker_scatter",
                    width = 7.2, height = 4)
    
    # Per-cell-type correlation matrix (for IPF cells only)
    ipf_cells <- expr_data[expr_data$condition == "IPF", ]
    ct_cor <- list()
    for (ct in unique(ipf_cells$celltype_predicted)) {
      sub <- ipf_cells[ipf_cells$celltype_predicted == ct, genes_sc]
      if (nrow(sub) > 30) {
        cor_mat <- cor(sub, method = "spearman", use = "complete.obs")
        ct_cor[[ct]] <- cor_mat[TARGET_GENE, ]
      }
    }
    
    if (length(ct_cor) > 0) {
      ct_cor_mat <- do.call(rbind, ct_cor)
      
      pdf(file.path(fig_dir, "ipf_tmem164_celltype_marker_cor.pdf"),
          width = 5, height = 4)
      pheatmap(ct_cor_mat,
               color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
               breaks = seq(-1, 1, length.out = 101),
               display_numbers = TRUE, number_format = "%.2f",
               fontsize = 6,
               main = paste0(TARGET_GENE,
                             " correlation with IPF markers per cell type"))
      dev.off()
    }
  }
  
  rm(seu); gc()
  message("✔ IPF scRNA-seq co-expression saved")
}

message("\n✔ IPF marker interaction analysis complete")
save_session_info(report_dir)
