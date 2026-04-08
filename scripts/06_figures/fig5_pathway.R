###############################################################################
# fig5_pathway.R — Figure 5: TMEM164-associated pathway landscape
#
# Panels:
#   A. GSEA enrichment plot — top KEGG pathways (COPD)
#   B. AUCell heatmap — cell-type × pathway
#   C. Ferroptosis/Autophagy scatter (IPF)
#   D. Epithelial fate radar/heatmap: TMEM164-high vs low
#
# Run: Rscript scripts/06_figures/fig5_pathway.R
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
# Panel A: GSEA dot plot
# ══════════════════════════════════════════════════════════════════════════════

gsea_path <- file.path(table_dir, "copd_gsea_kegg.csv")
if (file.exists(gsea_path)) {
  gsea_df <- fread(gsea_path)
  gsea_sig <- gsea_df[p.adjust < 0.05]
  
  if (nrow(gsea_sig) > 0) {
    top_paths <- head(gsea_sig[order(p.adjust)], 20)
    
    p_A <- ggplot(top_paths, aes(x = NES, y = reorder(Description, NES),
                                   size = setSize, colour = -log10(p.adjust))) +
      geom_point() +
      scale_color_viridis_c(option = "plasma") +
      labs(x = "NES", y = NULL, size = "Gene set size",
           colour = "-log10(padj)",
           title = "COPD: TMEM164-correlated KEGG pathways") +
      theme_nature()
    
    save_nature_fig(p_A, "fig5_A_gsea_dotplot", width = 5, height = 5)
  }
} else {
  message("⚠ GSEA results not found — run copd_pathway.R first")
}

# ══════════════════════════════════════════════════════════════════════════════
# Panel B: AUCell pathway × cell-type heatmap
# ══════════════════════════════════════════════════════════════════════════════

aucell_path <- file.path(table_dir, "copd_aucell_celltype_summary.csv")
if (file.exists(aucell_path)) {
  auc_df <- fread(aucell_path)
  
  # Pivot: cell type × pathway, value = mean AUC
  auc_wide <- dcast(auc_df, celltype ~ pathway, value.var = "mean_auc",
                    fun.aggregate = mean)
  auc_mat <- as.matrix(auc_wide[, -1, drop = FALSE])
  rownames(auc_mat) <- auc_wide$celltype
  
  pdf(file.path(fig_dir, "fig5_B_aucell_heatmap.pdf"), width = 5, height = 4)
  pheatmap(auc_mat, scale = "column",
           color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
           fontsize = 6,
           main = paste0(TARGET_GENE, ": pathway activity per cell type"))
  dev.off()
} else {
  message("⚠ AUCell summary not found — run copd_pathway.R first")
}

# ══════════════════════════════════════════════════════════════════════════════
# Panel C: Ferroptosis gene scatter (IPF)
# ══════════════════════════════════════════════════════════════════════════════

ipf_cor_path <- file.path(table_dir, "ipf_tmem164_pathway_correlations.csv")
if (file.exists(ipf_cor_path)) {
  ipf_cor <- fread(ipf_cor_path)
  ferro <- ipf_cor[pathway == "Ferroptosis"]
  
  if (nrow(ferro) > 0) {
    ferro$sig <- ifelse(ferro$padj < 0.05, "Significant", "NS")
    
    p_C <- ggplot(ferro, aes(x = rho, y = reorder(gene, rho),
                               colour = sig)) +
      geom_point(size = 2) +
      geom_segment(aes(xend = 0, yend = gene), linewidth = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_manual(values = c(Significant = "#E64B35", NS = "grey60")) +
      labs(x = paste0("Spearman ρ with ", TARGET_GENE),
           y = NULL, title = "Ferroptosis genes (IPF)") +
      theme_nature()
    
    save_nature_fig(p_C, "fig5_C_ferroptosis_lollipop", width = 3.5, height = 4)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# Panel D: Epithelial fate heatmap
# ══════════════════════════════════════════════════════════════════════════════

for (disease in c("copd", "ipf")) {
  fate_path <- file.path(table_dir,
                          paste0("epithelial_fate_", disease, "_sig_axes.csv"))
  if (file.exists(fate_path)) {
    sig_axes <- fread(fate_path)
    message(toupper(disease), " significant epithelial fate axes: ",
            nrow(sig_axes))
  }
}

message("\n✔ Figure 5 panels generated")
save_session_info(report_dir)
