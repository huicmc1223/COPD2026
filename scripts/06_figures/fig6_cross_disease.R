###############################################################################
# fig6_cross_disease.R — Figure 6: Cross-disease integration
#
# Panels:
#   A. Forest plot overlay — COPD + IPF meta-analysis
#   B. Cell-type barplot — shared/divergent regulation
#   C. WGCNA hub gene network (schematic)
#   D. CellChat differential interaction (if available)
#
# Run: Rscript scripts/06_figures/fig6_cross_disease.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# Panel A: Combined forest
# ══════════════════════════════════════════════════════════════════════════════

combined_path <- file.path(table_dir, "cross_disease_deg_combined.csv")
if (file.exists(combined_path)) {
  all_deg <- fread(combined_path)
  
  p_A <- ggplot(all_deg, aes(x = log2FC, y = reorder(geo_id, log2FC),
                               colour = disease)) +
    geom_errorbarh(aes(xmin = log2FC - 1.96 * se,
                       xmax = log2FC + 1.96 * se), height = 0.2) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c(COPD = "#4DBBD5", IPF = "#F39B7F")) +
    labs(x = paste0(TARGET_GENE, " log2FC (Disease vs Control)"),
         y = NULL, title = "Cross-disease meta-comparison") +
    theme_nature()
  
  save_nature_fig(p_A, "fig6_A_cross_forest", width = 5, height = 4)
}

# ══════════════════════════════════════════════════════════════════════════════
# Panel B: reference to cross_disease_celltype_barplot.pdf
# ══════════════════════════════════════════════════════════════════════════════

bar_path <- file.path(fig_dir, "cross_disease_celltype_barplot.pdf")
if (file.exists(bar_path)) {
  message("Panel B: → cross_disease_celltype_barplot.pdf")
} else {
  message("⚠ Generate cross_disease_comparison.R first")
}

# ══════════════════════════════════════════════════════════════════════════════
# Panel C: WGCNA hub gene network — top connections
# ══════════════════════════════════════════════════════════════════════════════

hub_files <- list.files(table_dir, pattern = "wgcna_.*module.*\\.csv",
                        full.names = TRUE)

if (length(hub_files) > 0) {
  hub_df <- fread(hub_files[1])
  top30 <- head(hub_df, 30)
  
  # Lollipop chart of kME values
  top30$is_tmem <- top30$gene == TARGET_GENE
  
  p_C <- ggplot(top30, aes(x = kME, y = reorder(gene, kME),
                             colour = is_tmem)) +
    geom_point(size = 2) +
    geom_segment(aes(xend = 0, yend = gene), linewidth = 0.4) +
    scale_color_manual(values = c("FALSE" = "#3C5488", "TRUE" = "#E64B35"),
                       guide = "none") +
    labs(x = "Module membership (kME)", y = NULL,
         title = paste0(TARGET_GENE, " module hub genes")) +
    theme_nature()
  
  save_nature_fig(p_C, "fig6_C_wgcna_hub", width = 3.5, height = 5)
}

# ══════════════════════════════════════════════════════════════════════════════
# Panel D: CellChat differential (reference)
# ══════════════════════════════════════════════════════════════════════════════

cc_path <- file.path(fig_dir, "cellchat_diff_network.pdf")
if (file.exists(cc_path)) {
  message("Panel D: → cellchat_diff_network.pdf")
}

# Composite
plots_available <- list()
if (exists("p_A")) plots_available$A <- p_A
if (exists("p_C")) plots_available$C <- p_C

if (length(plots_available) >= 2) {
  p_fig6 <- (plots_available$A | plots_available$C) +
    plot_annotation(tag_levels = "A")
  save_nature_fig(p_fig6, "fig6_composite", width = 7.2, height = 5)
}

message("\n✔ Figure 6 panels generated")
save_session_info(report_dir)
