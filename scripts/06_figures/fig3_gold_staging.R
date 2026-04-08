###############################################################################
# fig3_gold_staging.R — Figure 3: TMEM164 across GOLD 1-4 stages
#
# Panels:
#   A. Box/violin: TMEM164 per GOLD stage (bulk)
#   B. Trend test p-value annotation
#   C. ROC: COPD early detection
#   D. Comparison: single gene vs multi-gene panel
#
# Run: Rscript scripts/06_figures/fig3_gold_staging.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(pROC)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# Panel A-B: GOLD staging box/violin (re-use data from copd_gold_staging.R)
# ══════════════════════════════════════════════════════════════════════════════

gold_data_path <- file.path(table_dir, "copd_gold_staging_data.csv")
if (file.exists(gold_data_path)) {
  gold_df <- fread(gold_data_path)
  
  p_A <- ggplot(gold_df, aes(x = gold_stage, y = TMEM164_expr,
                               fill = gold_stage)) +
    geom_violin(linewidth = 0.3) +
    geom_boxplot(width = 0.15, linewidth = 0.3, outlier.size = 0.5) +
    scale_fill_manual(values = col_gold) +
    labs(x = "GOLD stage", y = paste0(TARGET_GENE, " expression"),
         title = paste0(TARGET_GENE, " across COPD severity")) +
    theme_nature() + NoLegend()
  
  # Per-dataset facet
  if ("dataset" %in% colnames(gold_df)) {
    p_A_facet <- p_A + facet_wrap(~ dataset, scales = "free_y")
    save_nature_fig(p_A_facet, "fig3_A_gold_facet", width = 7.2, height = 4)
  }
  
  save_nature_fig(p_A, "fig3_A_gold_violin", width = 3.5, height = 3.5)
} else {
  message("⚠ GOLD staging data not found — run copd_gold_staging.R first")
  p_A <- ggplot() + theme_void()
}

# ══════════════════════════════════════════════════════════════════════════════
# Panel C-D: ROC curves (re-use from diagnostic_model.R)
# ══════════════════════════════════════════════════════════════════════════════

roc_single_path <- file.path(fig_dir, "copd_tmem164_roc_single.pdf")
roc_panel_path  <- file.path(fig_dir, "copd_diagnostic_roc_comparison.pdf")

if (file.exists(roc_single_path)) {
  message("Panel C: → copd_tmem164_roc_single.pdf")
} else {
  message("⚠ ROC plot not found — run diagnostic_model.R first")
}

if (file.exists(roc_panel_path)) {
  message("Panel D: → copd_diagnostic_roc_comparison.pdf")
}

# ══════════════════════════════════════════════════════════════════════════════
# Composite (panels A + placeholder for C/D)
# ══════════════════════════════════════════════════════════════════════════════

if (exists("p_A")) {
  save_nature_fig(p_A, "fig3_composite_partial", width = 7.2, height = 3.5)
}

message("\n✔ Figure 3 panels generated")
save_session_info(report_dir)
