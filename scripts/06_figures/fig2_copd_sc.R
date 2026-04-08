###############################################################################
# fig2_copd_sc.R — Figure 2: COPD single-cell TMEM164 analysis
#
# Panels:
#   A. UMAP coloured by cell type
#   B. UMAP coloured by TMEM164 expression
#   C. Violin: TMEM164 per cell type (COPD vs Control)
#   D. Dot plot: TMEM164 + COPD markers per cell type
#   E. Proportion barplot: TMEM164-high cells per cell type
#
# Run: Rscript scripts/06_figures/fig2_copd_sc.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

seu_path <- file.path(proc_sc, "COPD_sc_annotated.rds")
if (!file.exists(seu_path)) {
  message("⚠ COPD annotated object not found")
  quit(save = "no")
}

seu <- readRDS(seu_path)

# Panel A: UMAP by cell type
p_A <- DimPlot(seu, group.by = "celltype_predicted", pt.size = 0.1,
               cols = col_celltype, label = TRUE, label.size = 2) +
  labs(title = "Cell types") +
  theme_nature() + NoLegend()

# Panel B: TMEM164 FeaturePlot
p_B <- FeaturePlot(seu, features = TARGET_GENE, pt.size = 0.1,
                   order = TRUE) +
  scale_color_viridis_c(option = "magma") +
  labs(title = TARGET_GENE) +
  theme_nature()

# Panel C: Violin
p_C <- VlnPlot(seu, features = TARGET_GENE,
               group.by = "celltype_predicted",
               split.by = "condition", pt.size = 0,
               cols = col_condition) +
  labs(title = paste0(TARGET_GENE, ": COPD vs Control")) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel D: DotPlot
copd_genes <- unique(c(TARGET_GENE, copd_known_markers))
copd_genes <- intersect(copd_genes, rownames(seu))

p_D <- DotPlot(seu, features = copd_genes,
               group.by = "celltype_predicted",
               split.by = "condition",
               cols = c("#4DBBD5", "#E64B35")) +
  RotatedAxis() +
  labs(title = "COPD markers") +
  theme_nature(base_size = 5)

# Panel E: Proportion of TMEM164-high cells
if (TARGET_GENE %in% rownames(seu)) {
  tmem_expr <- FetchData(seu, vars = TARGET_GENE)[, 1]
  seu$tmem_high <- tmem_expr > median(tmem_expr)
  
  prop_df <- as.data.frame(prop.table(
    table(seu$celltype_predicted, seu$tmem_high, seu$condition),
    margin = c(1, 3)
  ))
  colnames(prop_df) <- c("celltype", "tmem_high", "condition", "prop")
  prop_df <- prop_df[prop_df$tmem_high == TRUE, ]
  
  p_E <- ggplot(prop_df, aes(x = celltype, y = prop, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = col_condition) +
    labs(x = NULL, y = paste0("Proportion ", TARGET_GENE, "-high"),
         title = paste0(TARGET_GENE, "-high cell proportion")) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  p_E <- ggplot() + theme_void()
}

# Composite
p_top <- (p_A | p_B) + plot_layout(widths = c(1, 1))
p_mid <- p_C
p_bot <- (p_D | p_E)

p_fig2 <- p_top / p_mid / p_bot +
  plot_layout(heights = c(1, 0.8, 1)) +
  plot_annotation(tag_levels = "A")

save_nature_fig(p_fig2, "fig2_copd_sc", width = 7.2, height = 9)

message("\n✔ Figure 2 generated")
save_session_info(report_dir)
