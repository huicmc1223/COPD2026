###############################################################################
# fig4_ipf_trajectory.R — Figure 4: IPF epithelial-mesenchymal trajectory
#
# Panels:
#   A. UMAP pseudotime
#   B. TMEM164 along AT2→KRT8+→AT1 trajectory
#   C. EMT score vs pseudotime
#   D. KRT8+ transitional cells highlighted
#   E. TMEM164 in KRT8+ vs other epithelial
#
# Run: Rscript scripts/06_figures/fig4_ipf_trajectory.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# Load trajectory data
traj_path <- file.path(table_dir, "ipf_emt_trajectory_data.csv")
epi_path  <- file.path(proc_sc, "IPF_epithelial_trajectory.rds")

if (!file.exists(traj_path)) {
  message("⚠ Trajectory data not found — run ipf_emt_trajectory.R first")
  quit(save = "no")
}

traj_df <- fread(traj_path)

# Panel B: TMEM164 along pseudotime
p_B <- ggplot(traj_df, aes_string(x = "pseudotime", y = TARGET_GENE)) +
  geom_point(aes(colour = celltype_predicted), size = 0.2, alpha = 0.3) +
  geom_smooth(method = "loess", colour = "#E64B35", linewidth = 1) +
  scale_color_manual(values = col_celltype) +
  labs(x = "Pseudotime (AT2 → AT1)", y = paste0(TARGET_GENE, " expression"),
       title = paste0(TARGET_GENE, " trajectory")) +
  theme_nature()

# Panel C: EMT score
if ("EMT_score" %in% colnames(traj_df)) {
  p_C <- ggplot(traj_df, aes(x = pseudotime, y = EMT_score)) +
    geom_point(aes(colour = celltype_predicted), size = 0.2, alpha = 0.3) +
    geom_smooth(method = "loess", colour = "#3C5488", linewidth = 1) +
    scale_color_manual(values = col_celltype) +
    labs(x = "Pseudotime", y = "EMT score") +
    theme_nature()
} else {
  p_C <- ggplot() + theme_void()
}

# Panel D: IPF vs Control overlay
if ("condition" %in% colnames(traj_df)) {
  p_D <- ggplot(traj_df, aes_string(x = "pseudotime", y = TARGET_GENE,
                                      colour = "condition")) +
    geom_smooth(method = "loess", linewidth = 1) +
    scale_color_manual(values = col_condition) +
    labs(x = "Pseudotime", y = paste0(TARGET_GENE, " expression"),
         title = "IPF vs Control") +
    theme_nature()
} else {
  p_D <- ggplot() + theme_void()
}

# Additional panels from saved objects
if (file.exists(epi_path)) {
  epi <- readRDS(epi_path)
  
  # Panel A: UMAP coloured by pseudotime
  if ("pseudotime" %in% colnames(epi@meta.data)) {
    p_A <- FeaturePlot(epi, features = "pseudotime", pt.size = 0.3) +
      scale_color_viridis_c(option = "magma") +
      labs(title = "Pseudotime") +
      theme_nature()
  } else {
    p_A <- DimPlot(epi, group.by = "celltype_predicted",
                   cols = col_celltype, pt.size = 0.3) +
      theme_nature()
  }
  
  # Panel E: TMEM164 in KRT8+ vs non-KRT8
  if ("is_KRT8_transitional" %in% colnames(epi@meta.data) &&
      TARGET_GENE %in% rownames(epi)) {
    p_E <- VlnPlot(epi, features = TARGET_GENE,
                   group.by = "is_KRT8_transitional",
                   cols = c("FALSE" = "#4DBBD5", "TRUE" = "#E64B35"),
                   pt.size = 0) +
      labs(x = "KRT8+ transitional", title = paste0(TARGET_GENE,
                                                      " in KRT8+ state")) +
      theme_nature()
  } else {
    p_E <- ggplot() + theme_void()
  }
  
  rm(epi); gc()
} else {
  p_A <- ggplot() + annotate("text", x = 0.5, y = 0.5,
                              label = "Run ipf_emt_trajectory.R") + theme_void()
  p_E <- p_A
}

# Composite
p_fig4 <- (p_A | p_B) / (p_C | p_D) / p_E +
  plot_layout(heights = c(1, 1, 0.8)) +
  plot_annotation(tag_levels = "A")

save_nature_fig(p_fig4, "fig4_ipf_trajectory", width = 7.2, height = 9)

message("\n✔ Figure 4 generated")
save_session_info(report_dir)
