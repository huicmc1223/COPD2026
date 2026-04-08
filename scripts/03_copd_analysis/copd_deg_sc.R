###############################################################################
# copd_deg_sc.R — TMEM164 single-cell level differential expression in COPD
#
# Step 12: Cell-type–resolved TMEM164 expression (COPD vs Control)
#
# Input:  data/processed/scrnaseq/COPD_sc_annotated.rds
# Output: results/tables/copd_tmem164_sc_deg.csv
#         results/figures/copd_tmem164_violin_celltype.pdf
#         results/figures/copd_tmem164_dotplot.pdf
#
# Run: Rscript scripts/03_copd_analysis/copd_deg_sc.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# Load annotated COPD scRNA-seq object
# ══════════════════════════════════════════════════════════════════════════════

seu_path <- file.path(proc_sc, "COPD_sc_annotated.rds")
if (!file.exists(seu_path)) {
  message("⚠ COPD annotated object not found — run celltype_annotation.R first")
  quit(save = "no")
}

seu <- readRDS(seu_path)
message("Loaded: ", ncol(seu), " cells, ",
        length(unique(seu$celltype_predicted)), " cell types")

# Ensure condition and celltype columns exist
stopifnot("condition" %in% colnames(seu@meta.data))
stopifnot("celltype_predicted" %in% colnames(seu@meta.data))

# ══════════════════════════════════════════════════════════════════════════════
# Per-cell-type differential expression of TMEM164
# ══════════════════════════════════════════════════════════════════════════════

celltypes <- unique(seu$celltype_predicted)
deg_results <- list()

for (ct in celltypes) {
  sub <- subset(seu, celltype_predicted == ct)
  n_copd <- sum(sub$condition == "COPD")
  n_ctrl <- sum(sub$condition == "Control")
  
  if (n_copd < 10 || n_ctrl < 10) {
    message("  ", ct, ": skipping (COPD=", n_copd, ", Ctrl=", n_ctrl, ")")
    next
  }
  
  # Wilcoxon test for TMEM164
  if (TARGET_GENE %in% rownames(sub)) {
    Idents(sub) <- "condition"
    markers <- FindMarkers(sub, ident.1 = "COPD", ident.2 = "Control",
                           features = TARGET_GENE, logfc.threshold = 0,
                           min.pct = 0, test.use = "wilcox")
    
    if (nrow(markers) > 0) {
      markers$celltype <- ct
      markers$gene <- TARGET_GENE
      markers$n_copd <- n_copd
      markers$n_ctrl <- n_ctrl
      deg_results[[ct]] <- markers
    }
  }
}

# Combine
if (length(deg_results) > 0) {
  deg_df <- do.call(rbind, deg_results)
  deg_df <- deg_df[order(deg_df$p_val), ]
  write.csv(deg_df, file.path(table_dir, "copd_tmem164_sc_deg.csv"))
  message("\n── TMEM164 per cell type (COPD vs Control) ──")
  print(deg_df[, c("celltype", "avg_log2FC", "p_val_adj", "pct.1", "pct.2")])
}

# ══════════════════════════════════════════════════════════════════════════════
# Visualisation
# ══════════════════════════════════════════════════════════════════════════════

# 1. Violin plot: TMEM164 per cell type, split by condition
p_violin <- VlnPlot(seu, features = TARGET_GENE,
                    group.by = "celltype_predicted",
                    split.by = "condition",
                    pt.size = 0, cols = col_condition) +
  labs(title = paste0(TARGET_GENE, " expression in COPD vs Control"),
       x = "Cell type", y = "Expression") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_nature_fig(p_violin, "copd_tmem164_violin_celltype",
                width = 7.2, height = 4)

# 2. DotPlot: TMEM164 + key COPD markers per cell type per condition
genes_plot <- unique(c(TARGET_GENE, copd_known_markers))
genes_plot <- intersect(genes_plot, rownames(seu))

p_dot <- DotPlot(seu, features = genes_plot,
                 group.by = "celltype_predicted",
                 split.by = "condition",
                 cols = c("#4DBBD5", "#E64B35"),
                 dot.scale = 5) +
  RotatedAxis() +
  labs(title = "TMEM164 + COPD markers by cell type & condition") +
  theme_nature(base_size = 6)

save_nature_fig(p_dot, "copd_tmem164_dotplot",
                width = 7.2, height = 5)

# 3. FeaturePlot split by condition
p_feat <- FeaturePlot(seu, features = TARGET_GENE, split.by = "condition",
                      pt.size = 0.1, order = TRUE) &
  scale_color_viridis_c() &
  theme_nature()

save_nature_fig(p_feat, "copd_tmem164_feature_split",
                width = 7.2, height = 3.5)

# 4. Cell proportion barplot
prop_df <- as.data.frame(prop.table(
  table(seu$celltype_predicted, seu$condition), margin = 2
))
colnames(prop_df) <- c("CellType", "Condition", "Proportion")

p_prop <- ggplot(prop_df, aes(x = Condition, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = col_celltype) +
  labs(title = "Cell type proportions", y = "Proportion") +
  theme_nature()

save_nature_fig(p_prop, "copd_celltype_proportions", width = 3.5, height = 4)

message("\n✔ COPD single-cell DEG analysis complete")
save_session_info(report_dir)
