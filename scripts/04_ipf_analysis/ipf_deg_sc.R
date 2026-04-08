###############################################################################
# ipf_deg_sc.R — TMEM164 single-cell differential expression in IPF
#
# Step 16: Cell-type–resolved analysis + KRT8+ transitional state identification
#
# Focus: AT1, AT2, Basal (atypical), Myofibroblast, Profibrotic AM
#
# Run: Rscript scripts/04_ipf_analysis/ipf_deg_sc.R
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
# Load annotated IPF scRNA-seq object
# ══════════════════════════════════════════════════════════════════════════════

seu_path <- file.path(proc_sc, "IPF_sc_annotated.rds")
if (!file.exists(seu_path)) {
  message("⚠ IPF annotated object not found")
  quit(save = "no")
}

seu <- readRDS(seu_path)
message("Loaded: ", ncol(seu), " cells")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Per-cell-type TMEM164 differential expression
# ══════════════════════════════════════════════════════════════════════════════

celltypes <- unique(seu$celltype_predicted)
deg_results <- list()

for (ct in celltypes) {
  sub <- subset(seu, celltype_predicted == ct)
  n_ipf  <- sum(sub$condition == "IPF")
  n_ctrl <- sum(sub$condition == "Control")
  
  if (n_ipf < 10 || n_ctrl < 10) next
  if (!(TARGET_GENE %in% rownames(sub))) next
  
  Idents(sub) <- "condition"
  markers <- FindMarkers(sub, ident.1 = "IPF", ident.2 = "Control",
                         features = TARGET_GENE, logfc.threshold = 0,
                         min.pct = 0, test.use = "wilcox")
  
  if (nrow(markers) > 0) {
    markers$celltype <- ct
    markers$n_ipf <- n_ipf
    markers$n_ctrl <- n_ctrl
    deg_results[[ct]] <- markers
  }
}

if (length(deg_results) > 0) {
  deg_df <- rbindlist(deg_results, fill = TRUE)
  write.csv(deg_df, file.path(table_dir, "ipf_tmem164_sc_deg.csv"))
  message("TMEM164 per cell type (IPF vs Control):")
  print(deg_df[, .(celltype, avg_log2FC, p_val_adj)])
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. KRT8+ transitional state identification
# ══════════════════════════════════════════════════════════════════════════════

krt8_markers <- markers_ipf_specific$KRT8_transitional  # KRT8, CLDN4, LGALS3
krt8_present <- intersect(krt8_markers, rownames(seu))

if (length(krt8_present) >= 2) {
  message("\n── KRT8+ transitional state analysis ──")
  
  seu <- AddModuleScore(seu, features = list(krt8_present),
                        name = "KRT8_transitional", seed = 2026)
  
  # Threshold: cells in top 10% of KRT8 transitional score
  threshold <- quantile(seu$KRT8_transitional1, 0.90, na.rm = TRUE)
  seu$is_KRT8_transitional <- seu$KRT8_transitional1 > threshold
  
  n_krt8 <- sum(seu$is_KRT8_transitional)
  message("  KRT8+ transitional cells: ", n_krt8, " (",
          round(100 * n_krt8 / ncol(seu), 1), "%)")
  
  # TMEM164 in KRT8+ vs other epithelial
  epithelial_types <- c("AT1", "AT2", "Basal", "Club")
  epi <- subset(seu, celltype_predicted %in% epithelial_types)
  
  if (TARGET_GENE %in% rownames(epi) && ncol(epi) > 100) {
    Idents(epi) <- "is_KRT8_transitional"
    krt8_deg <- FindMarkers(epi, ident.1 = TRUE, ident.2 = FALSE,
                            features = TARGET_GENE, logfc.threshold = 0,
                            min.pct = 0)
    if (nrow(krt8_deg) > 0) {
      message("  TMEM164 in KRT8+ vs other epithelial: log2FC=",
              round(krt8_deg$avg_log2FC, 3),
              ", padj=", signif(krt8_deg$p_val_adj, 3))
    }
  }
  
  # Highlight on UMAP
  p_krt8 <- DimPlot(seu, cells.highlight = colnames(seu)[seu$is_KRT8_transitional],
                    cols.highlight = "#DC0000", cols = "grey80",
                    pt.size = 0.1) +
    labs(title = "KRT8+ transitional cells") +
    theme_nature() + NoLegend()
  
  p_krt8_score <- FeaturePlot(seu, features = "KRT8_transitional1",
                               pt.size = 0.1) +
    scale_color_viridis_c() +
    labs(title = "KRT8+ transitional score") +
    theme_nature()
  
  save_nature_fig(p_krt8 | p_krt8_score, "ipf_krt8_transitional",
                  width = 7.2, height = 3.5)
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Visualisation
# ══════════════════════════════════════════════════════════════════════════════

# Violin plot
p_violin <- VlnPlot(seu, features = TARGET_GENE,
                    group.by = "celltype_predicted",
                    split.by = "condition", pt.size = 0,
                    cols = col_condition) +
  labs(title = paste0(TARGET_GENE, " in IPF vs Control")) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_nature_fig(p_violin, "ipf_tmem164_violin_celltype",
                width = 7.2, height = 4)

# DotPlot with IPF markers
ipf_genes <- unique(c(TARGET_GENE, ipf_known_markers,
                       markers_ipf_specific$Atypical_Basal[1:4]))
ipf_genes <- intersect(ipf_genes, rownames(seu))

p_dot <- DotPlot(seu, features = ipf_genes,
                 group.by = "celltype_predicted",
                 split.by = "condition",
                 cols = c("#4DBBD5", "#F39B7F")) +
  RotatedAxis() +
  labs(title = "TMEM164 + IPF markers") +
  theme_nature(base_size = 6)

save_nature_fig(p_dot, "ipf_tmem164_dotplot", width = 7.2, height = 5)

# Save updated object
saveRDS(seu, file.path(proc_sc, "IPF_sc_annotated.rds"))

message("\n✔ IPF single-cell DEG analysis complete")
save_session_info(report_dir)
