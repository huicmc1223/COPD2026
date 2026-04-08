###############################################################################
# copd_pathway.R — GSEA and pathway enrichment for TMEM164 in COPD
#
# Step 14: GSEA (KEGG, Hallmark) + cell-type–specific pathway activity (AUCell)
#
# Run: Rscript scripts/03_copd_analysis/copd_pathway.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# 1. Bulk-level GSEA (using DESeq2 results)
# ══════════════════════════════════════════════════════════════════════════════

run_gsea_bulk <- function(dds_path, geo_id) {
  message("\n── GSEA: ", geo_id, " ──")
  
  dds <- readRDS(dds_path)
  res <- results(dds, contrast = c("condition", "COPD", "Control"))
  res <- res[!is.na(res$padj), ]
  
  # Ranked gene list by log2FC × -log10(pvalue)
  gene_list <- res$log2FoldChange
  names(gene_list) <- rownames(res)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Map to Entrez IDs
  gene_map <- bitr(names(gene_list), fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list_entrez <- gene_list[gene_map$SYMBOL]
  names(gene_list_entrez) <- gene_map$ENTREZID
  gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
  
  # KEGG GSEA
  gsea_kegg <- gseKEGG(geneList     = gene_list_entrez,
                        organism     = "hsa",
                        minGSSize    = 15,
                        maxGSSize    = 500,
                        pvalueCutoff = 0.1,
                        verbose      = FALSE,
                        seed         = TRUE)
  
  # Hallmark gene sets
  hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene)
  
  gsea_hallmark <- GSEA(gene_list_entrez,
                        TERM2GENE = hallmark,
                        minGSSize = 15,
                        pvalueCutoff = 0.1,
                        verbose = FALSE,
                        seed = TRUE)
  
  return(list(kegg = gsea_kegg, hallmark = gsea_hallmark, ranked = gene_list))
}

# Run on available datasets
gsea_results <- list()
for (geo_id in c("GSE57148")) {
  dds_path <- file.path(proc_bulk, geo_id, "dds.rds")
  if (file.exists(dds_path)) {
    gsea_results[[geo_id]] <- tryCatch(
      run_gsea_bulk(dds_path, geo_id),
      error = function(e) { message("  ⚠ ", e$message); NULL }
    )
  }
}

# Save GSEA results and plots
for (geo_id in names(gsea_results)) {
  res <- gsea_results[[geo_id]]
  if (is.null(res)) next
  
  out_dir <- file.path(table_dir, paste0("gsea_", geo_id))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # KEGG results
  if (!is.null(res$kegg) && nrow(as.data.frame(res$kegg)) > 0) {
    write.csv(as.data.frame(res$kegg),
              file.path(out_dir, "gsea_kegg.csv"), row.names = FALSE)
    
    # Dot plot
    p_kegg <- dotplot(res$kegg, showCategory = 20, font.size = 6) +
      labs(title = paste0(geo_id, " — KEGG GSEA")) +
      theme_nature(base_size = 6)
    save_nature_fig(p_kegg, paste0(geo_id, "_gsea_kegg_dotplot"),
                    width = 5, height = 5)
    
    # Check focus pathways (ferroptosis, cytokine, NFkB etc.)
    focus_ids <- unlist(kegg_focus)
    focus_res <- as.data.frame(res$kegg)
    focus_res <- focus_res[focus_res$ID %in% focus_ids, ]
    if (nrow(focus_res) > 0) {
      message("  Focus pathways enriched:")
      print(focus_res[, c("ID", "Description", "NES", "pvalue", "p.adjust")])
    }
  }
  
  # Hallmark results
  if (!is.null(res$hallmark) && nrow(as.data.frame(res$hallmark)) > 0) {
    write.csv(as.data.frame(res$hallmark),
              file.path(out_dir, "gsea_hallmark.csv"), row.names = FALSE)
    
    p_hall <- dotplot(res$hallmark, showCategory = 20, font.size = 6) +
      labs(title = paste0(geo_id, " — Hallmark GSEA")) +
      theme_nature(base_size = 6)
    save_nature_fig(p_hall, paste0(geo_id, "_gsea_hallmark_dotplot"),
                    width = 5, height = 5)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. Cell-type–specific pathway activity (AUCell on scRNA-seq)
# ══════════════════════════════════════════════════════════════════════════════

seu_path <- file.path(proc_sc, "COPD_sc_annotated.rds")
if (file.exists(seu_path) && requireNamespace("AUCell", quietly = TRUE)) {
  message("\n── AUCell pathway scoring ──")
  
  library(AUCell)
  seu <- readRDS(seu_path)
  
  # Build gene sets for focus pathways
  pathway_genesets <- list(
    Ferroptosis = ferroptosis_genes,
    Autophagy   = autophagy_genes,
    Hippo_YAP   = hippo_yap_genes,
    TGFb        = epithelial_fate$TGFb_axis,
    NFkB        = epithelial_fate$NFkB_axis,
    Glycolysis  = epithelial_fate$Glycolysis,
    Lipid_meta  = epithelial_fate$Lipid_meta
  )
  
  # Filter to genes present in dataset
  pathway_genesets <- lapply(pathway_genesets, function(gs) {
    intersect(gs, rownames(seu))
  })
  pathway_genesets <- pathway_genesets[sapply(pathway_genesets, length) >= 3]
  
  # AUCell scoring
  expr_mat <- GetAssayData(seu, layer = "counts")
  cells_rankings <- AUCell_buildRankings(expr_mat, plotStats = FALSE)
  
  geneSets <- lapply(names(pathway_genesets), function(nm) {
    GeneSet(pathway_genesets[[nm]], setName = nm)
  })
  names(geneSets) <- names(pathway_genesets)
  
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  
  # Add AUC scores to Seurat metadata
  auc_mat <- getAUC(cells_AUC)
  for (pw in rownames(auc_mat)) {
    seu[[paste0("AUC_", pw)]] <- auc_mat[pw, colnames(seu)]
  }
  
  # Heatmap: pathway AUC by cell type × condition
  auc_summary <- list()
  for (pw in rownames(auc_mat)) {
    col_name <- paste0("AUC_", pw)
    agg <- aggregate(seu@meta.data[[col_name]],
                     by = list(celltype = seu$celltype_predicted,
                               condition = seu$condition),
                     FUN = mean)
    agg$pathway <- pw
    colnames(agg)[3] <- "mean_AUC"
    auc_summary[[pw]] <- agg
  }
  auc_df <- rbindlist(auc_summary)
  write.csv(auc_df, file.path(table_dir, "copd_pathway_aucell_summary.csv"),
            row.names = FALSE)
  
  # Bubble plot
  p_auc <- ggplot(auc_df, aes(x = pathway, y = celltype,
                                size = mean_AUC, colour = condition)) +
    geom_point(alpha = 0.8) +
    scale_size_continuous(range = c(1, 6)) +
    scale_color_manual(values = col_condition) +
    labs(title = "Pathway activity (AUCell) by cell type & condition",
         x = "Pathway", y = "Cell type") +
    theme_nature(base_size = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  save_nature_fig(p_auc, "copd_pathway_aucell_bubble",
                  width = 6, height = 5)
  
  rm(seu); gc()
  message("  ✔ AUCell pathway analysis saved")
}

message("\n✔ COPD pathway analysis complete")
save_session_info(report_dir)
