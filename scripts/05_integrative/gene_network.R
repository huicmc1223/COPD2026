###############################################################################
# gene_network.R — TMEM164 gene regulatory network analysis
#
# Step 22: WGCNA hub analysis, CellChat ligand-receptor, TF enrichment
#
# Run: Rscript scripts/05_integrative/gene_network.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(data.table)
})

allowWGCNAThreads()

# ══════════════════════════════════════════════════════════════════════════════
# 1. WGCNA: identify TMEM164 module & hub genes
# ══════════════════════════════════════════════════════════════════════════════

run_wgcna_hub <- function(expr, meta, geo_id) {
  # Transpose for WGCNA: rows = samples, cols = genes
  expr_t <- t(expr)
  
  # Variance filter: top 5000 variable genes
  gene_var <- apply(expr_t, 2, var)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(5000, ncol(expr_t))]
  if (!(TARGET_GENE %in% top_genes)) {
    top_genes <- c(TARGET_GENE, top_genes[1:(length(top_genes) - 1)])
  }
  expr_t <- expr_t[, top_genes]
  
  # Pick soft threshold
  powers <- c(1:10, seq(12, 20, by = 2))
  sft <- pickSoftThreshold(expr_t, powerVector = powers,
                            networkType = "signed", verbose = 0)
  power_use <- sft$powerEstimate
  if (is.na(power_use) || power_use < 1) power_use <- 6
  message("  Soft threshold power: ", power_use)
  
  # Build network
  net <- blockwiseModules(expr_t, power = power_use,
                          networkType = "signed",
                          TOMType = "signed",
                          minModuleSize = 30,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          numericLabels = TRUE,
                          pamRespectsDendro = FALSE,
                          verbose = 0)
  
  # Which module contains TMEM164?
  module_labels <- net$colors
  names(module_labels) <- colnames(expr_t)
  tmem_module <- module_labels[TARGET_GENE]
  message("  TMEM164 module: ", tmem_module, " (",
          sum(module_labels == tmem_module), " genes)")
  
  # Module membership (kME)
  MEs <- moduleEigengenes(expr_t, module_labels)$eigengenes
  kME <- cor(expr_t, MEs, use = "pairwise.complete.obs")
  me_col <- paste0("ME", tmem_module)
  
  if (me_col %in% colnames(kME)) {
    module_genes <- names(module_labels[module_labels == tmem_module])
    kME_tmem <- kME[module_genes, me_col]
    hub_df <- data.frame(gene = module_genes,
                         kME = kME_tmem[module_genes])
    hub_df <- hub_df[order(-hub_df$kME), ]
    
    write.csv(hub_df,
              file.path(table_dir, paste0("wgcna_", geo_id, "_module",
                                           tmem_module, "_genes.csv")),
              row.names = FALSE)
    
    message("  Top hub genes: ",
            paste(head(hub_df$gene, 10), collapse = ", "))
    
    return(list(hub_df = hub_df, module = tmem_module, net = net,
                MEs = MEs, module_labels = module_labels))
  }
  NULL
}

# Run on best COPD dataset
for (geo_id in c("GSE76925", "GSE57148")) {
  paths <- c(file.path(proc_bulk, geo_id, "vst_matrix.csv"),
             file.path(proc_array, geo_id, "expr_quantile_gene.csv"))
  expr_path <- paths[file.exists(paths)][1]
  
  if (is.na(expr_path)) next
  
  expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
  if (!(TARGET_GENE %in% rownames(expr))) next
  
  message("\n── WGCNA: ", geo_id, " ──")
  wgcna_res <- tryCatch(run_wgcna_hub(expr, NULL, geo_id), error = function(e) {
    message("  WGCNA failed: ", e$message)
    NULL
  })
  
  if (!is.null(wgcna_res)) break
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. CellChat: ligand-receptor interactions involving TMEM164
# ══════════════════════════════════════════════════════════════════════════════

run_cellchat_analysis <- function() {
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    message("  ⚠ CellChat not installed — skipping")
    return(NULL)
  }
  
  library(CellChat)
  library(Seurat)
  
  seu_path <- file.path(proc_sc, "COPD_sc_annotated.rds")
  if (!file.exists(seu_path)) return(NULL)
  
  seu <- readRDS(seu_path)
  
  # Create CellChat objects for COPD and Control
  cellchat_list <- list()
  for (cond in c("COPD", "Control")) {
    sub <- subset(seu, condition == cond)
    if (ncol(sub) < 100) next
    
    data_input <- GetAssayData(sub, slot = "data")
    meta <- sub@meta.data[, c("celltype_predicted")]
    meta <- data.frame(labels = meta, row.names = colnames(sub))
    
    cc <- createCellChat(object = data_input, meta = meta,
                         group.by = "labels")
    CellChatDB <- CellChatDB.human
    cc@DB <- CellChatDB
    
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc)
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    
    cellchat_list[[cond]] <- cc
  }
  
  if (length(cellchat_list) == 2) {
    merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
    saveRDS(merged, file.path(proc_sc, "cellchat_merged.rds"))
    
    # Differential interactions
    pdf(file.path(fig_dir, "cellchat_diff_network.pdf"),
        width = 7, height = 5)
    netVisual_diffInteraction(merged, weight.scale = TRUE)
    dev.off()
    
    message("  ✔ CellChat analysis saved")
  }
  
  cellchat_list
}

cc_res <- tryCatch(run_cellchat_analysis(), error = function(e) {
  message("CellChat failed: ", e$message)
  NULL
})

# ══════════════════════════════════════════════════════════════════════════════
# 3. Transcription factor enrichment (TMEM164 co-expressed genes)
# ══════════════════════════════════════════════════════════════════════════════

if (!is.null(wgcna_res)) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  top_hub_genes <- head(wgcna_res$hub_df$gene, 200)
  
  # GO enrichment of hub genes
  ego <- enrichGO(gene = top_hub_genes, OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL", ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05)
  
  if (nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego),
              file.path(table_dir, "wgcna_hub_GO_enrichment.csv"),
              row.names = FALSE)
    
    p_go <- dotplot(ego, showCategory = 15) +
      labs(title = paste0(TARGET_GENE, " module hub genes: GO BP")) +
      theme_nature(base_size = 6)
    save_nature_fig(p_go, "wgcna_hub_go_dotplot", width = 5, height = 5)
    
    message("  ✔ Hub gene GO enrichment saved")
  }
}

message("\n✔ Gene network analysis complete")
save_session_info(report_dir)
