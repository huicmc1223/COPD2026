###############################################################################
# diagnostic_model.R — TMEM164 as COPD early diagnostic biomarker
#
# Step 21: ROC curves, logistic regression, GOLD stage discrimination
#
# Methods: pROC, glmnet (LASSO for multi-gene), cross-validation
#
# Run: Rscript scripts/05_integrative/diagnostic_model.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(pROC)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# 1. Single-gene ROC: TMEM164 for COPD vs Control
# ══════════════════════════════════════════════════════════════════════════════

roc_results <- list()

for (geo_id in c("GSE57148", "GSE76925", "GSE47460")) {
  paths <- c(file.path(proc_bulk, geo_id, "vst_matrix.csv"),
             file.path(proc_array, geo_id, "expr_quantile_gene.csv"),
             file.path(proc_array, geo_id, "expr_rma_gene.csv"))
  
  expr_path <- paths[file.exists(paths)][1]
  meta_path <- file.path(raw_array, geo_id, "phenodata_1.csv")
  if (is.na(expr_path) || !file.exists(meta_path)) next
  
  expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
  meta <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)
  common <- intersect(colnames(expr), rownames(meta))
  
  if (length(common) < 10 || !(TARGET_GENE %in% rownames(expr))) next
  
  # Binary label
  cond_col <- grep("condition|disease|group", colnames(meta),
                   ignore.case = TRUE, value = TRUE)[1]
  if (is.na(cond_col)) next
  
  meta <- meta[common, , drop = FALSE]
  label <- ifelse(grepl("COPD|copd", meta[[cond_col]], ignore.case = TRUE), 1, 0)
  tmem_val <- as.numeric(expr[TARGET_GENE, common])
  
  roc_obj <- roc(label, tmem_val, quiet = TRUE)
  roc_results[[geo_id]] <- list(roc = roc_obj, geo_id = geo_id,
                                 auc = auc(roc_obj))
  message(geo_id, ": AUC = ", round(auc(roc_obj), 3))
}

# Combined ROC plot
if (length(roc_results) > 0) {
  pdf(file.path(fig_dir, "copd_tmem164_roc_single.pdf"),
      width = 4, height = 4)
  
  colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")
  plot(roc_results[[1]]$roc, col = colors[1],
       main = paste0(TARGET_GENE, " ROC: COPD vs Control"),
       cex.main = 0.8)
  
  legend_text <- paste0(roc_results[[1]]$geo_id, " (AUC=",
                        round(roc_results[[1]]$auc, 3), ")")
  
  for (i in seq_along(roc_results)[-1]) {
    lines(roc_results[[i]]$roc, col = colors[i])
    legend_text <- c(legend_text,
                     paste0(roc_results[[i]]$geo_id, " (AUC=",
                            round(roc_results[[i]]$auc, 3), ")"))
  }
  
  legend("bottomright", legend = legend_text,
         col = colors[seq_along(roc_results)],
         lwd = 2, cex = 0.6, bty = "n")
  abline(a = 0, b = 1, lty = 2, col = "grey50")
  
  dev.off()
  message("✔ Single-gene ROC saved")
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. GOLD stage discrimination: early (GOLD 1) vs severe (GOLD 3-4)
# ══════════════════════════════════════════════════════════════════════════════

gold_roc <- function() {
  for (geo_id in c("GSE310058", "GSE76925")) {
    paths <- c(file.path(proc_bulk, geo_id, "vst_matrix.csv"),
               file.path(proc_array, geo_id, "expr_quantile_gene.csv"))
    expr_path <- paths[file.exists(paths)][1]
    meta_path <- file.path(raw_array, geo_id, "phenodata_1.csv")
    
    if (is.na(expr_path) || !file.exists(meta_path)) next
    
    expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
    meta <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)
    common <- intersect(colnames(expr), rownames(meta))
    
    gold_col <- grep("gold|stage|grade", colnames(meta),
                     ignore.case = TRUE, value = TRUE)[1]
    if (is.na(gold_col) || !(TARGET_GENE %in% rownames(expr))) next
    
    meta <- meta[common, , drop = FALSE]
    gold <- as.character(meta[[gold_col]])
    
    # Early (GOLD 1) vs Severe (GOLD 3-4)
    early   <- which(gold %in% c("1", "GOLD1", "GOLD 1"))
    severe  <- which(gold %in% c("3", "4", "GOLD3", "GOLD4",
                                   "GOLD 3", "GOLD 4"))
    
    if (length(early) >= 5 && length(severe) >= 5) {
      idx <- c(early, severe)
      label <- c(rep(0, length(early)), rep(1, length(severe)))
      tmem <- as.numeric(expr[TARGET_GENE, common[idx]])
      
      r <- roc(label, tmem, quiet = TRUE)
      message(geo_id, ": GOLD Early vs Severe AUC = ", round(auc(r), 3))
      return(r)
    }
  }
  NULL
}

gold_roc_obj <- gold_roc()

# ══════════════════════════════════════════════════════════════════════════════
# 3. Multi-gene panel (LASSO)
# ══════════════════════════════════════════════════════════════════════════════

build_lasso_panel <- function() {
  # Use largest COPD bulk dataset
  for (geo_id in c("GSE76925", "GSE57148")) {
    paths <- c(file.path(proc_bulk, geo_id, "vst_matrix.csv"),
               file.path(proc_array, geo_id, "expr_quantile_gene.csv"))
    expr_path <- paths[file.exists(paths)][1]
    meta_path <- file.path(raw_array, geo_id, "phenodata_1.csv")
    
    if (is.na(expr_path) || !file.exists(meta_path)) next
    
    expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
    meta <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)
    common <- intersect(colnames(expr), rownames(meta))
    
    cond_col <- grep("condition|disease|group", colnames(meta),
                     ignore.case = TRUE, value = TRUE)[1]
    if (is.na(cond_col)) next
    
    meta <- meta[common, , drop = FALSE]
    y <- ifelse(grepl("COPD|copd", meta[[cond_col]]), 1, 0)
    
    # Candidate panel: TMEM164 + known markers
    candidates <- intersect(c(TARGET_GENE, copd_known_markers),
                            rownames(expr))
    X <- t(expr[candidates, common])
    
    # LASSO with CV
    set.seed(2026)
    cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 5)
    
    # Selected features
    coefs <- coef(cv_fit, s = "lambda.min")
    selected <- rownames(coefs)[coefs[, 1] != 0]
    selected <- setdiff(selected, "(Intercept)")
    message("LASSO selected: ", paste(selected, collapse = ", "))
    
    # ROC of panel
    pred <- predict(cv_fit, X, s = "lambda.min", type = "response")
    panel_roc <- roc(y, as.numeric(pred), quiet = TRUE)
    message("Panel AUC: ", round(auc(panel_roc), 3))
    
    # Save coefficients
    coef_df <- data.frame(gene = rownames(coefs),
                          coefficient = as.numeric(coefs[, 1]))
    coef_df <- coef_df[coef_df$coefficient != 0, ]
    write.csv(coef_df,
              file.path(table_dir, "diagnostic_lasso_coefficients.csv"),
              row.names = FALSE)
    
    return(list(cv_fit = cv_fit, roc = panel_roc, selected = selected,
                dataset = geo_id))
  }
  NULL
}

lasso_res <- tryCatch(build_lasso_panel(), error = function(e) {
  message("LASSO panel failed: ", e$message)
  NULL
})

if (!is.null(lasso_res)) {
  # Single gene vs panel comparison
  if (length(roc_results) > 0) {
    best_single <- roc_results[[which.max(sapply(roc_results,
                                                   function(x) x$auc))]]
    
    pdf(file.path(fig_dir, "copd_diagnostic_roc_comparison.pdf"),
        width = 4, height = 4)
    plot(best_single$roc, col = "#4DBBD5",
         main = paste0(TARGET_GENE, ": single vs panel"),
         cex.main = 0.8)
    lines(lasso_res$roc, col = "#E64B35")
    legend("bottomright",
           legend = c(paste0("Single (AUC=",
                             round(auc(best_single$roc), 3), ")"),
                      paste0("Panel (AUC=",
                             round(auc(lasso_res$roc), 3), ")")),
           col = c("#4DBBD5", "#E64B35"), lwd = 2, cex = 0.7, bty = "n")
    abline(a = 0, b = 1, lty = 2, col = "grey50")
    dev.off()
    
    message("✔ Diagnostic model comparison saved")
  }
}

message("\n✔ Diagnostic model analysis complete")
save_session_info(report_dir)
