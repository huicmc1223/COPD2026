###############################################################################
# copd_gold_staging.R — TMEM164 expression across GOLD 1→2→3→4 stages
#
# Step 11: Stage-dependent expression trajectory
#
# Key datasets:
#   GSE310058 — snRNA-seq GOLD 1-4 (primary)
#   GSE38974  — Microarray GOLD 1,2,4
#   GSE76925  — Microarray GOLD 3-4
#   GSE57148  — Bulk RNA-seq GOLD 1-2
#
# Output: results/figures/copd_tmem164_gold_staging.pdf
#         results/tables/copd_tmem164_gold_trend.csv
#
# Run: Rscript scripts/03_copd_analysis/copd_gold_staging.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# Collect TMEM164 expression with GOLD stage annotation
# ══════════════════════════════════════════════════════════════════════════════

all_staging <- list()

# --- Bulk RNA-seq: GSE57148 (GOLD 1-2) ---
vst_path <- file.path(proc_bulk, "GSE57148", "vst_matrix.csv")
if (file.exists(vst_path)) {
  vst_mat <- read.csv(vst_path, row.names = 1, check.names = FALSE)
  meta <- read.csv(file.path(raw_bulk, "GSE57148", "phenodata_1.csv"),
                   row.names = 1, stringsAsFactors = FALSE)
  if (TARGET_GENE %in% rownames(vst_mat)) {
    common <- intersect(colnames(vst_mat), rownames(meta))
    df <- data.frame(
      sample    = common,
      TMEM164   = as.numeric(vst_mat[TARGET_GENE, common]),
      dataset   = "GSE57148",
      stringsAsFactors = FALSE
    )
    # GOLD stage needs to be extracted from metadata — adapt column names
    # Placeholder: assume a "gold_stage" column exists after curation
    if ("gold_stage" %in% colnames(meta)) {
      df$gold <- meta[common, "gold_stage"]
    } else {
      df$gold <- "1-2"  # Default from dataset description
    }
    all_staging[["GSE57148"]] <- df
  }
}

# --- Microarray datasets ---
for (geo_id in c("GSE76925", "GSE38974", "GSE47460")) {
  expr_path <- file.path(proc_array, geo_id, "expr_quantile_gene.csv")
  meta_path <- file.path(raw_array, geo_id, "phenodata_1.csv")
  
  if (file.exists(expr_path) && file.exists(meta_path)) {
    expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
    meta <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)
    
    if (TARGET_GENE %in% rownames(expr)) {
      common <- intersect(colnames(expr), rownames(meta))
      df <- data.frame(
        sample  = common,
        TMEM164 = as.numeric(expr[TARGET_GENE, common]),
        dataset = geo_id,
        stringsAsFactors = FALSE
      )
      if ("gold_stage" %in% colnames(meta)) {
        df$gold <- meta[common, "gold_stage"]
      }
      all_staging[[geo_id]] <- df
    }
  }
}

# --- snRNA-seq: GSE310058 (best coverage: GOLD 1-4) ---
sc_path <- file.path(proc_sc, "GSE310058")
rds_files <- list.files(sc_path, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) > 0) {
  seu <- readRDS(rds_files[1])
  if (TARGET_GENE %in% rownames(seu)) {
    expr_vals <- FetchData(seu, vars = c(TARGET_GENE, "condition"))
    # GSE310058 has GOLD 1-4 in metadata — needs curation
    if ("gold_stage" %in% colnames(seu@meta.data)) {
      expr_vals$gold <- seu@meta.data$gold_stage
    }
    # Pseudobulk: mean expression per sample per GOLD stage
    if ("orig.ident" %in% colnames(seu@meta.data)) {
      expr_vals$sample <- seu@meta.data$orig.ident
      pb <- aggregate(as.formula(paste0(TARGET_GENE, " ~ sample")),
                      data = expr_vals, FUN = mean)
      colnames(pb) <- c("sample", "TMEM164")
      pb$dataset <- "GSE310058"
      # Map sample → GOLD stage from metadata
      all_staging[["GSE310058"]] <- pb
    }
  }
  rm(seu); gc()
}

# ══════════════════════════════════════════════════════════════════════════════
# Combine & analyse
# ══════════════════════════════════════════════════════════════════════════════

staging_df <- rbindlist(all_staging, fill = TRUE)

if (nrow(staging_df) == 0) {
  message("⚠ No staging data available — ensure preprocessing is complete")
  quit(save = "no")
}

# Convert GOLD to ordered factor
staging_df$gold <- factor(staging_df$gold,
                          levels = c("Control", "1", "2", "3", "4",
                                     "1-2", "3-4"),
                          ordered = TRUE)

write.csv(staging_df, file.path(table_dir, "copd_tmem164_gold_staging.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Trend test: Jonckheere-Terpstra
# ══════════════════════════════════════════════════════════════════════════════

# Subset to individual GOLD stages only
staging_individual <- staging_df[gold %in% c("1", "2", "3", "4")]

if (nrow(staging_individual) > 0 && length(unique(staging_individual$gold)) >= 3) {
  if (requireNamespace("clinfun", quietly = TRUE)) {
    jt <- clinfun::jonckheere.test(
      staging_individual$TMEM164,
      as.numeric(staging_individual$gold),
      alternative = "increasing"
    )
    message("Jonckheere-Terpstra trend test:")
    message("  Statistic: ", round(jt$statistic, 2),
            "  P-value: ", signif(jt$p.value, 3))
    
    trend_result <- data.frame(
      test = "Jonckheere-Terpstra",
      statistic = jt$statistic,
      p_value = jt$p.value,
      direction = ifelse(jt$p.value < 0.05, "significant trend", "no trend")
    )
    write.csv(trend_result,
              file.path(table_dir, "copd_tmem164_gold_trend.csv"),
              row.names = FALSE)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# Visualisation
# ══════════════════════════════════════════════════════════════════════════════

# 1. Violin/box plot across GOLD stages
p1 <- ggplot(staging_df[gold %in% c("Control", "1", "2", "3", "4")],
             aes(x = gold, y = TMEM164, fill = gold)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.9) +
  scale_fill_manual(values = c(Control = "#4DBBD5",
                                `1` = "#FDD49E", `2` = "#FDBB84",
                                `3` = "#FC8D59", `4` = "#D7301F")) +
  labs(x = "GOLD Stage", y = paste0(TARGET_GENE, " expression"),
       title = paste0(TARGET_GENE, " across COPD GOLD stages")) +
  theme_nature() +
  theme(legend.position = "none")

# 2. Per-dataset faceted view
p2 <- ggplot(staging_df, aes(x = gold, y = TMEM164, fill = gold)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#4DBBD5",
                                `1` = "#FDD49E", `2` = "#FDBB84",
                                `3` = "#FC8D59", `4` = "#D7301F",
                                `1-2` = "#FDBB84", `3-4` = "#D7301F")) +
  labs(x = "GOLD Stage", y = paste0(TARGET_GENE, " expression"),
       title = "Per-dataset staging") +
  theme_nature(base_size = 6) +
  theme(legend.position = "none")

p_combined <- p1 / p2 + plot_layout(heights = c(1, 1))
save_nature_fig(p_combined, "copd_tmem164_gold_staging",
                width = 7.2, height = 6)

message("\n✔ COPD GOLD staging analysis complete")
save_session_info(report_dir)
