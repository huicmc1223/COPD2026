###############################################################################
# geo_helpers.R — GEO download / parsing helper functions
# Usage: source(here::here("scripts/utils/geo_helpers.R"))
###############################################################################

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(R.utils)
})

#' Download a GEO dataset (expression + metadata)
#'
#' @param geo_id GEO accession (e.g., "GSE57148").
#' @param dest_dir Destination directory for raw data.
#' @param type One of "microarray", "bulk_rnaseq", "scrnaseq", "spatial".
#' @return List with `eset` (ExpressionSet) and `supp_dir` (supplementary path).
download_geo <- function(geo_id, dest_dir, type = "microarray") {
  
  message("── Downloading ", geo_id, " ──")
  
  # Determine destination sub-directory
  sub_dir <- file.path(dest_dir, type, geo_id)
  dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
  
  result <- list(geo_id = geo_id, type = type, dir = sub_dir)
  
  # Download series matrix (metadata + expression for arrays)
  tryCatch({
    gse <- getGEO(geo_id, destdir = sub_dir, GSEMatrix = TRUE, getGPL = TRUE)
    if (is.list(gse) && length(gse) > 0) {
      result$eset <- gse
      # Save phenotype data
      for (i in seq_along(gse)) {
        pdata <- pData(gse[[i]])
        fwrite(pdata, file.path(sub_dir, paste0("phenodata_", i, ".csv")))
      }
    }
    message("  ✔ Series matrix downloaded")
  }, error = function(e) {
    message("  ⚠ Series matrix failed: ", e$message)
  })
  
  # Download supplementary files (raw counts, H5 etc.)
  tryCatch({
    supp <- getGEOSuppFiles(geo_id, baseDir = sub_dir, makeDirectory = FALSE)
    result$supp_dir <- sub_dir
    message("  ✔ Supplementary files downloaded")
  }, error = function(e) {
    message("  ⚠ Supplementary files failed: ", e$message)
  })
  
  return(result)
}

#' Download scRNA-seq data via wget (for large files from GEO FTP)
#'
#' @param geo_id GEO accession.
#' @param dest_dir Destination directory.
#' @param file_patterns Character vector of file patterns to download.
download_geo_ftp <- function(geo_id, dest_dir,
                             file_patterns = c("barcodes", "features", "matrix",
                                               "h5", "h5ad")) {
  sub_dir <- file.path(dest_dir, "scrnaseq", geo_id)
  dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Construct FTP URL
  prefix <- substring(geo_id, 1, nchar(geo_id) - 3)
  nnn <- paste0(prefix, "nnn")
  ftp_base <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", nnn, "/",
                     geo_id, "/suppl/")
  
  message("── FTP download ", geo_id, " from: ", ftp_base, " ──")
  
  # List and download files
  cmd <- paste0("wget -r -np -nd -A '",
                paste(paste0("*", file_patterns, "*"), collapse = ","),
                "' -P '", sub_dir, "' '", ftp_base, "'")
  
  system(cmd, intern = FALSE)
  
  # Decompress .gz files
  gz_files <- list.files(sub_dir, pattern = "\\.gz$", full.names = TRUE)
  for (f in gz_files) {
    tryCatch({
      R.utils::gunzip(f, remove = TRUE, overwrite = TRUE)
      message("  ✔ Decompressed: ", basename(f))
    }, error = function(e) {
      message("  ⚠ Decompress failed: ", basename(f))
    })
  }
  
  return(sub_dir)
}

#' Extract and standardise phenotype metadata from a GEO ExpressionSet
#'
#' @param eset An ExpressionSet from getGEO().
#' @param geo_id GEO accession (for labelling).
#' @param condition_col Name of the column containing condition info.
#' @return data.frame with columns: sample_id, condition, dataset, platform.
standardise_pheno <- function(eset, geo_id, condition_col = NULL) {
  pd <- pData(eset)
  result <- data.frame(
    sample_id = rownames(pd),
    geo_id    = geo_id,
    platform  = as.character(pd$platform_id[1]),
    stringsAsFactors = FALSE
  )
  
  # Try to auto-detect condition column
  if (is.null(condition_col)) {
    candidates <- grep("disease|condition|group|status|source|characteristic",
                       colnames(pd), ignore.case = TRUE, value = TRUE)
    if (length(candidates) > 0) {
      condition_col <- candidates[1]
      message("  Auto-detected condition column: ", condition_col)
    }
  }
  
  if (!is.null(condition_col) && condition_col %in% colnames(pd)) {
    result$condition <- as.character(pd[[condition_col]])
  } else {
    result$condition <- NA_character_
  }
  
  return(result)
}

#' Read 10x Genomics format (barcodes / features / matrix) into Seurat
#'
#' @param data_dir Directory containing the three files.
#' @param project_name Name for the Seurat object.
#' @return Seurat object.
read_10x_to_seurat <- function(data_dir, project_name = "scRNA") {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not installed")
  
  counts <- Seurat::Read10X(data_dir)
  seu <- Seurat::CreateSeuratObject(counts = counts,
                                    project = project_name,
                                    min.cells = 3,
                                    min.features = 200)
  message("  ✔ Seurat object: ", ncol(seu), " cells × ", nrow(seu), " genes")
  return(seu)
}

message("✔ geo_helpers.R loaded")
