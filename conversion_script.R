#### CONVERTING FROM SCANPY TO SEURAT ----

library(Seurat)
library(Matrix)
library(data.table)

SCALE_FACTOR <- 1000000
input_dir    <- "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/anndata_to_seurat_break/"
output_rds   <- "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/fishing_2lines_dualrun_5481genes_20pcs_epithelial_mesenchymal_2026_05_06.rds"

build_rds_lean <- function(input_dir, output_rds) {
  
  cat("Reading metadata with data.table...\n")
  obs_dt <- fread(file.path(input_dir, "obs.csv"), header = TRUE)
  obs <- as.data.frame(obs_dt[, -1, with = FALSE])
  rownames(obs) <- obs_dt[[1]]
  rm(obs_dt)
  cat("obs done\n")
  
  var_dt <- fread(file.path(input_dir, "var.csv"), header = TRUE)
  if (ncol(var_dt) == 1) {
    var <- data.frame(row.names = var_dt[[1]])
  } else {
    var <- as.data.frame(var_dt[, -1, with = FALSE])
    rownames(var) <- var_dt[[1]]
  }
  rm(var_dt)
  cat("var done\n")
  
  gene_names <- rownames(var)
  cell_names <- rownames(obs)
  cat("  Detected", length(cell_names), "cells and", length(gene_names), "genes.\n")
  
  cat("Loading raw counts...\n")
  raw_mat <- readMM(file.path(input_dir, "layers", "counts.mtx"))
  if (inherits(raw_mat, "dgTMatrix")) raw_mat <- as(raw_mat, "CsparseMatrix")
  cat("raw_mat dims:", nrow(raw_mat), "x", ncol(raw_mat), "\n")
  cat("expected:", length(gene_names), "x", length(cell_names), "\n")
  
  if (ncol(raw_mat) != length(cell_names) || nrow(raw_mat) != length(gene_names)) {
    stop(paste0("Dimension mismatch! Matrix: ", nrow(raw_mat), "x", ncol(raw_mat), 
                ". Metadata: ", length(gene_names), "x", length(cell_names)))
  }
  
  dimnames(raw_mat) <- list(gene_names, cell_names)
  cat("dimnames set\n")
  
  cat("Creating Seurat Object...\n")
  seurat_obj <- CreateSeuratObject(counts = raw_mat, meta.data = obs, assay = "RNA")
  cat("Seurat object created\n")
  rm(raw_mat, obs); gc()
  
  cat("Loading log1p data...\n")
  log1p_mat <- readMM(file.path(input_dir, "layers", "log1p.mtx"))
  log1p_mat <- log1p_mat / SCALE_FACTOR
  if (inherits(log1p_mat, "dgTMatrix")) log1p_mat <- as(log1p_mat, "CsparseMatrix")
  dimnames(log1p_mat) <- list(gene_names, cell_names)
  cat("log1p done\n")
  
  seurat_obj[["RNA"]]$data <- log1p_mat
  cat("log1p assigned\n")
  rm(log1p_mat); gc()
  
  seurat_obj@misc[["var"]] <- var
  cat("var stored in misc\n")
  rm(var); gc()
  
  cat("Loading embeddings...\n")
  embed_files <- list.files(file.path(input_dir, "embeddings"), pattern = "\\.csv$", full.names = TRUE)
  for (f in embed_files) {
    embed_name <- gsub("\\.csv$", "", basename(f))
    cat("  Adding reduction:", embed_name, "\n")
    emb_dt <- fread(f, header = TRUE)
    emb_mat <- as.matrix(emb_dt[, -1, with = FALSE])
    rownames(emb_mat) <- emb_dt[[1]]
    rm(emb_dt)
    cat("  emb_mat dims:", nrow(emb_mat), "x", ncol(emb_mat), "\n")
    
    # Skip if cell count doesn't match
    if (nrow(emb_mat) != length(cell_names)) {
      cat("  SKIPPING", embed_name, "- row mismatch:", nrow(emb_mat), "vs", length(cell_names), "\n")
      next
    }
    
    key <- tolower(gsub("^X_", "", embed_name))
    seurat_obj[[toupper(key)]] <- CreateDimReducObject(
      embeddings = emb_mat, 
      key = paste0(key, "_"), 
      assay = "RNA"
    )
    cat("  Reduction added:", embed_name, "\n")
    rm(emb_mat); gc()
  }
  
  cat("Saving RDS...\n")
  saveRDS(seurat_obj, file = output_rds, compress = FALSE) 
  cat("Done!\n")
}

build_rds_lean(input_dir, output_rds)
