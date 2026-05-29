#### LOAD PACKAGES AND SET ENVIRONMENT ----
## Packages
library(tidyverse)
library(stringr)
library(rlist)
library(Seurat)
library(SeuratObject)
library(data.table)
library(Matrix)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(arrow)
library(SeuratDisk)
library(data.table)
library(doParallel)                                                                                                   
library(RANN)                                             
library(Matrix)                                                                                                       
library(ggrepel)             

## Set seed & working directory
set.seed(0712)
work_dir <- setwd("/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/")

options(scipen=999)

## "not in" operand function
'%!in%' <- function(x,y)!('%in%'(x,y))

gene_overlaps <- c("CXCR4", "KCNK3", "XBP1", "AGR3", "ATF4", "CCT8", "CD28", "CD36", 
                   "CD68", "COL4A1", "COL4A3", "CSPG4", "CXCL6", "FAP", "FGF2", "FGFR4", 
                   "GJA5", "HRAS", "ID4", "IFT57", "LAMP3", "LEF1", "MSLN", "NPNT", 
                   "PECAM1", "PPARG", "RSPO2", "RSPO3", "RTKN2", "SECISBP2L", "SELL", 
                   "SLC25A4", "SOX2", "TRPC6", "TSPAN4", "WNT2", "ZEB1", "SNAI2", "AQP5", 
                   "CCL18", "CCL22", "CD247", "CD44", "CXCL12", "DNAH10", "FKBP11", "GCLM", 
                   "HDAC9", "HES1", "HEY1", "IGF1", "INMT", "KDR", "MEG3", "MMP12", "MS4A1", 
                   "MYC", "NKX2-1", "PDGFRB", "POSTN", "STAT1", "WNT3A", "ADGRL4", "SLC1A3", 
                   "AXL", "BCL2L11", "BMPR2", "CCN2", "CD1E", "CHIT1", "CP", "CXCL13", "EPAS1", 
                   "EREG", "FOXJ1", "GATA6", "HMOX1", "IFI16", "IFIT1", "ITGA3", "KLRB1", 
                   "MAD1L1", "MFAP5", "MUC5AC", "PDGFRA", "PDIA4", "PLIN2", "PLN", "RUNX1", 
                   "SHH", "SRSF1", "TREM2", "CD163", "CD8A", "EPCAM", "ETV5", "FCGR3A", "FOXF1", 
                   "FSTL1", "HIF1A", "HOPX", "ITGB6", "PI16", "PLVAP", "SOX9", "TCF21", "TGFBR2", 
                   "TMEM100", "ARF6", "BANK1", "CD3G", "CD47", "CH25H", "FGF10", "FGF7", "GZMB", 
                   "HSPG2", "IDH1", "KLRG1", "MARCO", "MZB1", "NUCB2", "NUTF2", "SEMA3C", "STAT6", 
                   "TP73", "CDK1", "CENPF", "EGFR", "FGFBP2", "ICAM1", "SEMA3A", "UPK3B", "WWTR1", 
                   "CD1C", "ACAA2", "AGER", "ATF3", "BCL2L1", "EDNRB", "HMGA1", "PTPN6", "SPATS2L", 
                   "ADAMTS1", "ITGAV", "LTA4H", "NRP1", "TCL1A", "TGFB2", "GBP1", "NFKB1", "LAMA3", 
                   "PTTG1", "SHANK3", "SNAI1", "ZBTB16", "BCL2", "IRF1", "VEGFA", "ATG7", "CDH1", 
                   "FCER1A", "TMEM33", "APLNR", "CA4", "CD79A", "CLDN5", "TNFRSF13C", "TP53", "BMP4", 
                   "DNAH12", "ELN", "CD3E", "CD2", "CFTR", "GDF15", "LGR5", "NTN4", "CFH", "FILIP1L", 
                   "MRC1", "PROX1", "CD34", "CES1", "FCGBP", "MUC5B", "SFRP1", "SMAD4", "TGFB1", "CD4", 
                   "SLIT3", "CD151", "CD19", "CD79B", "COL2A1", "ITGA6", "KLRC1", "MKI67", "SVEP1", 
                   "TBX4", "UNC119", "ASPN", "BAX", "FCN1", "GSR", "IFITM1", "IL4R", "KIT", "LTBP2", 
                   "TGFB3", "TMPRSS2", "WNT5A", "CTHRC1", "PROM1", "TOLLIP", "TOP2A", "NOTCH3", "FASN", 
                   "MLPH", "APLN", "ELANE", "PIM2", "LCK", "SMPD3", "TSPAN8", "CD14", "CD27", "LILRA4", 
                   "SFRP4", "SPRY2", "TRIM47", "TP63", "CD8B", "GZMA", "INHBA", "WT1", "DERL3", "PTGS1", "IRF7")

#### IMPORTING RAW TRANSCRIPT FILES FOR TRANSCRIPT-SPECIFIC CORRELATIONS ----

## Pulling the count and metadata files from Xenium output
transcript_filepaths <- list.files(pattern = "transcripts.parquet", recursive = TRUE)

## Processing files of choice and naming them accordingly
transcript_filepaths_subset <- transcript_filepaths[1:4]
transcript_files <- lapply(transcript_filepaths_subset, function(XX)
  read_parquet(XX, as_data_frame = TRUE))

ids <- c("prime_solo", "prime_v1", "v1_prime", "v1_solo")
names(transcript_files) <- ids

## For loop to aggregate transcripts per gene
transcript_count_files <- list()
for (i in seq_along(ids)) {
  print(i)
  
  df <- transcript_files[[i]]
  
  # Remove negative controls
  df <- df[!grepl("NegControl|Codeword|Intergenic", df$feature_name, ignore.case = TRUE), ]
  
  gene_counts_unfiltered <- df %>%
    group_by(feature_name) %>%
    summarise(
      gene_count_unfiltered = n(),
      .groups = "drop"
    )
  
  # Filtered counts (qv >= 20)
  gene_counts_filtered <- df %>%
    filter(qv >= 20) %>%
    group_by(feature_name) %>%
    summarise(
      gene_count_filtered = n(),
      .groups = "drop"
    )
  
  # Combine both totals
  gene_counts <- full_join(gene_counts_unfiltered, gene_counts_filtered, by = "feature_name") %>%
    replace_na(list(
      gene_count_unfiltered = 0,
      gene_count_filtered = 0
    ))
  
  # Set row names and remove feature_name column
  gene_counts <- as.data.frame(gene_counts)
  rownames(gene_counts) <- gene_counts$feature_name
  gene_counts$feature_name <- NULL
  
  # Save back to list
  transcript_count_files[[i]] <- gene_counts
}
names(transcript_count_files) <- ids

write_csv(transcript_count_files$prime_solo, 
          "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/prime_solo_transcripts_per_gene_qv20_2025_11_19.csv")
write_csv(transcript_count_files$prime_v1, 
          "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/prime_v1_transcripts_per_gene_qv20_2025_11_19.csv")
write_csv(transcript_count_files$v1_prime, 
          "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/v1_prime_transcripts_per_gene_qv20_2025_11_19.csv")
write_csv(transcript_count_files$v1_solo, 
          "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/v1_solo_transcripts_per_gene_qv20_2025_11_19.csv")

#### SEURAT OBJECT CREATION (DUAL + ALIGNED ONLY) ----

## Pulling the count and metadata files from Xenium output
transcript_filepaths <- list.files(pattern = "transcripts.parquet", 
                                   recursive = TRUE)
meta_filepaths <- list.files(pattern = "cells.csv.gz",
                             recursive = TRUE)

## Only v1_prime and segmentation
subset_transcripts <- transcript_filepaths[c(3, 5)]
subset_meta <- meta_filepaths[c(3, 5)]

## Parsing through all files and fixing the formats for both
transcript_files <- lapply(subset_transcripts, function(XX) {
  df <- arrow::read_parquet(XX, as_data_frame = TRUE)
  as.data.frame(df)
})
meta_files <- lapply(subset_meta, function(XX)
  read.delim(XX, sep = ","))

## ID list 
slide_type <- c("v1_prime", "prime_with_v1_seg")
names(transcript_files) <- slide_type
names(meta_files) <- slide_type

# Function to match v1_prime and prime + v1 segmentation
int_to_str26 <- function(x, len = 6) {
  x <- as.integer(x)
  n <- length(x)
  if (n == 0) return(character(0))
  
  m <- matrix(0L, nrow = n, ncol = len)
  for (i in seq_len(len)) {
    m[, len - i + 1] <- x %% 26
    x <- x %/% 26
  }
  apply(m, 1, function(r) paste0(letters[r + 1], collapse = ""))
}

# Create matched cell IDs for both metadata sets
meta_files$v1_prime$cell_index <- seq_len(nrow(meta_files$v1_prime))
meta_files$prime_with_v1_seg$cell_index <- seq_len(nrow(meta_files$prime_with_v1_seg))

meta_files$v1_prime$matched_cell_ids <- int_to_str26(meta_files$v1_prime$cell_index - 1, len = 6)
meta_files$prime_with_v1_seg$matched_cell_ids <- int_to_str26(meta_files$prime_with_v1_seg$cell_index - 1, len = 6)

# Processing transcript and metadata files
processed_transcripts <- vector("list", length(slide_type))
names(processed_transcripts) <- slide_type

for (i in seq_along(slide_type)) {
  message("Processing: ", slide_type[i])
  
  # Load transcript data
  df <- transcript_files[[i]]
  if (!is.data.frame(df)) df <- as.data.frame(df)
  setDT(df)
  
  # Basic filtering
  df <- df[overlaps_nucleus == "1"]
  df <- df[qv > 20]
  df <- df[!grepl("NegControl|Codeword|Intergenic", feature_name, ignore.case = TRUE)]
  
  # Append slide suffix to feature names
  df[, feature_name := paste0(feature_name, "_", slide_type[i])]
  
  # Select correct metadata object (ensure correct source)
  meta <- meta_files[[slide_type[i]]]
  if (!"matched_cell_ids" %in% colnames(meta)) {
    stop("Metadata for ", slide_type[i], " does not have matched_cell_ids.")
  }
  
  # Merge matched cell strings into transcript table
  meta_subset <- meta[, c("cell_id", "matched_cell_ids")]
  df <- merge(df, meta_subset, by = "cell_id", all.x = TRUE)
  
  # Check merge quality
  if (any(is.na(df$matched_cell_ids))) {
    warning(sum(is.na(df$matched_cell_ids)), " transcripts could not be matched for ", slide_type[i])
  }
  
  # Count feature occurrences per matched cell
  counts <- df[, .N, by = .(feature_name, matched_cell_ids)]
  
  # Make factor indices
  feature_levels <- sort(unique(counts$feature_name))
  cell_levels <- sort(unique(counts$matched_cell_ids))
  
  counts[, feature_index := match(feature_name, feature_levels)]
  counts[, cell_index := match(matched_cell_ids, cell_levels)]
  
  # Build sparse matrix (rows = features, cols = matched cell IDs)
  mat <- sparseMatrix(
    i = counts$feature_index,
    j = counts$cell_index,
    x = counts$N,
    dims = c(length(feature_levels), length(cell_levels)),
    dimnames = list(feature_levels, cell_levels)
  )
  
  processed_transcripts[[i]] <- mat
  
  # Update metadata 
  rownames(meta) <- meta$cell_id
  meta$cell_index <- seq_len(nrow(meta))
  meta_files[[slide_type[i]]] <- meta
}

## Bind transcript files
common_cells <- intersect(colnames(processed_transcripts$v1_prime), colnames(processed_transcripts$prime_with_v1_seg))

# Subset to shared columns
v1_prime_mat <- processed_transcripts$v1_prime[, common_cells, drop = FALSE]
prime_w_v1_seg_mat <- processed_transcripts$prime_with_v1_seg[, common_cells, drop = FALSE]

# Combine (stack rows = rbind, stack cells = cbind)
mat_combined <- rbind(v1_prime_mat, prime_w_v1_seg_mat)

### Seurat object 
seurat_obj <- CreateSeuratObject(counts = mat_combined, assay = "RNA")

# Adjusting metadata to append on
rownames(meta_files$v1_prime) <- meta_files$v1_prime$matched_cell_ids
rownames(meta_files$prime_with_v1_seg) <- meta_files$prime_with_v1_seg$matched_cell_ids

test <- meta_files$v1_prime %>%
  select(cell_id, x_centroid, y_centroid, 
         transcript_counts, cell_area, nucleus_area, 
         nucleus_count, matched_cell_ids) %>%
  relocate(matched_cell_ids, .after = cell_id) %>%
  rename(cell_id_original_v1_prime = cell_id,
         x_centroid_v1_prime = x_centroid, 
         y_centroid_v1_prime = y_centroid,
         transcript_counts_v1_prime = transcript_counts)

test1 <- meta_files$prime_with_v1_seg %>%
  select(cell_id, x_centroid, y_centroid, 
         transcript_counts, matched_cell_ids) %>%
  relocate(matched_cell_ids, .after = cell_id) %>%
  rename(cell_id_original_prime_v1_seg = cell_id,
         x_centroid_prime_v1_seg = x_centroid, 
         y_centroid_prime_v1_seg = y_centroid,
         transcript_counts_prime_v1_seg = transcript_counts)

metadata <- left_join(test, test1, by = "matched_cell_ids")
rownames(metadata) <- metadata$matched_cell_ids

seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)

### Prep for RAPIDS
seurat_obj_slot <- NormalizeData(seurat_obj)

rna <- seurat_obj_slot[["RNA"]]
rna_assay <- CreateAssayObject(counts = GetAssayData(rna, layer = "counts"))
seurat_obj_slot[["RNA"]] <- rna_assay

# pip install anndata in terminal

convertFormat(seurat_obj_slot, from="seurat", to="anndata",
              outFile='v1_prime_unfiltered_joined_2025_11_05.h5ad')

### QC and filtering

smoothScatter(seurat_obj_slot@meta.data$nCount_RNA,
              seurat_obj_slot@meta.data$nFeature_RNA,
              cex = 0.5, pch = 16)

smoothScatter(seurat_obj_slot@meta.data$cell_area,
              seurat_obj_slot@meta.data$nCount_RNA,
              cex = 0.5, pch = 16)

smoothScatter(seurat_obj_slot@meta.data$cell_area,
              seurat_obj_slot@meta.data$nFeature_RNA,
              cex = 0.5, pch = 16)

smoothScatter(seurat_obj_slot@meta.data$nucleus_area,
              seurat_obj_slot@meta.data$nCount_RNA,
              cex = 0.5, pch = 16)

smoothScatter(seurat_obj_slot@meta.data$nucleus_area,
              seurat_obj_slot@meta.data$nFeature_RNA,
              cex = 0.5, pch = 16)

seurat_obj_slot_filtered <- subset(seurat_obj_slot, subset = nCount_RNA >= 50
                     & nFeature_RNA >= 5
                     & cell_area >= 5
                     & cell_area <= 140
                     & nucleus_area >= 3)

#seurat_obj_slot[["RNA"]] <- as(seurat_obj_slot[["RNA"]], "Assay")

saveRDS(seurat_obj_slot_filtered, "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/v1_prime_filtered_joined_2025_11_25.rds")

### Adding raw counts back into split objects post-RAPIDS 

## Convert RAPIDS objects into Seurat
dual_obj_combined <- schard::h5ad2seurat("/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/v1_prime_comparison_allgenes_20pcs_filtered_split_2025_12_05.h5ad")
saveRDS(dual_obj_combined, "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/v1_prime_comparison_allgenes_20pcs_filtered_split_postrapids_2025_12_07.rds")

dual_obj_v1 <- schard::h5ad2seurat("/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/v1_prime_comparison_v1genes_20pcs_filtered_split_2025_12_05.h5ad")
saveRDS(dual_obj_v1, "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/v1_prime_comparison_v1genes_20pcs_filtered_split_postrapids_2025_12_07.rds")

dual_obj_prime <- schard::h5ad2seurat("/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/v1_prime_primegenes_combined_50pcs_filtered_split_2025_12_05.h5ad")
saveRDS(dual_obj_prime, "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/v1_prime_comparison_primegenes_50pcs_filtered_split_postrapids_2025_12_07.rds")

## Import counts file

raw_counts_matrix <- fread("/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/test_counts.csv",
                           sep = ",", header = TRUE)

raw_counts_matrix_adj <- raw_counts_matrix %>%
  column_to_rownames(var = "V1")

raw_counts_tibble <- read_delim(
  file = "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/test_counts.csv",
  delim = ",",
  col_names = TRUE,
  col_select = everything(), 
  id = "Gene_ID" 
)

raw_counts_df <- as.data.frame(raw_counts_tibble)

your_seurat_object[["RNA"]] <- SetAssayData(
  object = your_seurat_object[["RNA"]], 
  slot = "counts", 
  new.data = raw_counts_sparse
)

#### SEURAT OBJECT CREATION (V1 SOLO + V1_PRIME ONLY) ----

## Pulling the count and metadata files from Xenium output
transcript_filepaths <- list.files(pattern = "transcripts.parquet", 
                                   recursive = TRUE)
meta_filepaths <- list.files(pattern = "cells.csv.gz",
                             recursive = TRUE)

## Only v1_prime and v1_solo
subset_transcripts <- transcript_filepaths[c(3, 4)]
subset_meta <- meta_filepaths[c(3, 4)]

## Parsing through all files and fixing the formats for both
transcript_files <- lapply(subset_transcripts, function(XX) {
  df <- arrow::read_parquet(XX, as_data_frame = TRUE)
  as.data.frame(df)
})
meta_files <- lapply(subset_meta, function(XX)
  read.delim(XX, sep = ","))

## ID list 
slide_type <- c("v1_prime", "v1_solo")
names(transcript_files) <- slide_type
names(meta_files) <- slide_type

# Processing transcript and metadata files
processed_transcripts <- vector("list", length(slide_type))
names(processed_transcripts) <- slide_type

for (i in seq_along(slide_type)) {
  message("Processing: ", slide_type[i])
  
  # Load transcript data
  df <- transcript_files[[i]]
  if (!is.data.frame(df)) df <- as.data.frame(df)
  setDT(df)
  
  # Basic filtering
  df <- df[overlaps_nucleus == "1"]
  df <- df[qv > 20]
  df <- df[!grepl("NegControl|Codeword|Intergenic", feature_name, ignore.case = TRUE)]
  
  # Append slide suffix to feature names
  df[, cell_id := paste0(cell_id, "_", slide_type[i])]

  # Count feature occurrences per cell id
  counts <- df[, .N, by = .(feature_name, cell_id)]
  
  # Make factor indices
  feature_levels <- sort(unique(counts$feature_name))
  cell_levels <- sort(unique(counts$cell_id))
  
  counts[, feature_index := match(feature_name, feature_levels)]
  counts[, cell_index := match(cell_id, cell_levels)]
  
  # Build sparse matrix (rows = features, cols = matched cell IDs)
  mat <- sparseMatrix(
    i = counts$feature_index,
    j = counts$cell_index,
    x = counts$N,
    dims = c(length(feature_levels), length(cell_levels)),
    dimnames = list(feature_levels, cell_levels)
  )
  
  processed_transcripts[[i]] <- mat
  
  # Load and update the corresponding metadata file
  meta <- meta_files[[slide_type[i]]]
  setDT(meta)
  meta[, cell_id := paste0(cell_id, "_", slide_type[i])]
  rownames(meta) <- meta$cell_id 
  meta_files[[slide_type[i]]] <- meta
  
}

### Seurat objects into one
obj_list <- list()
for (i in seq_along(slide_type)) {
  
  seurat_obj <- CreateSeuratObject(counts = processed_transcripts[[i]], assay = "RNA")
  seurat_obj <- AddMetaData(seurat_obj, metadata = meta_files[[i]])
  
  # Add file name to metadata
  seurat_obj$slide_type <- slide_type[[i]]
  
  # Add each Seurat object into a list of objects
  obj_list[[i]] <- seurat_obj
  
}

## Merging all Seurat objects in obj_list together into one 
merged_obj_unfiltered_v1 <- merge(x = obj_list[[1]], y = obj_list[2:length(obj_list)])

### Adding sample-level information

## Importing cell stats for each slide
v1_solo_cell_stats <- list.files(path = "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/cell_stats", 
                                 pattern = "V1_solo_cells", 
                                 full.names = TRUE, 
                                 recursive = FALSE)
v1_prime_cell_stats <- list.files(path = "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/cell_stats", 
                                  pattern = "V1_prime_cells", 
                                  full.names = TRUE, 
                                  recursive = FALSE)

## Parsing through files and fixing the format
v1_solo_cell_stats_csv <- lapply(v1_solo_cell_stats, function(XX)
  read.delim(XX, sep = ",", skip = 2))
v1_prime_cell_stats_csv <- lapply(v1_prime_cell_stats, function(XX)
  read.delim(XX, sep = ",", skip = 2))

sample_ids <- c("PDL018D", "PDL026A", "PDL029A", "PDL031A",
                "PDL032D", "PDL033T", "PDL034D", "PDL040A", 
                "PDL042D", "PDL047T", "PDL055T", "PDL061",
                "PDL072", "PDL077", "PDL085A", "PDL095D", "PDL097")

names(v1_solo_cell_stats_csv) <- sample_ids
names(v1_prime_cell_stats_csv) <- sample_ids

## Removing slide_type from cell_id
merged_obj_unfiltered_v1@meta.data <- merged_obj_unfiltered_v1@meta.data %>%
  mutate(cell_id = str_remove(cell_id, "_v1_prime")) %>%
  mutate(cell_id = str_remove(cell_id, "_v1_solo"))
  
## Making a new column named "Sample" that assigns cells to their individual samples 
merged_obj_unfiltered_v1@meta.data <- merged_obj_unfiltered_v1@meta.data %>%
  mutate(sample = case_when(cell_id %in% v1_solo_cell_stats_csv[["PDL018D"]]$Cell.ID ~ "PDL018D",
                            cell_id %in% v1_solo_cell_stats_csv[["PDL026A"]]$Cell.ID ~ "PDL026A",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL029A"]]$Cell.ID ~ "PDL029A",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL031A"]]$Cell.ID ~ "PDL031A",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL032D"]]$Cell.ID ~ "PDL032D",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL033T"]]$Cell.ID ~ "PDL033T",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL034D"]]$Cell.ID ~ "PDL034D",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL040A"]]$Cell.ID ~ "PDL040A",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL042D"]]$Cell.ID ~ "PDL042D",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL047T"]]$Cell.ID ~ "PDL047T",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL055T"]]$Cell.ID ~ "PDL055T",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL061"]]$Cell.ID ~ "PDL061",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL072"]]$Cell.ID ~ "PDL072",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL077"]]$Cell.ID ~ "PDL077",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL085A"]]$Cell.ID ~ "PDL085A",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL095D"]]$Cell.ID ~ "PDL095D",
                             cell_id %in% v1_solo_cell_stats_csv[["PDL097"]]$Cell.ID ~ "PDL097",
                            
                             cell_id %in% v1_prime_cell_stats_csv[["PDL018D"]]$Cell.ID ~ "PDL018D",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL026A"]]$Cell.ID ~ "PDL026A",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL029A"]]$Cell.ID ~ "PDL029A",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL031A"]]$Cell.ID ~ "PDL031A",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL032D"]]$Cell.ID ~ "PDL032D",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL033T"]]$Cell.ID ~ "PDL033T",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL034D"]]$Cell.ID ~ "PDL034D",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL040A"]]$Cell.ID ~ "PDL040A",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL042D"]]$Cell.ID ~ "PDL042D",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL047T"]]$Cell.ID ~ "PDL047T",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL055T"]]$Cell.ID ~ "PDL055T",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL061"]]$Cell.ID ~ "PDL061",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL072"]]$Cell.ID ~ "PDL072",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL077"]]$Cell.ID ~ "PDL077",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL085A"]]$Cell.ID ~ "PDL085A",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL095D"]]$Cell.ID ~ "PDL095D",
                             cell_id %in% v1_prime_cell_stats_csv[["PDL097"]]$Cell.ID ~ "PDL097",
                            TRUE ~ "dropped"))

saveRDS(merged_obj_unfiltered_v1, 
        "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/v1_slides_2025_11_20.rds")

#### SEURAT OBJECT CREATION (PRIME SOLO + PRIME_V1 ONLY) ----

## Pulling the count and metadata files from Xenium output
transcript_filepaths <- list.files(pattern = "transcripts.parquet", 
                                   recursive = TRUE)
meta_filepaths <- list.files(pattern = "cells.csv.gz",
                             recursive = TRUE)

## Only prime_solo and prime_v1
subset_transcripts <- transcript_filepaths[c(1, 2)]
subset_meta <- meta_filepaths[c(1, 2)]

## Parsing through all files and fixing the formats for both
transcript_files <- lapply(subset_transcripts, function(XX) {
  df <- arrow::read_parquet(XX, as_data_frame = TRUE)
  as.data.frame(df)
})
meta_files <- lapply(subset_meta, function(XX)
  read.delim(XX, sep = ","))

## ID list 
slide_type <- c("prime_solo", "prime_v1")
names(transcript_files) <- slide_type
names(meta_files) <- slide_type

# Processing transcript and metadata files
processed_transcripts <- vector("list", length(slide_type))
names(processed_transcripts) <- slide_type

for (i in seq_along(slide_type)) {
  message("Processing: ", slide_type[i])
  
  # Load transcript data
  df <- transcript_files[[i]]
  if (!is.data.frame(df)) df <- as.data.frame(df)
  setDT(df)
  
  # Basic filtering
  df <- df[overlaps_nucleus == "1"]
  df <- df[qv > 20]
  df <- df[!grepl("NegControl|Codeword|Intergenic", feature_name, ignore.case = TRUE)]
  
  # Append slide suffix to feature names
  df[, cell_id := paste0(cell_id, "_", slide_type[i])]
  
  # Count feature occurrences per cell id
  counts <- df[, .N, by = .(feature_name, cell_id)]
  
  # Make factor indices
  feature_levels <- sort(unique(counts$feature_name))
  cell_levels <- sort(unique(counts$cell_id))
  
  counts[, feature_index := match(feature_name, feature_levels)]
  counts[, cell_index := match(cell_id, cell_levels)]
  
  # Build sparse matrix (rows = features, cols = matched cell IDs)
  mat <- sparseMatrix(
    i = counts$feature_index,
    j = counts$cell_index,
    x = counts$N,
    dims = c(length(feature_levels), length(cell_levels)),
    dimnames = list(feature_levels, cell_levels)
  )
  
  processed_transcripts[[i]] <- mat
  
  # Load and update the corresponding metadata file
  meta <- meta_files[[slide_type[i]]]
  setDT(meta)
  meta[, cell_id := paste0(cell_id, "_", slide_type[i])]
  rownames(meta) <- meta$cell_id 
  meta_files[[slide_type[i]]] <- meta
  
}

### Seurat objects into one
obj_list <- list()
for (i in seq_along(slide_type)) {
  
  seurat_obj <- CreateSeuratObject(counts = processed_transcripts[[i]], assay = "RNA")
  seurat_obj <- AddMetaData(seurat_obj, metadata = meta_files[[i]])
  
  # Add file name to metadata
  seurat_obj$slide_type <- slide_type[[i]]
  
  # Add each Seurat object into a list of objects
  obj_list[[i]] <- seurat_obj
  
}

## Merging all Seurat objects in obj_list together into one 
merged_obj_unfiltered_prime <- merge(x = obj_list[[1]], y = obj_list[2:length(obj_list)])

### Adding sample-level information

## Importing cell stats for each slide
prime_solo_cell_stats <- list.files(path = "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/cell_stats", 
                                 pattern = "prime_solo_cells", 
                                 full.names = TRUE, 
                                 recursive = FALSE)
prime_v1_cell_stats <- list.files(path = "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/cell_stats", 
                                  pattern = "prime_V1_cells", 
                                  full.names = TRUE, 
                                  recursive = FALSE)

## Parsing through files and fixing the format
prime_solo_cell_stats_csv <- lapply(prime_solo_cell_stats, function(XX)
  read.delim(XX, sep = ",", skip = 2))
prime_v1_cell_stats_csv <- lapply(prime_v1_cell_stats, function(XX)
  read.delim(XX, sep = ",", skip = 2))

sample_ids <- c("PDL018D", "PDL026A", "PDL029A", "PDL031A",
                "PDL032D", "PDL033T", "PDL034D", "PDL040A", 
                "PDL042D", "PDL047T", "PDL055T", "PDL061",
                "PDL072", "PDL077", "PDL085A", "PDL095D", "PDL097")

names(prime_solo_cell_stats_csv) <- sample_ids
names(prime_v1_cell_stats_csv) <- sample_ids

## Removing slide_type from cell_id
merged_obj_unfiltered_prime@meta.data <- merged_obj_unfiltered_prime@meta.data %>%
  mutate(cell_id = str_remove(cell_id, "_prime_v1")) %>%
  mutate(cell_id = str_remove(cell_id, "_prime_solo"))

## Making a new column named "Sample" that assigns cells to their individual samples 
merged_obj_unfiltered_prime@meta.data <- merged_obj_unfiltered_prime@meta.data %>%
  mutate(sample = case_when(cell_id %in% prime_solo_cell_stats_csv[["PDL018D"]]$Cell.ID ~ "PDL018D",
                            cell_id %in% prime_solo_cell_stats_csv[["PDL026A"]]$Cell.ID ~ "PDL026A",
                            cell_id %in%  prime_solo_cell_stats_csv[["PDL029A"]]$Cell.ID ~ "PDL029A",
                            cell_id %in%  prime_solo_cell_stats_csv[["PDL031A"]]$Cell.ID ~ "PDL031A",
                            cell_id %in%  prime_solo_cell_stats_csv[["PDL032D"]]$Cell.ID ~ "PDL032D",
                            cell_id %in%  prime_solo_cell_stats_csv[["PDL033T"]]$Cell.ID ~ "PDL033T",
                            cell_id %in%  prime_solo_cell_stats_csv[["PDL034D"]]$Cell.ID ~ "PDL034D",
                            cell_id %in%  prime_solo_cell_stats_csv[["PDL040A"]]$Cell.ID ~ "PDL040A",
                            cell_id %in%  prime_solo_cell_stats_csv[["PDL042D"]]$Cell.ID ~ "PDL042D",
                            cell_id %in%  prime_solo_cell_stats_csv[["PDL047T"]]$Cell.ID ~ "PDL047T",
                            cell_id %in% prime_solo_cell_stats_csv[["PDL055T"]]$Cell.ID ~ "PDL055T",
                            cell_id %in% prime_solo_cell_stats_csv[["PDL061"]]$Cell.ID ~ "PDL061",
                            cell_id %in% prime_solo_cell_stats_csv[["PDL072"]]$Cell.ID ~ "PDL072",
                            cell_id %in% prime_solo_cell_stats_csv[["PDL077"]]$Cell.ID ~ "PDL077",
                            cell_id %in% prime_solo_cell_stats_csv[["PDL085A"]]$Cell.ID ~ "PDL085A",
                            cell_id %in% prime_solo_cell_stats_csv[["PDL095D"]]$Cell.ID ~ "PDL095D",
                            cell_id %in% prime_solo_cell_stats_csv[["PDL097"]]$Cell.ID ~ "PDL097",
                            
                            cell_id %in% prime_v1_cell_stats_csv[["PDL018D"]]$Cell.ID ~ "PDL018D",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL026A"]]$Cell.ID ~ "PDL026A",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL029A"]]$Cell.ID ~ "PDL029A",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL031A"]]$Cell.ID ~ "PDL031A",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL032D"]]$Cell.ID ~ "PDL032D",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL033T"]]$Cell.ID ~ "PDL033T",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL034D"]]$Cell.ID ~ "PDL034D",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL040A"]]$Cell.ID ~ "PDL040A",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL042D"]]$Cell.ID ~ "PDL042D",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL047T"]]$Cell.ID ~ "PDL047T",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL055T"]]$Cell.ID ~ "PDL055T",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL061"]]$Cell.ID ~ "PDL061",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL072"]]$Cell.ID ~ "PDL072",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL077"]]$Cell.ID ~ "PDL077",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL085A"]]$Cell.ID ~ "PDL085A",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL095D"]]$Cell.ID ~ "PDL095D",
                            cell_id %in% prime_v1_cell_stats_csv[["PDL097"]]$Cell.ID ~ "PDL097",
                            TRUE ~ "dropped"))

saveRDS(merged_obj_unfiltered_prime, 
        "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/prime_slides_2025_11_20.rds")

#### Z-SCORING TO ASSESS COMPETITION ----

transcript_count_prime_df <- inner_join(
  tibble(feature_name = rownames(transcript_count_files$prime_v1), prime_v1 = transcript_count_files$prime_v1$gene_count_filtered),
  tibble(feature_name = rownames(transcript_count_files$prime_solo), prime_solo = transcript_count_files$prime_solo$gene_count_filtered),
  by = "feature_name"
) 

transcript_count_prime_df_zscore <- transcript_count_prime_df %>%
  mutate(
    # Column-wise mean and SD
    mean_prime_v1   = mean(prime_v1, na.rm = TRUE),
    sd_prime_v1     = sd(prime_v1, na.rm = TRUE),
    mean_prime_solo = mean(prime_solo, na.rm = TRUE),
    sd_prime_solo   = sd(prime_solo, na.rm = TRUE),
    
    # Z-scores per gene (row) relative to column
    zscore_prime_v1   = (prime_v1   - mean_prime_v1) / sd_prime_v1,
    zscore_prime_solo = (prime_solo - mean_prime_solo) / sd_prime_solo
  )

transcript_count_prime_df_zscore_overlaps <- transcript_count_prime_df_zscore %>%
  filter(feature_name %in% gene_overlaps) %>%
  rename_with(~ "overlapping_genes", .cols = "feature_name") %>%
  select(overlapping_genes, zscore_prime_v1, zscore_prime_solo) %>%
  rename_with(~ "overlapping_zscore_prime_v1", .cols = "zscore_prime_v1") %>%
  rename_with(~ "overlapping_zscore_prime_solo", .cols = "zscore_prime_solo")

transcript_count_prime_df_zscore_nonoverlaps <- transcript_count_prime_df_zscore %>%
  filter(feature_name %!in% gene_overlaps) %>%
  rename_with(~ "nonoverlapping_genes", .cols = "feature_name") %>%
  select(nonoverlapping_genes, zscore_prime_v1, zscore_prime_solo) %>%
  rename_with(~ "nonoverlapping_zscore_prime_v1", .cols = "zscore_prime_v1") %>%
  rename_with(~ "nonoverlapping_zscore_prime_solo", .cols = "zscore_prime_solo") %>%
  slice_head(n = 239)

test_zscores <- cbind(transcript_count_prime_df_zscore_overlaps, transcript_count_prime_df_zscore_nonoverlaps)

test_zscores <- cbind(
  transcript_count_prime_df_zscore_overlaps,
  transcript_count_prime_df_zscore_nonoverlaps
)

# Pivot both V1 and solo z-scores into long format
test_long_zscores <- test_zscores %>%
  select(
    overlapping_v1 = overlapping_zscore_prime_v1,
    overlapping_solo = overlapping_zscore_prime_solo,
    nonoverlapping_v1 = nonoverlapping_zscore_prime_v1,
    nonoverlapping_solo = nonoverlapping_zscore_prime_solo
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("category", "type"),
    names_sep = "_",
    values_to = "zscore"
  ) %>%
  mutate(
    box = case_when(
      category == "overlapping" & type == "v1" ~ "Overlapping_prime_v1",
      category == "overlapping" & type == "solo" ~ "Overlapping_prime_solo",
      category == "nonoverlapping" & type == "v1" ~ "Non-overlapping_prime_v1",
      category == "nonoverlapping" & type == "solo" ~ "Non-overlapping_prime_solo"
    )
  ) %>%
  mutate(box = factor(box, levels = unique(box)))

ggplot(test_long_zscores, aes(x = box, y = zscore, fill = category)) +
  geom_boxplot(outlier.color = "black", alpha = 0.7) +
  labs(
    x = NULL,
    y = "Z-score",
  ) +
  theme(legend.position = "none") +
  ylim(-1, 10)

#### CELL PROXIMITY ----

filter <- dplyr::filter                                   
select <- dplyr::select
pull   <- dplyr::pull                                                                                                 

outDir <- '/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/proximity_results_epi_mes/'                  
obj <- readRDS('/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/fishing_2lines_dualrun_5481genes_20pcs_epithelial_mesenchymal_2026_05_06.rds')                                                                        
x   <- obj@meta.data                                      

sample_ids <- x %>% pull(sample) %>% as.character() %>% unique()                                                      

r <- 10                                                                                                               

nCores <- 24                                                                                                          
cl <- makeCluster(nCores)
registerDoParallel(cores = nCores)                                                                                    

clusterEvalQ(cl, {                                                                                                    
  library(dplyr)
  library(RANN)                                                                                                       
})                                                        

proximal_epi <- foreach(i = 1:length(sample_ids)) %dopar% {                                                           
  sid     <- sample_ids[i]
  obj.sid <- x %>% filter(sample == sid)                                                                              
  
  epi <- obj.sid %>% filter(Lineage == 'Epithelial')                                                                  
  mes <- obj.sid %>% filter(Lineage == 'Mesenchymal')                                                                 
  
  if (nrow(epi) == 0 || nrow(mes) == 0) return(NULL)                                                                  
  
  epi_coords <- epi %>% select(x_centroid_v1_prime, y_centroid_v1_prime)                                              
  mes_coords <- mes %>% select(x_centroid_v1_prime, y_centroid_v1_prime)
  
  nn <- RANN::nn2(data = mes_coords, query = epi_coords, k = 1)                                                       
  
  epi$nearest_mes_id       <- rownames(mes)[nn$nn.idx[, 1]]                                                           
  epi$nearest_mes_distance <- nn$nn.dists[, 1]            
  epi$sample               <- sid                                                                                     
  
  epi %>% filter(nearest_mes_distance <= r)
}                                                                                                                     

proximal_epi_clean <- proximal_epi[                                                                                   
  sapply(proximal_epi, function(XX) is.data.frame(XX) && nrow(XX) > 0)
]                                                                                                                     

proximal_epi_clean <- do.call('rbind', proximal_epi_clean)                                                            
stopImplicitCluster()                                     

### Add prox results to object

obj$nearest_mes_id       <- proximal_epi_clean$nearest_mes_id[match(rownames(obj@meta.data),                          
                                                                    rownames(proximal_epi_clean))]
obj$nearest_mes_distance <- proximal_epi_clean$nearest_mes_distance[match(rownames(obj@meta.data),                    
                                                                          rownames(proximal_epi_clean))]              
obj$proximal_to_mes      <- rownames(obj@meta.data) %in% rownames(proximal_epi_clean)

saveRDS(obj, '/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/objects/fishing_2lines_dualrun_5481genes_20pcs_epithelial_mesenchymal_proximity_2026_05_07.rds')                                                                


#### LIGAND - RECEPTOR PAIR IDENTIFICATION ----

.libPaths(c("/home/avannan/R/rstudio-4.3.0-4-with_modules.sif",
            "/usr/local/lib/R/site-library",
            "/usr/local/lib/R/library"))
library(rlang, lib.loc = "/home/avannan/R/rstudio-4.3.0-4-with_modules.sif")
library(tidyverse)
library(biomaRt,   lib.loc = "/home/avannan/R/rstudio-4.3.0-4-with_modules.sif")

set.seed(823)
options(scipen = 99999)
filter <- dplyr::filter
select <- dplyr::select

## Paths
LR_PATH       <- "/home/avannan/references/human_lr_pair.txt"
GENES_V1_PATH    <- "/scratch/avannan/projects/fishing2lines/results/genes_v1.txt"
GENES_PRIME_PATH <- "/scratch/avannan/projects/fishing2lines/results/genes_prime.txt"
GENES_DUAL_PATH  <- "/scratch/avannan/projects/fishing2lines/results/genes_dual.txt"
OUT_DIR       <- "/scratch/avannan/projects/fishing2lines/results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Load LR pairs
stopifnot(file.exists(LR_PATH))
lr_df <- read_tsv(LR_PATH, show_col_types = FALSE)

cat("Loaded LR pair list:\n")
cat("  Pairs:", nrow(lr_df), "\n")
cat("  Unique ligands:", n_distinct(lr_df$ligand_gene_symbol), "\n")
cat("  Unique receptors:", n_distinct(lr_df$receptor_gene_symbol), "\n")

lr_genes_all <- unique(c(lr_df$ligand_gene_symbol, lr_df$receptor_gene_symbol))
cat("  Unique LR genes combined:", length(lr_genes_all), "\n\n")

## Load Xenium Gene panel lists

# Gene lists pre-extracted from h5ad /var/_index via Python (h5py)
stopifnot(file.exists(GENES_V1_PATH), file.exists(GENES_PRIME_PATH), file.exists(GENES_DUAL_PATH))

genes_v1    <- readLines(GENES_V1_PATH)
genes_prime <- readLines(GENES_PRIME_PATH)
genes_dual  <- readLines(GENES_DUAL_PATH)

cat("Panel gene counts:\n")
cat("  V1    (expected  480):", length(genes_v1),    "\n")
cat("  Prime (expected 5001):", length(genes_prime), "\n")
cat("  Dual  (expected 5242):", length(genes_dual),  "\n\n")

## Direct match
# Check which LR genes appear by name in each panel as-is
v1_direct    <- intersect(lr_genes_all, genes_v1)
prime_direct <- intersect(lr_genes_all, genes_prime)
dual_direct  <- intersect(lr_genes_all, genes_dual)

cat("Direct matches (symbol identical):\n")
cat("  V1:    ", length(v1_direct),    "/", length(lr_genes_all), "LR genes\n")
cat("  Prime: ", length(prime_direct), "/", length(lr_genes_all), "LR genes\n")
cat("  Dual:  ", length(dual_direct),  "/", length(lr_genes_all), "LR genes\n\n")

## BIOMART alias check
# Some genes may be in the panel under an old/alternate HGNC symbol.
# Two-pass query:
#   Pass 1: for each LR gene (by canonical symbol), get all known synonyms.
#   Pass 2: for each panel gene NOT already a direct match, check if biomaRt
#           recognises it as a synonym for any LR gene.
# Wrapped in tryCatch — if Ensembl is unreachable, falls back to direct matches only.
cat("Querying biomaRt for gene aliases ...\n")

panel_alias_df <- tryCatch({
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host = "https://www.ensembl.org")
  
  ### Pass 1: LR gene -> synonyms 
  lr_synonym_df <- getBM(
    attributes = c("hgnc_symbol", "external_synonym"),
    filters    = "hgnc_symbol",
    values     = lr_genes_all,
    mart       = mart
  ) %>%
    filter(external_synonym != "")   # drop rows where no synonym was returned
  
  cat("  LR genes with >=1 synonym:", n_distinct(lr_synonym_df$hgnc_symbol), "\n")
  
  ### Pass 2: panel genes not yet matched -> do they map to an LR gene? 
  panel_unmatched_v1    <- setdiff(genes_v1,    lr_genes_all)
  panel_unmatched_prime <- setdiff(genes_prime, lr_genes_all)
  panel_unmatched_dual  <- setdiff(genes_dual,  lr_genes_all)
  panel_unmatched_all   <- unique(c(panel_unmatched_v1, panel_unmatched_prime, panel_unmatched_dual))
  
  # Query: treat each unmatched panel gene as a potential synonym
  alias_df <- getBM(
    attributes = c("hgnc_symbol", "external_synonym"),
    filters    = "external_synonym",
    values     = panel_unmatched_all,
    mart       = mart
  ) %>%
    filter(hgnc_symbol %in% lr_genes_all) %>%
    dplyr::rename(panel_alias = external_synonym, lr_canonical = hgnc_symbol)
  
  cat("  Panel genes that are aliases for an LR gene:",
      n_distinct(alias_df$panel_alias), "\n")
  
  if (nrow(alias_df) > 0) {
    cat("\n  Alias matches found:\n")
    print(as.data.frame(alias_df %>% arrange(lr_canonical)), row.names = FALSE)
  } else {
    cat("  No alias matches found — all coverage is by direct symbol match.\n")
  }
  
  alias_df
  
}, error = function(e) {
  cat("  WARNING: biomaRt query failed:", conditionMessage(e), "\n")
  cat("  Proceeding with direct symbol matches only.\n")
  # Return an empty data frame with the expected columns
  tibble(panel_alias = character(), lr_canonical = character())
})

## Build coverage sets
# For each panel, the set of LR genes covered = direct matches UNION alias-resolved matches

# Panel genes that resolve to an LR gene via alias
v1_via_alias <- panel_alias_df %>%
  filter(panel_alias %in% genes_v1) %>%
  pull(lr_canonical) %>%
  unique()

prime_via_alias <- panel_alias_df %>%
  filter(panel_alias %in% genes_prime) %>%
  pull(lr_canonical) %>%
  unique()

dual_via_alias <- panel_alias_df %>%
  filter(panel_alias %in% genes_dual) %>%
  pull(lr_canonical) %>%
  unique()

covered_v1    <- unique(c(v1_direct,    v1_via_alias))
covered_prime <- unique(c(prime_direct, prime_via_alias))
covered_dual  <- unique(c(dual_direct,  dual_via_alias))

cat("\nTotal LR genes covered after alias resolution:\n")
cat("  V1:    ", length(covered_v1),    "/", length(lr_genes_all), "\n")
cat("  Prime: ", length(covered_prime), "/", length(lr_genes_all), "\n")
cat("  Dual:  ", length(covered_dual),  "/", length(lr_genes_all), "\n\n")

## Count complete LR pairs
# A pair is COMPLETE when BOTH the ligand AND receptor are covered by the panel.
coverage_df <- lr_df %>%
  mutate(
    ligand_in_v1      = ligand_gene_symbol   %in% covered_v1,
    receptor_in_v1    = receptor_gene_symbol %in% covered_v1,
    ligand_in_prime   = ligand_gene_symbol   %in% covered_prime,
    receptor_in_prime = receptor_gene_symbol %in% covered_prime,
    ligand_in_dual    = ligand_gene_symbol   %in% covered_dual,
    receptor_in_dual  = receptor_gene_symbol %in% covered_dual,
    complete_v1       = ligand_in_v1    & receptor_in_v1,
    complete_prime    = ligand_in_prime & receptor_in_prime,
    complete_dual     = ligand_in_dual  & receptor_in_dual,
    complete_any      = complete_v1 | complete_prime | complete_dual,
    complete_all      = complete_v1 & complete_prime & complete_dual,
  )

cat("============================================================\n")
cat("COMPLETE LR PAIRS (both ligand and receptor on-panel):\n")
cat("  V1   :", sum(coverage_df$complete_v1),    "/", nrow(lr_df), "pairs\n")
cat("  Prime:", sum(coverage_df$complete_prime), "/", nrow(lr_df), "pairs\n")
cat("  Dual :", sum(coverage_df$complete_dual),  "/", nrow(lr_df), "pairs\n")
cat("  Any panel :", sum(coverage_df$complete_any), "\n")
cat("  All panels:", sum(coverage_df$complete_all), "\n")
cat("  Prime only (not V1):", sum(coverage_df$complete_prime & !coverage_df$complete_v1), "\n")
cat("  Dual only  (not V1):", sum(coverage_df$complete_dual  & !coverage_df$complete_v1), "\n")
cat("============================================================\n\n")

## Pairs gained by larger panels over V1
dual_only_pairs <- coverage_df %>%
  filter(complete_dual, !complete_v1) %>%
  select(lr_pair, ligand_gene_symbol, receptor_gene_symbol,
         ligand_in_v1, receptor_in_v1)

cat("Pairs gained by Dual over V1 (", nrow(dual_only_pairs), " pairs):\n", sep = "")
if (nrow(dual_only_pairs) > 0) {
  dual_only_pairs %>%
    arrange(ligand_gene_symbol) %>%
    print(n = 30)
  if (nrow(dual_only_pairs) > 30) cat("  ... (see lr_coverage.tsv for full list)\n")
}

out_path <- file.path(OUT_DIR, "lr_coverage.tsv")
coverage_df %>%
  write_tsv(out_path)
cat("\nFull coverage table saved -> ", out_path, "\n")

#### LIGAND - RECEPTOR PAIR ANALYSIS ----                        

lr_pairs <- read.table(
  '/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/human_lr_pair.txt',
  header = TRUE, sep = '\t'
)

counts    <- GetAssayData(obj, layer = 'counts')                                                                      
all_genes <- rownames(counts)

genes_v1    <- all_genes[endsWith(all_genes, '-v1-prime')]                                                            
genes_prime <- all_genes[endsWith(all_genes, '-prime-with-v1-seg')]

make_gene_expr <- function(gene_name, panel_genes) {                                                                  
  matches <- panel_genes[startsWith(panel_genes, paste0(gene_name, '-'))]
  if (length(matches) == 0) return(NULL)                                                                              
  if (length(matches) == 1) return(setNames(as.numeric(counts[matches, ]), colnames(counts)))
  return(Matrix::colSums(counts[matches, , drop = FALSE]))                                                            
}                                                         

run_lr <- function(panel_genes, panel_name) {                                                                         
  cat(sprintf("\nBuilding expression vectors for %s panel...\n", panel_name))
  lr_genes  <- union(lr_pairs$ligand_gene_symbol, lr_pairs$receptor_gene_symbol)                                      
  gene_expr <- lapply(setNames(lr_genes, lr_genes), make_gene_expr, panel_genes = panel_genes)                        
  gene_expr <- gene_expr[!sapply(gene_expr, is.null)]                                                                 
  
  epi_ids <- rownames(proximal_epi_clean)                                                                             
  mes_ids <- proximal_epi_clean$nearest_mes_id            
  
  cat(sprintf("Running LR matching for %s panel...\n", panel_name))                                                   
  lr_hits <- lapply(1:nrow(lr_pairs), function(j) {
    if (j %% 500 == 0) cat(sprintf("  LR pair %d / %d\n", j, nrow(lr_pairs)))                                         
    
    lig <- lr_pairs$ligand_gene_symbol[j]                                                                             
    rec <- lr_pairs$receptor_gene_symbol[j]                                                                           
    
    if (!lig %in% names(gene_expr) || !rec %in% names(gene_expr)) return(NULL)
    
    lig_expr <- gene_expr[[lig]][epi_ids]                                                                             
    rec_expr <- gene_expr[[rec]][mes_ids]
    
    hit_idx <- which(!is.na(lig_expr) & !is.na(rec_expr) & lig_expr > 0 & rec_expr > 0)                               
    if (length(hit_idx) == 0) return(NULL)
    
    data.frame(                                           
      epi_cell      = epi_ids[hit_idx],                                                                               
      mes_cell      = mes_ids[hit_idx],                   
      ligand        = lig,                                                                                            
      receptor      = rec,
      lr_pair       = lr_pairs$lr_pair[j],                                                                            
      ligand_expr   = lig_expr[hit_idx],                  
      receptor_expr = rec_expr[hit_idx]                                                                               
    )
  })                                                                                                                  
  
  hits <- do.call('rbind', lr_hits[!sapply(lr_hits, is.null)])                                                        
  saveRDS(hits, file.path(outDir, sprintf('epi_mes_lr_hits_r10_%s_2026_05_07.rds', panel_name)))
  write.csv(hits, file.path(outDir, sprintf('epi_mes_lr_hits_r10_%s_2026_05_07.csv', panel_name)), row.names = FALSE) 
  cat(sprintf("Saved %s results: %d hits, %d unique LR pairs\n", panel_name, nrow(hits),                              
              length(unique(hits$lr_pair))))                                                                          
  return(hits)                                                                                                        
} 

hits_v1    <- run_lr(genes_v1,    'v1')                                                                               
hits_prime <- run_lr(genes_prime, 'prime')

cat("\nDone!\n")                                                                                                      
cat(sprintf("V1    — unique LR pairs: %d\n", length(unique(hits_v1$lr_pair))))
cat(sprintf("Prime — unique LR pairs: %d\n", length(unique(hits_prime$lr_pair))))

cat(sprintf("V1    — total cell pairs with LR hits: %d\n", nrow(hits_v1)))
cat(sprintf("Prime — total cell pairs with LR hits: %d\n", nrow(hits_prime)))
cat(sprintf("Overlapping LR pairs: %d\n", length(intersect(unique(hits_v1$lr_pair), unique(hits_prime$lr_pair)))))

# Summary df and saving it
all_pairs <- data.frame(
  lr_pair = union(unique(hits_v1$lr_pair), unique(hits_prime$lr_pair))
)

all_pairs$V1      <- ifelse(all_pairs$lr_pair %in% unique(hits_v1$lr_pair),    "yes", "no")
all_pairs$Prime   <- ifelse(all_pairs$lr_pair %in% unique(hits_prime$lr_pair), "yes", "no")
all_pairs$Overlap <- ifelse(all_pairs$V1 == "yes" & all_pairs$Prime == "yes",  "yes", "no")

print(all_pairs)

write_csv(all_pairs, 
          "/scratch/smallapragada/v1_prime5k_comparison_finished_2026_01/ligand_receptor_pairs_all_panels.csv")

