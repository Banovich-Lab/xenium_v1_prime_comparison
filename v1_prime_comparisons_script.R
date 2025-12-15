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

## Set seed & working directory
set.seed(0712)
work_dir <- setwd("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/")

options(scipen=999)

## "not in" operand function
'%!in%' <- function(x,y)!('%in%'(x,y))

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
          "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/prime_solo_transcripts_per_gene_qv20_2025_11_19.csv")
write_csv(transcript_count_files$prime_v1, 
          "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/prime_v1_transcripts_per_gene_qv20_2025_11_19.csv")
write_csv(transcript_count_files$v1_prime, 
          "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/v1_prime_transcripts_per_gene_qv20_2025_11_19.csv")
write_csv(transcript_count_files$v1_solo, 
          "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/v1_solo_transcripts_per_gene_qv20_2025_11_19.csv")

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

saveRDS(seurat_obj_slot_filtered, "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_filtered_joined_2025_11_25.rds")

### Adding raw counts back into split objects post-RAPIDS 

## Convert RAPIDS objects into Seurat
dual_obj_combined <- schard::h5ad2seurat("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_comparison_allgenes_20pcs_filtered_split_2025_12_05.h5ad")
saveRDS(dual_obj_combined, "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_comparison_allgenes_20pcs_filtered_split_postrapids_2025_12_07.rds")

dual_obj_v1 <- schard::h5ad2seurat("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_comparison_v1genes_20pcs_filtered_split_2025_12_05.h5ad")
saveRDS(dual_obj_v1, "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_comparison_v1genes_20pcs_filtered_split_postrapids_2025_12_07.rds")

dual_obj_prime <- schard::h5ad2seurat("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_primegenes_combined_50pcs_filtered_split_2025_12_05.h5ad")
saveRDS(dual_obj_prime, "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_comparison_primegenes_50pcs_filtered_split_postrapids_2025_12_07.rds")

## Import counts file

raw_counts_matrix <- fread("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/test_counts.csv",
                           sep = ",", header = TRUE)

raw_counts_matrix_adj <- raw_counts_matrix %>%
  column_to_rownames(var = "V1")

raw_counts_tibble <- read_delim(
  file = "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/test_counts.csv",
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
v1_solo_cell_stats <- list.files(path = "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/cell_stats", 
                                 pattern = "V1_solo_cells", 
                                 full.names = TRUE, 
                                 recursive = FALSE)
v1_prime_cell_stats <- list.files(path = "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/cell_stats", 
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
        "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_slides_2025_11_20.rds")

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
prime_solo_cell_stats <- list.files(path = "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/cell_stats", 
                                 pattern = "prime_solo_cells", 
                                 full.names = TRUE, 
                                 recursive = FALSE)
prime_v1_cell_stats <- list.files(path = "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/cell_stats", 
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
        "/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/prime_slides_2025_11_20.rds")

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

