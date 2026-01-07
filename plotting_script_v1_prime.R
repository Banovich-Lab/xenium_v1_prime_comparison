#### LOAD PACKAGES AND SET ENVIRONMENT ----
## Packages

my_lib_path <- .libPaths()[2]

library(tidyverse)
library(stringr)
library(rlist)
library(Seurat)
library(SeuratObject)
library(data.table)
library(Matrix)
library(ggplot2, lib.loc = my_lib_path)
library(reshape2)
library(scales)
library(forcats)
library(stats)
library(circlize)
library(viridisLite)
library(RColorBrewer)
library(ComplexHeatmap)
library(cowplot)
library(patchwork)
library(scCustomize)
library(ggalluvial)
library(ggraph, lib.loc = my_lib_path)

## Set seed & working directory
set.seed(0712)
work_dir <- setwd("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/")

options(scipen=999)

## "not in" operand function
'%!in%' <- function(x,y)!('%in%'(x,y))

#### FIGURE 1 ----

# Identifying genes as overlapping between V1 and prime vs not
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

### Figure 1B - V1 correlation scatterplot

# V1_solo vs. V1_Prime
transcript_count_files$v1_prime$gene_count_filtered <- as.numeric(transcript_count_files$v1_prime$gene_count_filtered)
transcript_count_files$v1_solo$gene_count_filtered <- as.numeric(transcript_count_files$v1_solo$gene_count_filtered)

merged_df_v1 <- inner_join(
  tibble(feature_name = rownames(transcript_count_files$v1_prime), x = transcript_count_files$v1_prime$gene_count_filtered),
  tibble(feature_name = rownames(transcript_count_files$v1_solo), y = transcript_count_files$v1_solo$gene_count_filtered),
  by = "feature_name"
)

# correlation stats
cor_test_v1 <- cor.test(merged_df_v1$x, merged_df_v1$y, method = "pearson")

merged_df_v1 <- merged_df_v1 %>%
  mutate(gene_type = case_when(merged_df_v1$feature_name %in% gene_overlaps ~ "Overlapping gene",
                               TRUE ~ "V1 panel gene"))

v1_plot <- ggplot(merged_df_v1, aes(x = x, y = y, color = gene_type, shape = gene_type)) +
  geom_point(alpha = 1, size = 1.0) +
  theme_classic() +
  labs(
    x = "Total transcripts (V1 + Prime)",
    y = "Total transcripts (V1 only)"
  ) +
  scale_x_continuous(
    limits = c(0, 800000), 
    labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
  ) +
  scale_y_continuous(
    limits = c(0, 800000),
    labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(
    values = c("Overlapping gene" = "#FFB000", "V1 panel gene" = "#DC267F")
  ) +
  scale_shape_manual(
    values = c("Overlapping gene" = 17, "V1 panel gene" = 15) 
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 1, shape = c(17, 15))), 
    shape = "none"
  ) + 
  theme(
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(color = "black", face = "bold", size = 7),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = element_text(color = "black", size = 7),
    legend.position = "bottom",
    legend.title = element_blank()
  )

pdf("/home/smallapragada/v1_5K_panel_comparison_project/v1_panel_scatterplot.pdf", width = 2.3, height = 2.3)
v1_plot
dev.off()

# Prime_solo vs. Prime_v1
transcript_count_files$prime_v1$gene_count_filtered <- as.numeric(transcript_count_files$prime_v1$gene_count_filtered)
transcript_count_files$prime_solo$gene_count_filtered <- as.numeric(transcript_count_files$prime_solo$gene_count_filtered)

merged_df_prime <- inner_join(
  tibble(feature_name = rownames(transcript_count_files$prime_v1), x = transcript_count_files$prime_v1$gene_count_filtered),
  tibble(feature_name = rownames(transcript_count_files$prime_solo), y = transcript_count_files$prime_solo$gene_count_filtered),
  by = "feature_name"
)

# correlation stats
cor_test_prime <- cor.test(merged_df_prime$x, merged_df_prime$y, method = "pearson")

merged_df_prime <- merged_df_prime %>%
  mutate(gene_type = case_when(merged_df_prime$feature_name %in% gene_overlaps ~ "Overlapping gene",
                               TRUE ~ "Prime panel gene"))

prime_plot <-  ggplot(merged_df_prime, aes(x = x, y = y, color = gene_type, shape = gene_type)) +
  geom_point(alpha = 1, size = 1.0) +
  theme_classic() +
  labs(
    x = "Total transcripts (Dual)",
    y = "Total transcripts (Prime only)"
  ) +
  scale_x_continuous(
    limits = c(0, 800000), 
    labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
  ) +
  scale_y_continuous(
    limits = c(0, 800000),
    labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(
    values = c("Overlapping gene" = "#FFB000", "Prime panel gene" = "#648FFF")
  ) +
  scale_shape_manual(
    values = c("Overlapping gene" = 17, "Prime panel gene" = 19) 
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 1, shape = c(17, 19))), 
    shape = "none"
  ) + 
  theme(
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(color = "black", face = "bold", size = 7),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 6),
    axis.title.x = element_text(color = "black", face = "bold", size = 7),
    legend.text = element_text(color = "black", size = 7),
    legend.position = "top",
    legend.title = element_blank()
  )

pdf("/home/smallapragada/v1_5K_panel_comparison_project/prime_panel_scatterplot.pdf", width = 2.3, height = 2.5)
prime_plot
dev.off()

### Figure 1C - sample correlation plots

# Split objects by sample and by slide type
merged_obj_unfiltered_v1_solo <- subset(merged_obj_unfiltered_v1, subset = slide_type == "v1_solo")
merged_obj_unfiltered_v1_dual <- subset(merged_obj_unfiltered_v1, subset = slide_type == "v1_prime")
merged_obj_unfiltered_v1_solo_split <- SplitObject(merged_obj_unfiltered_v1_solo, split.by = "sample")
merged_obj_unfiltered_v1_dual_split <- SplitObject(merged_obj_unfiltered_v1_dual, split.by = "sample")

merged_obj_unfiltered_prime_solo <- subset(merged_obj_unfiltered_prime, subset = slide_type == "prime_solo")
merged_obj_unfiltered_prime_dual <- subset(merged_obj_unfiltered_prime, subset = slide_type == "prime_v1")
merged_obj_unfiltered_prime_solo_split <- SplitObject(merged_obj_unfiltered_prime_solo, split.by = "sample")
merged_obj_unfiltered_prime_dual_split <- SplitObject(merged_obj_unfiltered_prime_dual, split.by = "sample")

# Extracting the gene counts and pulling them into dfs

listed_list <- list(
  v1_solo = merged_obj_unfiltered_v1_solo_split,
  v1_prime = merged_obj_unfiltered_v1_dual_split,
  prime_solo = merged_obj_unfiltered_prime_solo_split,
  prime_v1 = merged_obj_unfiltered_prime_dual_split
)

all_results_list <- list()

for (l in names(listed_list)) {
  
  current_split_list <- listed_list[[l]]
  sample_ids <- names(current_split_list)
  
  for (i in seq_along(current_split_list)) {
    
    current_sample_id <- sample_ids[i]
    seurat_object <- current_split_list[[i]]
    
    counts_matrix <- GetAssayData(object = seurat_object, slot = "counts")
    gene_counts_vector <- rowSums(counts_matrix)
    
    temp_df <- tibble(
      feature_name = names(gene_counts_vector),
      gene_count = gene_counts_vector,
      sample = current_sample_id,
      data_source = l
    )
    
    combined_key <- paste(l, current_sample_id, sep = "_")
    all_results_list[[combined_key]] <- temp_df
  }
}

gene_counts_df_long <- bind_rows(all_results_list)

gene_counts_df_wide <- gene_counts_df_long %>%
  pivot_wider(
    id_cols = c(feature_name, sample),
    names_from = data_source,
    values_from = gene_count
  )

# Group by sample and calculate the two required correlations (R and P-value)
final_correlation_df <- gene_counts_df_wide %>%
  group_by(sample) %>%
  do({
    df_sample <- . 
    
    # v1_solo against v1_prime
    cor_v1 <- cor.test(df_sample$v1_solo, df_sample$v1_prime, method = "pearson")
    
    # prime_solo against prime_v1
    cor_prime <- cor.test(df_sample$prime_solo, df_sample$prime_v1, method = "pearson")
    
    tibble(
      R_v1_solo_v1_prime = cor_v1$estimate,
      P_v1_solo_v1_prime = cor_v1$p.value,
      
      R_prime_solo_prime_v1 = cor_prime$estimate,
      P_prime_solo_prime_v1 = cor_prime$p.value
    )
  }) %>%
  ungroup()

## Heatmap of R-values

final_correlation_mat <- final_correlation_df %>%
  select(sample, R_prime_solo_prime_v1, R_v1_solo_v1_prime) %>% 
  filter(sample %!in% "dropped") %>%
  rename_with(~ "V1", .cols = "R_v1_solo_v1_prime") %>%
  rename_with(~ "Prime", .cols = "R_prime_solo_prime_v1") %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()

write_csv(final_correlation_mat, "/home/smallapragada/v1_5K_panel_comparison_project/correlation_heatmap_r_vals.csv")

magma_colors <- rev(magma(256))
magma_colors <- magma_colors[-(200:256)]
color_mapping <- colorRamp2(
  breaks = seq(0.90, 1.0, length.out = length(magma_colors)),
  colors = magma_colors
)

heatmap <- Heatmap(final_correlation_mat, 
                   col = color_mapping,
                   name = "R value", 
                   show_column_names = TRUE,
                   show_row_names = TRUE,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_names_gp = gpar(fontsize = 7),
                   row_names_side = "left",
                   column_names_side = "top",
                   column_names_gp = gpar(fontsize = 7, 
                                          fontface = "bold"),
                   column_names_rot = 0,
                   heatmap_legend_param = list(
                     title = "R-value",
                     legend_width = unit(1.5, "cm"),
                     title_gp = gpar(fontsize = 6, fontface = "bold"),
                     labels_gp = gpar(fontsize = 5)
                   ))

draw(heatmap)

pdf("/home/smallapragada/v1_5K_panel_comparison_project/heatmap_label.pdf", width = 1.7, height = 3.1)
draw(heatmap)
dev.off()

heatmap <- Heatmap(final_correlation_mat, 
                   col = color_mapping,
                   name = "R value", 
                   show_column_names = TRUE,
                   show_row_names = TRUE,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_names_gp = gpar(fontsize = 7),
                   row_names_side = "left",
                   column_names_side = "top",
                   column_names_gp = gpar(fontsize = 7, 
                                          fontface = "bold"),
                   column_names_rot = 0,
                   show_heatmap_legend = FALSE)

draw(heatmap)

pdf("/home/smallapragada/v1_5K_panel_comparison_project/heatmap_plot.pdf", width = 1.7, height = 3.1)
draw(heatmap)
dev.off()

### Figure 1D - Z-score boxplots V1 & Prime

## V1 boxplots

transcript_count_v1_df <- inner_join(
  tibble(feature_name = rownames(transcript_count_files$v1_prime), v1_prime = transcript_count_files$v1_prime$gene_count_filtered),
  tibble(feature_name = rownames(transcript_count_files$v1_solo), v1_solo = transcript_count_files$v1_solo$gene_count_filtered),
  by = "feature_name"
) 

transcript_count_v1_df_zscore <- transcript_count_v1_df %>%
  mutate(
    # Column-wise mean and SD
    mean_v1_prime = mean(v1_prime, na.rm = TRUE),
    sd_v1_prime = sd(v1_prime, na.rm = TRUE),
    mean_v1_solo = mean(v1_solo, na.rm = TRUE),
    sd_v1_solo = sd(v1_solo, na.rm = TRUE),
    # Z-scores per gene (row) relative to column
    zscore_v1_prime   = (v1_prime   - mean_v1_prime) / sd_v1_prime,
    zscore_v1_solo = (v1_solo - mean_v1_solo) / sd_v1_solo
  )

transcript_count_v1_df_zscore_overlaps <- transcript_count_v1_df_zscore %>%
  filter(feature_name %in% gene_overlaps) %>%
  rename_with(~ "overlapping_genes", .cols = "feature_name") %>%
  select(overlapping_genes, zscore_v1_prime, zscore_v1_solo) %>%
  rename_with(~ "overlapping_zscore_v1_prime", .cols = "zscore_v1_prime") %>%
  rename_with(~ "overlapping_zscore_v1_solo", .cols = "zscore_v1_solo")

# Pulling a random subset of 239 genes that don't overlap
transcript_count_v1_df_zscore_nonoverlaps <- transcript_count_v1_df_zscore %>%
  filter(feature_name %!in% gene_overlaps) %>%
  rename_with(~ "nonoverlapping_genes", .cols = "feature_name") %>%
  select(nonoverlapping_genes, zscore_v1_prime, zscore_v1_solo) %>%
  rename_with(~ "nonoverlapping_zscore_v1_prime", .cols = "zscore_v1_prime") %>%
  rename_with(~ "nonoverlapping_zscore_v1_solo", .cols = "zscore_v1_solo") %>%
  slice_head(n = 239)

df_zscores_v1 <- cbind(transcript_count_v1_df_zscore_overlaps, transcript_count_v1_df_zscore_nonoverlaps)

# Pivot both V1 and solo z-scores into long format
df_long_zscores_v1 <- df_zscores_v1 %>%
  select(
    overlapping_dual = overlapping_zscore_v1_prime,
    overlapping_solo = overlapping_zscore_v1_solo,
    nonoverlapping_dual = nonoverlapping_zscore_v1_prime,
    nonoverlapping_solo = nonoverlapping_zscore_v1_solo
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("category", "type"),
    names_sep = "_",
    values_to = "zscore"
  ) %>%
  mutate(
    box = case_when(
      category == "overlapping" & type == "dual" ~ "Overlap_dual",
      category == "overlapping" & type == "solo" ~ "Overlap_prime",
      category == "nonoverlapping" & type == "dual" ~ "Nonoverlap_dual",
      category == "nonoverlapping" & type == "solo" ~ "Nonoverlap_prime"
    )
  ) %>%
  mutate(
    box = factor(box, levels = c(
      "Nonoverlap_prime",
      "Nonoverlap_dual",
      "Overlap_prime",
      "Overlap_dual"
    ))
  )

# Plot
axis_labels <- c("V1",
                 "Dual",
                 "V1",
                 "Dual")

v1_zscore_plot <- ggplot(df_long_zscores_v1, aes(x = box, y = zscore, fill = category)) +
  geom_boxplot(width = 0.7, outlier.color = "black", alpha = 1, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  labs(
    x = NULL,
    y = "Z-score"
  ) +
  scale_x_discrete(labels = axis_labels) +
  ylim(-1, 2) +
  scale_fill_manual(
    labels = c("Overlapping genes", "Non-overlapping genes"),
    values = c("overlapping" = "#FFB000", "nonoverlapping" = "#DC267F")  
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 7, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(color = "black", size = 7, face = "bold"),
        strip.text = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.text = element_text(color = "black", size = 7),
        legend.title = element_blank())

pdf("/home/smallapragada/v1_5K_panel_comparison_project/v1_boxplots_zscore.pdf", width = 2.75, height = 2.38)
v1_zscore_plot
dev.off()

## Prime boxplots

transcript_count_prime_df <- inner_join(
  tibble(feature_name = rownames(transcript_count_files$prime_v1), prime_v1 = transcript_count_files$prime_v1$gene_count_filtered),
  tibble(feature_name = rownames(transcript_count_files$prime_solo), prime_solo = transcript_count_files$prime_solo$gene_count_filtered),
  by = "feature_name"
) 

transcript_count_prime_df_zscore <- transcript_count_prime_df %>%
  mutate(
    # Column-wise mean and SD
    mean_prime_v1 = mean(prime_v1, na.rm = TRUE),
    sd_prime_v1 = sd(prime_v1, na.rm = TRUE),
    mean_prime_solo = mean(prime_solo, na.rm = TRUE),
    sd_prime_solo = sd(prime_solo, na.rm = TRUE),
    # Z-scores per gene (row) relative to column
    zscore_prime_v1 = (prime_v1   - mean_prime_v1) / sd_prime_v1,
    zscore_prime_solo = (prime_solo - mean_prime_solo) / sd_prime_solo
  )

transcript_count_prime_df_zscore_overlaps <- transcript_count_prime_df_zscore %>%
  filter(feature_name %in% gene_overlaps) %>%
  rename_with(~ "overlapping_genes", .cols = "feature_name") %>%
  select(overlapping_genes, zscore_prime_v1, zscore_prime_solo) %>%
  rename_with(~ "overlapping_zscore_prime_v1", .cols = "zscore_prime_v1") %>%
  rename_with(~ "overlapping_zscore_prime_solo", .cols = "zscore_prime_solo")

# Pulling a random subset of 239 genes that don't overlap
transcript_count_prime_df_zscore_nonoverlaps <- transcript_count_prime_df_zscore %>%
  filter(feature_name %!in% gene_overlaps) %>%
  rename_with(~ "nonoverlapping_genes", .cols = "feature_name") %>%
  select(nonoverlapping_genes, zscore_prime_v1, zscore_prime_solo) %>%
  rename_with(~ "nonoverlapping_zscore_prime_v1", .cols = "zscore_prime_v1") %>%
  rename_with(~ "nonoverlapping_zscore_prime_solo", .cols = "zscore_prime_solo") %>%
  slice_head(n = 239)

df_zscores_prime <- cbind(transcript_count_prime_df_zscore_overlaps, transcript_count_prime_df_zscore_nonoverlaps)

# Pivot both V1 and solo z-scores into long format
df_long_zscores_prime <- df_zscores_prime %>%
  select(
    overlapping_dual = overlapping_zscore_prime_v1,
    overlapping_solo = overlapping_zscore_prime_solo,
    nonoverlapping_dual = nonoverlapping_zscore_prime_v1,
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
      category == "overlapping" & type == "dual" ~ "Overlap_dual",
      category == "overlapping" & type == "solo" ~ "Overlap_prime",
      category == "nonoverlapping" & type == "dual" ~ "Nonoverlap_dual",
      category == "nonoverlapping" & type == "solo" ~ "Nonoverlap_prime"
    )
  ) %>%
  mutate(
    box = factor(box, levels = c(
      "Nonoverlap_prime",
      "Nonoverlap_dual",
      "Overlap_prime",
      "Overlap_dual"
    ))
  )

# Plot
axis_labels <- c("Prime",
                 "Dual",
                 "Prime",
                 "Dual")

prime_zscore_plot <- ggplot(df_long_zscores_prime, aes(x = box, y = zscore, fill = category)) +
  geom_boxplot(width = 0.7, outlier.color = "black", alpha = 1, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  labs(
    x = NULL,
    y = "Z-score"
  ) +
  ylim(-1, 4) +
  scale_x_discrete(labels = axis_labels) +
  scale_fill_manual(
    labels = c("Overlapping genes", "Non-overlapping genes"),
    values = c("overlapping" = "#FFB000", "nonoverlapping" = "#648FFF")  
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 7, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(color = "black", size = 7, face = "bold"),
        strip.text = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.text = element_text(color = "black", size = 7),
        legend.title = element_blank())

pdf("/home/smallapragada/v1_5K_panel_comparison_project/prime_boxplots_zscore.pdf", width = 2.75, height = 2.4)
prime_zscore_plot
dev.off()

### Figure 1E - Violin plots

## V1 only plot

merged_obj_unfiltered_v1 <- readRDS("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_slides_2025_11_20.rds")

plot_data_v1 <- FetchData(
  object = merged_obj_unfiltered_v1, 
  vars = c("slide_type", "nCount_RNA", "nFeature_RNA"),
  slot = "counts" 
)

plot_data_v1$slide_type <- factor(
  plot_data_v1$slide_type, 
  level = c("v1_solo", "v1_prime"),
  labels = c("V1", "Dual")
)

count_plot_v1 <- ggplot(plot_data_v1, aes(x = slide_type, y = nCount_RNA, fill = slide_type)) +
  geom_violin(position = position_dodge(1), alpha = 1, scale = "width") +
  theme_classic() + 
  scale_y_continuous(breaks = pretty_breaks(n = 3)) + 
  scale_fill_manual(values= c("#FE6100", "#662D91")) +
  labs(y = "Total transcripts", title = "V1 panel") +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(color = "black", size = 7, face = "bold"),
        strip.text = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 8, color = "black", hjust = 0.5),
        legend.position = "none")

pdf("/home/smallapragada/v1_5K_panel_comparison_project/v1_violin_ncount.pdf", width = 1.7, height = 1.2)
count_plot_v1
dev.off()

feature_plot_v1 <- ggplot(plot_data_v1, aes(x = slide_type, y = nFeature_RNA, fill = slide_type)) +
  geom_violin(position = position_dodge(1), alpha = 1, scale = "width") +
  theme_classic() + 
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  scale_fill_manual(values= c("#FE6100", "#662D91")) +
  labs(y = "Total unique transcripts", title = " ") +
  theme(axis.text.x = element_text(color = "black", size = 7, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(color = "black", size = 7, face = "bold"),
        strip.text = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none")

pdf("/home/smallapragada/v1_5K_panel_comparison_project/v1_violin_nfeature.pdf", width = 1.7, height = 1.4)
feature_plot_v1
dev.off()

## Prime only plot

merged_obj_unfiltered_prime <- readRDS("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/prime_slides_2025_11_20.rds")

plot_data_prime <- FetchData(
  object = merged_obj_unfiltered_prime, 
  vars = c("slide_type", "nCount_RNA", "nFeature_RNA"),
  slot = "counts" 
)

plot_data_prime$slide_type <- factor(
  plot_data_prime$slide_type, 
  level = c("prime_solo", "prime_v1"),
  labels = c("Prime", "Dual")
)

count_plot_prime <- ggplot(plot_data_prime, aes(x = slide_type, y = nCount_RNA, fill = slide_type)) +
  geom_violin(position = position_dodge(1), alpha = 1, scale = "width") +
  theme_classic() + 
  scale_y_continuous(breaks = pretty_breaks(n = 3)) + 
  scale_fill_manual(values= c("#2E6F40", "#662D91")) +
  labs(y = " ", title = "Prime panel") +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(color = "black", size = 7, face = "bold"),
        strip.text = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 8, color = "black", hjust = 0.5),
        legend.position = "none")

pdf("/home/smallapragada/v1_5K_panel_comparison_project/prime_violin_ncount.pdf", width = 1.7, height = 1.2)
count_plot_prime
dev.off()

feature_plot_prime <- ggplot(plot_data_prime, aes(x = slide_type, y = nFeature_RNA, fill = slide_type)) +
  geom_violin(position = position_dodge(1), alpha = 1, scale = "width") +
  theme_classic() + 
  #scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  scale_fill_manual(values= c("#2E6F40", "#662D91")) +
  labs(y = "", title = " ") +
  theme(axis.text.x = element_text(color = "black", size = 7, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(color = "black", size = 7, face = "bold"),
        strip.text = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none")

pdf("/home/smallapragada/v1_5K_panel_comparison_project/prime_violin_nfeature.pdf", width = 1.7, height = 1.4)
feature_plot_prime
dev.off()

#### FIGURE 2 ----

## Figure 2B - Cell count bar plots pre and post filtering

# CSV from python
cell_counts <- read.delim("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/cell_counts_prepostfiltering.csv", sep = ",")

# Setting up data for plot
cell_counts$Filtering_Stage <- factor(cell_counts$Filtering_Stage,
                                      level = c("Unfiltered", "Filtered (combined)",
                                                "Filtered (V1)", "Filtered (Prime)"))
cell_counts <- cell_counts %>%
  mutate(groupings = case_when(cell_counts$Filtering_Stage %in% c("Unfiltered") ~ "Pre-filtering",
                               TRUE ~ "Post-filtering")) %>%
  mutate(groupings = factor(groupings, levels = c("Pre-filtering", "Post-filtering")))


colors_stage <- c("Unfiltered" = "#97c1a9", 
                  "Filtered (combined)" = "#FDDD5C",
                  "Filtered (V1)" = "#E26D5C", 
                  "Filtered (Prime)" = "#723D46")

plot_cell_counts <- ggplot(cell_counts, aes(x = Filtering_Stage, y = Count, fill = Filtering_Stage)) +
  geom_bar(stat = "identity", width = 0.5) + 
  geom_text(aes(label = Count), vjust = -0.5, size = 1.5) +
  scale_y_continuous(
    limits = c(0, 1100000),
    labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
  ) +
  scale_fill_manual(values = colors_stage) + 
  facet_grid(. ~ groupings, space = "free", scales = "free_x", switch = "x") +
  theme_classic() + 
  labs(y = "Cell counts", title = " ", fill = " ", x = " ") +
  guides(fill = guide_legend(override.aes = list(size = 2),
                             title.position = "right", nrow = 1,
                             label.theme = element_text(size = 5),
                             by_row = TRUE)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 7),
        axis.title.y = element_text(color = "black", face = "bold", size = 8),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 8, face = "bold"),
        axis.ticks.y = element_line(color = "black"),
        plot.title = element_blank(),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(color = "black", size = 6),
        legend.position = "top")

pdf("/home/smallapragada/v1_5K_panel_comparison_project/bar_chart_pre_post_filtering.pdf", width = 4, height = 3.1)
plot_cell_counts
dev.off()

## Figure 2D - Heatmap of overlapping genes

seurat_obj_combined <- readRDS("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_filtered_joined_2025_11_25.rds")

## Keeping all cells

# Split seurat objects based on panel
genes <- rownames(seurat_obj_combined@assays$RNA@counts)

string_a <- "-v1-prime" 
string_b <- "-prime-with-v1-seg" 

genes_v1_prime <- grep(pattern = string_a, x = genes, value = TRUE)
genes_prime_with_v1_seg <- grep(pattern = string_b, x = genes, value = TRUE)

seurat_obj_v1_prime <- seurat_obj_combined[genes_v1_prime, ]
seurat_obj_prime_with_v1_seg <- seurat_obj_combined[genes_prime_with_v1_seg, ]

# Ensure that V1 gene names don't have that string 
v1_gene_names <- rownames(seurat_obj_v1_prime@assays$RNA@counts)

updated_v1_gene_names <- gsub(
  pattern = string_a,
  replacement = "",
  x = v1_gene_names
)

rownames(seurat_obj_v1_prime@assays$RNA@counts) <- updated_v1_gene_names

if (nrow(seurat_obj_v1_prime@assays$RNA@data) > 0) {
  rownames(seurat_obj_v1_prime@assays$RNA@data) <- updated_v1_gene_names
}

if (nrow(seurat_obj_v1_prime@assays$RNA@scale.data) > 0) {
  rownames(seurat_obj_v1_prime@assays$RNA@scale.data) <- updated_v1_gene_names
}

# Ensure that prime gene names don't have that string 
prime_gene_names <- rownames(seurat_obj_prime_with_v1_seg@assays$RNA@counts)

updated_prime_gene_names <- gsub(
  pattern = string_b,
  replacement = "",
  x = prime_gene_names
)

rownames(seurat_obj_prime_with_v1_seg@assays$RNA@counts) <- updated_prime_gene_names

if (nrow(seurat_obj_prime_with_v1_seg@assays$RNA@data) > 0) {
  rownames(seurat_obj_prime_with_v1_seg@assays$RNA@data) <- updated_prime_gene_names
}

if (nrow(seurat_obj_prime_with_v1_seg@assays$RNA@scale.data) > 0) {
  rownames(seurat_obj_prime_with_v1_seg@assays$RNA@scale.data) <- updated_prime_gene_names
}

# Subset both objects to overlapping genes
seurat_obj_v1_prime_overlap <- seurat_obj_v1_prime[gene_overlaps, ]
seurat_obj_prime_with_v1_seg_overlap <- seurat_obj_prime_with_v1_seg[gene_overlaps, ]

# Create your matrices for the heatmap
cell_gene_matrix_v1 <- as.matrix(seurat_obj_v1_prime_overlap@assays$RNA@counts)
cell_gene_matrix_prime <- as.matrix(seurat_obj_prime_with_v1_seg_overlap@assays$RNA@counts)

avg_expression_v1 <- rowMeans(cell_gene_matrix_v1)
feature_names_v1 <- rownames(cell_gene_matrix_v1)
avg_expression_prime <- rowMeans(cell_gene_matrix_prime)
feature_names_prime <- rownames(cell_gene_matrix_prime)

cell_gene_matrix_df_exp <- data.frame(
  feature_names_v1 = feature_names_v1,
  V1 = avg_expression_v1,
  feature_names_prime = feature_names_prime,
  Prime = avg_expression_prime
) %>%
  select(feature_names_v1, V1, Prime) %>%
  rename(genes = feature_names_v1)

# Function to calculate the correlation between each gene across both panels
run_cor_test_row <- function(row_cell_gene_matrix_v1, row_cell_gene_matrix_prime) {
  cor_result <- cor.test(row_cell_gene_matrix_v1, row_cell_gene_matrix_prime, method = "spearman")
  return(c(
    r_value = cor_result$estimate,
    p_value = cor_result$p.value
  ))
}

results_matrix <- mapply(
  FUN = run_cor_test_row, 
  as.data.frame(t(cell_gene_matrix_v1)), 
  as.data.frame(t(cell_gene_matrix_prime))
)

cor_mat <- as.data.frame(t(results_matrix)) %>%
  rownames_to_column() %>%
  rename(genes = rowname) %>%
  rename(`Rho-value` = r_value.rho) %>%
  inner_join(., cell_gene_matrix_df_exp) %>%
  column_to_rownames("genes") %>%
  select(-p_value) %>%
  as.matrix()

## Retaining cells that have at least one count or more in both chemistries

pattern <- paste0("\\b(", paste(gene_overlaps, collapse = "|"), ")\\b")
matched_genes <- grep(pattern, rownames(seurat_obj_combined), value = TRUE)

matched_genes <- matched_genes[!grepl("SOX2-OT", matched_genes)]

# Subset the object by these features
seurat_obj_combined_overlaps_both <- seurat_obj_combined[matched_genes, ]
seurat_obj_combined_overlaps_both_counts <- GetAssayData(seurat_obj_combined_overlaps_both, slot = "counts")

names_a <- paste0(gene_overlaps, string_a)
names_b <- paste0(gene_overlaps, string_b)

mat_a <- seurat_obj_combined_overlaps_both_counts[names_a, ]
mat_b <- seurat_obj_combined_overlaps_both_counts[names_b, ]

keep_mask <- (mat_a != 0 & mat_b != 0)

cells_kept_count <- rowSums(keep_mask)

gene_keep_summary <- data.frame(
  gene = gene_overlaps,
  cells_kept = cells_kept_count,
  total_cells = ncol(seurat_obj_combined_overlaps_both_counts),
  stringsAsFactors = FALSE
)

## Plot

# R matrix
r_mat <- cor_mat[, "R-value", drop = FALSE]

# Expression matrix
exp_mat <- cor_mat_all_cells[, c("V1", "Prime")]

# Color schemes
magma_colors <- rev(magma(256))
r_colors <- colorRamp2(
  breaks = seq(0, 1.0, length.out = length(magma_colors)),
  colors = magma_colors
)

gnbu_colors <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
gnbu_colors <- gnbu_colors[-(1:20)]

exp_colors <- colorRamp2(
  breaks = seq(0, 1, length.out = length(gnbu_colors)),
  colors = gnbu_colors
)

ht_r <- Heatmap(
  r_mat,
  name = "Rho value",      
  col = r_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_side = "top",
  show_row_names = FALSE,
  column_names_rot = 0,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 7, fontface = "bold"),
  show_heatmap_legend = FALSE
)

ht_exp <- Heatmap(
  exp_mat,
  name = "Expression",          
  col = exp_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,            
  column_names_rot = 0,
  column_names_side = "top",
  column_names_gp = gpar(fontsize = 7, fontface = "bold"),
  show_heatmap_legend = FALSE
)

final_heatmap <- ht_exp + ht_r

pdf("/home/smallapragada/v1_5K_panel_comparison_project/heatmap_r_val_avg_exp.pdf", width = 3, height = 3.6)
draw(final_heatmap)
dev.off()

ht_r <- Heatmap(
  r_mat,
  name = "R value",      
  col = r_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_side = "top",
  show_row_names = TRUE,
  row_names_side = "left",
  column_names_rot = 0,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 7, fontface = "bold"),
  heatmap_legend_param = list(
    title = "R-value",
    legend_width = unit(1.5, "cm"),
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5)
  )
)

ht_exp <- Heatmap(
  exp_mat,
  name = "Expression",          
  col = exp_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,            
  column_names_rot = 0,
  column_names_side = "top",
  column_names_gp = gpar(fontsize = 7, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Average expression",
    legend_width = unit(1.5, "cm"),
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5)
  )
)

final_heatmap <- ht_r + ht_exp

pdf("/home/smallapragada/v1_5K_panel_comparison_project/heatmap_r_val_avg_exp_label.pdf", width = 3, height = 5)
draw(final_heatmap)
dev.off()

### Figure 2E - secreted genes histogram

## Importing dual slide and splitting based on panel genes
dual_obj_combined <- readRDS("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_comparison_allgenes_20pcs_2025_12_02.rds")

dual_obj_combined_v1 <- dual_obj_combined[genes_v1_prime, ]
dual_obj_combined_prime <- dual_obj_combined[genes_prime_with_v1_seg, ]

# Importing secreted gene lists
v1_secreted_genes <- read.delim("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/saved_csvs/Secreted_v1.csv", 
                                sep = ",", header = FALSE)
prime_secreted_genes <- read.delim("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/saved_csvs/Secreted_5K.csv", 
                                   sep = ",", header = FALSE)

## Stats about secreted genes
summary(dual_obj_combined_v1$nCount_RNA)
summary(dual_obj_combined_v1$nCount_RNA)

## Subsetting the combined seurat object into V1 vs Prime secreted genes
v1_secreted_genes_full_names <- paste0(v1_secreted_genes$V1, "-v1-prime")
matching_genes <- v1_secreted_genes_full_names[v1_secreted_genes_full_names %in% rownames(dual_obj_combined_v1)]
seurat_obj_v1_secreted_genes <- dual_obj_combined_v1[matching_genes, ]

prime_secreted_genes_full_names <- paste0(prime_secreted_genes$V1, "-prime-with-v1-seg")
matching_genes <- prime_secreted_genes_full_names[prime_secreted_genes_full_names %in% rownames(dual_obj_combined_prime)]
seurat_obj_prime_secreted_genes <- dual_obj_combined_prime[matching_genes, ]

## Calculating nFeature per each panel secreted genes
v1_counts_matrix <- GetAssayData(
  object = seurat_obj_v1_secreted_genes, 
  assay = "RNA", 
  slot = "counts"
)

v1_n_feature_per_cell <- colSums(v1_counts_matrix > 0)

v1_secreted_seurat_obj_final <- AddMetaData(
  object = seurat_obj_v1_secreted_genes,
  metadata = v1_n_feature_per_cell,
  col.name = "nFeature_RNA_v1_secreted" 
)

prime_counts_matrix <- GetAssayData(
  object = seurat_obj_prime_secreted_genes, 
  assay = "RNA", 
  slot = "counts"
)

prime_n_feature_per_cell <- colSums(prime_counts_matrix > 0)

prime_secreted_seurat_obj_final <- AddMetaData(
  object = seurat_obj_prime_secreted_genes,
  metadata = prime_n_feature_per_cell,
  col.name = "nFeature_RNA_prime_secreted" 
)

v1_data <- data.frame(
  nFeature_RNA = v1_secreted_seurat_obj_final$nFeature_RNA_v1_secreted,
  Sample = "V1 Secreted"
)

prime_data <- data.frame(
  nFeature_RNA = prime_secreted_seurat_obj_final$nFeature_RNA_prime_secreted,
  Sample = "Prime Secreted"
)

# Combine the data frames
combined_data <- bind_rows(v1_data, prime_data)

density_plot <- ggplot(combined_data, aes(x = nFeature_RNA, fill = Sample)) +
  geom_density(
    aes(y = after_stat(count)),
    alpha = 0.5,          
    color = "black", 
    bw = 0.75
  ) +
  xlim(0, 40) + 
  scale_fill_manual(values = c("V1 Secreted" = "purple4", "Prime Secreted" = "darkorange")) +
  labs(
    x = "Unique genes per cell",
    y = "Cell count", 
    fill = " "
  ) +
  theme_classic() +
  theme(
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(color = "black", size = 7, face = "bold"),
    axis.ticks.x = element_line(color = "black"), 
    axis.text.x = element_text(color = "black", size = 6), 
    axis.title.x = element_text(color = "black", size = 7, face = "bold"),
    legend.key.size = unit(2.5, "mm"),
    legend.text = element_text(color = "black", size = 6)
  )

pdf("/home/smallapragada/v1_5K_panel_comparison_project/density_plot.pdf", width = 3.7, height = 2.5)
density_plot
dev.off()

#### FIGURE 3 ----

### Figure 3A - Triple UMAP across all panels

## Mapping cell ids to each object

dual_obj_combined <- readRDS("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_comparison_allgenes_20pcs_2025_12_02.rds")
dual_obj_v1 <- readRDS("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_comparison_v1genes_20pcs_2025_12_02.rds")
dual_obj_prime <- readRDS("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_prime_primegenes_combined_20pcs_2025_12_02.rds")

# Metadata
dual_obj_combined_meta <- dual_obj_combined@meta.data
dual_obj_v1_meta <- dual_obj_v1@meta.data
dual_obj_prime_meta <- dual_obj_prime@meta.data

# Pulling cell clusters from each object
combined_df_labels <- dual_obj_combined@meta.data %>%
  rownames_to_column(var = "cell_ids_rownames") %>%
  select(cell_ids_rownames, matched_cell_ids, leiden_0.5_combined, leiden_1.0_combined) 
  
v1_df_labels <- dual_obj_v1@meta.data %>% 
  rownames_to_column(var = "cell_ids_rownames") %>%
  select(cell_ids_rownames, matched_cell_ids, leiden_0.5_v1, leiden_1.0_v1) 
  
prime_df_labels <- dual_obj_prime@meta.data %>%
  rownames_to_column(var = "cell_ids_rownames") %>%
  select(cell_ids_rownames, matched_cell_ids, leiden_0.5_prime, leiden_1.0_prime) 

## V1 metadata 
dual_obj_v1_meta_start <- dual_obj_v1@meta.data %>%
  rownames_to_column(var = "cell_ids_rownames_initial")

dual_obj_v1_meta_prime <- merge(
  x = dual_obj_v1_meta_start,
  y = prime_df_labels,
  by = "matched_cell_ids",
  all.x = TRUE
)

dual_obj_v1_meta_prime_combo <- merge(
  x = dual_obj_v1_meta_prime,
  y = combined_df_labels,
  by = "matched_cell_ids",
  all.x = TRUE
) %>%
  column_to_rownames(var = "cell_ids_rownames_initial")

# Remove the extraneous cell_ids_rownames columns generated by the merges
cols_to_remove <- c("cell_ids_rownames.x", "cell_ids_rownames.y", "cell_ids_rownames") 
cols_to_remove_final <- intersect(cols_to_remove, colnames(dual_obj_v1_meta_prime_combo))

dual_obj_v1_meta_prime_combo <- dual_obj_v1_meta_prime_combo %>%
  select(-one_of(cols_to_remove_final))

# Assign back to the Seurat object
dual_obj_v1@meta.data <- dual_obj_v1_meta_prime_combo

## Combined metadata 
dual_obj_combined_meta_start <- dual_obj_combined@meta.data %>%
  rownames_to_column(var = "cell_ids_rownames_initial")

dual_obj_combined_meta_v1 <- merge(
  x = dual_obj_combined_meta_start,
  y = v1_df_labels,
  by = "matched_cell_ids",
  all.x = TRUE
)

dual_obj_combined_meta_v1_prime <- merge(
  x = dual_obj_combined_meta_v1,
  y = prime_df_labels,
  by = "matched_cell_ids",
  all.x = TRUE
) %>%
  column_to_rownames(var = "cell_ids_rownames_initial")

# Remove the extraneous cell_ids_rownames columns generated by the merges
cols_to_remove <- c("cell_ids_rownames.x", "cell_ids_rownames.y", "cell_ids_rownames") 
cols_to_remove_final <- intersect(cols_to_remove, colnames(dual_obj_combined_meta_v1_prime))

dual_obj_combined_meta_v1_prime <- dual_obj_combined_meta_v1_prime %>%
  select(-one_of(cols_to_remove_final))

# Assign back to the Seurat object
dual_obj_combined@meta.data <- dual_obj_combined_meta_v1_prime

## Prime metadata 
dual_obj_prime_meta_start <- dual_obj_prime@meta.data %>%
  rownames_to_column(var = "cell_ids_rownames_initial")

dual_obj_prime_meta_v1 <- merge(
  x = dual_obj_prime_meta_start,
  y = v1_df_labels,
  by = "matched_cell_ids",
  all.x = TRUE
)

dual_obj_prime_meta_v1_combined <- merge(
  x = dual_obj_prime_meta_v1,
  y = combined_df_labels,
  by = "matched_cell_ids",
  all.x = TRUE
) %>%
  column_to_rownames(var = "cell_ids_rownames_initial")

# Remove the extraneous cell_ids_rownames columns generated by the merges
cols_to_remove <- c("cell_ids_rownames.x", "cell_ids_rownames.y", "cell_ids_rownames") 
cols_to_remove_final <- intersect(cols_to_remove, colnames(dual_obj_prime_meta_v1_combined))

dual_obj_prime_meta_v1_combined <- dual_obj_prime_meta_v1_combined %>%
  select(-one_of(cols_to_remove_final))

# Assign back to the Seurat object
dual_obj_prime@meta.data <- dual_obj_prime_meta_v1_combined

### Figure 4A - Triple UMAP

# Color palettes
get_clean_palette <- function(n) {
  cols <- randomcoloR::distinctColorPalette(n + 20)
  is_grey <- function(hex) {
    rgb_val <- col2rgb(hex)
    return(sd(rgb_val) < 15) 
  }
  clean_cols <- cols[!sapply(cols, is_grey)]
  return(clean_cols[1:n])
}

# Generate palettes
v1_colors <- get_clean_palette(10)
prime_colors <- get_clean_palette(10)
combo_colors <- get_clean_palette(14)

# Set the NA color for UMAPS
na_grey <- "grey40" 

# V1 plot
plot_v1_v1 <- DimPlot(dual_obj_v1, 
                      reduction = "X_umap", 
                      group.by = "leiden_0.5_v1", 
                      raster = T, 
                      cols = v1_colors) +
  coord_equal() + 
  NoLegend() + 
  theme_void() +
  theme(plot.title = element_blank()) +
  scale_color_manual(values = v1_colors, na.value = na_grey)

plot_v1_prime <- DimPlot(dual_obj_v1, 
                         reduction = "X_umap", 
                         group.by = "leiden_0.5_prime", 
                         raster = T, 
                         cols = prime_colors) +
  coord_equal() + 
  NoLegend() + 
  theme_void() +
  theme(plot.title = element_blank()) +
  scale_color_manual(values = prime_colors, na.value = na_grey)

plot_v1_combo <- DimPlot(dual_obj_v1, 
                         reduction = "X_umap", 
                         group.by = "leiden_0.5_combined", 
                         raster = T, 
                         cols = combo_colors) +
  coord_equal() + 
  NoLegend() + 
  theme_void() +
  theme(plot.title = element_blank()) +
  scale_color_manual(values = combo_colors, na.value = na_grey)

plot_v1 <- plot_v1_v1 / plot_v1_prime / plot_v1_combo

pdf("/home/smallapragada/v1_5K_panel_comparison_project/umap_v1_all_projections.pdf", width = 2.5, height = 6)
print(plot_v1)
dev.off()

# Prime plots 
plot_prime_v1 <- DimPlot(dual_obj_prime, 
                         reduction = "X_umap", 
                         group.by = "leiden_0.5_v1", 
                         raster = T, 
                         cols = v1_colors) +
  coord_equal() + 
  NoLegend() + 
  theme_void() + 
  theme(plot.title = element_blank()) +
  scale_color_manual(values = v1_colors, na.value = na_grey)

plot_prime_prime <- DimPlot(dual_obj_prime, 
                            reduction = "X_umap", 
                            group.by = "leiden_0.5_prime", 
                            raster = T, 
                            cols = prime_colors) +
  coord_equal() + 
  NoLegend() + 
  theme_void() + 
  theme(plot.title = element_blank()) +
  scale_color_manual(values = prime_colors, na.value = na_grey)

plot_prime_combo <- DimPlot(dual_obj_prime, 
                            reduction = "X_umap", 
                            group.by = "leiden_0.5_combined", 
                            raster = T, 
                            cols = combo_colors) +
  coord_equal() + 
  NoLegend() + 
  theme_void() + 
  theme(plot.title = element_blank()) +
  scale_color_manual(values = combo_colors, na.value = na_grey)

plot_prime <- plot_prime_v1 / plot_prime_prime / plot_prime_combo

pdf("/home/smallapragada/v1_5K_panel_comparison_project/umap_prime_all_projections.pdf", width = 2.5, height = 6)
print(plot_prime)
dev.off()

# Combo plots 
plot_combo_v1 <- DimPlot(dual_obj_combined, 
                         reduction = "X_umap", 
                         group.by = "leiden_0.5_v1", 
                         raster = T, 
                         cols = v1_colors) +
  coord_equal() + 
  NoLegend() + 
  theme_void() + 
  theme(plot.title = element_blank()) +
  scale_color_manual(values = v1_colors, na.value = na_grey)

plot_combo_prime <- DimPlot(dual_obj_combined, 
                            reduction = "X_umap", 
                            group.by = "leiden_0.5_prime", 
                            raster = T, 
                            cols = prime_colors) +
  coord_equal() + 
  NoLegend() + 
  theme_void() + 
  theme(plot.title = element_blank()) +
  scale_color_manual(values = prime_colors, na.value = na_grey)

plot_combo_combo <- DimPlot(dual_obj_combined, 
                            reduction = "X_umap", 
                            group.by = "leiden_0.5_combined", 
                            raster = T, 
                            cols = combo_colors) +
  coord_equal() + 
  NoLegend() + 
  theme_void() + 
  theme(plot.title = element_blank()) +
  scale_color_manual(values = combo_colors, na.value = na_grey)

plot_combo <- plot_combo_v1 / plot_combo_prime / plot_combo_combo

pdf("/home/smallapragada/v1_5K_panel_comparison_project/umap_combo_all_projections.pdf", width = 2.5, height = 6)
print(plot_combo)
dev.off()

### Figure 3B - alluvial 

combo_origin_leiden <- dual_obj_combined_meta_v1_prime %>%
  select(leiden_0.5_combined, leiden_0.5_v1, leiden_0.5_prime)

clustree_df <- combo_origin_leiden %>%
  mutate(
    across(
      .cols = everything(),
      .fns = ~factor(replace_na(as.character(.), "NA"))
    )
  )

colnames(clustree_df) <- c("leiden_1", "leiden_2", "leiden_3")

custom_strata_order <- c(as.character(0:13), "NA") 

clustree_df <- clustree_df %>%
  mutate(
    leiden_1 = factor(leiden_1, levels = custom_strata_order),
    leiden_2 = factor(leiden_2, levels = custom_strata_order),
    leiden_3 = factor(leiden_3, levels = custom_strata_order)
  )

alluvial_data <- clustree_df %>%
  group_by(leiden_1, leiden_2, leiden_3) %>%
  summarise(Freq = n(), .groups = "drop")

alluvial <- ggplot(alluvial_data,
                   aes(axis1 = leiden_1, axis2 = leiden_2, axis3 = leiden_3, y = Freq)) +
  geom_alluvium(aes(fill = leiden_1),
                width = 0.3,
                curve_type = "sigmoid",
                alpha = 1) +
  geom_stratum(width = 0.3, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
  scale_fill_manual(
    values = combo_colors,
    name = "Cluster label"
  ) +
  theme_classic() +
  labs(y = "Number of cells", x = "") +
  theme(
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 9, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

pdf("/home/smallapragada/v1_5K_panel_comparison_project/combined_alluvial.pdf", width = 7, height = 4)
alluvial
dev.off()

#### SUPPLEMENTAL MATERIAL ----

### Figure S1 - Correlation scatterplots per sample (V1_solo vs. V1_Prime)

merged_obj_unfiltered_v1 <- readRDS("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/v1_slides_2025_11_20.rds")
merged_obj_unfiltered_prime <- readRDS("/scratch/smallapragada/bbl_project/v1_prime5k_comparison/objects/prime_slides_2025_11_20.rds")

## V1 solo and dual counts
merged_obj_unfiltered_v1_split <- SplitObject(merged_obj_unfiltered_v1, split.by = "slide_type")

v1_solo_counts <- AggregateExpression(
  merged_obj_unfiltered_v1_split$v1_solo,
  assays = "RNA",
  slot = "counts",
  group.by = "sample",
  return.seurat = FALSE
)$RNA
v1_solo_counts <- v1_solo_counts[, colnames(v1_solo_counts) != "dropped"]

v1_dual_counts <- AggregateExpression(
  merged_obj_unfiltered_v1_split$v1_prime,
  assays = "RNA",
  slot = "counts",
  group.by = "sample",
  return.seurat = FALSE
)$RNA
v1_dual_counts <- v1_dual_counts[, colnames(v1_dual_counts) != "dropped"]

sample_names <- unique(colnames(v1_dual_counts))

# Merging dfs together
v1_solo_counts_renamed <- as.data.frame(v1_solo_counts) %>%
  rownames_to_column(var = "gene") %>%
  rename_with(~ paste0(., "_solo"), .cols = -gene)

v1_dual_counts_renamed <- as.data.frame(v1_dual_counts) %>%
  rownames_to_column(var = "gene") %>%
  rename_with(~ paste0(., "_dual"), .cols = -gene)

merged_counts_v1 <- full_join(
  v1_solo_counts_renamed,
  v1_dual_counts_renamed,
  by = "gene"
)

# Creating list of dataframes per sample
v1_r_df <- data.frame(sample = character(), R = numeric(), stringsAsFactors = FALSE)
sample_list_of_dfs <- list()

for (sample in sample_names) {
  
  solo_col <- paste0(sample, "_solo")
  dual_col <- paste0(sample, "_dual")
  
  sample_df <- merged_counts_v1 %>%
    select(gene, !!solo_col, !!dual_col) %>%
    rename(solo_count = !!solo_col, dual_count = !!dual_col)
  
  # stats
  cor_test_result <- cor.test(sample_df$solo_count, sample_df$dual_count, method = "pearson")
  
  # Append results to the data frame
  v1_r_df <- v1_r_df %>%
    add_row(
      sample = sample,
      R = cor_test_result$estimate
    )
  
  sample_df <- sample_df %>%
    mutate(gene_type = case_when(gene %in% gene_overlaps ~ "Overlapping gene",
                                 TRUE ~ "V1 panel gene"))
  
  sample_list_of_dfs[[sample]] <- sample_df
  
  v1_plot <- ggplot(sample_df, aes(x = dual_count, y = solo_count, color = gene_type, shape = gene_type)) +
    geom_point(alpha = 1, size = 1.5) +
    theme_classic() +
    labs(
      title = paste(sample), 
      x = "Total transcripts (V1 + Prime)",
      y = "Total transcripts (V1 only)"
    ) +
    scale_x_continuous(
      limits = c(0, 650000),
      labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
    ) +
    scale_y_continuous(
      limits = c(0, 650000),
      labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(
      values = c("Overlapping gene" = "#FFB000", "V1 panel gene" = "#DC267F")
    ) +
    scale_shape_manual(
      values = c("Overlapping gene" = 17, "V1 panel gene" = 15)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 1, shape = c(17, 15))),
      shape = "none"
    ) +
    theme(
      plot.title = element_text(color = "black", face = "bold", size = 8, hjust = 0.5),
      axis.ticks.y = element_line(color = "black"),
      axis.text.y = element_text(color = "black", size = 6),
      axis.title.y = element_blank(),
      axis.ticks.x = element_line(color = "black"), 
      axis.text.x = element_text(color = "black", size = 6), 
      axis.title.x = element_blank(), 
      legend.text = element_text(color = "black", size = 7),
      legend.position = "none",
      legend.title = element_blank()
    )
  
  # Save the Plot to a PDF File
  pdf_filename <- paste0("/home/smallapragada/v1_5K_panel_comparison_project/v1_correlation_all_samples_supps/", 
                         sample, "_v1_panel_scatterplot_supp.pdf")
  
  pdf(pdf_filename, width = 1.75, height = 1.75)
  print(v1_plot)
  dev.off()
}

### Figure S2 - Correlation scatterplots per sample (Prime_solo vs. Prime_V1)
merged_obj_unfiltered_prime_split <- SplitObject(merged_obj_unfiltered_prime, split.by = "slide_type")

prime_solo_counts <- AggregateExpression(
  merged_obj_unfiltered_prime_split$prime_solo,
  assays = "RNA",
  slot = "counts",
  group.by = "sample",
  return.seurat = FALSE
)$RNA
prime_solo_counts <- prime_solo_counts[, colnames(prime_solo_counts) != "dropped"]

prime_dual_counts <- AggregateExpression(
  merged_obj_unfiltered_prime_split$prime_v1,
  assays = "RNA",
  slot = "counts",
  group.by = "sample",
  return.seurat = FALSE
)$RNA
prime_dual_counts <- prime_dual_counts[, colnames(prime_dual_counts) != "dropped"]

# Merging dfs together
prime_solo_counts_renamed <- as.data.frame(prime_solo_counts) %>%
  rownames_to_column(var = "gene") %>%
  rename_with(~ paste0(., "_solo"), .cols = -gene)

prime_dual_counts_renamed <- as.data.frame(prime_dual_counts) %>%
  rownames_to_column(var = "gene") %>%
  rename_with(~ paste0(., "_dual"), .cols = -gene)

merged_counts_prime <- full_join(
  prime_solo_counts_renamed,
  prime_dual_counts_renamed,
  by = "gene"
)

# Creating list of dataframes per sample
prime_r_df <- data.frame(sample = character(), R = numeric(), stringsAsFactors = FALSE)
sample_list_of_dfs <- list()

for (sample in sample_names) {
  
  solo_col <- paste0(sample, "_solo")
  dual_col <- paste0(sample, "_dual")
  
  sample_df <- merged_counts_prime %>%
    select(gene, !!solo_col, !!dual_col) %>%
    rename(solo_count = !!solo_col, dual_count = !!dual_col)
  
  # stats
  cor_test_result <- cor.test(sample_df$solo_count, sample_df$dual_count, method = "pearson")
  
  # Append results to the data frame
  prime_r_df <- prime_r_df %>%
    add_row(
      sample = sample,
      R = cor_test_result$estimate
    )
  
  sample_df <- sample_df %>%
    mutate(gene_type = case_when(gene %in% gene_overlaps ~ "Overlapping gene",
                                 TRUE ~ "Prime panel gene"))
  
  sample_list_of_dfs[[sample]] <- sample_df
  
  prime_plot <- ggplot(sample_df, aes(x = dual_count, y = solo_count, color = gene_type, shape = gene_type)) +
    geom_point(alpha = 1, size = 1.5) +
    theme_classic() +
    labs(
      title = paste(sample), 
      x = "Total transcripts (V1 + Prime)",
      y = "Total transcripts (V1 only)"
    ) +
    scale_x_continuous(
      limits = c(0, 100000),
      labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
    ) +
    scale_y_continuous(
      limits = c(0, 100000),
      labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(
      values = c("Overlapping gene" = "#FFB000", "Prime panel gene" = "#648FFF")
    ) +
    scale_shape_manual(
      values = c("Overlapping gene" = 17, "Prime panel gene" = 15)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 1, shape = c(17, 15))),
      shape = "none"
    ) +
    theme(
      plot.title = element_text(color = "black", face = "bold", size = 8, hjust = 0.5),
      axis.ticks.y = element_line(color = "black"),
      axis.text.y = element_text(color = "black", size = 6),
      axis.title.y = element_blank(),
      axis.ticks.x = element_line(color = "black"), 
      axis.text.x = element_text(color = "black", size = 6), 
      axis.title.x = element_blank(), 
      legend.text = element_text(color = "black", size = 7),
      legend.position = "none",
      legend.title = element_blank()
    )
  
  # Save the Plot to a PDF File
  pdf_filename <- paste0("/home/smallapragada/v1_5K_panel_comparison_project/prime_correlation_all_samples_supps/", 
                         sample, "_prime_panel_scatterplot_supp.pdf")
  
  pdf(pdf_filename, width = 1.75, height = 1.75)
  print(prime_plot)
  dev.off()
}

### Table 1 - R value + p value per sample correlation stats

final_r_df <- bind_rows(
  v1_r_df %>% mutate(Panel = "V1"),
  prime_r_df %>% mutate(Panel = "Prime")
) %>%
  rename(Sample = sample) %>%
  rename(`R-value` = R)

write_csv(final_r_df, '/home/smallapragada/v1_5K_panel_comparison_project/r_correlation_all_samples_supp.csv')

### Figure S3 - violin plots of nCount and nFeature split by sample

plot_data_v1 <- FetchData(
  object = merged_obj_unfiltered_v1, 
  vars = c("slide_type", "sample", "nCount_RNA", "nFeature_RNA"),
  slot = "counts" 
) %>%
  filter(sample != "dropped")

plot_data_v1$slide_type <- factor(
  plot_data_v1$slide_type, 
  level = c("v1_solo", "v1_prime"),
  labels = c("V1", "Dual")
)

count_plot_v1 <- ggplot(plot_data_v1, aes(x = sample, y = nCount_RNA, fill = slide_type)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = 1, scale = "width") +
  theme_classic() + 
  scale_y_continuous(breaks = pretty_breaks(n = 3)) + 
  scale_fill_manual(values = c("V1" = "#FE6100", "Dual" = "#662D91")) +
  labs(
    y = "Total transcripts",
    title = "V1 panel",
    fill = " "
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black"), 
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(color = "black", size = 7, face = "bold"),
    strip.text = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 8, color = "black", hjust = 0.5),
    legend.position = "bottom", 
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6)
  )

pdf("/home/smallapragada/v1_5K_panel_comparison_project/v1_violin_ncount_sample_supp.pdf", width = 9, height = 4)
count_plot_v1
dev.off()

feature_plot_v1 <- ggplot(plot_data_v1, aes(x = sample, y = nFeature_RNA, fill = slide_type)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = 1, scale = "width") +
  theme_classic() + 
  scale_y_continuous(breaks = pretty_breaks(n = 3)) + 
  scale_fill_manual(values = c("V1" = "#FE6100", "Dual" = "#662D91")) +
  labs(
    y = "Total unique transcripts",
    title = "V1 panel",
    fill = " "
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black"), 
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(color = "black", size = 7, face = "bold"),
    strip.text = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 8, color = "black", hjust = 0.5),
    legend.position = "bottom", 
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6)
  )

pdf("/home/smallapragada/v1_5K_panel_comparison_project/v1_violin_nfeature_sample_supp.pdf", width = 9, height = 4)
feature_plot_v1
dev.off()

## Figure S4 - nCount and nFeature split by sample (prime)

plot_data_prime <- FetchData(
  object = merged_obj_unfiltered_prime, 
  vars = c("slide_type", "sample", "nCount_RNA", "nFeature_RNA"),
  slot = "counts" 
) %>%
  filter(sample != "dropped")

plot_data_prime$slide_type <- factor(
  plot_data_prime$slide_type, 
  level = c("prime_solo", "prime_v1"),
  labels = c("Prime", "Dual")
)

count_plot_prime <- ggplot(plot_data_prime, aes(x = sample, y = nCount_RNA, fill = slide_type)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = 1, scale = "width") +
  theme_classic() + 
  scale_y_continuous(breaks = pretty_breaks(n = 3)) + 
  scale_fill_manual(values = c("Prime" = "#2E6F40", "Dual" = "#662D91")) +
  labs(
    y = "Total transcripts",
    title = "Prime panel",
    fill = " "
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black"), 
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(color = "black", size = 7, face = "bold"),
    strip.text = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 8, color = "black", hjust = 0.5),
    legend.position = "bottom", 
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6)
  )

pdf("/home/smallapragada/v1_5K_panel_comparison_project/prime_violin_ncount_sample_supp.pdf", width = 9, height = 4)
count_plot_prime
dev.off()

feature_plot_prime <- ggplot(plot_data_prime, aes(x = sample, y = nFeature_RNA, fill = slide_type)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = 1, scale = "width") +
  theme_classic() + 
  scale_y_continuous(breaks = pretty_breaks(n = 3)) + 
  scale_fill_manual(values = c("Prime" = "#2E6F40", "Dual" = "#662D91")) +
  labs(
    y = "Total unique transcripts",
    title = "Prime panel",
    fill = " "
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black"), 
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(color = "black", size = 7, face = "bold"),
    strip.text = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 8, color = "black", hjust = 0.5),
    legend.position = "bottom", 
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6)
  )

pdf("/home/smallapragada/v1_5K_panel_comparison_project/prime_violin_nfeature_sample_supp.pdf", width = 9, height = 4)
feature_plot_prime
dev.off()

### Figure S6 - Correlation scatterplots of V1 solo vs prime solo & V1 Dual vs prime dual for overlapped genes

# V1_solo vs. Prime_solo
transcript_count_files$prime_solo$gene_count_filtered <- as.numeric(transcript_count_files$prime_solo$gene_count_filtered)
transcript_count_files$v1_solo$gene_count_filtered <- as.numeric(transcript_count_files$v1_solo$gene_count_filtered)

merged_df_solo <- inner_join(
  tibble(feature_name = rownames(transcript_count_files$prime_solo), x = transcript_count_files$prime_solo$gene_count_filtered),
  tibble(feature_name = rownames(transcript_count_files$v1_solo), y = transcript_count_files$v1_solo$gene_count_filtered),
  by = "feature_name"
)

solo_plot <- ggplot(merged_df_solo, aes(x = x, y = y)) +
  geom_point(alpha = 1, size = 1.0, color = "blue") +
  theme_classic() +
  labs(
    x = "Total transcripts (Prime solo)",
    y = "Total transcripts (V1 solo)"
  ) +
  scale_x_continuous(
    limits = c(0, 900000), 
    labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
  ) +
  scale_y_continuous(
    limits = c(0, 900000),
    labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") 

# V1_prime vs. Prime_V1
transcript_count_files$prime_v1$gene_count_filtered <- as.numeric(transcript_count_files$prime_v1$gene_count_filtered)
transcript_count_files$v1_prime$gene_count_filtered <- as.numeric(transcript_count_files$v1_prime$gene_count_filtered)

merged_df_dual <- inner_join(
  tibble(feature_name = rownames(transcript_count_files$prime_v1), x = transcript_count_files$prime_v1$gene_count_filtered),
  tibble(feature_name = rownames(transcript_count_files$v1_prime), y = transcript_count_files$v1_prime$gene_count_filtered),
  by = "feature_name"
)

dual_plot <- ggplot(merged_df_dual, aes(x = x, y = y)) +
  geom_point(alpha = 1, size = 1.0, color = "red") +
  theme_classic() +
  labs(
    x = "Total transcripts (Prime dual)",
    y = "Total transcripts (V1 dual)"
  ) +
  scale_x_continuous(
    limits = c(0, 900000), 
    labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
  ) +
  scale_y_continuous(
    limits = c(0, 900000),
    labels = label_number(suffix = "K", scale = 1e-3, big.mark = "")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") 

plot <- solo_plot + dual_plot

pdf("/home/smallapragada/v1_5K_panel_comparison_project/supp_scatterplot_overlap_solo_dual.pdf", width = 7.5, height = 3.5)
plot
dev.off()

### Figure S7 - V1 - Prime - Combined + Prime - V1 - Combined alluvial

prime_origin_leiden <- dual_obj_prime_meta_v1_combined %>%
  select(leiden_0.5_prime, leiden_0.5_v1, leiden_0.5_combined)

clustree_df <- prime_origin_leiden %>%
  mutate(
    across(
      .cols = everything(),
      .fns = ~factor(replace_na(as.character(.), "NA"))
    )
  )

colnames(clustree_df) <- c("leiden_1", "leiden_2", "leiden_3")

custom_strata_order <- c(as.character(0:13), "NA") 

clustree_df <- clustree_df %>%
  mutate(
    leiden_1 = factor(leiden_1, levels = custom_strata_order),
    leiden_2 = factor(leiden_2, levels = custom_strata_order),
    leiden_3 = factor(leiden_3, levels = custom_strata_order)
  )

cluster_colors <- c(
  "#E6194B",  
  "#3CB44B", 
  "#FFE119",  
  "#4363D8", 
  "#F58231", 
  "#911EB4", 
  "#46F0F0", 
  "#F032E6", 
  "#BCF60C",  
  "#FABEBE", 
  "#008080", 
  "#E6BEFF",  
  "#AA6E28", 
  "#800000" 
)

alluvial_data <- clustree_df %>%
  group_by(leiden_1, leiden_2, leiden_3) %>%
  summarise(Freq = n(), .groups = "drop")

alluvial_prime <- ggplot(alluvial_data,
                   aes(axis1 = leiden_1, axis2 = leiden_2, axis3 = leiden_3, y = Freq)) +
  geom_alluvium(aes(fill = leiden_1),
                width = 0.3,
                curve_type = "sigmoid") +
  geom_stratum(width = 0.3, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
  scale_fill_manual(
    values = cluster_colors,
    name = "Cluster label"
  ) +
  theme_classic() +
  labs(y = "Number of cells", x = "") +
  theme(
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 9, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

v1_origin_leiden <- dual_obj_v1_meta_prime_combo %>%
  select(leiden_0.5_v1, leiden_0.5_prime, leiden_0.5_combined)

clustree_df <- v1_origin_leiden %>%
  mutate(
    across(
      .cols = everything(),
      .fns = ~factor(replace_na(as.character(.), "NA"))
    )
  )

colnames(clustree_df) <- c("leiden_1", "leiden_2", "leiden_3")

custom_strata_order <- c(as.character(0:13), "NA") 

clustree_df <- clustree_df %>%
  mutate(
    leiden_1 = factor(leiden_1, levels = custom_strata_order),
    leiden_2 = factor(leiden_2, levels = custom_strata_order),
    leiden_3 = factor(leiden_3, levels = custom_strata_order)
  )

cluster_colors <- c(
  "#E6194B",  
  "#3CB44B", 
  "#FFE119",  
  "#4363D8", 
  "#F58231", 
  "#911EB4", 
  "#46F0F0", 
  "#F032E6", 
  "#BCF60C",  
  "#FABEBE", 
  "#008080", 
  "#E6BEFF",  
  "#AA6E28", 
  "#800000" 
)

alluvial_data <- clustree_df %>%
  group_by(leiden_1, leiden_2, leiden_3) %>%
  summarise(Freq = n(), .groups = "drop")

alluvial_v1 <- ggplot(alluvial_data,
                         aes(axis1 = leiden_1, axis2 = leiden_2, axis3 = leiden_3, y = Freq)) +
  geom_alluvium(aes(fill = leiden_1),
                width = 0.3,
                curve_type = "sigmoid") +
  geom_stratum(width = 0.3, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
  scale_fill_manual(
    values = cluster_colors,
    name = "Cluster label"
  ) +
  theme_classic() +
  labs(y = "Number of cells", x = "") +
  theme(
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 9, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

alluvial <- alluvial_v1 / alluvial_prime

pdf("/home/smallapragada/v1_5K_panel_comparison_project/prime_v1_alluvial_supp.pdf", width = 7, height = 8)
alluvial
dev.off()
