# Date last modified: May 31, 2024 by Zhining Sui
# Program: clustering_preprocess.R
# Purpose: Preprocess and filter spatial transcriptomics (ST) expression data to prepare for clustering.
#-------------------------------------------------------------------------------
# Data Inputs:
# - `raw_data_dir`: Path to the raw ST data directory (where the count matrices locate).
# - `processed_data_dir`: Path to the processed data directory (where the Seurat objects will be saved).
# - `figure_dir`: Path to the directory where figures will be saved.
# - metadata file for the sample (includes information for each sample)
# - metadata file for spatial transcriptomics data (includes the file names of count matrices).
#
# Data Outputs:
# - Processed and filtered Seurat objects saved in `processed_data_dir`.
# - Quality control (QC) plots saved in `figure_dir`.
#-------------------------------------------------------------------------------
# Notes:
# 1. This script includes spot filtering to remove spots with abnormal gene detection. (Section `Filter spots`)
#   - Filter out spots that have very few or a lot features detected, < 500 or > 4000.
#   - Only filter out genes that are not expressed in all the spots. 
#   - Output Seurat object: `merged_st_obj_filt_spots.rds`
# 2. After spot filtering, three ways of gene filtering is performed to remove genes that are rarely expressed:
#   - Section `Filter genes 1`: 
#     - For each sample, after filtering out poor spots, filter out genes that are not expressed in more than 90% of spots. 
#     - Output Seurat object: `merged_st_obj_filt_zeros_per_sample.rds`
#   - Section `Filter genes 2`: 
#     - Filter out any genes that are not expressed in more than 90% of spots in at least one sample.
#     - Identify common highly expressed genes that are expressed in more than 10% of spots in all samples. Filter out all other genes. 
#     - Output Seurat object: `merged_st_obj_filt_common_high_expr_genes.rds`
#   - Section `Filter genes 3`: 
#     - First merge all the samples to get a pooled expression matrix, and then filter out genes that are not expressed in more than 90% of spots. 
#     - Output Seurat object: `merged_st_obj_filt_zeros_all.rds`
#-------------------------------------------------------------------------------



filtered_zeros_all_90: f

library(Seurat) 
library(tidyverse)
library(Matrix) 
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot()) 
options(future.globals.maxSize = 1e9) 

# Define directory paths
raw_data_dir = "../../../st2image_data/He_et_al_2020/data/"
processed_data_dir = "../../data/He_2020/"
figure_dir = "../../figure/He_2020/"

# Load processed metadata for spots with information, proportions and annotations.
load(file.path(processed_data_dir, "prop_ann_meta.RData"))
# Read ST data metadata
st_meta <- read.csv(file.path(raw_data_dir, "metadata.csv"))

# Get unique IDs for spots in ST data
spots <- unique(prop_ann_meta$X) 

# Function to process the count matrix of each sample and create Seurat object
process_sample <- function(sample_meta) {
  # Get sample ID
  s1 <- paste0(sample_meta$patient, "_", sample_meta$replicate) 
  print(s1)
  
  # Get the cancer subtype
  type1 <- sample_meta$type
  # Simplify subtype
  type1 <- ifelse(type1 %in% c('Luminal_A', 'Luminal_B'), "ER+", 
                  ifelse(type1 %in% c('HER2_luminal', 'HER2_non_luminal'), "HER2+", type1)) 
  
  # Read and preprocess count matrix
  st_counts_s <- read.csv(file.path(raw_data_dir, sample_meta$count_matrix), sep = "\t") %>%
    remove_rownames() %>%
    column_to_rownames(var="X") %>%
    t() %>%
    as("CsparseMatrix") %>%
    .[startsWith(rownames(.), "ENSG"),] # Keep only rows starting with "ENSG"
  
  # Select spots that are in ST data and belongs to this sample and subset the corresponding count matrix 
  spots_s <- spots[startsWith(spots, s1)]
  spots_s <- sapply(strsplit(spots_s, "_"), function(x) x[3])
  st_counts_s <- st_counts_s[, spots_s]
  
  # Create Seurat object with spots in ST data
  CreateSeuratObject(counts = st_counts_s, min.features = 100, project = s1)
}

# Function to create QC plots
create_qc_plots <- function(obj, prefix) {
  pdf(file.path(figure_dir, paste0(prefix, "_QC.pdf")), width = 24, height = 8)
  
  features <- c("nFeature_RNA", "nCount_RNA") # Features for QC plots
  Idents(obj) <- "sid"
  for (feature in features) {
    p1 <- VlnPlot(obj, features = feature, ncol = 1, split.by = "subtype", cols = scales::hue_pal()(5)) + # Violin plots split by subtype
      VlnPlot(obj, features = feature, ncol = 1, split.by = "lab") # Violin plots split by lab
    
    print(p1)
  }
  
  # Scatter plots for QC
  p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "subtype", split.by = "sid", ncol = 9, cols = scales::hue_pal()(5)) +
    FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "lab", split.by = "sid", ncol = 9)
  
  print(p2)
  dev.off()
}

# Create a merged Seurat object for all samples ---------------------------

# Process all samples and create Seurat object for each sample
all_st_obj <- lapply(1:nrow(st_meta), function(i) process_sample(st_meta[i,]))
# Merge Seurat objects into one over
merged_st_obj <- merge(x=all_st_obj[[1]],y = c(all_st_obj)[-1],
                       add.cell.ids = paste0(st_meta$patient, "_", st_meta$replicate))
# Add metadata to merged Seurat object
metadata <- merged_st_obj@meta.data %>%
  dplyr::rename(sid = orig.ident) %>%
  mutate(X = rownames(.),
         pid = sapply(str_split(sid, "_"), `[[`, 1),
         rid = sapply(str_split(sid, "_"), `[[`, 2),
         subtype = plyr::mapvalues(pid, from = unique(st_meta[,c("type", "patient")])$patient, to = unique(st_meta[,c("type", "patient")])$type),
         label = plyr::mapvalues(X, from = unique(prop_ann_meta[,c("X", "label")])$X, to = unique(prop_ann_meta[,c("X", "label")])$label),
         lab = ifelse(rid %in% c("C1", "C2", "D1"), "lab_C1_C2_D1", "lab_D2_E1_E2"))
# Assign metadata to Seurat object
merged_st_obj@meta.data <- metadata 
# Save the Seurat object
saveRDS(merged_st_obj, file = file.path(processed_data_dir, "merged_st_obj_unfilt.rds")) # Save the Seurat object
# Generate QC plots for unfiltered data
create_qc_plots(merged_st_obj, "unfilt") 


# Filter spots ------------------------------------------------------------

# Gene detection filtering (filter out spots with abnormal gene detection)
ggplot() +
  geom_density(aes(x=log10(merged_st_obj$nFeature_RNA)))

tmp <- merged_st_obj$nFeature_RNA
mean(tmp>4000)
mean(tmp<500)
# Filter spots that have unique feature counts over 4000 or less than 500
high_det <- WhichCells(merged_st_obj, expression = nFeature_RNA > 4000)
low_det <- WhichCells(merged_st_obj, expression = nFeature_RNA < 500)
filtered_cells <- setdiff(WhichCells(merged_st_obj), c(low_det, high_det))
merged_st_obj_tmp <- subset(merged_st_obj, cells = filtered_cells)

# Number of cells after filtering
cbind(table(Idents(merged_st_obj)), table(Idents(merged_st_obj_tmp)))

# Remove genes that are no longer expressed in any cell
all_st_obj_filt <- list()
for(l1 in Layers(merged_st_obj_tmp)){
  s1 = gsub("counts.", "", l1)
  print(s1)
  filt_matrix <- GetAssayData(merged_st_obj_tmp, "RNA", layer = l1)
  filtered_obj <- CreateSeuratObject(counts = filt_matrix, project = s1, min.cells = 3)
  filtered_obj$orig.ident <- s1
  all_st_obj_filt[[s1]] <- filtered_obj
}
merged_st_obj_filt <- merge(x=all_st_obj_filt[[1]], y = c(all_st_obj_filt)[-1])
# Update metadata for filtered object
merged_st_obj_filt@meta.data <- cbind(select(merged_st_obj_filt@meta.data, -c("orig.ident")), 
                                      select(merged_st_obj@meta.data, -c("nCount_RNA", "nFeature_RNA"))[rownames(merged_st_obj_filt@meta.data), ])
# Save the Seurat object 
saveRDS(merged_st_obj_filt, file=file.path(processed_data_dir, "merged_st_obj_filt_spots.rds"))
# Generate QC plots for filtered data
create_qc_plots(merged_st_obj_filt, "filt_spots")


# Get the list of genes that are expressed in more than 10% of spots in each sample. 
pdf(file.path(figure_dir, "zero_proportions.pdf"), width = 8, height = 5)
high_expr_genes <- list()
low_expr_genes <- list()
for(l1 in Layers(merged_st_obj_tmp)){
  s1 = gsub("counts.", "", l1)
  print(s1)
  filt_matrix <- GetAssayData(merged_st_obj_tmp, "RNA", layer = l1)
  zero_prop <- rowSums(filt_matrix == 0)/ncol(filt_matrix)
  p <- ggplot() +
    geom_histogram(aes(x=zero_prop), binwidth = 0.01) +
    geom_vline(xintercept = 0.9, col = "red") +
    geom_vline(xintercept = 0.8, col = "blue") +
    labs(title = s1) +
    annotate("text", x= 0.3, y=2500, col = "red",
             label= paste0(sum(zero_prop < 0.9), " genes are expressed in more than 10% spots.")) + 
    annotate("text", x = 0.3, y=2300, col = "blue",
             label = paste0(sum(zero_prop < 0.8), " genes are expressed in more than 20% spots."))
  print(p)
  high_expr_genes[[l1]] <- rownames(filt_matrix)[zero_prop < 0.9]
  low_expr_genes[[l1]] <- rownames(filt_matrix)[zero_prop >= 0.9]
}
dev.off()

# Get intersection of genes that are expressed in more than 10% of spots in each sample.
common_high_expr_genes <- Reduce(intersect, high_expr_genes)
length(common_high_expr_genes)
# Get intersection of genes that are not expressed in more than 90% of spots in each sample.
common_low_expr_genes <- Reduce(intersect, low_expr_genes)
length(common_low_expr_genes)


# Filter genes 1 ----------------------------------------------------------

# 1. Remove genes that are not expressed in more than 90% of spots for each sample.
all_st_obj_filt <- list()
for(l1 in Layers(merged_st_obj_tmp)){
  s1 = gsub("counts.", "", l1)
  print(s1)
  filt_matrix <- GetAssayData(merged_st_obj_tmp, "RNA", layer = l1)
  filtered_obj <- CreateSeuratObject(counts = filt_matrix, project = s1, min.cells = 0.1*ncol(filt_matrix))
  filtered_obj$orig.ident <- s1
  all_st_obj_filt[[s1]] <- filtered_obj
}
merged_st_obj_filt <- merge(x=all_st_obj_filt[[1]], y = c(all_st_obj_filt)[-1])
# Update metadata for filtered object
merged_st_obj_filt@meta.data <- cbind(select(merged_st_obj_filt@meta.data, -c("orig.ident")), 
                                      select(merged_st_obj@meta.data, -c("nCount_RNA", "nFeature_RNA"))[rownames(merged_st_obj_filt@meta.data), ])
# Save the Seurat object 
saveRDS(merged_st_obj_filt, file=file.path(processed_data_dir, "merged_st_obj_filt_zeros_per_sample.rds"))
# Generate QC plots for filtered data
create_qc_plots(merged_st_obj_filt, "filt_zeros_per_sample")

# Filter genes 2 ----------------------------------------------------------

# 2. Remove genes that are not expressed in more than 90% of spots in at least one sample.
all_st_obj_filt <- list()
for(l1 in Layers(merged_st_obj_tmp)){
  s1 = gsub("counts.", "", l1)
  print(s1)
  filt_matrix <- GetAssayData(merged_st_obj_tmp, "RNA", layer = l1)
  filt_matrix <- filt_matrix[common_high_expr_genes,]
  filtered_obj <- CreateSeuratObject(counts = filt_matrix, project = s1, min.cells = 0)
  filtered_obj$orig.ident <- s1
  all_st_obj_filt[[s1]] <- filtered_obj
}
merged_st_obj_filt <- merge(x=all_st_obj_filt[[1]], y = c(all_st_obj_filt)[-1])
# Update metadata for filtered object
merged_st_obj_filt@meta.data <- cbind(select(merged_st_obj_filt@meta.data, -c("orig.ident")), 
                                      select(merged_st_obj@meta.data, -c("nCount_RNA", "nFeature_RNA"))[rownames(merged_st_obj_filt@meta.data), ])
# Save the Seurat object 
saveRDS(merged_st_obj_filt, file=file.path(processed_data_dir, "merged_st_obj_filt_common_high_expr_genes.rds"))
# Generate QC plots for filtered data
create_qc_plots(merged_st_obj_filt, "filt_common_high_expr_genes")

# Filter genes 3 ----------------------------------------------------------

# 3. Remove genes that are not expressed in more than 90% of all pooled spots
merged_st_obj_tmp <- JoinLayers(merged_st_obj_tmp)
filt_matrix <- GetAssayData(merged_st_obj_tmp, "RNA", layer = "counts")
zero_prop <- rowSums(filt_matrix == 0)/ncol(filt_matrix)

pdf(file.path(figure_dir, "zero_proportions_all.pdf"), width = 8, height = 5)
ggplot() +
  geom_histogram(aes(x=zero_prop), binwidth = 0.01) +
  geom_vline(xintercept = 0.9, col = "red") +
  geom_vline(xintercept = 0.8, col = "blue") +
  labs(title = "All Samples") +
  annotate("text", x= 0.3, y=9000, col = "red",
           label= paste0(sum(zero_prop < 0.9), " genes are expressed in more than 10% spots.")) + 
  annotate("text", x = 0.3, y=8500, col = "blue",
           label = paste0(sum(zero_prop < 0.8), " genes are expressed in more than 20% spots."))

dev.off()

filtered_obj <- CreateSeuratObject(counts = filt_matrix, project = "all", min.cells = 0.1*ncol(filt_matrix))
# Update metadata for filtered object
filtered_obj@meta.data <- cbind(select(filtered_obj@meta.data, -c("orig.ident")), 
                                select(merged_st_obj@meta.data, -c("nCount_RNA", "nFeature_RNA"))[rownames(merged_st_obj_filt@meta.data), ])
filtered_obj[["RNA"]] <- split(filtered_obj[["RNA"]], f = filtered_obj$sid)
merged_st_obj_filt <- filtered_obj
# Save the Seurat object 
saveRDS(merged_st_obj_filt, file=file.path(processed_data_dir, "merged_st_obj_filt_zeros_all.rds"))
# Generate QC plots for filtered data
create_qc_plots(merged_st_obj_filt, "filt_zeros_all")

