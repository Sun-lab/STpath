# Date last modified: May 31, 2024 by Zhining Sui
# Program: clustering.R
# Purpose: Normalize, integrate, and cluster spatial transcriptomics (ST) data and identify differentially expressed (DE) markers.
#-------------------------------------------------------------------------------
# Data Inputs:
# - `processed_data_dir`: Path to the processed data directory (where the Seurat objects to be clustered are saved).
# - `cluster_rslt_dir`: Path to the directory where clustering results will be saved.
# - `merged_st_obj`: Seurat object containing filtered ST data (should be in `processed_data_dir`).
# - `id_symbol_mapping`: Data frame for mapping gene IDs to gene symbols (optional).
#
# Data Outputs:
# - Normalized, integrated, and clustered Seurat objects saved in `processed_data_dir` and `cluster_rslt_dir`.
# - Identified DE markers and top 10 DE markers for each cluster saved in `cluster_rslt_dir`.
#-------------------------------------------------------------------------------
# Notes:
# 1. **Normalization** (function `normalize_data`):
#    - SCTransform normalization followed by PCA, neighbor finding, clustering, and UMAP embedding.
#    - Log-normalization followed by PCA, neighbor finding, clustering, and UMAP embedding.
#    - Output Seurat object: `merged_st_obj_normalized_<suffix>.rds`
# 2. **Integration** (function `integrate_data`):
#    - Integrate data using CCA, RPCA, and Harmony methods with both log-normalized and SCTransform data.
#    - Output Seurat object: `integrated_st_obj_<suffix>.rds`
# 3. **Clustering** (function `perform_clustering`):
#    - Perform clustering using integrated reductions at multiple resolutions.
#    - Output Seurat object: `clustered_st_obj_<suffix>.rds`
# 4. **Identification of DE Markers** (function `identify_DE_markers`):
#    - Identify DE markers for each cluster.
#    - Output list of DE markers: `DE_markers_identified_<suffix>.rds`
# 5. **Identification of Top 10 DE Markers** (function `identify_top10_markers`):
#    - Extract top 10 DE markers for each cluster based on adjusted p-value and log fold change.
#    - Output list of top 10 DE markers: `top10_DE_markers_identified_<suffix>.rds`
#-------------------------------------------------------------------------------


library(Seurat) 
library(ggplot2) 
library(cowplot) 
library(tidyverse)
theme_set(theme_cowplot()) 
options(future.globals.maxSize = 1e9) 

# Define functions --------------------------------------------------------

# Normalization function
normalize_data <- function(obj) {
  obj <- SCTransform(obj, return.only.var.genes = F) # SCTransform normalization
  obj <- RunPCA(obj, reduction.name = "pca.unintegrated.sct") # PCA
  obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca.unintegrated.sct") # Find neighbors
  obj <- FindClusters(obj, resolution = 0.1, cluster.name = "unintegrated.sct_clusters_0.1") # Find clusters
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca.unintegrated.sct", reduction.name = "umap.unintegrated.sct") # Run UMAP
  
  DefaultAssay(obj) <- "RNA" # Switch to RNA assay
  obj <- NormalizeData(obj) # Log-normalization
  obj <- FindVariableFeatures(obj, nfeatures = 2000) # Identify variable features
  obj <- ScaleData(obj, assay="RNA", features = rownames(obj)) # Scale data
  obj <- RunPCA(obj, reduction.name = "pca.unintegrated.lognorm") # PCA
  obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca.unintegrated.lognorm") # Find neighbors
  obj <- FindClusters(obj, resolution = 0.1, cluster.name = "unintegrated.lognorm_clusters_0.1") # Find clusters
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca.unintegrated.lognorm", reduction.name = "umap.unintegrated.lognorm") # Run UMAP
  obj
}

# Integration function
integrate_data <- function(obj) {
  obj <- IntegrateLayers(object = obj, method = CCAIntegration, normalization.method = "LogNormalize", orig.reduction = "pca.unintegrated.lognorm", new.reduction = "integrated.cca.lognorm")
  obj <- IntegrateLayers(object = obj, method = RPCAIntegration, normalization.method = "LogNormalize", orig.reduction = "pca.unintegrated.lognorm", new.reduction = "integrated.rpca.lognorm", k.weight = 80)
  obj <- IntegrateLayers(object = obj, method = HarmonyIntegration, normalization.method = "LogNormalize", orig.reduction = "pca.unintegrated.lognorm", new.reduction = "integrated.harmony.lognorm", verbose = TRUE)
  DefaultAssay(obj) <- "SCT"
  obj <- IntegrateLayers(object = obj, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "pca.unintegrated.sct", new.reduction = "integrated.cca.sct")
  obj <- IntegrateLayers(object = obj, method = RPCAIntegration, normalization.method = "SCT", orig.reduction = "pca.unintegrated.sct", new.reduction = "integrated.rpca.sct", k.weight = 80)
  obj <- IntegrateLayers(object = obj, method = HarmonyIntegration, normalization.method = "SCT", orig.reduction = "pca.unintegrated.sct", new.reduction = "integrated.harmony.sct", verbose = TRUE)
  obj
}

# Clustering function
perform_clustering <- function(obj) {
  reductions <- names(obj@reductions)[startsWith(names(obj@reductions), "integrated")]
  DefaultAssay(obj) <- "RNA"
  for (red in reductions[grepl("lognorm", reductions)]) {
    obj <- FindNeighbors(obj, reduction = red, dims = 1:30)
    obj <- FindClusters(obj, resolution = c(0.3, 0.2, 0.1, 0.05), cluster.name = paste0(red, "_clusters_", c(0.3, 0.2, 0.1, 0.05)))
    obj <- RunUMAP(obj, reduction = red, dims = 1:30, reduction.name = paste0("umap.", red))
  }
  DefaultAssay(obj) <- "SCT"
  for (red in reductions[grepl("sct", reductions)]) {
    obj <- FindNeighbors(obj, reduction = red, dims = 1:30)
    obj <- FindClusters(obj, resolution = c(0.3, 0.2, 0.1, 0.05), cluster.name = paste0(red, "_clusters_", c(0.3, 0.2, 0.1, 0.05)))
    obj <- RunUMAP(obj, reduction = red, dims = 1:30, reduction.name = paste0("umap.", red))
  }
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)
  DefaultAssay(obj) <- "SCT"
  obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)
  obj
}

# Identify differentially expressed genes
identify_DE_markers <- function(obj) {
  cluster_methods <- sort(colnames(obj@meta.data)[startsWith(colnames(obj@meta.data), "integrated")])
  markers_identified <- list()
  for (cm in cluster_methods) {
    m <- strsplit(cm, "_")[[1]][1]
    norm <- strsplit(m, ".", fixed = TRUE)[[1]][3]
    DefaultAssay(obj) <- ifelse(norm == "lognorm", "RNA", "SCT")
    Idents(obj) <- cm
    markers <- FindAllMarkers(obj, verbose = TRUE)
    markers_identified[[cm]] <- markers
  }
  markers_identified
}

# Identify top 10 differentially expressed markers
identify_top10_markers <- function(markers_identified) {
  top10_DE_markers <- list()
  for (cm in names(markers_identified)) {
    markers_cm <- markers_identified[[cm]]
    if (nrow(markers_cm) == 0) {
      top10_DE_markers[[cm]] <- data.frame()
    } else {
      top10 <- markers_cm %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.001 & avg_log2FC > 0) %>%
        slice_max(avg_log2FC, n = 10) %>%
        ungroup()
      top10$gene_symbol <- plyr::mapvalues(top10$gene, from = id_symbol_mapping$id, to = id_symbol_mapping$symbol)
      top10_DE_markers[[cm]] <- top10
    }
  }
  top10_DE_markers
}


# Start run  --------------------------------------------------------------

# Define directory paths
processed_data_dir = "../../data/He_2020/"
cluster_rslt_dir = "../../output/He_2020/clustering/"

# Load Seurat object to be clustered
suffix <- "filt_zeros_all"
merged_st_obj <- readRDS(file.path(processed_data_dir, paste0("merged_st_obj_", suffix, ".rds")))
# Load symbol mapping data
load(file.path(processed_data_dir, "id_symbol_mapping.RData"))

# Normalize data
merged_st_obj <- normalize_data(merged_st_obj)
saveRDS(merged_st_obj, file = file.path(processed_data_dir, paste0("merged_st_obj_normalized_", suffix, ".rds")))

# Integrate data
integrated_st_obj <- integrate_data(merged_st_obj) 
saveRDS(integrated_st_obj, file = file.path(processed_data_dir, paste0("integrated_st_obj_", suffix, ".rds")))

# Perform clustering 
clustered_st_obj <- perform_clustering(integrated_st_obj) 
saveRDS(clustered_st_obj, file = file.path(cluster_rslt_dir, paste0("clustered_st_obj_", suffix, ".rds"))) 

# Identify all markers
markers_identified <- identify_DE_markers(clustered_st_obj)
saveRDS(markers_identified, file = file.path(cluster_rslt_dir, paste0("DE_markers_identified_", suffix, ".rds")))

# Extract the top10 markers
top10_DE_markers <- identify_top10_markers(markers_identified)
saveRDS(top10_DE_markers, file = file.path(cluster_rslt_dir, paste0("top10_DE_markers_identified_", suffix, ".rds")))
