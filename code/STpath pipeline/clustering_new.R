library(dplyr) 
library(tidyr) 
library(stringr) 
library(ggplot2) 
library(Seurat) 
library(Matrix) 
library(gridExtra) 
library(cowplot) 
theme_set(theme_cowplot()) 
options(future.globals.maxSize = 1e9) 

# Define directory paths
raw_data_dir = "../../../st2image_data/He_et_al_2020/data/"
processed_data_dir = "../../data/He_2020/"
figure_dir = "../../figure/He_2020/"
# sc_dir <- "../../st2image_data/Wu_2021/data/scRNASeq/"
# st_dir <- "../../st2image_data/He_et_al_2020/data/"
# st_output_dir <- "../../st2image_data/He_et_al_2020/output/"
# decon_dir <- "../output/He_et_al_2020/deconvolution/"
# cluster_dir <- "../output/He_et_al_2020/clustering/"

# Load symbol mapping data
load(file.path(processed_data_dir, "id_symbol_mapping.RData"))
# Load Seurat object for all samples 
load(file.path(processed_data_dir, "merged_st_obj_unfilt.RData"))

# spots <- unique(prop_ann_meta$X) # Get unique spot identifiers
# annotations <- unique(prop_ann_meta[, c("X", "label")]) # Get unique annotations
# rownames(annotations) <- annotations$X # Set row names for annotations

# Normalization function
normalize_data <- function(obj) {
  obj <- SCTransform(obj, return.only.var.genes = F) # SCTransform normalization
  obj <- RunPCA(obj, reduction.name = "pca.unintegrated.sct") # PCA
  obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca.unintegrated.sct") # Find neighbors
  obj <- FindClusters(obj, resolution = 0.1, cluster.name = "unintegrated.sct_clusters_0.1") # Find clusters
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca.unintegrated.sct", reduction.name = "umap.unintegrated.sct") # Run UMAP
  DefaultAssay(obj) <- "RNA" # Switch to RNA assay
  obj <- NormalizeData(obj) # Log-normalization
  obj <- FindVariableFeatures(obj) # Identify variable features
  obj <- ScaleData(obj, assay="RNA", features = rownames(obj)) # Scale data
  obj <- RunPCA(obj, reduction.name = "pca.unintegrated.lognorm") # PCA
  obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca.unintegrated.lognorm") # Find neighbors
  obj <- FindClusters(obj, resolution = 0.1, cluster.name = "unintegrated.lognorm_clusters_0.1") # Find clusters
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca.unintegrated.lognorm", reduction.name = "umap.unintegrated.lognorm") # Run UMAP
  obj
}

merged_st_obj <- normalize_data(merged_st_obj) # Normalize data
save(merged_st_obj, file = file.path(data_dir, "merged_st_obj_normalized_unfilt.RData")) # Save normalized object

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

integrated_st_obj <- integrate_data(merged_st_obj) # Integrate data
save(integrated_st_obj, file = file.path(st_dir, "integrated_st_obj_unfilt.RData")) # Save integrated object

# Clustering function
cluster_data <- function(obj, prefix) {
  reductions <- names(obj@reductions)[startsWith(names(obj@reductions), "integrated")] # Get integrated reductions
  for (red in reductions) {
    obj <- FindNeighbors(obj, reduction = red, dims = 1:30) # Find neighbors
    obj <- FindClusters(obj, resolution = c(0.3, 0.2, 0.1, 0.05), cluster.name = paste0(red, "_clusters")) # Find clusters at multiple resolutions
    obj <- RunUMAP(obj, reduction = red, dims = 1:30, reduction.name = paste0("umap.", red)) # Run UMAP
  }
  obj <- JoinLayers(obj) # Join layers
  DefaultAssay(obj) <- "SCT"
  obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE) # Prepare for finding markers
  obj
}

integrated_st_obj_clusters <- cluster_data(integrated_st_obj, "unfilt") # Cluster data
save(integrated_st_obj_clusters, file = file.path(st_dir, "integrated_st_obj_clusters_unfilt.RData")) # Save clustered object




