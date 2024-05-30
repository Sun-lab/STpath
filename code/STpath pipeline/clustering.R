library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
library(ggplot2)
library(Seurat)
library(Matrix)
library(tidyverse)
library(gridExtra)
library(cowplot)
theme_set(theme_cowplot())
options(future.globals.maxSize = 1e9)

### directories ----
sc_dir = "../../st2image_data/Wu_2021/data/scRNASeq/"
st_dir = "../../st2image_data/BRCA/data/"
st_output_dir = "../../st2image_data/BRCA/output/"
decon_dir = "../output/BRCA/deconvolution/"
cluster_dir = "../output/BRCA/clustering/"

load(file.path(st_dir, "prop_ann_meta.RData"))
spots <- unique(prop_ann_meta$X)
annotations <- unique(prop_ann_meta[, c("X", "label")])
rownames(annotations) <- annotations$X

matched_symbols <- read.csv(file.path(st_output_dir, "id_symbol_conversion/matched_symbols_all.csv"))[,-1]
id_symbol_mapping <- unique(matched_symbols[,1:5])
for (i in 1:nrow(id_symbol_mapping)) {
  symbol <- unique(unlist(id_symbol_mapping[i,2:5]))
  if(sum(!is.na(symbol))==0) next
  if(length(symbol)>1) symbol = paste(na.omit(symbol), collapse = "/")
  id_symbol_mapping[i,"symbol"] <- symbol
}
id_symbol_mapping <- na.omit(id_symbol_mapping[,c("id", "symbol")])

save(id_symbol_mapping, file = "id_symbol_mapping.RData")

st_meta = read.csv(file.path(st_dir, "metadata.csv"))
# 1. Read in samples ------------------------------------------------------

all_st_obj <- list()
for(i in 1:nrow(st_meta)){
  s1 = paste0(st_meta[i,"patient"], "_", st_meta[i,"replicate"])
  print(s1)
  type1 = st_meta[i,"type"]
  type1 = ifelse(type1 %in% c('Luminal_A', 'Luminal_B'), "ER+", 
                 ifelse(type1 %in% c('HER2_luminal', 'HER2_non_luminal'), "HER2+", type1))
  
  ## Read in 10X ST data for a single sample
  st_counts_s <- read.csv(file = file.path(st_dir, st_meta[i, "count_matrix"]), sep = "\t")
  st_counts_s <- st_counts_s %>% 
    remove_rownames %>% 
    column_to_rownames(var="X") %>%
    t()
  st_counts_s <- as(st_counts_s, "sparseMatrix")
  st_counts_s <- as(st_counts_s, "CsparseMatrix")
  
  st_counts_s <- st_counts_s[startsWith(rownames(st_counts_s), "ENSG"),]
  
  # rownames(st_counts_s) <- plyr::mapvalues(rownames(st_counts_s), 
  #                                          from = id_symbol_mapping$id, 
  #                                          to = id_symbol_mapping$symbol)
  spots_s <- spots[startsWith(spots, s1)]
  spots_s <- sapply(strsplit(spots_s, "_"), "[[", 3)
  st_counts_s <- st_counts_s[,spots_s]
  
  gene_ids <- rownames(st_counts_s)
  cell_ids <- colnames(st_counts_s)
  
  # gene_ids[duplicated(gene_ids)]
  ## Turn count matrix into a Seurat object (output is a Seurat object)
  st_s <- CreateSeuratObject(counts = st_counts_s,
                             min.features = 100,
                             project = s1)
  all_st_obj[[i]] <- st_s
}

## Create a merged Seurat object
merged_st_obj <- merge(x=all_st_obj[[1]],y = c(all_st_obj)[-1],
                       add.cell.ids = paste0(st_meta[1:length(all_st_obj),"patient"], 
                                             "_", 
                                             st_meta[1:length(all_st_obj),"replicate"])
                       )

# Add number of genes per UMI for each cell to metadata
merged_st_obj$log10GenesPerUMI <- log10(merged_st_obj$nFeature_RNA) / log10(merged_st_obj$nCount_RNA)

# Create metadata dataframe
metadata <- merged_st_obj@meta.data

# Rename columns
metadata <- metadata %>%
  dplyr::rename(sid = orig.ident)

# Add additional information to metadata
metadata$X <- rownames(metadata)
metadata$pid <- sapply(str_split(metadata$sid, "_"), "[[", 1)
metadata$rid <- sapply(str_split(metadata$sid, "_"), "[[", 2)
metadata$subtype <- plyr::mapvalues(metadata$pid, 
                                 from = unique(st_meta[,c("type", "patient")])$patient, 
                                 to = unique(st_meta[,c("type", "patient")])$type)

metadata$label <- plyr::mapvalues(metadata$X, 
                                  from = unique(prop_ann_meta[,c("X", "label")])$X, 
                                  to = unique(prop_ann_meta[,c("X", "label")])$label)

metadata$lab <- ifelse(metadata$rid %in% c("C1", "C2", "D1"), "lab_C1_C2_D1", "lab_D2_E1_E2")

# Add metadata back to Seurat object
merged_st_obj@meta.data <- metadata
save(merged_st_obj, file=file.path(st_dir, "merged_st_obj_unfilt.RData"))

# ========== UNFILTERED ========== ----------------------------------------

# load(file.path(st_dir, "merged_st_obj_unfilt.RData"))
merged_st_obj_unfilt <- merged_st_obj
rm(merged_st_obj)

# 2. Quality Control ---------------------------------------------------------

table(Idents(merged_st_obj_unfilt))
colors_subtype <- scales::hue_pal()(5)
colors_subtype <- setNames(colors_subtype, unique(metadata$subtype))

pdf("QC.pdf", width = 16, height = 8)

# Visualize QC metrics as a violin plot
VlnPlot(merged_st_obj_unfilt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_unfilt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_unfilt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_unfilt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_unfilt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_unfilt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "lab")

FeatureScatter(merged_st_obj_unfilt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "subtype",
               split.by = "sid", ncol = 9,
               cols = colors_subtype)

FeatureScatter(merged_st_obj_unfilt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "lab",
               split.by = "sid", ncol = 9)


dev.off()


# 3. Run normalization ----------------------------------------------------

# SCTransform-normalization
merged_st_obj_unfilt <- SCTransform(merged_st_obj_unfilt, return.only.var.genes = F)
merged_st_obj_unfilt <- RunPCA(merged_st_obj_unfilt, reduction.name = "pca.unintegrated.sct")
merged_st_obj_unfilt <- FindNeighbors(merged_st_obj_unfilt, dims = 1:30, reduction = "pca.unintegrated.sct")
merged_st_obj_unfilt <- FindClusters(merged_st_obj_unfilt, resolution = 0.1, cluster.name = "unintegrated.sct_clusters_0.1")
merged_st_obj_unfilt <- RunUMAP(merged_st_obj_unfilt, dims = 1:30, reduction = "pca.unintegrated.sct", reduction.name = "umap.unintegrated.sct")

# Log-normalization
DefaultAssay(merged_st_obj_unfilt) <- "RNA"
merged_st_obj_unfilt <- NormalizeData(merged_st_obj_unfilt)
merged_st_obj_unfilt <- FindVariableFeatures(merged_st_obj_unfilt)

# top10 <- head(VariableFeatures(merged_st_obj_unfilt), 10) # the 10 most highly variable genes
# plot variable features without labels
# VariableFeaturePlot(merged_st_obj_unfilt)

merged_st_obj_unfilt <- ScaleData(merged_st_obj_unfilt,  
                                  assay="RNA", 
                                  features = rownames(merged_st_obj_unfilt))

merged_st_obj_unfilt <- RunPCA(merged_st_obj_unfilt, reduction.name = "pca.unintegrated.lognorm")
# print(merged_st_obj_unfilt[["pca.unintegrated.lognorm"]], dims = 1:5, nfeatures = 10)

# VizDimLoadings(merged_st_obj_unfilt, dims = 1:2, reduction = "pca.unintegrated.lognorm")
# DimPlot(merged_st_obj_unfilt, reduction = "pca.unintegrated.lognorm")
# DimHeatmap(merged_st_obj_unfilt, reduction = "pca.unintegrated.lognorm", dims = 1:15, cells = 500, balanced = TRUE)

merged_st_obj_unfilt <- FindNeighbors(merged_st_obj_unfilt, dims = 1:30, reduction = "pca.unintegrated.lognorm")
merged_st_obj_unfilt <- FindClusters(merged_st_obj_unfilt, resolution = 0.1, cluster.name = "unintegrated.lognorm_clusters_0.1")

merged_st_obj_unfilt <- RunUMAP(merged_st_obj_unfilt, dims = 1:30, reduction = "pca.unintegrated.lognorm", 
                         reduction.name = "umap.unintegrated.lognorm")
# DimPlot(merged_st_obj_unfilt, reduction = "umap.unintegrated.lognorm", group.by = c("sid"))

save(merged_st_obj_unfilt, file = file.path(st_dir, "merged_st_obj_normalized_unfilt.RData"))

# 4. Run integration ------------------------------------------------------
integrated_st_obj_unfilt <- merged_st_obj_unfilt
rm(merged_st_obj_unfilt)

DefaultAssay(integrated_st_obj_unfilt) <- "RNA"
integrated_st_obj_unfilt <- IntegrateLayers(
  object = integrated_st_obj_unfilt,
  method = CCAIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm",
  new.reduction = "integrated.cca.lognorm"
)

integrated_st_obj_unfilt <- IntegrateLayers(
  object = integrated_st_obj_unfilt, 
  method = RPCAIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.rpca.lognorm",
  k.weight = 80
)

integrated_st_obj_unfilt <- IntegrateLayers(
  object = integrated_st_obj_unfilt, 
  method = HarmonyIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.harmony.lognorm",
  verbose = T
)

?IntegrateLayers
DefaultAssay(integrated_st_obj_unfilt) <- "SCT"
integrated_st_obj_unfilt <- IntegrateLayers(
  object = integrated_st_obj_unfilt,
  method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct",
  new.reduction = "integrated.cca.sct"
)

integrated_st_obj_unfilt <- IntegrateLayers(
  object = integrated_st_obj_unfilt, 
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.rpca.sct",
  k.weight = 80
)

integrated_st_obj_unfilt <- IntegrateLayers(
  object = integrated_st_obj_unfilt, 
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.harmony.sct",
  verbose = T
)

save(integrated_st_obj_unfilt, file = file.path(st_dir, "integrated_st_obj_unfilt.RData"))


# 5. Clustering -----------------------------------------------------------

integrated_st_obj_clusters_unfilt <- integrated_st_obj_unfilt
rm(integrated_st_obj_unfilt)

reductions <- names(integrated_st_obj_clusters_unfilt@reductions)[startsWith(names(integrated_st_obj_clusters_unfilt@reductions), "integrated")]

DefaultAssay(integrated_st_obj_clusters_unfilt) <- "RNA"
for(red in reductions[1:3]){
  integrated_st_obj_clusters_unfilt <- FindNeighbors(integrated_st_obj_clusters_unfilt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_unfilt <- FindClusters(integrated_st_obj_clusters_unfilt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_unfilt <- FindClusters(integrated_st_obj_clusters_unfilt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_unfilt <- FindClusters(integrated_st_obj_clusters_unfilt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_unfilt <- FindClusters(integrated_st_obj_clusters_unfilt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_unfilt <- RunUMAP(integrated_st_obj_clusters_unfilt, reduction = red, 
                               dims = 1:30, reduction.name = paste0("umap.", red))
}

DefaultAssay(integrated_st_obj_clusters_unfilt) <- "SCT"
for(red in reductions[4:6]){
  integrated_st_obj_clusters_unfilt <- FindNeighbors(integrated_st_obj_clusters_unfilt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_unfilt <- FindClusters(integrated_st_obj_clusters_unfilt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_unfilt <- FindClusters(integrated_st_obj_clusters_unfilt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_unfilt <- FindClusters(integrated_st_obj_clusters_unfilt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_unfilt <- FindClusters(integrated_st_obj_clusters_unfilt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_unfilt <- RunUMAP(integrated_st_obj_clusters_unfilt, reduction = red, 
                               dims = 1:30, reduction.name = paste0("umap.", red))
}

# cluster_methods <- colnames(integrated_st_obj_clusters_unfilt@meta.data)[startsWith(colnames(integrated_st_obj_clusters_unfilt@meta.data), "integrated")]
# for (m in cluster_methods) {
#   print(m)
#   print(table(integrated_st_obj_clusters_unfilt@meta.data[, m]))
# }

DefaultAssay(integrated_st_obj_clusters_unfilt) <- "RNA"
integrated_st_obj_clusters_unfilt <- JoinLayers(integrated_st_obj_clusters_unfilt)

DefaultAssay(integrated_st_obj_clusters_unfilt) <- "SCT"
integrated_st_obj_clusters_unfilt <- PrepSCTFindMarkers(integrated_st_obj_clusters_unfilt,
                                                        assay = "SCT", verbose = TRUE)

save(integrated_st_obj_clusters_unfilt, file = file.path(st_dir, "integrated_st_obj_clusters_unfilt.RData"))

# 6. Identify differential expressed genes across conditions --------------

# load("id_symbol_mapping.RData")
# cell_type_markers <- data.frame(celltype = c("Epithelial cells", 
#                                              "Cyclying cells",
#                                              "T cells",
#                                              "Myeloid cells",
#                                              "B cells",
#                                              "Plasmablasts",
#                                              "Endothelial cells",
#                                              "Mesenchymal cells"),
#                                 marker_symbol = c("EPCAM", 
#                                                   "MKI67", 
#                                                   "CD3D",
#                                                   "CD68",
#                                                   "MS4A1",
#                                                   "JCHAIN",
#                                                   "PECAM1",
#                                                   "PDGFRB"))
# 
# matched_symbols <- read.csv("/Users/zhiningsui/GitHub/st2image_data/BRCA/output/id_symbol_conversion/matched_symbols_all.csv")
# marker_ids <- matched_symbols[, 2:6] %>% 
#   filter_at(vars(-1), any_vars(. %in% cell_type_markers$marker_symbol)) %>%
#   unique() %>%
#   column_to_rownames("biomart")
# cell_type_markers <- cbind(cell_type_markers, marker_ids[cell_type_markers$marker_symbol,1])
# colnames(cell_type_markers)[3] <- "marker_ID"


cluster_methods <- sort(colnames(integrated_st_obj_clusters_unfilt@meta.data[startsWith(colnames(integrated_st_obj_clusters_unfilt@meta.data), "integrated")]))
seurat_integrated <- integrated_st_obj_clusters_unfilt
rm(integrated_st_obj_clusters_unfilt)

markers_identified_unfilt <- list()

for(cm in cluster_methods){
  m <- strsplit(cm, "_")[[1]][1]
  norm <- strsplit(m, ".", fixed = TRUE)[[1]][3]
  
  if(norm == "lognorm"){
    DefaultAssay(seurat_integrated) <- "RNA"
  }else if(norm == "sct"){
    DefaultAssay(seurat_integrated) <- "SCT"
  }else{
    break
  }
  
  Idents(seurat_integrated) <- cm
  markers <- FindAllMarkers(seurat_integrated, verbose = T)
  markers_identified_unfilt[[cm]] <- markers
}

save(markers_identified_unfilt, file = file.path(st_dir, "DE_markers_identified_unfilt.RData"))


top10_DE_markers_unfilt <- list()

for (cm in cluster_methods) {
  markers_cm <- markers_identified_unfilt[[cm]] 
  
  if(nrow(markers_cm)==0){
    top10_DE_markers_unfilt[[cm]] <- data.frame()
  }else{
    top10 <- markers_cm %>%
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < 0.001 & avg_log2FC > 0) %>%
      slice_max(avg_log2FC, n = 10) %>%
      ungroup()
    top10$gene_symbol <- plyr::mapvalues(top10$gene, 
                                         from = id_symbol_mapping$id,
                                         to = id_symbol_mapping$symbol)
    top10_DE_markers_unfilt[[cm]] <- top10
  }
}

save(top10_DE_markers_unfilt, file = file.path(st_dir, "top10_DE_markers_identified_unfilt.RData"))


rm(seurat_integrated)

# ========== FILTERED (spots) ========== ----------------------------------

load(file.path(st_dir, "merged_st_obj_unfilt.RData"))

# Gene detection filtering
ggplot() +
  geom_density(aes(x=log10(merged_st_obj$nFeature_RNA)))

tmp <- merged_st_obj$nFeature_RNA
mean(tmp>4000)
mean(tmp<500)

# filter spots that have unique feature counts over 4000 or less than 500
high_det <- WhichCells(merged_st_obj, expression = nFeature_RNA > 4000)
low_det <- WhichCells(merged_st_obj, expression = nFeature_RNA < 500)

# remove these cells
merged_st_obj_tmp <- subset(merged_st_obj, 
                            cells=setdiff(WhichCells(merged_st_obj),
                                          c(low_det, high_det)))
# number of cells after filtering
cbind(table(Idents(merged_st_obj)), table(Idents(merged_st_obj_tmp)))

# filter out genes that no longer show any expression in the filtered dataset. 
all_st_obj_filt <- list()
for(l1 in Layers(merged_st_obj_tmp)){
  s1 = gsub("counts.", "", l1)
  print(s1)
  
  filt_matrix <- GetAssayData(merged_st_obj_tmp, "RNA", layer = l1)
  filtered_obj <- CreateSeuratObject(counts = filt_matrix,
                                     project = s1,
                                     min.cells = 3)
  filtered_obj$orig.ident <- s1
  all_st_obj_filt[[s1]] <- filtered_obj
}

## Create a merged Seurat object
merged_st_obj_filt <- merge(x=all_st_obj_filt[[1]], y = c(all_st_obj_filt)[-1])

rm(merged_st_obj)
rm(merged_st_obj_tmp)
rm(filtered_obj, filt_matrix)
rm(all_st_obj_filt)

# Add number of genes per UMI for each cell to metadata
merged_st_obj_filt$log10GenesPerUMI <- log10(merged_st_obj_filt$nFeature_RNA) / log10(merged_st_obj_filt$nCount_RNA)

# Create metadata dataframe
metadata <- merged_st_obj_filt@meta.data

# Rename columns
metadata <- metadata %>%
  dplyr::rename(sid = orig.ident)

# Add additional information to metadata
metadata$X <- rownames(metadata)
metadata$pid <- sapply(str_split(metadata$sid, "_"), "[[", 1)
metadata$rid <- sapply(str_split(metadata$sid, "_"), "[[", 2)
metadata$subtype <- plyr::mapvalues(metadata$pid, 
                                    from = unique(st_meta[,c("type", "patient")])$patient, 
                                    to = unique(st_meta[,c("type", "patient")])$type)
metadata$label <- plyr::mapvalues(metadata$X, 
                                  from = unique(prop_ann_meta[,c("X", "label")])$X, 
                                  to = unique(prop_ann_meta[,c("X", "label")])$label)
metadata$lab <- ifelse(metadata$rid %in% c("C1", "C2", "D1"), "lab_C1_C2_D1", "lab_D2_E1_E2")

# Add metadata back to Seurat object
merged_st_obj_filt@meta.data <- metadata
save(merged_st_obj_filt, file=file.path(st_dir, "merged_st_obj_filt.RData"))


# 2. Quality Control ---------------------------------------------------------
# QC
Idents(merged_st_obj_filt) <- "sid"
table(Idents(merged_st_obj_filt))
colors_subtype <- scales::hue_pal()(5)
colors_subtype <- setNames(colors_subtype, unique(metadata$subtype))

pdf("QC_filt.pdf", width = 16, height = 8)

# Visualize QC metrics as a violin plot
VlnPlot(merged_st_obj_filt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_filt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_filt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "lab")

FeatureScatter(merged_st_obj_filt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "subtype",
               split.by = "sid", ncol = 9,
               cols = colors_subtype)

FeatureScatter(merged_st_obj_filt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "lab",
               split.by = "sid", ncol = 9)
dev.off()

# 3. Run normalization ----------------------------------------------------

# SCTransform-normalization
merged_st_obj_filt <- SCTransform(merged_st_obj_filt, return.only.var.genes = F)
merged_st_obj_filt <- RunPCA(merged_st_obj_filt, reduction.name = "pca.unintegrated.sct")
merged_st_obj_filt <- FindNeighbors(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.sct")
merged_st_obj_filt <- FindClusters(merged_st_obj_filt, resolution = 0.1, cluster.name = "unintegrated.sct_clusters_0.1")
merged_st_obj_filt <- RunUMAP(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.sct", reduction.name = "umap.unintegrated.sct")

# Log-normalization
DefaultAssay(merged_st_obj_filt) <- "RNA"
merged_st_obj_filt <- NormalizeData(merged_st_obj_filt)
merged_st_obj_filt <- FindVariableFeatures(merged_st_obj_filt)
# top10 <- head(VariableFeatures(merged_st_obj_filt), 10) # the 10 most highly variable genes
# plot variable features without labels
# VariableFeaturePlot(merged_st_obj_filt)

merged_st_obj_filt <- ScaleData(merged_st_obj_filt, 
                                assay="RNA", 
                                features = rownames(merged_st_obj_filt))
merged_st_obj_filt <- RunPCA(merged_st_obj_filt, reduction.name = "pca.unintegrated.lognorm")
# print(merged_st_obj_filt[["pca.unintegrated.lognorm"]], dims = 1:5, nfeatures = 10)

# VizDimLoadings(merged_st_obj_filt, dims = 1:2, reduction = "pca.unintegrated.lognorm")
# DimPlot(merged_st_obj_filt, reduction = "pca.unintegrated.lognorm")
# DimHeatmap(merged_st_obj_filt, reduction = "pca.unintegrated.lognorm", dims = 1:15, cells = 500, balanced = TRUE)

merged_st_obj_filt <- FindNeighbors(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.lognorm")
merged_st_obj_filt <- FindClusters(merged_st_obj_filt, resolution = 0.1, cluster.name = "unintegrated.lognorm_clusters_0.1")

merged_st_obj_filt <- RunUMAP(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.lognorm", 
                         reduction.name = "umap.unintegrated.lognorm")

# DimPlot(merged_st_obj_filt, reduction = "umap.unintegrated.lognorm", group.by = c("sid"))

save(merged_st_obj_filt, file = file.path(st_dir, "merged_st_obj_normalized_filt.RData"))


# 4. Run integration ------------------------------------------------------
integrated_st_obj_filt <- merged_st_obj_filt
rm(merged_st_obj_filt)

DefaultAssay(integrated_st_obj_filt) <- "RNA"
integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt,
  method = CCAIntegration,
  orig.reduction = "pca.unintegrated.lognorm",
  new.reduction = "integrated.cca.lognorm"
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = RPCAIntegration,
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.rpca.lognorm",
  k.weight = 80
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = HarmonyIntegration,
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.harmony.lognorm",
  verbose = T
)

DefaultAssay(integrated_st_obj_filt) <- "SCT"
integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt,
  method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct",
  new.reduction = "integrated.cca.sct"
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.rpca.sct",
  k.weight = 75
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.harmony.sct",
  verbose = T
)

save(integrated_st_obj_filt, file = file.path(st_dir, "integrated_st_obj_filt.RData"))

# 5. Clustering -----------------------------------------------------------

integrated_st_obj_clusters_filt <- integrated_st_obj_filt
rm(integrated_st_obj_filt)

reductions <- names(integrated_st_obj_clusters_filt@reductions)[startsWith(names(integrated_st_obj_clusters_filt@reductions), "integrated")]

DefaultAssay(integrated_st_obj_clusters_filt) <- "RNA"
for(red in reductions[1:3]){
  integrated_st_obj_clusters_filt <- FindNeighbors(integrated_st_obj_clusters_filt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_filt <- RunUMAP(integrated_st_obj_clusters_filt, reduction = red, 
                               dims = 1:30, reduction.name = paste0("umap.", red))
}

DefaultAssay(integrated_st_obj_clusters_filt) <- "SCT"
for(red in reductions[4:6]){
  integrated_st_obj_clusters_filt <- FindNeighbors(integrated_st_obj_clusters_filt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_filt <- RunUMAP(integrated_st_obj_clusters_filt, reduction = red, 
                               dims = 1:30, reduction.name = paste0("umap.", red))
}

# cluster_methods <- colnames(integrated_st_obj_clusters_filt@meta.data)[startsWith(colnames(integrated_st_obj_clusters_filt@meta.data), "integrated")]
# 
# for (m in cluster_methods) {
#   print(m)
#   print(table(integrated_st_obj_clusters_filt@meta.data[, m]))
# }

DefaultAssay(integrated_st_obj_clusters_filt) <- "RNA"
integrated_st_obj_clusters_filt <- JoinLayers(integrated_st_obj_clusters_filt)

DefaultAssay(integrated_st_obj_clusters_filt) <- "SCT"
integrated_st_obj_clusters_filt <- PrepSCTFindMarkers(integrated_st_obj_clusters_filt,
                                                        assay = "SCT", verbose = TRUE)

save(integrated_st_obj_clusters_filt, file = file.path(st_dir, "integrated_st_obj_clusters_filt.RData"))


# 6. Identify differential expressed genes across conditions --------------

cluster_methods <- sort(colnames(integrated_st_obj_clusters_filt@meta.data[startsWith(colnames(integrated_st_obj_clusters_filt@meta.data), "integrated")]))
seurat_integrated <- integrated_st_obj_clusters_filt
rm(integrated_st_obj_clusters_filt)

markers_identified_filt <- list()

for(cm in cluster_methods){
  m <- strsplit(cm, "_")[[1]][1]
  norm <- strsplit(m, ".", fixed = TRUE)[[1]][3]
  
  if(norm == "lognorm"){
    DefaultAssay(seurat_integrated) <- "RNA"
  }else if(norm == "sct"){
    DefaultAssay(seurat_integrated) <- "SCT"
  }else{
    break
  }
  
  Idents(seurat_integrated) <- cm
  markers <- FindAllMarkers(seurat_integrated, verbose = T)
  
  markers_identified_filt[[cm]] <- markers
}

save(markers_identified_filt, file = file.path(st_dir, "DE_markers_identified_filt.RData"))

top10_DE_markers_filt <- list()
for (cm in cluster_methods) {
  markers_cm <- markers_identified_filt[[cm]] 
  
  if(nrow(markers_cm)==0){
    top10_DE_markers_filt[[cm]] <- data.frame()
  }else{
    top10 <- markers_cm %>%
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < 0.001 & avg_log2FC > 0) %>%
      slice_max(avg_log2FC, n = 10) %>%
      ungroup()
    
    top10$gene_symbol <- plyr::mapvalues(top10$gene, 
                                         from = id_symbol_mapping$id,
                                         to = id_symbol_mapping$symbol)
    
    top10_DE_markers_filt[[cm]] <- top10
  }
}
save(top10_DE_markers_filt, file = file.path(st_dir, "top10_DE_markers_identified_filt.RData"))


# ========== FILTERED (zeros per sample, union) ========== ----------------

load(file.path(st_dir, "merged_st_obj_unfilt.RData"))

# Gene detection filtering
ggplot() +
  geom_density(aes(x=log10(merged_st_obj$nFeature_RNA)))

tmp <- merged_st_obj$nFeature_RNA
mean(tmp>4000)
mean(tmp<500)

# filter spots that have unique feature counts over 4000 or less than 500
high_det <- WhichCells(merged_st_obj, expression = nFeature_RNA > 4000)
low_det <- WhichCells(merged_st_obj, expression = nFeature_RNA < 500)

# remove these cells
merged_st_obj_tmp <- subset(merged_st_obj, 
                            cells=setdiff(WhichCells(merged_st_obj),
                                          c(low_det, high_det)))
# number of cells after filtering
cbind(table(Idents(merged_st_obj)), table(Idents(merged_st_obj_tmp)))


# filter out genes that no longer show any expression in the filtered dataset. 
all_st_obj_filt <- list()
pdf("zero_proportions.pdf", width = 8, height = 5)
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
  filtered_obj <- CreateSeuratObject(counts = filt_matrix,
                                     project = s1,
                                     min.cells = 0.1*ncol(filt_matrix))
  filtered_obj$orig.ident <- s1
  all_st_obj_filt[[s1]] <- filtered_obj
}
dev.off()

## Create a merged Seurat object
merged_st_obj_filt <- merge(x=all_st_obj_filt[[1]], y = c(all_st_obj_filt)[-1])


rm(merged_st_obj)
rm(merged_st_obj_tmp)
rm(filtered_obj, filt_matrix)
rm(all_st_obj_filt)

# Add number of genes per UMI for each cell to metadata
merged_st_obj_filt$log10GenesPerUMI <- log10(merged_st_obj_filt$nFeature_RNA) / log10(merged_st_obj_filt$nCount_RNA)

# Create metadata dataframe
metadata <- merged_st_obj_filt@meta.data

# Rename columns
metadata <- metadata %>%
  dplyr::rename(sid = orig.ident)

# Add additional information to metadata
metadata$X <- rownames(metadata)
metadata$pid <- sapply(str_split(metadata$sid, "_"), "[[", 1)
metadata$rid <- sapply(str_split(metadata$sid, "_"), "[[", 2)
metadata$subtype <- plyr::mapvalues(metadata$pid, 
                                    from = unique(st_meta[,c("type", "patient")])$patient, 
                                    to = unique(st_meta[,c("type", "patient")])$type)
metadata$label <- plyr::mapvalues(metadata$X, 
                                  from = unique(prop_ann_meta[,c("X", "label")])$X, 
                                  to = unique(prop_ann_meta[,c("X", "label")])$label)
metadata$lab <- ifelse(metadata$rid %in% c("C1", "C2", "D1"), "lab_C1_C2_D1", "lab_D2_E1_E2")

# Add metadata back to Seurat object
merged_st_obj_filt@meta.data <- metadata
save(merged_st_obj_filt, file=file.path(st_dir, "merged_st_obj_filt_zeros.RData"))


# 2. Quality Control ---------------------------------------------------------
# QC
Idents(merged_st_obj_filt) <- "sid"
table(Idents(merged_st_obj_filt))
colors_subtype <- scales::hue_pal()(5)
colors_subtype <- setNames(colors_subtype, unique(metadata$subtype))

pdf("QC_filt_zeros.pdf", width = 16, height = 8)

# Visualize QC metrics as a violin plot
VlnPlot(merged_st_obj_filt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_filt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_filt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "lab")

FeatureScatter(merged_st_obj_filt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "subtype",
               split.by = "sid", ncol = 9,
               cols = colors_subtype)

FeatureScatter(merged_st_obj_filt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "lab",
               split.by = "sid", ncol = 9)
dev.off()

# 3. Run normalization ----------------------------------------------------

# SCTransform-normalization
merged_st_obj_filt <- SCTransform(merged_st_obj_filt, return.only.var.genes = FALSE)
merged_st_obj_filt <- RunPCA(merged_st_obj_filt, reduction.name = "pca.unintegrated.sct")
merged_st_obj_filt <- FindNeighbors(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.sct")
merged_st_obj_filt <- FindClusters(merged_st_obj_filt, resolution = 0.1, cluster.name = "unintegrated.sct_clusters_0.1")
merged_st_obj_filt <- RunUMAP(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.sct", reduction.name = "umap.unintegrated.sct")

# Log-normalization
DefaultAssay(merged_st_obj_filt) <- "RNA"
merged_st_obj_filt <- NormalizeData(merged_st_obj_filt)
merged_st_obj_filt <- FindVariableFeatures(merged_st_obj_filt, nfeatures = 2000)
# top10 <- head(VariableFeatures(merged_st_obj_filt), 10) # the 10 most highly variable genes
# plot variable features without labels
# VariableFeaturePlot(merged_st_obj_filt)

merged_st_obj_filt <- ScaleData(merged_st_obj_filt, 
                                assay="RNA", 
                                features = rownames(merged_st_obj_filt))
merged_st_obj_filt <- RunPCA(merged_st_obj_filt, reduction.name = "pca.unintegrated.lognorm")
# print(merged_st_obj_filt[["pca.unintegrated.lognorm"]], dims = 1:5, nfeatures = 10)

# VizDimLoadings(merged_st_obj_filt, dims = 1:2, reduction = "pca.unintegrated.lognorm")
# DimPlot(merged_st_obj_filt, reduction = "pca.unintegrated.lognorm")
# DimHeatmap(merged_st_obj_filt, reduction = "pca.unintegrated.lognorm", dims = 1:15, cells = 500, balanced = TRUE)

merged_st_obj_filt <- FindNeighbors(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.lognorm")
merged_st_obj_filt <- FindClusters(merged_st_obj_filt, resolution = 0.1, cluster.name = "unintegrated.lognorm_clusters_0.1")

merged_st_obj_filt <- RunUMAP(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.lognorm", 
                              reduction.name = "umap.unintegrated.lognorm")

# DimPlot(merged_st_obj_filt, reduction = "umap.unintegrated.lognorm", group.by = c("sid"))

save(merged_st_obj_filt, file = file.path(st_dir, "merged_st_obj_normalized_filt_zeros.RData"))


# 4. Run integration ------------------------------------------------------

integrated_st_obj_filt <- merged_st_obj_filt
rm(merged_st_obj_filt)

DefaultAssay(integrated_st_obj_filt) <- "RNA"
integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt,
  method = CCAIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm",
  new.reduction = "integrated.cca.lognorm"
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = RPCAIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.rpca.lognorm",
  k.weight = 80
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = HarmonyIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.harmony.lognorm",
  verbose = T
)

DefaultAssay(integrated_st_obj_filt) <- "SCT"
integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt,
  method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct",
  new.reduction = "integrated.cca.sct"
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.rpca.sct",
  k.weight = 75
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.harmony.sct",
  verbose = T
)

save(integrated_st_obj_filt, file = file.path(st_dir, "integrated_st_obj_filt_zeros.RData"))

# 5. Clustering -----------------------------------------------------------

integrated_st_obj_clusters_filt <- integrated_st_obj_filt
rm(integrated_st_obj_filt)

reductions <- names(integrated_st_obj_clusters_filt@reductions)[startsWith(names(integrated_st_obj_clusters_filt@reductions), "integrated")]

DefaultAssay(integrated_st_obj_clusters_filt) <- "RNA"
for(red in reductions[1:3]){
  integrated_st_obj_clusters_filt <- FindNeighbors(integrated_st_obj_clusters_filt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_filt <- RunUMAP(integrated_st_obj_clusters_filt, reduction = red, 
                                             dims = 1:30, reduction.name = paste0("umap.", red))
}

DefaultAssay(integrated_st_obj_clusters_filt) <- "SCT"
for(red in reductions[4:6]){
  integrated_st_obj_clusters_filt <- FindNeighbors(integrated_st_obj_clusters_filt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_filt <- RunUMAP(integrated_st_obj_clusters_filt, reduction = red, 
                                             dims = 1:30, reduction.name = paste0("umap.", red))
}

# cluster_methods <- colnames(integrated_st_obj_clusters_filt@meta.data)[startsWith(colnames(integrated_st_obj_clusters_filt@meta.data), "integrated")]
# 
# for (m in cluster_methods) {
#   print(m)
#   print(table(integrated_st_obj_clusters_filt@meta.data[, m]))
# }

DefaultAssay(integrated_st_obj_clusters_filt) <- "RNA"
integrated_st_obj_clusters_filt <- JoinLayers(integrated_st_obj_clusters_filt)

DefaultAssay(integrated_st_obj_clusters_filt) <- "SCT"
integrated_st_obj_clusters_filt <- PrepSCTFindMarkers(integrated_st_obj_clusters_filt,
                                                      assay = "SCT", verbose = TRUE)

save(integrated_st_obj_clusters_filt, file = file.path(st_dir, "integrated_st_obj_clusters_filt_zeros.RData"))


# 6. Identify differential expressed genes across conditions --------------

cluster_methods <- sort(colnames(integrated_st_obj_clusters_filt@meta.data[startsWith(colnames(integrated_st_obj_clusters_filt@meta.data), "integrated")]))
seurat_integrated <- integrated_st_obj_clusters_filt
rm(integrated_st_obj_clusters_filt)

markers_identified_filt <- list()

for(cm in cluster_methods){
  m <- strsplit(cm, "_")[[1]][1]
  norm <- strsplit(m, ".", fixed = TRUE)[[1]][3]
  
  if(norm == "lognorm"){
    DefaultAssay(seurat_integrated) <- "RNA"
  }else if(norm == "sct"){
    DefaultAssay(seurat_integrated) <- "SCT"
  }else{
    break
  }
  
  Idents(seurat_integrated) <- cm
  markers <- FindAllMarkers(seurat_integrated, verbose = T)
  
  markers_identified_filt[[cm]] <- markers
}

save(markers_identified_filt, file = file.path(st_dir, "DE_markers_identified_filt_zeros.RData"))

top10_DE_markers_filt <- list()
for (cm in cluster_methods) {
  markers_cm <- markers_identified_filt[[cm]] 
  
  if(nrow(markers_cm)==0){
    top10_DE_markers_filt[[cm]] <- data.frame()
  }else{
    top10 <- markers_cm %>%
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < 0.001 & avg_log2FC > 0) %>%
      slice_max(avg_log2FC, n = 10) %>%
      ungroup()
    
    top10$gene_symbol <- plyr::mapvalues(top10$gene, 
                                         from = id_symbol_mapping$id,
                                         to = id_symbol_mapping$symbol)
    
    top10_DE_markers_filt[[cm]] <- top10
  }
}
save(top10_DE_markers_filt, file = file.path(st_dir, "top10_DE_markers_identified_filt_zeros.RData"))

# ========== FILTERED (zeros per sample, intersection, 300 hvgs) ========== ----------------

load(file.path(st_dir, "merged_st_obj_unfilt.RData"))

# Gene detection filtering
ggplot() +
  geom_density(aes(x=log10(merged_st_obj$nFeature_RNA)))

tmp <- merged_st_obj$nFeature_RNA
mean(tmp>4000)
mean(tmp<500)

# filter spots that have unique feature counts over 4000 or less than 500
high_det <- WhichCells(merged_st_obj, expression = nFeature_RNA > 4000)
low_det <- WhichCells(merged_st_obj, expression = nFeature_RNA < 500)

# remove these cells
merged_st_obj_tmp <- subset(merged_st_obj, 
                            cells=setdiff(WhichCells(merged_st_obj),
                                          c(low_det, high_det)))
# number of cells after filtering
cbind(table(Idents(merged_st_obj)), table(Idents(merged_st_obj_tmp)))

# intersection of genes that are expressed in more than 10% of spots in each sample.
high_expr_genes <- list()
for(l1 in Layers(merged_st_obj_tmp)){
  s1 = gsub("counts.", "", l1)
  print(s1)
  
  filt_matrix <- GetAssayData(merged_st_obj_tmp, "RNA", layer = l1)
  zero_prop <- rowSums(filt_matrix == 0)/ncol(filt_matrix)
  filt_matrix <- GetAssayData(merged_st_obj_tmp, "RNA", layer = l1)
  high_expr_genes[[l1]] <- rownames(filt_matrix)[zero_prop < 0.9]
}

common_high_expr_genes <- Reduce(intersect, high_expr_genes)

all_st_obj_filt <- list()
for(l1 in Layers(merged_st_obj_tmp)){
  s1 = gsub("counts.", "", l1)
  print(s1)
  
  filt_matrix <- GetAssayData(merged_st_obj_tmp, "RNA", layer = l1)
  
  filt_matrix <- filt_matrix[common_high_expr_genes,]
  zero_prop <- rowSums(filt_matrix == 0)/ncol(filt_matrix)
  
  p <- ggplot() +
    geom_histogram(aes(x=zero_prop), binwidth = 0.01) +
    geom_vline(xintercept = 0.9, col = "red") +
    geom_vline(xintercept = 0.8, col = "blue") +
    labs(title = s1) 
  print(p)
  
  filtered_obj <- CreateSeuratObject(counts = filt_matrix,
                                     project = s1,
                                     min.cells = 0)
  filtered_obj$orig.ident <- s1
  all_st_obj_filt[[s1]] <- filtered_obj
}

## Create a merged Seurat object
merged_st_obj_filt <- merge(x=all_st_obj_filt[[1]], y = c(all_st_obj_filt)[-1])

rm(merged_st_obj)
rm(merged_st_obj_tmp)
rm(filtered_obj, filt_matrix)
rm(all_st_obj_filt)

# Add number of genes per UMI for each cell to metadata
merged_st_obj_filt$log10GenesPerUMI <- log10(merged_st_obj_filt$nFeature_RNA) / log10(merged_st_obj_filt$nCount_RNA)

# Create metadata dataframe
metadata <- merged_st_obj_filt@meta.data

# Rename columns
metadata <- metadata %>%
  dplyr::rename(sid = orig.ident)

# Add additional information to metadata
metadata$X <- rownames(metadata)
metadata$pid <- sapply(str_split(metadata$sid, "_"), "[[", 1)
metadata$rid <- sapply(str_split(metadata$sid, "_"), "[[", 2)
metadata$subtype <- plyr::mapvalues(metadata$pid, 
                                    from = unique(st_meta[,c("type", "patient")])$patient, 
                                    to = unique(st_meta[,c("type", "patient")])$type)
metadata$label <- plyr::mapvalues(metadata$X, 
                                  from = unique(prop_ann_meta[,c("X", "label")])$X, 
                                  to = unique(prop_ann_meta[,c("X", "label")])$label)
metadata$lab <- ifelse(metadata$rid %in% c("C1", "C2", "D1"), "lab_C1_C2_D1", "lab_D2_E1_E2")

# Add metadata back to Seurat object
merged_st_obj_filt@meta.data <- metadata
merged_st_obj_filt_zeros_intersect <- merged_st_obj_filt
save(merged_st_obj_filt_zeros_intersect, file=file.path(st_dir, "merged_st_obj_filt_zeros_intersect.RData"))
rm(merged_st_obj_filt_zeros_intersect)

# 2. Quality Control ---------------------------------------------------------
# QC
Idents(merged_st_obj_filt) <- "sid"
table(Idents(merged_st_obj_filt))
colors_subtype <- scales::hue_pal()(5)
colors_subtype <- setNames(colors_subtype, unique(metadata$subtype))

pdf("QC_filt_zeros_intersect.pdf", width = 16, height = 8)

# Visualize QC metrics as a violin plot
VlnPlot(merged_st_obj_filt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_filt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_filt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "lab")

FeatureScatter(merged_st_obj_filt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "subtype",
               split.by = "sid", ncol = 9,
               cols = colors_subtype)

FeatureScatter(merged_st_obj_filt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "lab",
               split.by = "sid", ncol = 9)
dev.off()

# 3. Run normalization ----------------------------------------------------

# SCTransform-normalization
merged_st_obj_filt <- SCTransform(merged_st_obj_filt, return.only.var.genes = FALSE)
merged_st_obj_filt <- RunPCA(merged_st_obj_filt, reduction.name = "pca.unintegrated.sct")
merged_st_obj_filt <- FindNeighbors(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.sct")
merged_st_obj_filt <- FindClusters(merged_st_obj_filt, resolution = 0.1, cluster.name = "unintegrated.sct_clusters_0.1")
merged_st_obj_filt <- RunUMAP(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.sct", reduction.name = "umap.unintegrated.sct")

# Log-normalization
DefaultAssay(merged_st_obj_filt) <- "RNA"
merged_st_obj_filt <- NormalizeData(merged_st_obj_filt)
merged_st_obj_filt <- FindVariableFeatures(merged_st_obj_filt, nfeatures = 300)
# top10 <- head(VariableFeatures(merged_st_obj_filt), 10) # the 10 most highly variable genes
# plot variable features without labels
# VariableFeaturePlot(merged_st_obj_filt)

merged_st_obj_filt <- ScaleData(merged_st_obj_filt, 
                                assay="RNA", 
                                features = rownames(merged_st_obj_filt))
merged_st_obj_filt <- RunPCA(merged_st_obj_filt, reduction.name = "pca.unintegrated.lognorm")
# print(merged_st_obj_filt[["pca.unintegrated.lognorm"]], dims = 1:5, nfeatures = 10)

# VizDimLoadings(merged_st_obj_filt, dims = 1:2, reduction = "pca.unintegrated.lognorm")
# DimPlot(merged_st_obj_filt, reduction = "pca.unintegrated.lognorm")
# DimHeatmap(merged_st_obj_filt, reduction = "pca.unintegrated.lognorm", dims = 1:15, cells = 500, balanced = TRUE)

merged_st_obj_filt <- FindNeighbors(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.lognorm")
merged_st_obj_filt <- FindClusters(merged_st_obj_filt, resolution = 0.1, cluster.name = "unintegrated.lognorm_clusters_0.1")

merged_st_obj_filt <- RunUMAP(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.lognorm", 
                              reduction.name = "umap.unintegrated.lognorm")

# DimPlot(merged_st_obj_filt, reduction = "umap.unintegrated.lognorm", group.by = c("sid"))

merged_st_obj_filt_zeros_intersect <- merged_st_obj_filt
save(merged_st_obj_filt_zeros_intersect, file = file.path(st_dir, "merged_st_obj_normalized_filt_zeros_intersect_300hvg.RData"))
rm(merged_st_obj_filt_zeros_intersect)

# 4. Run integration ------------------------------------------------------

integrated_st_obj_filt <- merged_st_obj_filt
rm(merged_st_obj_filt)

DefaultAssay(integrated_st_obj_filt) <- "RNA"
integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt,
  method = CCAIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm",
  new.reduction = "integrated.cca.lognorm"
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = RPCAIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.rpca.lognorm",
  k.weight = 80
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = HarmonyIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.harmony.lognorm",
  verbose = T
)

DefaultAssay(integrated_st_obj_filt) <- "SCT"
integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt,
  method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct",
  new.reduction = "integrated.cca.sct"
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.rpca.sct",
  k.weight = 75
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.harmony.sct",
  verbose = T
)

integrated_st_obj_filt_zeros_intersect <- integrated_st_obj_filt
save(integrated_st_obj_filt_zeros_intersect, file = file.path(st_dir, "integrated_st_obj_filt_zeros_intersect.RData"))
rm(integrated_st_obj_filt_zeros_intersect)

# 5. Clustering -----------------------------------------------------------

integrated_st_obj_clusters_filt <- integrated_st_obj_filt
rm(integrated_st_obj_filt)

reductions <- names(integrated_st_obj_clusters_filt@reductions)[startsWith(names(integrated_st_obj_clusters_filt@reductions), "integrated")]

DefaultAssay(integrated_st_obj_clusters_filt) <- "RNA"
for(red in reductions[1:3]){
  integrated_st_obj_clusters_filt <- FindNeighbors(integrated_st_obj_clusters_filt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_filt <- RunUMAP(integrated_st_obj_clusters_filt, reduction = red, 
                                             dims = 1:30, reduction.name = paste0("umap.", red))
}

DefaultAssay(integrated_st_obj_clusters_filt) <- "SCT"
for(red in reductions[4:6]){
  integrated_st_obj_clusters_filt <- FindNeighbors(integrated_st_obj_clusters_filt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_filt <- RunUMAP(integrated_st_obj_clusters_filt, reduction = red, 
                                             dims = 1:30, reduction.name = paste0("umap.", red))
}

# cluster_methods <- colnames(integrated_st_obj_clusters_filt@meta.data)[startsWith(colnames(integrated_st_obj_clusters_filt@meta.data), "integrated")]
# 
# for (m in cluster_methods) {
#   print(m)
#   print(table(integrated_st_obj_clusters_filt@meta.data[, m]))
# }

DefaultAssay(integrated_st_obj_clusters_filt) <- "RNA"
integrated_st_obj_clusters_filt <- JoinLayers(integrated_st_obj_clusters_filt)

DefaultAssay(integrated_st_obj_clusters_filt) <- "SCT"
integrated_st_obj_clusters_filt <- PrepSCTFindMarkers(integrated_st_obj_clusters_filt,
                                                      assay = "SCT", verbose = TRUE)

integrated_st_obj_clusters_filt_zeros_intersect <- integrated_st_obj_clusters_filt
save(integrated_st_obj_clusters_filt_zeros_intersect, file = file.path(st_dir, "integrated_st_obj_clusters_filt_zeros_intersect.RData"))
rm(integrated_st_obj_clusters_filt_zeros_intersect)

# 6. Identify differential expressed genes across conditions --------------

cluster_methods <- sort(colnames(integrated_st_obj_clusters_filt@meta.data[startsWith(colnames(integrated_st_obj_clusters_filt@meta.data), "integrated")]))
seurat_integrated <- integrated_st_obj_clusters_filt
rm(integrated_st_obj_clusters_filt)

markers_identified_filt <- list()

for(cm in cluster_methods){
  m <- strsplit(cm, "_")[[1]][1]
  norm <- strsplit(m, ".", fixed = TRUE)[[1]][3]
  
  if(norm == "lognorm"){
    DefaultAssay(seurat_integrated) <- "RNA"
  }else if(norm == "sct"){
    DefaultAssay(seurat_integrated) <- "SCT"
  }else{
    break
  }
  
  Idents(seurat_integrated) <- cm
  markers <- FindAllMarkers(seurat_integrated, verbose = T)
  
  markers_identified_filt[[cm]] <- markers
}

markers_identified_filt_zeros_intersect <- markers_identified_filt
save(markers_identified_filt_zeros_intersect, file = file.path(st_dir, "DE_markers_identified_filt_zeros_intersect.RData"))
rm(markers_identified_filt_zeros_intersect)

top10_DE_markers_filt <- list()
for (cm in cluster_methods) {
  markers_cm <- markers_identified_filt[[cm]] 
  
  if(nrow(markers_cm)==0){
    top10_DE_markers_filt[[cm]] <- data.frame()
  }else{
    top10 <- markers_cm %>%
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < 0.001 & avg_log2FC > 0) %>%
      slice_max(avg_log2FC, n = 10) %>%
      ungroup()
    
    top10$gene_symbol <- plyr::mapvalues(top10$gene, 
                                         from = id_symbol_mapping$id,
                                         to = id_symbol_mapping$symbol)
    
    top10_DE_markers_filt[[cm]] <- top10
  }
}
top10_DE_markers_filt_zeros_intersect <- top10_DE_markers_filt
save(top10_DE_markers_filt_zeros_intersect, file = file.path(st_dir, "top10_DE_markers_identified_filt_zeros_intersect.RData"))
rm(top10_DE_markers_filt_zeros_intersect)


# ========== FILTERED (zeros all samples, 90%, 2000 hvgs) ========== ----------------------

load(file.path(st_dir, "merged_st_obj_unfilt.RData"))

# Gene detection filtering
ggplot() +
  geom_density(aes(x=log10(merged_st_obj$nFeature_RNA)))

tmp <- merged_st_obj$nFeature_RNA
mean(tmp>4000)
mean(tmp<500)

# filter spots that have unique feature counts over 4000 or less than 500
high_det <- WhichCells(merged_st_obj, expression = nFeature_RNA > 4000)
low_det <- WhichCells(merged_st_obj, expression = nFeature_RNA < 500)

# remove these cells
merged_st_obj_tmp <- subset(merged_st_obj, 
                            cells=setdiff(WhichCells(merged_st_obj),
                                          c(low_det, high_det)))
# number of cells after filtering
cbind(table(Idents(merged_st_obj)), table(Idents(merged_st_obj_tmp)))


# filter out genes that no longer show any expression in the filtered dataset. 
merged_st_obj_tmp <- JoinLayers(merged_st_obj_tmp)
filt_matrix <- GetAssayData(merged_st_obj_tmp, "RNA", layer = "counts")
zero_prop <- rowSums(filt_matrix == 0)/ncol(filt_matrix)

pdf("zero_proportions_all.pdf", width = 8, height = 5)
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

filtered_obj <- CreateSeuratObject(counts = filt_matrix,
                                   project = "all",
                                   min.cells = 0.1*ncol(filt_matrix))

meta <- merged_st_obj_tmp@meta.data
meta_f <- filtered_obj@meta.data[,-1]
meta_f <- cbind(meta_f, meta[rownames(meta_f), !colnames(meta) %in% c("nCount_RNA", "nFeature_RNA")])
meta_f$log10GenesPerUMI <- log10(meta_f$nFeature_RNA) / log10(meta_f$nCount_RNA)

filtered_obj@meta.data <- meta_f

filtered_obj[["RNA"]] <- split(filtered_obj[["RNA"]], f = filtered_obj$sid)
merged_st_obj_filt <- filtered_obj

rm(merged_st_obj)
rm(merged_st_obj_tmp)
rm(filtered_obj, filt_matrix)

save(merged_st_obj_filt, file=file.path(st_dir, "merged_st_obj_filt_zeros_all.RData"))


# 2. Quality Control ---------------------------------------------------------
# QC
Idents(merged_st_obj_filt) <- "sid"
table(Idents(merged_st_obj_filt))
colors_subtype <- scales::hue_pal()(5)
colors_subtype <- setNames(colors_subtype, unique(metadata$subtype))

pdf("QC_filt_zeros_all.pdf", width = 16, height = 8)

# Visualize QC metrics as a violin plot
VlnPlot(merged_st_obj_filt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_filt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_filt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "lab")

FeatureScatter(merged_st_obj_filt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "subtype",
               split.by = "sid", ncol = 9,
               cols = colors_subtype)

FeatureScatter(merged_st_obj_filt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "lab",
               split.by = "sid", ncol = 9)
dev.off()

# 3. Run normalization ----------------------------------------------------

# SCTransform-normalization
merged_st_obj_filt <- SCTransform(merged_st_obj_filt, return.only.var.genes = FALSE)
merged_st_obj_filt <- RunPCA(merged_st_obj_filt, reduction.name = "pca.unintegrated.sct")
merged_st_obj_filt <- FindNeighbors(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.sct")
merged_st_obj_filt <- FindClusters(merged_st_obj_filt, resolution = 0.1, cluster.name = "unintegrated.sct_clusters_0.1")
merged_st_obj_filt <- RunUMAP(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.sct", reduction.name = "umap.unintegrated.sct")

# Log-normalization
DefaultAssay(merged_st_obj_filt) <- "RNA"
merged_st_obj_filt <- NormalizeData(merged_st_obj_filt)
merged_st_obj_filt <- FindVariableFeatures(merged_st_obj_filt, nfeatures = 2000)
# top10 <- head(VariableFeatures(merged_st_obj_filt), 10) # the 10 most highly variable genes
# plot variable features without labels
# VariableFeaturePlot(merged_st_obj_filt)

merged_st_obj_filt <- ScaleData(merged_st_obj_filt, 
                                assay="RNA", 
                                features = rownames(merged_st_obj_filt))
merged_st_obj_filt <- RunPCA(merged_st_obj_filt, reduction.name = "pca.unintegrated.lognorm")
# print(merged_st_obj_filt[["pca.unintegrated.lognorm"]], dims = 1:5, nfeatures = 10)

# VizDimLoadings(merged_st_obj_filt, dims = 1:2, reduction = "pca.unintegrated.lognorm")
# DimPlot(merged_st_obj_filt, reduction = "pca.unintegrated.lognorm")
# DimHeatmap(merged_st_obj_filt, reduction = "pca.unintegrated.lognorm", dims = 1:15, cells = 500, balanced = TRUE)

merged_st_obj_filt <- FindNeighbors(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.lognorm")
merged_st_obj_filt <- FindClusters(merged_st_obj_filt, resolution = 0.1, cluster.name = "unintegrated.lognorm_clusters_0.1")

merged_st_obj_filt <- RunUMAP(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.lognorm", 
                              reduction.name = "umap.unintegrated.lognorm")

# DimPlot(merged_st_obj_filt, reduction = "umap.unintegrated.lognorm", group.by = c("sid"))

save(merged_st_obj_filt, file = file.path(st_dir, "merged_st_obj_normalized_filt_zeros_all.RData"))


# 4. Run integration ------------------------------------------------------

integrated_st_obj_filt <- merged_st_obj_filt
rm(merged_st_obj_filt)

DefaultAssay(integrated_st_obj_filt) <- "RNA"
integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt,
  method = CCAIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm",
  new.reduction = "integrated.cca.lognorm"
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = RPCAIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.rpca.lognorm",
  k.weight = 80
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = HarmonyIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.harmony.lognorm",
  verbose = T
)

DefaultAssay(integrated_st_obj_filt) <- "SCT"
integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt,
  method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct",
  new.reduction = "integrated.cca.sct"
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.rpca.sct",
  k.weight = 75
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.harmony.sct",
  verbose = T
)

save(integrated_st_obj_filt, file = file.path(st_dir, "integrated_st_obj_filt_zeros_all.RData"))

# 5. Clustering -----------------------------------------------------------

integrated_st_obj_clusters_filt <- integrated_st_obj_filt
rm(integrated_st_obj_filt)

reductions <- names(integrated_st_obj_clusters_filt@reductions)[startsWith(names(integrated_st_obj_clusters_filt@reductions), "integrated")]

DefaultAssay(integrated_st_obj_clusters_filt) <- "RNA"
for(red in reductions[1:3]){
  integrated_st_obj_clusters_filt <- FindNeighbors(integrated_st_obj_clusters_filt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_filt <- RunUMAP(integrated_st_obj_clusters_filt, reduction = red, 
                                             dims = 1:30, reduction.name = paste0("umap.", red))
}

DefaultAssay(integrated_st_obj_clusters_filt) <- "SCT"
for(red in reductions[4:6]){
  integrated_st_obj_clusters_filt <- FindNeighbors(integrated_st_obj_clusters_filt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_filt <- RunUMAP(integrated_st_obj_clusters_filt, reduction = red, 
                                             dims = 1:30, reduction.name = paste0("umap.", red))
}

# cluster_methods <- colnames(integrated_st_obj_clusters_filt@meta.data)[startsWith(colnames(integrated_st_obj_clusters_filt@meta.data), "integrated")]
# 
# for (m in cluster_methods) {
#   print(m)
#   print(table(integrated_st_obj_clusters_filt@meta.data[, m]))
# }

DefaultAssay(integrated_st_obj_clusters_filt) <- "RNA"
integrated_st_obj_clusters_filt <- JoinLayers(integrated_st_obj_clusters_filt)

DefaultAssay(integrated_st_obj_clusters_filt) <- "SCT"
integrated_st_obj_clusters_filt <- PrepSCTFindMarkers(integrated_st_obj_clusters_filt,
                                                      assay = "SCT", verbose = TRUE)

save(integrated_st_obj_clusters_filt, file = file.path(st_dir, "integrated_st_obj_clusters_filt_zeros_all.RData"))


# 6. Identify differential expressed genes across conditions --------------

cluster_methods <- sort(colnames(integrated_st_obj_clusters_filt@meta.data[startsWith(colnames(integrated_st_obj_clusters_filt@meta.data), "integrated")]))
seurat_integrated <- integrated_st_obj_clusters_filt
rm(integrated_st_obj_clusters_filt)

markers_identified_filt <- list()

for(cm in cluster_methods){
  m <- strsplit(cm, "_")[[1]][1]
  norm <- strsplit(m, ".", fixed = TRUE)[[1]][3]
  
  if(norm == "lognorm"){
    DefaultAssay(seurat_integrated) <- "RNA"
  }else if(norm == "sct"){
    DefaultAssay(seurat_integrated) <- "SCT"
  }else{
    break
  }
  
  Idents(seurat_integrated) <- cm
  markers <- FindAllMarkers(seurat_integrated, verbose = T)
  
  markers_identified_filt[[cm]] <- markers
}

save(markers_identified_filt, file = file.path(st_dir, "DE_markers_identified_filt_zeros_all.RData"))

top10_DE_markers_filt <- list()
for (cm in cluster_methods) {
  markers_cm <- markers_identified_filt[[cm]] 
  
  if(nrow(markers_cm)==0){
    top10_DE_markers_filt[[cm]] <- data.frame()
  }else{
    top10 <- markers_cm %>%
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < 0.001 & avg_log2FC > 0) %>%
      slice_max(avg_log2FC, n = 10) %>%
      ungroup()
    
    top10$gene_symbol <- plyr::mapvalues(top10$gene, 
                                         from = id_symbol_mapping$id,
                                         to = id_symbol_mapping$symbol)
    
    top10_DE_markers_filt[[cm]] <- top10
  }
}
save(top10_DE_markers_filt, file = file.path(st_dir, "top10_DE_markers_identified_filt_zeros_all.RData"))



# ========== FILTERED (zeros all samples, 80%, 500 hvgs) ========== ----------------------

load(file.path(st_dir, "merged_st_obj_unfilt.RData"))

# Gene detection filtering
# filter spots that have unique feature counts over 4000 or less than 500
high_det <- WhichCells(merged_st_obj, expression = nFeature_RNA > 4000)
low_det <- WhichCells(merged_st_obj, expression = nFeature_RNA < 500)

# remove these cells
merged_st_obj_tmp <- subset(merged_st_obj, 
                            cells=setdiff(WhichCells(merged_st_obj),
                                          c(low_det, high_det)))
# number of cells after filtering
cbind(table(Idents(merged_st_obj)), table(Idents(merged_st_obj_tmp)))


# filter out genes that no longer show any expression in the filtered dataset. 
merged_st_obj_tmp <- JoinLayers(merged_st_obj_tmp)
filt_matrix <- GetAssayData(merged_st_obj_tmp, "RNA", layer = "counts")
zero_prop <- rowSums(filt_matrix == 0)/ncol(filt_matrix)

filtered_obj <- CreateSeuratObject(counts = filt_matrix,
                                   project = "all",
                                   min.cells = 0.2*ncol(filt_matrix))

meta <- merged_st_obj_tmp@meta.data
meta_f <- filtered_obj@meta.data[,-1]
meta_f <- cbind(meta_f, meta[rownames(meta_f), !colnames(meta) %in% c("nCount_RNA", "nFeature_RNA")])
meta_f$log10GenesPerUMI <- log10(meta_f$nFeature_RNA) / log10(meta_f$nCount_RNA)

filtered_obj@meta.data <- meta_f

filtered_obj[["RNA"]] <- split(filtered_obj[["RNA"]], f = filtered_obj$sid)
merged_st_obj_filt <- filtered_obj

rm(merged_st_obj)
rm(merged_st_obj_tmp)
rm(filtered_obj, filt_matrix)

save(merged_st_obj_filt, file=file.path(st_dir, "merged_st_obj_filt_zeros_all_80.RData"))


# 2. Quality Control ---------------------------------------------------------
# QC
Idents(merged_st_obj_filt) <- "sid"
table(Idents(merged_st_obj_filt))
colors_subtype <- scales::hue_pal()(5)
colors_subtype <- setNames(colors_subtype, unique(meta_f$subtype))

pdf("QC_filt_zeros_all_80.pdf", width = 16, height = 8)

# Visualize QC metrics as a violin plot
VlnPlot(merged_st_obj_filt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("nFeature_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_filt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("nCount_RNA"), 
        ncol = 1,
        split.by = "lab")

VlnPlot(merged_st_obj_filt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "subtype",
        cols = colors_subtype)

VlnPlot(merged_st_obj_filt, 
        features = c("log10GenesPerUMI"), 
        ncol = 1,
        split.by = "lab")

FeatureScatter(merged_st_obj_filt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "subtype",
               split.by = "sid", ncol = 9,
               cols = colors_subtype)

FeatureScatter(merged_st_obj_filt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "lab",
               split.by = "sid", ncol = 9)
dev.off()

# 3. Run normalization ----------------------------------------------------

# SCTransform-normalization
merged_st_obj_filt <- SCTransform(merged_st_obj_filt, return.only.var.genes = FALSE)
merged_st_obj_filt <- RunPCA(merged_st_obj_filt, reduction.name = "pca.unintegrated.sct")
merged_st_obj_filt <- FindNeighbors(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.sct")
merged_st_obj_filt <- FindClusters(merged_st_obj_filt, resolution = 0.1, cluster.name = "unintegrated.sct_clusters_0.1")
merged_st_obj_filt <- RunUMAP(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.sct", reduction.name = "umap.unintegrated.sct")

# Log-normalization
DefaultAssay(merged_st_obj_filt) <- "RNA"
merged_st_obj_filt <- NormalizeData(merged_st_obj_filt)
merged_st_obj_filt <- FindVariableFeatures(merged_st_obj_filt, nfeatures = 500)
# top10 <- head(VariableFeatures(merged_st_obj_filt), 10) # the 10 most highly variable genes
# plot variable features without labels
# VariableFeaturePlot(merged_st_obj_filt)

merged_st_obj_filt <- ScaleData(merged_st_obj_filt, 
                                assay="RNA", 
                                features = rownames(merged_st_obj_filt))
merged_st_obj_filt <- RunPCA(merged_st_obj_filt, reduction.name = "pca.unintegrated.lognorm")
# print(merged_st_obj_filt[["pca.unintegrated.lognorm"]], dims = 1:5, nfeatures = 10)

# VizDimLoadings(merged_st_obj_filt, dims = 1:2, reduction = "pca.unintegrated.lognorm")
# DimPlot(merged_st_obj_filt, reduction = "pca.unintegrated.lognorm")
# DimHeatmap(merged_st_obj_filt, reduction = "pca.unintegrated.lognorm", dims = 1:15, cells = 500, balanced = TRUE)

merged_st_obj_filt <- FindNeighbors(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.lognorm")
merged_st_obj_filt <- FindClusters(merged_st_obj_filt, resolution = 0.1, cluster.name = "unintegrated.lognorm_clusters_0.1")

merged_st_obj_filt <- RunUMAP(merged_st_obj_filt, dims = 1:30, reduction = "pca.unintegrated.lognorm", 
                              reduction.name = "umap.unintegrated.lognorm")

# DimPlot(merged_st_obj_filt, reduction = "umap.unintegrated.lognorm", group.by = c("sid"))

save(merged_st_obj_filt, file = file.path(st_dir, "merged_st_obj_normalized_filt_zeros_all_80_500hvg.RData"))


# 4. Run integration ------------------------------------------------------

integrated_st_obj_filt <- merged_st_obj_filt
rm(merged_st_obj_filt)

DefaultAssay(integrated_st_obj_filt) <- "RNA"
integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt,
  method = CCAIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm",
  new.reduction = "integrated.cca.lognorm"
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = RPCAIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.rpca.lognorm",
  k.weight = 80
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = HarmonyIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "pca.unintegrated.lognorm", 
  new.reduction = "integrated.harmony.lognorm",
  verbose = T
)

DefaultAssay(integrated_st_obj_filt) <- "SCT"
integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt,
  method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct",
  new.reduction = "integrated.cca.sct"
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.rpca.sct",
  k.weight = 75
)

integrated_st_obj_filt <- IntegrateLayers(
  object = integrated_st_obj_filt, 
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca.unintegrated.sct", 
  new.reduction = "integrated.harmony.sct",
  verbose = T
)

save(integrated_st_obj_filt, file = file.path(st_dir, "integrated_st_obj_filt_zeros_all_80.RData"))

# 5. Clustering -----------------------------------------------------------

integrated_st_obj_clusters_filt <- integrated_st_obj_filt
rm(integrated_st_obj_filt)

reductions <- names(integrated_st_obj_clusters_filt@reductions)[startsWith(names(integrated_st_obj_clusters_filt@reductions), "integrated")]

DefaultAssay(integrated_st_obj_clusters_filt) <- "RNA"
for(red in reductions[1:3]){
  integrated_st_obj_clusters_filt <- FindNeighbors(integrated_st_obj_clusters_filt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_filt <- RunUMAP(integrated_st_obj_clusters_filt, reduction = red, 
                                             dims = 1:30, reduction.name = paste0("umap.", red))
}

DefaultAssay(integrated_st_obj_clusters_filt) <- "SCT"
for(red in reductions[4:6]){
  integrated_st_obj_clusters_filt <- FindNeighbors(integrated_st_obj_clusters_filt, reduction = red, dims = 1:30)
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.3, cluster.name = paste0(red,"_clusters_0.3"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.2, cluster.name = paste0(red,"_clusters_0.2"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.1, cluster.name = paste0(red,"_clusters_0.1"))
  integrated_st_obj_clusters_filt <- FindClusters(integrated_st_obj_clusters_filt, resolution = 0.05, cluster.name = paste0(red,"_clusters_0.05"))
  
  integrated_st_obj_clusters_filt <- RunUMAP(integrated_st_obj_clusters_filt, reduction = red, 
                                             dims = 1:30, reduction.name = paste0("umap.", red))
}

# cluster_methods <- colnames(integrated_st_obj_clusters_filt@meta.data)[startsWith(colnames(integrated_st_obj_clusters_filt@meta.data), "integrated")]
# 
# for (m in cluster_methods) {
#   print(m)
#   print(table(integrated_st_obj_clusters_filt@meta.data[, m]))
# }

DefaultAssay(integrated_st_obj_clusters_filt) <- "RNA"
integrated_st_obj_clusters_filt <- JoinLayers(integrated_st_obj_clusters_filt)

DefaultAssay(integrated_st_obj_clusters_filt) <- "SCT"
integrated_st_obj_clusters_filt <- PrepSCTFindMarkers(integrated_st_obj_clusters_filt,
                                                      assay = "SCT", verbose = TRUE)

save(integrated_st_obj_clusters_filt, file = file.path(st_dir, "integrated_st_obj_clusters_filt_zeros_all_80.RData"))


# 6. Identify differential expressed genes across conditions --------------

cluster_methods <- sort(colnames(integrated_st_obj_clusters_filt@meta.data[startsWith(colnames(integrated_st_obj_clusters_filt@meta.data), "integrated")]))
seurat_integrated <- integrated_st_obj_clusters_filt
rm(integrated_st_obj_clusters_filt)

markers_identified_filt <- list()

for(cm in cluster_methods){
  m <- strsplit(cm, "_")[[1]][1]
  norm <- strsplit(m, ".", fixed = TRUE)[[1]][3]
  
  if(norm == "lognorm"){
    DefaultAssay(seurat_integrated) <- "RNA"
  }else if(norm == "sct"){
    DefaultAssay(seurat_integrated) <- "SCT"
  }else{
    break
  }
  
  Idents(seurat_integrated) <- cm
  markers <- FindAllMarkers(seurat_integrated, verbose = T)
  
  markers_identified_filt[[cm]] <- markers
}

save(markers_identified_filt, file = file.path(st_dir, "DE_markers_identified_filt_zeros_all_80.RData"))

top10_DE_markers_filt <- list()
for (cm in cluster_methods) {
  markers_cm <- markers_identified_filt[[cm]] 
  
  if(nrow(markers_cm)==0){
    top10_DE_markers_filt[[cm]] <- data.frame()
  }else{
    top10 <- markers_cm %>%
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < 0.001 & avg_log2FC > 0) %>%
      slice_max(avg_log2FC, n = 10) %>%
      ungroup()
    
    top10$gene_symbol <- plyr::mapvalues(top10$gene, 
                                         from = id_symbol_mapping$id,
                                         to = id_symbol_mapping$symbol)
    
    top10_DE_markers_filt[[cm]] <- top10
  }
}
save(top10_DE_markers_filt, file = file.path(st_dir, "top10_DE_markers_identified_filt_zeros_all_80.RData"))
