# Load libraries

library(ggplot2)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(gtable)
library(grid)
library(patchwork)
library(reshape2)
library(cowplot)
theme_set(theme_cowplot())

rm(list = ls())
gc()


source("visualization_helper.R")


# ========== FILTERED (zeros all samples, 90%) ========== ----------------------

### directories ----
data_dir = "../../data/He_2020/"
cluster_dir = "../../output/He_2020/clustering/"

annotations <- read.csv(file.path(data_dir, "He_clusters_filt_zeros_all_90.csv")) %>%
  remove_rownames() %>%
  column_to_rownames("X") %>%
  mutate(across(where(is.integer), as.factor))

tab1 <- as.data.frame.matrix(table(annotations$pid, annotations$label))
tab2 <- as.data.frame.matrix(table(annotations$pid, annotations$integrated.cca.lognorm_clusters_0.1))

tab <- bind_cols(tab1, tab2)

library(kableExtra)
kbl(tab, col.names = c("Non", "Tumor", "Cluster 0", "Cluster 1", "Cluster 2"),
    format = "latex") %>%
  kable_classic() %>%
  add_header_above(c(" " = 1, "He et al." = 2, "Seurat v5" = 3))
  

load(file.path(data_dir, "integrated_st_obj_clusters_filt_zeros_all_90.RData"))
seurat_integrated <- integrated_st_obj_clusters_filt
rm(integrated_st_obj_clusters_filt)

spots <- rownames(annotations)
sids <- sort(unique(annotations$sid))
pids <- sort(unique(annotations$pid))
subtypes <- sort(unique(annotations$subtype))

# Seurat clustering results
cluster_methods <- sort(colnames(annotations)[startsWith(colnames(annotations), "integrated")])
clusters <- annotations[, cluster_methods]

for (cm in cluster_methods) {
  print(table(annotations[, c("label", cm)]))
}

# x and y axes scales
annotations$row <- sapply(strsplit(annotations$loc, "x"),"[[",1)
annotations$col <- sapply(strsplit(annotations$loc, "x"),"[[",2)
range_x <- range(as.integer(annotations$row) - min(as.integer(annotations$row)))
range_y <- range(max(as.integer(annotations$col)) - as.integer(annotations$col)) 

# colors for tumor and non labels
colors_ann <- c("#E41A1C", "#377EB8")
colors_ann <- setNames(colors_ann, c("tumor", "non"))

# 1. Plot spots: Tumor vs Non-tumor (Author's annotation) ------------------

annotations$label <- factor(annotations$label, levels = c("tumor", "non"))
annotations_tumor_non <- unique(annotations[, c("sid", "pid", "loc", "label", "row", "col")])
visualize_clusters_comparison_spotplot(annotation_df_list = list(annotations_tumor_non),
                                       colnameList = list("label"),
                                       titles = c(""),
                                       pdf_width = 12,
                                       pdf_height = 5,
                                       fn_names = "tumor_non",
                                       fn_suffix = "label",
                                       color_list = list(colors_ann))

# seurat_integrated@meta.data <- cbind(seurat_integrated@meta.data, 
#                                      annotations[rownames(seurat_integrated@meta.data), 1:7])

view(seurat_integrated@meta.data)
seurat_integrated$sid_label<-paste(seurat_integrated$sid, seurat_integrated$label, sep = "\n")
seurat_integrated$pid_label<-paste(seurat_integrated$pid, seurat_integrated$label, sep = "\n")
seurat_integrated$subtype_label<-paste(seurat_integrated$subtype, seurat_integrated$label, sep = "\n")

# 2. Plot UMAP: Integrated over samples ------------------------

visualize_umap(cluster_methods = "integrated.cca.lognorm_clusters_0.1",
               seurat_obj = seurat_integrated,
               fn_suffix = "integrated.cca.lognorm_clusters_0.1",
               split.by.2 = "label",
               n_col.1 = 1,
               n_col.2 = 2,
               strip.text.size = 14,
               axis.text.size = 12,
               legend.text.size = 12,
               axis.title.size = 14,
               legend.point.size = 4,
               tag.size = 14,
               plot_titles = "All",
               pdf_height = 6,
               pdf_width = 20, 
               includes_number = T,
               includes_tag = T)


visualize_umap(cluster_methods = "integrated.cca.lognorm_clusters_0.1",
               seurat_obj = seurat_integrated,
               fn_suffix = "stratified_by_subtype_integrated.cca.lognorm_clusters_0.1",
               split.by.1 = "subtype",
               split.by.2 = "subtype_label",
               n_col.1 = 5,
               n_col.2 = 10,
               pdf_width = 13, 
               pdf_height = 6,
               strip.text.size = 10,
               axis.text.size = 10,
               legend.text.size = 11,
               axis.title.size = 14,
               plot_titles = "")

visualize_umap(cluster_methods = "integrated.cca.lognorm_clusters_0.1",
               seurat_obj = seurat_integrated,
               fn_suffix = "stratified_by_patient_integrated.cca.lognorm_clusters_0.1",
               split.by.1 = "pid",
               n_col.1 = 6,
               pdf_width = 10, 
               pdf_height = 6,
               strip.text.size = 10,
               axis.text.size = 10,
               legend.text.size = 11,
               axis.title.size = 14,
               plot_titles = "")

# 3. Plot spots: Integrated over samples ---------------------------------------

visualize_clusters_comparison_spotplot(annotation_df_list = list(annotations,
                                                                 annotations),
                                       colnameList = list("label",
                                                          "integrated.cca.lognorm_clusters_0.1"),
                                       fn_names = "integrated.cca.lognorm_clusters_0.1",
                                       titles = c("He et al. Cluster", "Seurat Cluster"),
                                       pdf_width = 12,
                                       pdf_height = 8,
                                       fn_suffix = "filt_all_zeros_90",
                                       color_list = list(colors_ann, NULL))


# 4. Plot top 10 DE markers -----------------------------------------------

load(file.path(data_dir, "top10_DE_markers_identified_filt_zeros_all.RData"))

visualize_markers_heatmap(cluster_methods = cluster_methods, 
                          seurat_obj = seurat_integrated, 
                          markers_list = top10_DE_markers_filt, 
                          fn_suffix = "filt_zeros_all",
                          pdf_width = 45,
                          pdf_height = 30)


visualize_markers_vlnplot(cluster_methods = cluster_methods, 
                          seurat_obj = seurat_integrated, 
                          markers_list = top10_DE_markers_filt, 
                          annotations_df = annotations_cluster,
                          fn_suffix = "filt_zeros_all")

visualize_markers_umap(cluster_methods = cluster_methods, 
                       seurat_obj = seurat_integrated, 
                       markers_list = top10_DE_markers_filt, 
                       annotations_df = annotations_cluster,
                       fn_suffix = "filt_zeros_all")

# 5. Compare with unfiltered -----------------------------------------------

merged_ann_clusters <- merge(annotations_cluster_unfilt, 
                             annotations_cluster, 
                             by = c(31, 1:6), 
                             all = T,
                             suffixes = c("_unfilt", "_filt")) 

for (col in colnames(merged_ann_clusters)[startsWith(colnames(merged_ann_clusters), "integrated")]) {
  merged_ann_clusters[, col] <- factor(merged_ann_clusters[, col],
                                       levels = sort(unique(merged_ann_clusters[, col]), 
                                                     na.last = TRUE),
                                       exclude=NULL)
}

colnameList.1 <- paste0(cluster_methods, "_unfilt")
colnameList.2 <- paste0(cluster_methods, "_filt")
visualize_clusters_comparison_heatmap(annotation_df = merged_ann_clusters, 
                              colnameList.1 = colnameList.1,
                              colnameList.2 = colnameList.2,
                              name.1 = "Unfiltered", 
                              name.2 = "Filtered (Zeros for All)",
                              fn_suffix = "unfilt_vs_filt_zeros_all")


common <- paste0(rep(c("lognorm", "sct"), each = 4), 
                 "_clusters_", 
                 rep(c("0.05", "0.1", "0.2", "0.3"),2), 
                 "_filt")
colnameList.1 <- paste0("integrated.cca.", common)
colnameList.2 <- paste0("integrated.rpca.", common)
visualize_clusters_comparison_heatmap(annotation_df = merged_ann_clusters, 
                              colnameList.1 = colnameList.1,
                              colnameList.2 = colnameList.2,
                              name.1 = "CCA", 
                              name.2 = "RPCA",
                              fn_suffix = "rpca_vs_cca_filt_zeros_all")


common <- paste0(rep(c("lognorm", "sct"), each = 4), 
                 "_clusters_", 
                 rep(c("0.05", "0.1", "0.2", "0.3"),2), 
                 "_filt")
colnameList.1 <- paste0("integrated.harmony.", common)
colnameList.2 <- paste0("integrated.rpca.", common)
visualize_clusters_comparison_heatmap(annotation_df = merged_ann_clusters, 
                              colnameList.1 = colnameList.1,
                              colnameList.2 = colnameList.2,
                              name.1 = "Harmony", 
                              name.2 = "RPCA",
                              fn_suffix = "rpca_vs_harmony_filt_zeros_all")

common_intg <- paste0("integrated.", c("cca", "rpca", "harmony"))
common_res <- paste0("_clusters_", 
                     c("0.05", "0.1", "0.2", "0.3"), 
                     "_filt")
colnameList.1 <- paste0(rep(common_intg, each = length(common_res)) , ".lognorm", common_res)
colnameList.2 <- paste0(rep(common_intg, each = length(common_res)) , ".sct", common_res)
visualize_clusters_comparison_heatmap(annotation_df = merged_ann_clusters, 
                              colnameList.1 = colnameList.1,
                              colnameList.2 = colnameList.2,
                              name.1 = "Log-Normalization", 
                              name.2 = "SCTransform",
                              fn_suffix = "lognorm_vs_sct_filt_zeros_all")


visualize_clusters_comparison_spotplot(annotation_df_list = list(annotations_cluster_unfilt, 
                                                                 annotations_cluster),
                                       colnameList = list(cluster_methods,
                                                          cluster_methods),
                                       titles = c("Unfiltered", "Filtered (Zeros by All)"),
                                       pdf_width=15,
                                       pdf_height=10,
                                       fn_names = cluster_methods,
                                       fn_suffix = "unfilt_vs_filt_zeros_all")


  











p_list <- visualize_umap(cluster_methods = "integrated.cca.lognorm_clusters_0.1",
                         seurat_obj = seurat_integrated,
                         fn_suffix = "integrated.cca.lognorm_clusters_0.1",
                         split.by.1 = "label",
                         n_col.1 = 1,
                         strip.text.size = 14,
                         axis.text.size = 12,
                         legend.text.size = 12,
                         axis.title.size = 14,
                         legend.point.size = 4,
                         tag.size = 14,
                         pdf_height = 6,
                         pdf_width = 7, 
                         plot_titles = "",
                         includes_number = T,
                         returns_plot = 1)

p_umap <- p_list[[1]]

scatterplot_AUC <- scatterplot_AUC + 
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16)) 
p_AUC_He

fig6 <- free(scatterplot_AUC) / ((p_umap | free(p_AUC_He)) & 
                                   theme(axis.text = element_text(size = 14),
                                         axis.title = element_text(size = 16),
                                         legend.text = element_text(size = 14),
                                         legend.title = element_text(size = 16)))+ 
  plot_annotation(tag_levels = 'A',
                  tag_prefix = "(",
                  tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 16, face = "bold")) 

ggsave(paste0("Fig6.png"), fig6, width = 11.5, height = 10.5)





