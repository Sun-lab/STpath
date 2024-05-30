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
library(sva)
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






load(file.path(st_dir, "merged_st_obj_normalized_unfilt.RData"))
merged_st_obj_unfilt <- JoinLayers(merged_st_obj_unfilt)

# setting up the expression data from the normalized seurat 

edata <- GetAssayData(merged_st_obj_unfilt)
info <- merged_st_obj_unfilt@meta.data
batch = info$sid
mod_combat = model.matrix(~1, data=info)
combat.edata = ComBat(dat = edata, 
                      batch = batch, 
                      mod = mod_combat, 
                      par.prior = TRUE, 
                      prior.plots = FALSE)




