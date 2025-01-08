rm(list = ls())
setwd("~/GitHub/STpath/code/STpath pipeline")
# Load necessary libraries
library(ggplot2)
library(data.table)
library(MuSiC)
library(CARD)
library(stringr)
library(Matrix)
library(readr)
library(tidyverse)
library(biomaRt)
library(gprofiler2)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v79)
library(dplyr)
library(grid)

source("deconvolution.R")

# Set directories 
sc_dir <- "/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/Wu_2021/data/scRNASeq/"
st_dir <- "/Users/zhiningsui/GitHub/STimage-1K4M"
output_dir <- "/Users/zhiningsui/GitHub/STpath/output/STimage-1K4M/deconvolution/"

# Load metadata for the ST data
st_meta <- read.csv("/Users/zhiningsui/GitHub/STpath/data/STimage-1K4M/create_patches_input_brca.csv")
# Create spatial location data
spatial <- data.frame()
for (i in 1:nrow(st_meta)) {
  sp <- read.csv(st_meta[i, "spatial_path"], sep = ',')
  sp$spot_id <- sp$X
  sp$sid <- st_meta[i, "sid"]
  # sp$loc <- gsub(paste0(unique(sp$sid), "_"), "", sp$spot_id)
  # sp$x <- sapply(strsplit(sp$loc, "x"), "[[", 1)
  # sp$y <- sapply(strsplit(sp$loc, "x"), "[[", 2)
  sp$x <- sp$xaxis
  sp$y <- sp$yaxis
  spatial <- rbind(spatial, sp)
}

# Load reference scRNA-seq data 
sc_count <- readMM(file.path(sc_dir, "count_matrix_sparse.mtx"))
sc_meta <- read.csv(file = file.path(sc_dir, "metadata.csv"))
rownames(sc_meta) <- sc_meta$X
sc_gene <- read.csv(file = file.path(sc_dir, "count_matrix_genes.tsv"), sep = "\t", header = FALSE)
colnames(sc_gene) <- "symbol.sc"
sc_barcode <- read.csv(file = file.path(sc_dir, "count_matrix_barcodes.tsv"), sep = "\t", header = FALSE)
colnames(sc_count) <- sc_barcode$V1
rownames(sc_count) <- sc_gene$symbol.sc
sc_count <- as(sc_count, "sparseMatrix")
sc_count <- as(sc_count, "CsparseMatrix")
sc_count.log2trans <- log2(sc_count+1)
sc_count.log2trans <- as(sc_count.log2trans, "sparseMatrix")
sc_count.log2trans <- as(sc_count.log2trans, "CsparseMatrix")

# Run deconvolution --------------------------------------
run_deconvolution(sc_count = sc_count, 
                  sc_meta = sc_meta,
                  st_meta = st_meta,
                  spatial_location = spatial,
                  output_dir = output_dir,
                  save_CARD_objs = F)

# Visualize proportions for all datasets ----------------------------------

colors = c(`B cells` = "#2a5fbd",
           CAFs = "#7ddffa",
           `Cancer Epithelial` = "#117d30",
           Endothelial = "#71f25a", 
           Myeloid = "#ebc857",
           `Normal Epithelial` ="#D39200",
           Plasmablasts = "#F8766D",
           PVL = "#DB72FB",
           `T cells` = "#bd2a84")
visualize_proportions(st_meta = st_meta, 
                      output_dir = output_dir,
                      ct.visualize = c("B-cells",
                                       "CAFs",
                                       "Cancer Epithelial",
                                       "Endothelial", 
                                       "Myeloid" ,
                                       "Normal Epithelial",
                                       "Plasmablasts",
                                       "PVL",
                                       "T-cells"),
                      width = 10,
                      point_size = 1,
                      pie_radius = 0.5,
                      y_reverse = T,
                      x_reverse = F,
                      xy_flip = T)

# Create input csv file for STpath pipeline -------------------------------
STpath_input_dir <- "/Users/zhiningsui/GitHub/STpath/data/STimage-1K4M/"
prop_files = list.files(output_dir, pattern=".csv")
fn <- gsub("Proportion_|_celltype_major.csv", "", prop_files)
lst <- lapply(prop_files, function(x){ read.csv(file.path(output_dir, x), header=TRUE, stringsAsFactors=FALSE) })
lst <- lapply(lst, function(x) {
  x$invasive.cancer = x$Cancer.Epithelial
  x$stroma = x$Endothelial + x$PVL + x$CAFs
  x$lymphocyte = x$T.cells + x$B.cells + x$Plasmablasts
  x$others = x$Myeloid + x$Normal.Epithelial
  return(x)
})
lapply(seq_along(lst), function(i) {
  write.csv(lst[[i]], file.path(STpath_input_dir, 
                                paste0("Proportion_", fn[i], ".csv")))
})
prop <- do.call("bind_rows", lst)
write.csv(prop, file.path(STpath_input_dir, "Proportion_STimage-1K4M.csv"))


# Predict expression from proportions -------------------------------------

avg_sc_ref <- list()
for (subtype in unique(sc_meta$subtype)) {
  cells2use <- which(sc_meta$subtype == subtype)
  sc_count_i <- sc_count[, cells2use]
  sc_meta_i <- sc_meta[cells2use, ]
  avg_expr_type <- matrix(nrow = nrow(sc_count_i),
                          ncol = length(unique(sc_meta_i$celltype_major)))
  colnames(avg_expr_type) <- unique(sc_meta_i$celltype_major)
  rownames(avg_expr_type) <- rownames(sc_count_i)
  for (ct in unique(sc_meta_i$celltype_major)) {
    cells <- sc_meta_i[sc_meta_i$celltype_major == ct,]$X
    expr <- sc_count_i[, cells]
    avg_expr_type[, ct] <- rowMeans(expr)
  }
  avg_sc_ref[[subtype]] <- avg_expr_type
}

avg_expr_type <- matrix(nrow = nrow(sc_count),
                        ncol = length(unique(sc_meta$celltype_major)))
colnames(avg_expr_type) <- unique(sc_meta$celltype_major)
rownames(avg_expr_type) <- rownames(sc_count)
for (ct in unique(sc_meta$celltype_major)) {
  cells <- sc_meta[sc_meta$celltype_major == ct,]$X
  expr <- sc_count[, cells]
  avg_expr_type[, ct] <- rowMeans(expr)
}
avg_sc_ref[["All"]] <- avg_expr_type

# log2 transformation
for (subtype in unique(sc_meta$subtype)) {
  cells2use <- which(sc_meta$subtype == subtype)
  sc_count_i <- sc_count.log2trans[, cells2use]
  sc_meta_i <- sc_meta[cells2use, ]
  avg_expr_type <- matrix(nrow = nrow(sc_count_i),
                          ncol = length(unique(sc_meta_i$celltype_major)))
  colnames(avg_expr_type) <- unique(sc_meta_i$celltype_major)
  rownames(avg_expr_type) <- rownames(sc_count_i)
  for (ct in unique(sc_meta_i$celltype_major)) {
    cells <- sc_meta_i[sc_meta_i$celltype_major == ct,]$X
    expr <- sc_count_i[, cells]
    avg_expr_type[, ct] <- rowMeans(expr)
  }
  avg_sc_ref[[paste0(subtype, ".log2trans")]] <- avg_expr_type
}

avg_expr_type <- matrix(nrow = nrow(sc_count.log2trans),
                        ncol = length(unique(sc_meta$celltype_major)))
colnames(avg_expr_type) <- unique(sc_meta$celltype_major)
rownames(avg_expr_type) <- rownames(sc_count.log2trans)
for (ct in unique(sc_meta$celltype_major)) {
  cells <- sc_meta[sc_meta$celltype_major == ct,]$X
  expr <- sc_count.log2trans[, cells]
  avg_expr_type[, ct] <- rowMeans(expr)
}
avg_sc_ref[["All.log2trans"]] <- avg_expr_type

save(avg_sc_ref, file = file.path(output_dir, "sc_ref_avg_ct_expr.rda"))

load(file.path(output_dir, "sc_ref_avg_ct_expr.rda"))

expr_pred_all <- list()
expr_pred_all.log2trans <- list()
expr_obs_all <- list()
expr_obs_all.log2trans <- list()

slides <- c("Human_Breast_Wu_06052021_Visium_1142243F",
            "Human_Breast_Wu_06052021_Visium_1160920F",
            "Human_Breast_Wu_06052021_Visium_CID4465",
            "Human_Breast_Wu_06052021_Visium_CID44971",
            "Human_Breast_Wu_06052021_Visium_CID4535")

for (i in 1:length(slides)) {
  s1 <- slides[i]
  type1 <- st_meta[st_meta$sid == s1, "subtype"]
  # s1 <- st_meta[i, "sid"]
  prop <- read.csv(file.path(output_dir, sprintf("/Proportion_%s_celltype_major.csv", s1)),
                   check.names = F)
  # Set the first column as rownames and clean up the proportion data
  colnames(prop)[1] <- "X"
  prop <- prop %>%
    column_to_rownames("X") %>%
    dplyr::select(-sid) %>%
    t()
  # Retrieve the reference average expression matrix
  sc_ct_avg <- avg_sc_ref[[type1]]
  sc_ct_avg <- sc_ct_avg[, rownames(prop)]
  # Calculate predicted gene expression
  predicted_expr <- sc_ct_avg %*% prop
  # Store the predicted expression in a list
  expr_pred_all[[s1]] <- predicted_expr
  
  # Retrieve the reference average expression matrix
  sc_ct_avg <- avg_sc_ref[[paste0(type1, ".log2trans")]]
  sc_ct_avg <- sc_ct_avg[, rownames(prop)]
  # Calculate predicted gene expression
  predicted_expr <- sc_ct_avg %*% prop
  # Store the predicted expression in a list
  expr_pred_all.log2trans[[s1]] <- predicted_expr
  
  gc()
  
  # Read the spatial transcriptomics count matrix for the sample
  stdata_s <- read.csv(st_meta[st_meta$sid == s1, "count_matrix_dir"]) %>%
    remove_rownames %>% column_to_rownames(var="X") %>% t()
  stdata_s <- as(stdata_s, "sparseMatrix")
  stdata_s <- as(stdata_s, "CsparseMatrix")
  stdata_s.log2trans <- log2(stdata_s+1)
  stdata_s.log2trans <- as(stdata_s.log2trans, "sparseMatrix")
  stdata_s.log2trans <- as(stdata_s.log2trans, "CsparseMatrix")
  # Extract the count matrix for columns (samples) matching the predicted expression matrix
  st_count <- stdata_s[, colnames(predicted_expr)]
  expr_obs_all[[s1]] <- st_count
  st_count <- stdata_s.log2trans[, colnames(predicted_expr)]
  expr_obs_all.log2trans[[s1]] <- st_count

  gc()
}

save(expr_pred_all, expr_obs_all, expr_pred_all.log2trans, expr_obs_all.log2trans, 
     file = file.path(output_dir, "expr_pred_obs.rda"))
load(file.path(output_dir, "expr_pred_obs.rda"))

sample_wise_corr <- list()
gene_wise_corr <- list()
for (i in 1:3) {
  s1 <- slides[i]
  expr_pred <- expr_pred_all[[s1]]
  expr_obs <- expr_obs_all[[s1]]
  # Find common genes between the predicted and observed matrices
  common_genes <- intersect(rownames(expr_obs), rownames(expr_pred))
  expr_obs <- expr_obs[common_genes,]
  expr_pred <- expr_pred[common_genes, ]
  expr_obs <- expr_obs[, colnames(expr_pred)]
  # Compute sample-wise correlations
  colwise_corr <- sapply(seq_len(ncol(expr_obs)), function(i) {
    cor(expr_obs[, i], expr_pred[, i], method = "pearson")
  })
  names(colwise_corr) <- colnames(expr_obs)
  sample_wise_corr[[s1]] <- colwise_corr
  # Compute gene-wise correlations
  rowwise_corr <- sapply(seq_len(nrow(expr_obs)), function(i) {
    cor(expr_obs[i, ], expr_pred[i, ], method = "pearson")
  })
  names(rowwise_corr) <- rownames(expr_obs)
  gene_wise_corr[[s1]] <- rowwise_corr
  
  gc()
}

sample_wise_corr.log2trans <- list()
gene_wise_corr.log2trans <- list()
for (i in 151:189) {
  s1 <- st_meta[i, "sid"]
  expr_pred <- expr_pred_all.log2trans[[s1]]
  expr_obs <- expr_obs_all.log2trans[[s1]]
  # Find common genes between the predicted and observed matrices
  common_genes <- intersect(rownames(expr_obs), rownames(expr_pred))
  expr_obs <- expr_obs[common_genes,]
  expr_pred <- expr_pred[common_genes, ]
  expr_obs <- expr_obs[, colnames(expr_pred)]
  # Compute sample-wise correlations
  colwise_corr <- sapply(seq_len(ncol(expr_obs)), function(i) {
    cor(expr_obs[, i], expr_pred[, i], method = "pearson")
  })
  names(colwise_corr) <- colnames(expr_obs)
  sample_wise_corr.log2trans[[s1]] <- colwise_corr
  # Compute gene-wise correlations
  rowwise_corr <- sapply(seq_len(nrow(expr_obs)), function(i) {
    cor(expr_obs[i, ], expr_pred[i, ], method = "pearson")
  })
  names(rowwise_corr) <- rownames(expr_obs)
  gene_wise_corr.log2trans[[s1]] <- rowwise_corr
}


sample_wise_corr1 <- sample_wise_corr
gene_wise_corr1 <- gene_wise_corr
sample_wise_corr.log2trans1 <- sample_wise_corr.log2trans
gene_wise_corr.log2trans1 <- gene_wise_corr.log2trans


load("corr_pred_obs.rda")
sample_wise_corr <- c(sample_wise_corr1, sample_wise_corr)
gene_wise_corr <- c(gene_wise_corr1, gene_wise_corr)
sample_wise_corr.log2trans <- c(sample_wise_corr.log2trans1, sample_wise_corr.log2trans)
gene_wise_corr.log2trans <- c(gene_wise_corr.log2trans1, gene_wise_corr.log2trans)

save(sample_wise_corr, gene_wise_corr, 
     sample_wise_corr.log2trans, gene_wise_corr.log2trans, 
     file = "corr_pred_obs.rda")

cat("Summary of sample-wise correlations:\n")
sapply(sample_wise_corr, summary)

pdf("sample_wise_corr.pdf", width = 15, height = 15)
par(mfrow = c(3,3))
for (i in 1:length(sample_wise_corr)) {
  hist(sample_wise_corr[[i]], main = paste0("Sample-wise Correlation for ", names(sample_wise_corr)[i]), 
       xlab = "Correlation",  breaks = 30, xlim = c(0,1))
}
dev.off()

i = 1
hist(sample_wise_corr[[i]], main = paste0("Spot-wise Correlation for ", names(sample_wise_corr)[i]), 
     xlab = "Correlation",  breaks = 30, xlim = c(0,1))
hist(gene_wise_corr[[i]], main = paste0("Gene-wise Correlation for ", names(sample_wise_corr)[i]), 
     xlab = "Correlation",  breaks = 30, xlim = c(0,1))

spot <- names(sort(sample_wise_cor$`1142243F`, decreasing = T)[1])

spot_df <- data.frame(cbind( expr_pred_all$`1142243F`[, spot],
                             expr_obs_all$`1142243F`[, spot]))
colnames(spot_df) <- c("Predicted", "Observed")
spot_df$Predicted <- round(spot_df$Predicted)
ggplot(spot_df, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.5)




