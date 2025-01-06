# -----------------------------------------------------------------------------#
# Date last modified: June 16, 2024 by Zhining Sui
# Program: run_deconvolution.R
# Purpose: Perform cell type deconvolution for spatial transcriptomics data using reference scRNA-seq data.

# Data Inputs:
# - scRNA-seq reference count data: A sparse count matrix, metadata for cells, and gene names.
# - Spatial transcriptomics (ST) count data: Filtered sparse count matrix, metadata for spots, and spatial location information.

# Data Outputs:
# - CARD_obj_<sampleID>.rds: Deconvoluted object for each sample.
# - CARD_alg_matrix_<sampleID>.rds: Algorithm matrix for each sample.
# - Proportion_<sampleID>.csv: Cell type proportions for each spot in the sample.
# - Plots showing the cell type proportions and spatial distributions saved as JPG files in the specified output directory.
# -----------------------------------------------------------------------------#

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

source("func_deconvolution.R")

# Define paths for scRNA-seq and spatial transcriptomics (ST) data --------

sc_dir <- "/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/Wu_2021/data/scRNASeq/"
st_dirs <- list(
  He = "/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/He_2020/data/",
  `10x` = "/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/10X/data/",
  Wu = "/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/Wu_2021/data/",
  Janesick = "/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/Janesick_2023/data/"
)
output_dirs <- list(
  He = "../../output/He_2020/deconvolution/",
  `10x` = "../../output/10X/deconvolution/",
  Wu = "../../output/Wu_2021/deconvolution/",
  Janesick = "../../output/Janesick_2023/deconvolution/"
)


# Prepare st_meta and spatial_location for each ST dataset ----------------

st_meta_all <- list()
spatial_all <- list()

###### Wu_2021 
# Create metadata for the spatial transcriptomics data
st_dir <- st_dirs$Wu
metadata_files = list.files(file.path(st_dir, "metadata"), pattern="_metadata.csv")
metadata <- data.frame()
for(f1 in metadata_files){
  d1 <- fread(file.path(st_dir, "metadata", f1))
  metadata <- rbind(metadata, d1)
}
metadata$sid <- metadata$patientid
st_meta <- unique(metadata[, c("sid", "subtype")])
st_meta$subtype <- ifelse(st_meta$subtype %in% c("ER", "HER2"), paste0(st_meta$subtype, "+"), st_meta$subtype)
st_meta$count_matrix_dir <- paste0(st_dir, "filtered_count_matrices/", st_meta$sid, "_filtered_count_matrix/matrix.mtx.gz")
st_meta$feature_dir <- paste0(st_dir, "filtered_count_matrices/", st_meta$sid, "_filtered_count_matrix/features.tsv.gz")
st_meta$barcode_dir <- paste0(st_dir, "filtered_count_matrices/", st_meta$sid, "_filtered_count_matrix/barcodes.tsv.gz")
st_meta <- as.data.frame(st_meta)

# Create spatial location data
spatial_dir = paste0(st_dir, "spatial/", st_meta$sid, "_spatial/")
spatial <- data.frame()
for (i in 1:nrow(st_meta)) {
  sp <- read.csv(file.path(st_dir, "spatial/", paste0(st_meta[i, "sid"], "_spatial/tissue_positions_list.csv")), 
                 header = FALSE, sep = ',')
  sp <- sp[c(1,3,4)]
  colnames(sp)[2:3] <- c("x", "y")
  sp$sid <- st_meta[i, "sid"]
  sp$spot_id <- sp$V1
  # sp$spot_id <- paste0(sp$sid, "_", sp$V1)
  spatial <- rbind(spatial, sp)
}
st_meta_all[["Wu"]] <- st_meta
spatial_all[["Wu"]] <- spatial

###### 10x Visium 
# Create metadata for the spatial transcriptomics data
st_dir <- st_dirs$`10x`
sample_ids = c("Visium_FFPE_Human_Breast_Cancer", "Visium_Human_Breast")
st_meta <- data.frame(sid = sample_ids,
                      subtype = NA,
                      count_matrix_dir = file.path(st_dir, sample_ids, "filtered_feature_bc_matrix/matrix.mtx.gz"),
                      feature_dir = file.path(st_dir, sample_ids, "filtered_feature_bc_matrix/features.tsv.gz"),
                      barcode_dir = file.path(st_dir, sample_ids, "filtered_feature_bc_matrix/barcodes.tsv.gz"))

# Create spatial location data
spatial <- data.frame()
for (i in 1:nrow(st_meta)) {
  sp <- read.csv(file.path(st_dir, st_meta[i, "sid"], "spatial", "tissue_positions_list.csv"), 
                 header = FALSE, sep = ',')
  sp <- sp[c(1,3,4)]
  colnames(sp)[2:3] <- c("x", "y")
  sp$sid <- st_meta[i, "sid"]
  sp$spot_id <- sp$V1
  spatial <- rbind(spatial, sp)
}
st_meta_all[["10x"]] <- st_meta
spatial_all[["10x"]] <- spatial

###### He_2020 
# Create metadata for the spatial transcriptomics data
st_dir <- st_dirs$He
st_meta <- read.csv(file.path(st_dir, "metadata.csv"))
st_meta$sid <- paste0(st_meta$patient, "_", st_meta$replicate)
st_meta$subtype <- st_meta$type
st_meta$subtype <- ifelse(st_meta$subtype %in% c('Luminal_A', 'Luminal_B'), "ER+", 
                          ifelse(st_meta$subtype %in% c('HER2_luminal', 'HER2_non_luminal'), "HER2+", 
                                 st_meta$subtype))
st_meta$count_matrix_dir <- file.path(st_dir, st_meta$count_matrix)

# Create spatial location data
spatial <- data.frame()
for (i in 1:nrow(st_meta)) {
  sp <- read.csv(file.path(st_dir, st_meta[i, "spot_coordinates"]))
  sp$sid <- st_meta[i, "sid"]
  # sp$spot_id <- paste0(sp$sid, "_", sp$X.1)
  sp$spot_id <- sp$X.1
  sp$x <- as.integer(sapply(strsplit(sp$X.1, "x"), "[[", 1))
  sp$y <- as.integer(sapply(strsplit(sp$X.1, "x"), "[[", 2))
  spatial <- rbind(spatial, sp)
}
st_meta_all[["He"]] <- st_meta
spatial_all[["He"]] <- spatial

###### Janesick_2023
# Create metadata for the spatial transcriptomics data
st_dir <- st_dirs$Janesick
sample_ids = c("VisiumS1")
st_meta <- data.frame(sid = sample_ids,
                      subtype = NA,
                      count_matrix_dir = file.path(st_dir, sample_ids, "filtered_feature_bc_matrix/matrix.mtx.gz"),
                      feature_dir = file.path(st_dir, sample_ids, "filtered_feature_bc_matrix/features.tsv.gz"),
                      barcode_dir = file.path(st_dir, sample_ids, "filtered_feature_bc_matrix/barcodes.tsv.gz"))

# Create spatial location data
spatial <- data.frame()
for (i in 1:nrow(st_meta)) {
  sp <- read.csv(file.path(st_dir, st_meta[i, "sid"], "spatial", "tissue_positions.csv"), 
                 header = T, sep = ',')
  sp <- sp[c(1,3,4)]
  colnames(sp)[2:3] <- c("x", "y")
  sp$sid <- st_meta[i, "sid"]
  sp$spot_id <- sp$barcode
  spatial <- rbind(spatial, sp)
}
st_meta_all[["Janesick"]] <- st_meta
spatial_all[["Janesick"]] <- spatial


# Load reference scRNA-seq data -------------------------------------------

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

# sc_dir <- "/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/Janesick_2023/data/FRP/sample_filtered_feature_bc_matrix"
# sc_count <- readMM(file.path(sc_dir, "matrix.mtx.gz"))
# sc_meta <- read.csv(file = file.path(sc_dir, "metadata.csv"))
# rownames(sc_meta) <- sc_meta$Barcode
# sc_gene <- read.csv(file = file.path(sc_dir, "features.tsv.gz"), sep = "\t", header = FALSE)
# sc_barcode <- read.csv(file = file.path(sc_dir, "barcodes.tsv.gz"), sep = "\t", header = FALSE)
# colnames(sc_count) <- sc_barcode$V1
# rownames(sc_count) <- sc_gene$V1
# sc_count <- as(sc_count, "sparseMatrix")
# sc_count <- as(sc_count, "CsparseMatrix")
# sc_count <- sc_count[, intersect(sc_meta$Barcode, sc_barcode$V1)]
# sc_meta <- sc_meta[colnames(sc_count), ]
# sc_meta$sid <- "S1"

# Run deconvolution for each dataset --------------------------------------
for (dataset in names(st_dirs)) {
  run_deconvolution(st_meta = st_meta_all[[dataset]],
                    spatial_location = spatial_all[[dataset]],
                    output_dir = output_dirs[[dataset]])
}

run_deconvolution(st_meta = st_meta_all$Wu,
                  spatial_location = spatial_all$Wu,
                  output_dir = output_dirs$Wu)

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

widths <- list(He = 9, `10x` = 14, Wu = 12)
pie_radii <- list(He = 0.5, `10x` = 0.7, Wu = 0.7)
point_sizes <- list(He = 3.5, `10x` = 1, Wu = 1)
y_reverse_bool <- list(He = T, `10x` = F, Wu = F)
x_reverse_bool <- list(He = F, `10x` = T, Wu = T)
xy_flip_bool <- list(He = F, `10x` = T, Wu = T)

for (dataset in names(st_dirs)) {
  visualize_proportions(st_meta = st_meta_all[[dataset]], 
                        output_dir = output_dirs[[dataset]],
                        ct.visualize = c("B-cells",
                                         "CAFs",
                                         "Cancer Epithelial",
                                         "Endothelial", 
                                         "Myeloid" ,
                                         "Normal Epithelial",
                                         "Plasmablasts",
                                         "PVL",
                                         "T-cells"),
                        width = widths[[dataset]],
                        point_size = point_sizes[[dataset]],
                        pie_radius = pie_radii[[dataset]],
                        y_reverse = y_reverse_bool[[dataset]],
                        x_reverse = x_reverse_bool[[dataset]],
                        xy_flip = xy_flip_bool[[dataset]])
}



# Create input csv file for STpath pipeline -------------------------------

STpath_input_dirs <- list(
  He = "../../data/He_2020/",
  `10x` = "../../data/10X/",
  Wu = "../../data/Wu_2021/"
)

for (dataset in names(output_dirs)) {
  prop_files = list.files(output_dirs[[dataset]], pattern=".csv")
  fn <- gsub("Proportion_|_celltype_major.csv", "", prop_files)
  
  lst <- lapply(prop_files, function(x){ read.csv(file.path(output_dirs[[dataset]], x), header=TRUE, stringsAsFactors=FALSE) })
  lst <- lapply(lst, function(x) {
    x$X1 <- x$X
    x$X <- paste0(x$sid, "_", x$X1)
    x$invasive.cancer = x$Cancer.Epithelial
    x$stroma = x$Endothelial + x$PVL + x$CAFs
    x$lymphocyte = x$T.cells + x$B.cells + x$Plasmablasts
    x$others = x$Myeloid + x$Normal.Epithelial
    return(x)
  })
  
  lapply(seq_along(lst), function(i) {
    write.csv(lst[[i]], file.path(STpath_input_dirs[[dataset]], 
                                    paste0("Proportion_", fn[i], ".csv")))
  })
  
  prop <- do.call("bind_rows", lst)
  write.csv(prop, file.path(STpath_input_dirs[[dataset]], paste0("Proportion_", dataset, ".csv")))
}




# Predict expression from proportions -------------------------------------
# Wu et data
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

st_meta <- st_meta_all$Wu
output_dir <- output_dirs$Wu
predicted_expr_list <- list()
observed_st_expr_list <- list()
sample_wise_cor <- list()
gene_wise_cor <- list()
for (i in 1:nrow(st_meta)) {
  # Retrieve the reference average expression matrix for the specific cell subtype
  avg_type_ref <- avg_sc_ref[[st_meta[i, "subtype"]]]
  # Get the sample ID for the current iteration
  s1 <- st_meta[i, "sid"]
  # Read the cell type proportion data for the current sample
  prop <- read.csv(file.path(output_dir, sprintf("/Proportion_%s_celltype_major.csv", s1)),
                   check.names = F)
  # Set the first column as rownames and clean up the proportion data
  colnames(prop)[1] <- "X"
  prop <- prop %>%
    column_to_rownames("X") %>%
    dplyr::select(-sid) %>%
    t()
  # Calculate predicted gene expression using matrix multiplication
  predicted_expr <- avg_type_ref %*% prop
  # Store the predicted expression in a list
  predicted_expr_list[[s1]] <- predicted_expr
  
  # Read the spatial transcriptomics count matrix for the sample
  stdata_s <- readMM(st_meta[i, "count_matrix_dir"])
  # Handle missing column or row names for the count matrix
  if(is.null(colnames(stdata_s)) || is.null(rownames(stdata_s))){
    barcode <- read.table(file = gzfile(st_meta[i, "barcode_dir"]),
                          sep = '\t', header = FALSE)
    feature <- read.table(file = gzfile(st_meta[i, "feature_dir"]),
                          sep = '\t', header = FALSE)
    barcode_names <- barcode$V1
    # feature_names <- feature$V1
    feature_names <- feature[, sapply(feature, function(column) !("Gene Expression" %in% column) && !any(grepl("^ENSG", column)))]
    
    # Assign to the dimension with same length.
    if (ncol(stdata_s) == length(barcode_names) && nrow(stdata_s) == length(feature_names)) {
      colnames(stdata_s) <- barcode_names
      rownames(stdata_s) <- feature_names
    } else if (nrow(stdata_s) == length(barcode_names) && ncol(stdata_s) == length(feature_names)) {
      rownames(stdata_s) <- barcode_names
      colnames(stdata_s) <- feature_names
    } else {
      cat("The lengths of the barcode or feature names do not match the dimensions of the matrix.\n")
    }
  }
  # Convert the matrix to a sparse format for efficient computation
  stdata_s <- as(stdata_s, "sparseMatrix")
  stdata_s <- as(stdata_s, "CsparseMatrix")
  # Extract the count matrix for columns (samples) matching the predicted expression matrix
  st_count <- stdata_s[, colnames(predicted_expr)]
  observed_st_expr_list[[s1]] <- st_count
  # # Find common genes between the predicted and observed matrices
  # common_genes <- intersect(rownames(st_count), rownames(predicted_expr))
  # st_count <- st_count[common_genes,]
  # predicted_expr <- predicted_expr[common_genes, ]
  # 
  # # Compute sample-wise correlations directly
  # column_wise_correlation <- sapply(seq_len(ncol(st_count)), function(i) {
  #   cor(st_count[, i], predicted_expr[, i], method = "pearson")
  # })
  # names(column_wise_correlation) <- colnames(st_count)
  # sample_wise_cor[[s1]] <- column_wise_correlation
  # 
  # # Compute gene-wise correlations directly
  # row_wise_correlation <- sapply(seq_len(nrow(st_count)), function(i) {
  #   cor(st_count[i, ], predicted_expr[i, ], method = "pearson")
  # })
  # names(row_wise_correlation) <- rownames(st_count)
  # gene_wise_cor[[s1]] <- row_wise_correlation
}

save(predicted_expr_list, observed_st_expr_list, gene_wise_cor, sample_wise_cor, avg_sc_ref, file = "card_predict_expression.rda")

cat("Summary of sample-wise correlations:\n")
sapply(sample_wise_cor, summary)
par(mfrow = c(2,3))
hist(sample_wise_cor$`1142243F`, main = "Sample-wise Correlation for 1142243F", 
     xlab = "Correlation",  breaks = 30, xlim = c(0,1))
hist(sample_wise_cor$`1160920F`, main = "Sample-wise Correlation for 1160920F", 
     xlab = "Correlation", breaks = 30, xlim = c(0,1))
hist(sample_wise_cor$CID4290, main = "Sample-wise Correlation for CID4290", 
     xlab = "Correlation", breaks = 30, xlim = c(0,1))
hist(sample_wise_cor$CID4465, main = "Sample-wise Correlation for CID4465", 
     xlab = "Correlation", breaks = 30, xlim = c(0,1))
hist(sample_wise_cor$CID44971, main = "Sample-wise Correlation for CID44971", 
     xlab = "Correlation", breaks = 30, xlim = c(0,1))
hist(sample_wise_cor$CID4535, main = "Sample-wise Correlation for CID4535", 
     xlab = "Correlation", breaks = 30, xlim = c(0,1))


spot <- names(sort(sample_wise_cor$`1142243F`, decreasing = T)[1])

spot_df <- data.frame(cbind( predicted_expr_list$`1142243F`[, spot],
                  observed_st_expr_list$`1142243F`[, spot]))
colnames(spot_df) <- c("Predicted", "Observed")
spot_df$Predicted <- round(spot_df$Predicted)
ggplot(spot_df, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.5)




