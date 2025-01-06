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
st_meta$sid <- st_meta$slide

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

# Run deconvolution --------------------------------------
run_deconvolution(st_meta = st_meta,
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

widths <- list(He = 9, `10x` = 14, Wu = 12)
pie_radii <- list(He = 0.5, `10x` = 0.7, Wu = 0.7)
point_sizes <- list(He = 3.5, `10x` = 1, Wu = 1)
y_reverse_bool <- list(He = T, `10x` = F, Wu = F)
x_reverse_bool <- list(He = F, `10x` = T, Wu = T)
xy_flip_bool <- list(He = F, `10x` = T, Wu = T)


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




