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


# Define functions --------------------------------------------------------

# -----------------------------------------------------------------------------#
# Convert gene IDs to gene symbols using multiple sources.
# Args:
#     gene_ids (vector): A vector of gene IDs to convert.
# Returns:
#     all_matches (data.frame): A data frame with the original gene IDs and the corresponding gene symbols from different sources.
# -----------------------------------------------------------------------------#
get_gene_symbols <- function(gene_ids) {
  # biomaRt
  ensembl <- useEnsembl(biomart = "genes")
  ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
  gene_symbol_1 <- getBM(mart = ensembl, attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                         filters = "ensembl_gene_id", values = gene_ids)
  colnames(gene_symbol_1) <- c("id", "symbol.biomaRt")
  gene_symbol_1$symbol.biomaRt <- ifelse(gene_symbol_1$symbol.biomaRt == "", NA, gene_symbol_1$symbol.biomaRt)
  
  # gprofiler2
  gene_symbol_2 <- gconvert(gene_ids,
                            organism="hsapiens",
                            target="ENTREZGENE",
                            filter_na = FALSE)[, c("input", "target")]
  colnames(gene_symbol_2) <- c("id", "symbol.gprofiler2")
  
  # org.Hs.eg.db
  symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, keytype = "ENSEMBL", column = "SYMBOL")
  gene_symbol_3 <- data.frame(id = names(symbols), symbol.org.Hs.eg.db = symbols, stringsAsFactors = FALSE)
  
  # EnsDb.Hsapiens.v79
  gene_symbol_4 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= gene_ids, keytype = "GENEID", columns = c("GENEID","SYMBOL"))
  colnames(gene_symbol_4) <- c("id", "symbol.EnsDb.Hsapiens.v79")
  
  # Initialize data frame to store results
  all_matches <- data.frame(id = gene_ids)
  
  # Combine results from all sources
  all_matches <- all_matches %>%
    left_join(gene_symbol_1, by = "id") %>%
    left_join(gene_symbol_2, by = "id") %>%
    left_join(gene_symbol_3, by = "id") %>%
    left_join(gene_symbol_4, by = "id") %>%
    pivot_longer(2:5, values_to = "symbol")
  
  all_matches <- unique(all_matches[, c(1,3)])
  
  return(all_matches)
}


# -----------------------------------------------------------------------------#
# Perform cell type deconvolution on spatial transcriptomics data using reference scRNA-seq data.
# Args:
#     st_meta (data.frame): Metadata for spatial transcriptomics samples. Contains sample_id, subtype, count_matrix_dir, feature_dir, and barcode_dir.
#     spatial_location (data.frame): Spatial location information for ST spots, with row names being barcodes, having one column of row position and one of column position.
#     output_dir (str): Directory to save the deconvolution results.
# Outputs:
#     Saves the following files to the output directory:
#     - CARD_obj_<sampleID>.rds: Deconvoluted object for each sample.
#     - CARD_alg_matrix_<sampleID>.rds: Algorithm matrix for each sample.
#     - Proportion_<sampleID>.csv: Cell type proportions for each spot in the sample.
# -----------------------------------------------------------------------------#

run_deconvolution <- function(st_meta, spatial_location, output_dir) {
  
  for (i in 1:nrow(st_meta)) {
    s1 <- st_meta[i, "sid"]
    type1 <- st_meta[i, "subtype"]
    
    # Load and prepare spatial transcriptomics count data
    stdata_s <- tryCatch({
      readMM(st_meta[i, "count_matrix_dir"])
    }, error = function(e) {
      # If an error occurs, try to read the matrix using read.csv
      tryCatch({
        read.csv(file = st_meta[i, "count_matrix_dir"], sep = "\t") %>%
          remove_rownames %>% column_to_rownames(var="X")
      }, error = function(e) {
        # If both attempts fail, print an error message
        cat("Both readMM and read.csv failed to read the file:", st_meta[i, "count_matrix_dir"], "\n")
      })
    })
    # Check if the matrix was successfully read
    if (is.null(stdata_s)) {
      cat("Failed to read the matrix from the file:", st_meta[i, "count_matrix_dir"], "\n")
    } else {
      cat("Count matrix successfully read.\n")
    }

    # Check if the rowname and colname exist
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
    # Verify the names are assigned correctly
    if (is.null(colnames(stdata_s)) || is.null(rownames(stdata_s))) {
      cat("Names are not assigned to dimensions.\n")
    } else {
      cat("Names are assigned to dimensions.\n")
    }

    stdata_s <- as(stdata_s, "sparseMatrix")
    stdata_s <- as(stdata_s, "CsparseMatrix")
    if(any(rownames(stdata_s) %in% spatial_location[,"spot_id"])){
      stdata_s <- t(stdata_s)
    }

    st_count <- stdata_s[!startsWith(rownames(stdata_s), "X__ambiguous"), ]

    # Check if gene ID of ST data need to be converted to gene symbol to match with the reference
    if (any(startsWith(rownames(st_count), "ENSG"))) {
      # Retrieve the gene symbols for ST data
      symbols_st <- get_gene_symbols(rownames(st_count))
      colnames(symbols_st) <- c("id.st", "symbol.st")

      # Match symbols between scRNA-seq and ST data
      symbols_matched <- symbols_st[symbols_st$symbol.st %in% sc_gene$symbol.sc, ]
      colnames(symbols_matched) <- c("id.st", "symbol.sc")
      # ST genes matched with more than one reference gene
      id_multi <- symbols_matched$id.st[duplicated(symbols_matched$id.st)]
      # Reference genes matched with more than one gene in ST
      symbol_multi <- symbols_matched$symbol.sc[duplicated(symbols_matched$symbol.sc)]
      # Extract ST genes that match with a common reference gene
      id_common_match <- unlist(symbols_matched[symbols_matched$symbol.sc %in% symbol_multi, "id.st"])
      # Remove ST genes that share the same gene symbol in reference and ST genes that match with more than one reference gene
      genes_bijective <- symbols_matched[!symbols_matched$id.st %in% unique(c(id_multi, id_common_match)), ]

      # Convert ST gene ID to gene symbol
      st_count <- st_count[genes_bijective$id.st, ]
      rownames(st_count) <- genes_bijective$symbol.sc
    }

    # Select cells of the same subtype
    if(!is.na(type1)){
      cells2use <- which(sc_meta$subtype == type1)
    }else{
      cells2use <- 1:nrow(sc_meta)
    }
    sc_count_i <- sc_count[, cells2use]
    sc_meta_i <- sc_meta[cells2use, ]
    
    # Prepare spatial location
    spatial_location_i <- spatial_location %>%
      dplyr::filter(spot_id %in% colnames(st_count) & sid == s1) %>%
      column_to_rownames("spot_id") %>%
      dplyr::select(x,y)

    if(nrow(spatial_location_i) > ncol(st_count)){
      spatial_location_i <- spatial_location_i[colnames(st_count), ]
    }else{
      st_count <- st_count[, rownames(spatial_location_i)]
    }

    # Create CARD object
    cat(s1, "Create CARD object.\n")
    cat(sprintf("dim(sc_count) for subtype %s:", type1), dim(sc_count_i), "\n")
    cat("dim(st_count): ", dim(st_count), "\n")
    # cell_type_res <- "Annotation"
    cell_type_res <- "celltype_major"

    CARD_obj <- createCARDObject(
      sc_count = sc_count_i,
      sc_meta = sc_meta_i,
      spatial_count = st_count,
      spatial_location = spatial_location_i,
      ct.varname = cell_type_res,
      ct.select = unique(sc_meta_i[[cell_type_res]]),
      sample.varname = "orig.ident",
      minCountGene = 100,
      minCountSpot = 5
    )
    gc()
    
    # Run deconvolution
    cat("start CARD deconvolution.\n")
    CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)
    cat("finish CARD deconvolution.\n")

    dir.create(output_dir, showWarnings = T)

    # Save CARD object and results
    saveRDS(CARD_obj, file = file.path(output_dir, sprintf("/CARD_obj_%s_%s.rds", s1, cell_type_res)))
    amat <- CARD_obj@algorithm_matrix
    saveRDS(amat, file = file.path(output_dir, sprintf("/CARD_alg_matrix_%s_%s.rds", s1, cell_type_res)))

    # Get proportions of each spot and save results
    prop <- as.data.frame(CARD_obj@Proportion_CARD)
    prop$sid <- s1
    write.csv(prop, file.path(output_dir, sprintf("/Proportion_%s_%s.csv", s1, cell_type_res)))
  }
}


# -----------------------------------------------------------------------------#
# Visualize the cell type proportions and spatial distribution for each sample.
# Args:
#     st_meta (data.frame): Metadata for spatial transcriptomics samples.
#     output_dir (str): Directory to save the visualizations.
#     ct.visualize (vector): Cell types of interest to visualize.
#     width (numeric): Width of the output plots.
#     pie_radius (numeric): Radius of the pie charts in the plots.
#     point_size (numeric): Size of the points in the proportion plots.
#     x_reverse (bool): Whether to reverse the x-axis.
#     y_reverse (bool): Whether to reverse the y-axis.
#     xy_flip (bool): Whether to flip the x and y axes.
# Outputs:
#     Saves the following plots to the output directory:
#     - cell_type_proportion_<sampleID>.jpg: Visualization of the cell type proportions for all spots.
#     - specific_type_proportion_<sampleID>.jpg: Visualization of the spatial distribution of the cell type proportions.
# -----------------------------------------------------------------------------#
visualize_proportions <- function(st_meta, output_dir, ct.visualize = NULL, 
                                  width = 10, point_size = 2, pie_radius = 1,
                                  x_reverse = F, y_reverse = F, xy_flip = F) {
  colors = c(`B-cells` = "#2a5fbd",
             CAFs = "#7ddffa",
             `Cancer Epithelial` = "#117d30",
             Endothelial = "#71f25a", 
             Myeloid = "#ebc857",
             `Normal Epithelial` ="#D39200",
             Plasmablasts = "#F8766D",
             PVL = "#DB72FB",
             `T-cells` = "#bd2a84")
  cell_type_res = "celltype_major"
  
  for(i in 1:nrow(st_meta)){
    s1 = st_meta[i, "sid"]
    # Load deconvoluted object
    CARD_obj <- readRDS(file.path(output_dir, sprintf("/CARD_obj_%s_%s.rds", s1, cell_type_res)))
    
    h_w_ratio <-  (max(CARD_obj@spatial_location$y)-min(CARD_obj@spatial_location$y)) / (max(CARD_obj@spatial_location$x)-min(CARD_obj@spatial_location$x))
    # Get the proportion of cell types of each cell
    prop <- as.data.frame(CARD_obj@Proportion_CARD)
    
    # Visualization of the proportions of all spots
    p1 <- CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD, 
                             spatial_location = CARD_obj@spatial_location, 
                             colors = colors,
                             radius = pie_radius)
    
    # Visualization of the spatial distribution of the cell type proportion
    p2 <- CARD.visualize.prop(
      proportion = CARD_obj@Proportion_CARD,        
      spatial_location = CARD_obj@spatial_location, 
      ct.visualize = ct.visualize,                 
      colors = c("lightblue","lightyellow","red"), 
      NumCols = 3,
      pointSize = point_size) 
    
    if(x_reverse){
      p1 <- p1 + scale_x_reverse()
      p2 <- p2 + scale_x_reverse()
    }
    if(y_reverse){
      p1 <- p1 + scale_y_reverse()
      p2 <- p2 + scale_y_reverse()
    }
    if(xy_flip){
      p1 <- p1 + coord_flip()
      p2 <- p2 + coord_flip()
      h_w_ratio <- 1/h_w_ratio
    }
    
    # Save the figures
   
    jpeg(file.path(output_dir, paste0("cell_type_proportion_", s1, ".jpg")), 
         width = width, height = width*h_w_ratio, units = 'in', res = 300)
    print(p1)
    dev.off()
    
    jpeg(file.path(output_dir, paste0("specific_type_proportion_", s1, ".jpg")), 
         width = width, height = width*h_w_ratio, units = 'in', res = 300)
    print(p2)
    dev.off()
    
    gc()
  }
}



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


setwd("/Users/zhiningsui/GitHub/STpath/data/STImage-1K4M/Visium")
input <- read.csv("create_patches_input.csv")[,-1]
fns <- list.files("/Users/zhiningsui/GitHub/STImage-1K4M/Visium/coord")

ids <- gsub("_coord.csv", "", fns)
info <- input[input$sample_id %in% ids, ]
info$spatial_path <- paste0("/Users/zhiningsui/GitHub/STImage-1K4M/Visium/coord/", info$sample_id, "_coord.csv")
all(file.exists(info$spatial_path))
info$image_path <- paste0("/Users/zhiningsui/GitHub/STImage-1K4M/Visium/image/", info$sample_id, ".png")
all(file.exists(info$image_path))

for (id in info$sample_id ) {
  sp <- read.csv(info[info$sample_id == id, "spatial_path"])
  info[info$sample_id == id, "diameter"] <- unique(sp$r) * 2
}

info$pxl_x_col <- 3
info$pxl_y_col <- 2
info$barcode_col <- 1
info$output_dir <- "../../output/STImage-1K4M/Visium/patch/"
write.csv(info, "create_patches_input.csv", row.names = F)
 


###############
###############
setwd("/Users/zhiningsui/GitHub/STpath/data/STImage-1K4M/ST")
fns <- list.files("/Users/zhiningsui/GitHub/STImage-1K4M/ST/coord")
ids <- gsub("_coord.csv", "", fns)
info <- data.frame(sample_id = ids)
info$spatial_path <- paste0("/Users/zhiningsui/GitHub/STImage-1K4M/ST/coord/", info$sample_id, "_coord.csv")
all(file.exists(info$spatial_path))
info$image_path <- paste0("/Users/zhiningsui/GitHub/STImage-1K4M/ST/image/", info$sample_id, ".png")
all(file.exists(info$image_path))

for (id in info$sample_id ) {
  sp <- read.csv(info[info$sample_id == id, "spatial_path"])
  info[info$sample_id == id, "diameter"] <- unique(sp$r) * 2
}

info$pxl_x_col <- 3
info$pxl_y_col <- 2
info$barcode_col <- 1
info$output_dir <- "../../output/STImage-1K4M/Visium/patch/"
write.csv(info, "create_patches_input.csv", row.names = F)

