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

run_deconvolution <- function(st_meta, spatial_location, output_dir, save_CARD_objs = T) {
  
  for (i in 1:nrow(st_meta)) {
    s1 <- st_meta[i, "sid"]
    if(sum(colnames(st_meta) %in% "subtype")==0){
      type1 <- NA
    } else {
      type1 <- st_meta[i, "subtype"]
    }
    
    # Load and prepare spatial transcriptomics count data
    stdata_s <- tryCatch({
      readMM(st_meta[i, "count_matrix_dir"])
    }, error = function(e) {
      # If an error occurs, try to read the matrix using read.csv
      tryCatch({
        read.csv(file = st_meta[i, "count_matrix_dir"], sep = "\t") %>%
          remove_rownames %>% column_to_rownames(var="X")
      }, error = function(e) {
        tryCatch({
          read.csv(file = st_meta[i, "count_matrix_dir"]) %>%
            remove_rownames %>% column_to_rownames(var="X")
        }, error = function(e){
          # If both attempts fail, print an error message
          cat("Both readMM and read.csv failed to read the file:", st_meta[i, "count_matrix_dir"], "\n")
        })
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
    if(save_CARD_objs){
      saveRDS(CARD_obj, file = file.path(output_dir, sprintf("/CARD_obj_%s_%s.rds", s1, cell_type_res)))
      amat <- CARD_obj@algorithm_matrix
      saveRDS(amat, file = file.path(output_dir, sprintf("/CARD_alg_matrix_%s_%s.rds", s1, cell_type_res)))
    }
    
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


