
st_dir <- "/Users/zhiningsui/GitHub/STimage-1K4M"

# Create input files for create_patches.py --------------------------------

metadata_all <- read.csv(file.path(st_dir, "meta/meta_all_gene.csv"))
sp_info <- metadata_all[, c("slide", "tech", "species", "tissue", "involve_cancer")]
sp_info$dir <- file.path(st_dir, sp_info$tech)
sp_info$spatial_path <- file.path(sp_info$dir, "coord", paste0(sp_info$slide, "_coord.csv"))
all(file.exists(sp_info$spatial_path))
sp_info <- sp_info[file.exists(sp_info$spatial_path), ]
sp_info$image_path <- file.path(sp_info$dir, "image", paste0(sp_info$slide, ".png"))
all(file.exists(sp_info$image_path))
sp_info <- sp_info[file.exists(sp_info$image_path), ]
sp_info$output_dir <- file.path(sp_info$dir, "patch") 
sp_info$pxl_x_col <- 3
sp_info$pxl_y_col <- 2
sp_info$barcode_col <- 1
for (id in sp_info$slide ) {
  sp <- read.csv(sp_info[sp_info$slide == id, "spatial_path"])
  sp_info[sp_info$slide == id, "diameter"] <- unique(sp$r) * 2
}
sp_info$count_matrix_dir <- file.path(sp_info$dir, "gene_exp", paste0(sp_info$slide, "_count.csv"))
write.csv(sp_info, 
          "/Users/zhiningsui/GitHub/STpath/data/STimage-1K4M/create_patches_input.csv",
          row.names = F)

sp_info <- sp_info[sp_info$involve_cancer == "True" &
                     sp_info$species == "human" &
                     sp_info$tissue == "breast", ]
write.csv(sp_info, 
          "/Users/zhiningsui/GitHub/STpath/data/STimage-1K4M/create_patches_input_brca.csv",
          row.names = F)

directory <- "/Users/zhiningsui/GitHub/STimage-1K4M/Visium/image"
all_files <- list.files(directory)
files_to_delete <- all_files[!all_files %in% paste0(sp_info$slide, ".png")]
file_paths <- file.path(directory, files_to_delete)
file.remove(file_paths)

# Create metadata for deconvolution ---------------------------------------

st_dir <- "/Users/zhiningsui/GitHub/STimage-1K4M/"

metadata_all <- read.csv("/Users/zhiningsui/GitHub/STpath/data/STImage-1K4M/meta/meta_all_gene.csv")
metadata_all$sid <- metadata_all$slide
st_meta <- metadata_all[metadata_all$sid %in% sp_info$slide, ]
st_meta$count_matrix_dir <- paste0(st_dir, st_meta$tech, "/gene_exp/", st_meta$sid, "_count.csv")
st_meta <- as.data.frame(st_meta)
write.csv(st_meta, "/Users/zhiningsui/GitHub/STpath/data/STImage-1K4M/metadata.csv", row.names = F)

sc_dir <- "/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/Wu_2021/data/scRNASeq/"
output_dir <- "../../output/STimage-1K4M/deconvolution/ST"

tmp <- read.csv(st_meta$count_matrix_dir[170])


# Prepare st_meta and spatial_location for each ST dataset ----------------
st_meta_all <- list()
spatial_all <- list()


# Create metadata for the ST data
metadata <- metadata_all[metadata_all$tech == "ST", ]
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










