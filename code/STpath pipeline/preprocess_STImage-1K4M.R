
st_dir <- "/Users/zhiningsui/GitHub/STimage-1K4M"
github_dir <- "/Users/zhiningsui/GitHub/STpath/data/STimage-1K4M"
sc_dir <- "/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/Wu_2021/data/scRNASeq/"
output_dir <- "../../output/STimage-1K4M/deconvolution/ST"

# Create input files for create_patches.py and deconvolution -------------------

metadata <- read.csv(file.path(github_dir, "metadata_brca.csv"))
metadata_visium <- metadata[metadata$tech == "Visium", ]
metadata_visium$subtype <- ifelse(metadata_visium$subtype %in% c("ER", "HER2"),
                                  paste0(metadata_visium$subtype, "+"), 
                                  metadata_visium$subtype)
metadata_visium$subtype <- ifelse(metadata_visium$subtype == "Luminal",
                                  "ER+", 
                                  metadata_visium$subtype)
write.csv(metadata_visium, file = file.path(github_dir, "metadata_brca_visium.csv"),
          row.names = F)


mtx <- readMM("/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/Wu_2021/data/filtered_count_matrices/CID4290_filtered_count_matrix/matrix.mtx.gz")
mtx <- as.matrix(mtx)
barcodes <- read.csv("/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/Wu_2021/data/filtered_count_matrices/CID4290_filtered_count_matrix/barcodes.tsv.gz",
                     header = F)
features <- read.csv("/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/Wu_2021/data/filtered_count_matrices/CID4290_filtered_count_matrix/features.tsv.gz",
                     header = F)
colnames(mtx) <- barcodes$V1
rownames(mtx) <- features$V1
mtx <- t(mtx)
rownames(mtx) <- paste0("Human_Breast_Wu_06052021_Visium_CID4290_", rownames(mtx))
write.csv(mtx, "/Users/zhiningsui/GitHub/STimage-1K4M/Visium/gene_exp/Human_Breast_Wu_06052021_Visium_CID4290_count.csv")


tmp <- read.csv("/Users/zhiningsui/Library/CloudStorage/Dropbox/st2image_data/Wu_2021/data/spatial/CID4290_spatial/tissue_positions_list.csv", 
                header = F)

tmp <- tmp[, c(1, 5, 6)]
tmp$V1 <- paste0("Human_Breast_Wu_06052021_Visium_CID4290_", tmp$V1)
tmp <- tmp %>%
  column_to_rownames("V1")
colnames(tmp) <- c("yaxis", "xaxis")
  
metadata_visium$sid <- metadata_visium$slide
sp_info <- metadata_visium[, c("sid", "subtype", "tech")]
sp_info$dir <- file.path(st_dir, sp_info$tech)
sp_info$spatial_path <- file.path(sp_info$dir, "coord", paste0(sp_info$sid, "_coord.csv"))
all(file.exists(sp_info$spatial_path))
sp_info <- sp_info[file.exists(sp_info$spatial_path), ]
sp_info$image_path <- file.path(sp_info$dir, "image", paste0(sp_info$sid, ".png"))
all(file.exists(sp_info$image_path))
sp_info <- sp_info[file.exists(sp_info$image_path), ]
sp_info$output_dir <- file.path(sp_info$dir, "patch") 
sp_info$pxl_x_col <- 3
sp_info$pxl_y_col <- 2
sp_info$barcode_col <- 1
for (id in sp_info$sid ) {
  sp <- read.csv(sp_info[sp_info$sid == id, "spatial_path"])
  sp_info[sp_info$sid == id, "diameter"] <- unique(sp$r) * 2
}
sp_info$count_matrix_dir <- file.path(sp_info$dir, "gene_exp", paste0(sp_info$sid, "_count.csv"))
write.csv(sp_info, 
          file.path(github_dir, "create_patches_input_brca.csv"),
          row.names = F)


directory <- "/Users/zhiningsui/GitHub/STimage-1K4M/Visium/image"
all_files <- list.files(directory)
files_to_delete <- all_files[!all_files %in% paste0(sp_info$slide, ".png")]
file_paths <- file.path(directory, files_to_delete)
file.remove(file_paths)






