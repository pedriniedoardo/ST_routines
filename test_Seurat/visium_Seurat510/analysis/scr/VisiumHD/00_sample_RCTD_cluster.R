# ======================================================================
# == load libraries ==
# ======================================================================

library(Seurat)
library(spacexr)
library(rhdf5)
library(tidyverse)
library(schard)

# ======================================================================
# == custom functions ==
# ======================================================================

# this can be swapped with schard implementation
# just make suret to change the slot accordingly: sobj@assays$Spatial or sobj@assays$RNA
# h5ad_to_seurat = function(file_in) {

#   stopifnot(require(rhdf5))
#   stopifnot(require(Matrix))
#   stopifnot(require(Seurat))

#   cell_counts = h5read(file = file_in, 'X')
#   cell_obs = h5read(file = file_in, 'obs')
#   cell_var = h5read(file = file_in, 'var')

#   n_rows <- length(cell_counts$indptr) - 1  # Number of rows in the matrix
#   n_cols <- max(cell_counts$indices) + 1    # Number of columns (assuming 0-based indices)

#   i = rep(1:n_rows, times = diff(cell_counts$indptr))
#   j = cell_counts$indices + 1
#   x = as.numeric(cell_counts$data)

#   dcg_mat = sparseMatrix(
#     i = i, j = j, x = as.numeric(x), dims = c(n_rows, n_cols),
#     dimnames = list(
#       as.character(cell_obs$`_index`),
#       as.character(cell_var$`_index`)
#     ))

#   CreateSeuratObject(CreateAssayObject(counts = t(dcg_mat)), meta.data = cell_obs)

# }

classify_type = function(RCTD_results_nuclei, RCTD_results_cells) {
  cg = intersect(rownames(RCTD_results_nuclei),
                 rownames(RCTD_results_cells))
  df_nuclei = RCTD_results_nuclei[cg,]
  df_cells = RCTD_results_cells[cg,]
  feature_type = ifelse(df_nuclei$first_type == as.character(df_cells$first_type), "cells", "nuclei")
  feature_type[df_nuclei$spot_class != "singlet" & df_cells$spot_class != "singlet"] = "rejected"
  feature_type[df_nuclei$spot_class == "singlet" & df_cells$spot_class != "singlet"] = "nuclei"
  feature_type[df_nuclei$spot_class != "singlet" & df_cells$spot_class == "singlet"] = "cells"
  feature_class = ifelse(feature_type=="cells", as.character(df_cells$first_type), as.character(df_nuclei$first_type))
  feature_class[feature_type=="rejected"] = "rejected"
  data.frame(type=feature_type, class=feature_class, row.names = cg)
}


# ======================================================================
# == Snakemake integation ==
# ======================================================================

# # Access inputs/outputs/params from Snakemake object
# reference_file <- snakemake@input$ref
# print(reference_file)
# nuclei_grouped_file <- snakemake@input$nuclei_grouped_adata
# print(nuclei_grouped_file)
# nuclei_expanded_file <- snakemake@input$nuclei_expanded_adata
# print(nuclei_expanded_file)
# output_csv <- snakemake@output$csv_filters
# print(output_csv)
# 
# message("Reference file: ", reference_file)
# message("Nuclei grouped anndata: ", nuclei_grouped_file)
# message("Nuclei expanded anndata: ", nuclei_expanded_file)
# message("Output CSV: ", output_csv)

reference_file <- "../../../../test_data/references/reference_tabula_muris_tiny.rds"
print(reference_file)
nuclei_grouped_file <- "../../../../test_data/Mouse_Embryo/Mouse_Embryo_nuclei_grouped_geometry.h5ad"
print(nuclei_grouped_file)
nuclei_expanded_file <- "../../../../test_data/Mouse_Embryo/Mouse_Embryo_nuclei_expanded_geometry.h5ad"
print(nuclei_expanded_file)
output_csv <- "../../../../test_data/Mouse_Embryo/Mouse_Embryo_RCTDedo_filters.csv"
print(output_csv)


# ======================================================================
# == load input ==
# ======================================================================

message("Load the inputs.")
# seu_obj <- h5ad_to_seurat("../../../../data/Mouse_Embryo/Mouse_Embryo_nuclei_grouped.h5ad")

reference <- readRDS(reference_file)
message("Loadaded ref.")

# seu_obj_nuclei <- h5ad_to_seurat(nuclei_grouped_file)
seu_obj_nuclei <- schard::h5ad2seurat_spatial(nuclei_grouped_file)
# seu_obj_nuclei <- readRDS(nuclei_grouped_file)
# meta_obs_nuclei <- seu_obj_nuclei@meta.data[,c("array_row", "array_col")]
message("Loadaded nuclei.")

# seu_obj_cell <- h5ad_to_seurat(nuclei_expanded_file)
seu_obj_cell <- schard::h5ad2seurat_spatial(nuclei_expanded_file)
# seu_obj_cell <- readRDS(nuclei_expanded_file)
# meta_obs_cell <- seu_obj_cell@meta.data[,c("array_row", "array_col")]
message("Loadaded cells.")

# ======================================================================
# == RCTD processing ==
# ======================================================================

message("Starting RCTD analysis script.")

# ADATA: NUCLEI
# spatial dataset
counts_hd_nuclei <- seu_obj_nuclei[["Spatial"]]$counts
obj_br_hd_nuclei <- colnames(seu_obj_nuclei[["Spatial"]])
coords_nuclei <- GetTissueCoordinates(seu_obj_nuclei)[obj_br_hd_nuclei, 1:2]

# Stefano explained this is legal as the stardist processing do not produce a whole counting per cell.
query_nuclei <- SpatialRNA(coords_nuclei, round(counts_hd_nuclei), colSums(round(counts_hd_nuclei)))

# run RCTD
RCTD_nuclei <- create.RCTD(query_nuclei, reference, max_cores = 20)
message("RCTD nuclei created.")
# use "doublet" in this case instead of "full", which is more appropriate for VisiumHD dataset with higher resolution
RCTD_nuclei <- run.RCTD(RCTD_nuclei, doublet_mode = "doublet")
message("RCTD nuclei run clompleted.")

RCTD_nuclei_result <- RCTD_nuclei@results
# write.csv(RCTD_nuclei$results_df,"../../../../data/Mouse_Embryo/Mouse_Embryo_RCTD_nuclei_tiny.csv")


# ADATA: CELLS
counts_hd_cells <- seu_obj_cell[["Spatial"]]$counts
obj_br_hd_cells <- colnames(seu_obj_cell[["Spatial"]])
coords_cells <- GetTissueCoordinates(seu_obj_cell)[obj_br_hd_cells, 1:2]

# Stefano explained this is legal as the stardist processing do not produce a whole counting per cell.
query_cells <- SpatialRNA(coords_cells, round(counts_hd_cells), colSums(round(counts_hd_cells)))

# run RCTD
RCTD_cells <- create.RCTD(query_cells, reference, max_cores = 20)
message("RCTD cells created.")
# use "doublet" in this case instead of "full", which is more appropriate for VisiumHD dataset with higher resolution
RCTD_cells <- run.RCTD(RCTD_cells, doublet_mode = "doublet")
message("RCTD cells run clompleted.")

RCTD_cells_result <- RCTD_cells@results
# write.csv(RCTD_cells$results_df,"../../../../data/Mouse_Embryo/Mouse_Embryo_RCTD_cells_tiny.csv")


# PARTE 2: FILTRI ON RCTD RESULTS

# INPUT: 
# RCTD_cells = read.csv(paste0(path, "RCTD_cells.csv"), row.names = "X")
# RCTD_nuclei = read.csv(paste0(path, "RCTD_nuclei.csv"), row.names = "X")

df_class_RCTD <- classify_type(RCTD_results_nuclei = RCTD_nuclei_result$results_df, RCTD_results_cells = RCTD_cells_result$results_df) 
table(df_class_RCTD$type)    # vedo quante ok e quante rejected da filtri RCTD
table(df_class_RCTD$class)   # number cell types detected

write.csv(df_class_RCTD, output_csv)
message("RCTD filters saved.")
