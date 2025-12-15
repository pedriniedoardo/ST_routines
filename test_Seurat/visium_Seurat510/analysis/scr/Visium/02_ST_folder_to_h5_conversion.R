# AIM ---------------------------------------------------------------------
# sample snippet to convert 10x genomics data from folder to .h5 file

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DropletUtils)
library(Seurat)
library(patchwork)

# convert folder to .h5 file ----------------------------------------------
# define the folder 
folder <- dir("../data/Schirmer_ST", full.names = TRUE)
sample_name <- str_remove(folder, "../data/Schirmer_ST/")

# run the conversion
# x <- "../../data/Schirmer_ST/GSM8563697"
lapply(folder,function(x){
  # check the progress
  print(x)
  
  # read in the folder data
  filter_matrix <- Read10X(file.path(x,"misc/filtered_feature_bc_matrix"))

  # write the .h5 file
  write10xCounts(file.path(x,"filtered_feature_bc_matrix.h5"), filter_matrix, type = "HDF5",
                 genome = "GRCh38", version = "3", overwrite = TRUE,
                 gene.id = rownames(filter_matrix),
                 gene.symbol = rownames(filter_matrix))
})
