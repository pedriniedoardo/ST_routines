# AIM ---------------------------------------------------------------------
# load sample dataset from image and count matrices

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# try to load a full object -----------------------------------------------
# identify the location of the output matrix
matrix_dir = 'filtered_feature_bc_matrix/' 
counts = Seurat::Read10X(data.dir = matrix_dir)  

data = Seurat::CreateSeuratObject(
  counts = counts , 
  project = 'test', 
  assay = 'Spatial')

data$slice = 1 
data$region = 'test' 

imgpath = "spatial"
img = Seurat::Read10X_Image(image.dir = imgpath)  

Seurat::DefaultAssay(object = img) <- 'Spatial'  

img = img[colnames(x = data)]  
data[['image']] = img  

SpatialFeaturePlot(data, features = "nCount_Spatial")


# try to load the output from BayesPrism ----------------------------------



