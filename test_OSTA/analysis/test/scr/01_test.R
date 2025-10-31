# python version ----------------------------------------------------------
Sys.setenv(RETICULATE_PYTHON = "/home/edo/micromamba/envs/env_sc/bin/python")
library(reticulate)
reticulate::use_condaenv(condaenv = "/home/edo/micromamba/envs/env_sc")

py_config()

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(sceasy)
library(reticulate)
library(zellkonverter)

# convert the file --------------------------------------------------------
h5ad_path = "cortex.h5ad"
cortex.seurat <- sceasy::convertFormat(h5ad_path, from="anndata", to="seurat",outFile='filename.rds')

cortex.sce <- as.SingleCellExperiment(cortex.seurat)

(cortex.ad <- SCE2AnnData(cortex.sce))
unique(cortex.ad$obs$cell_type)

cortex.sce2 <- AnnData2SCE(cortex.ad)

scanpy <- import("scanpy")

