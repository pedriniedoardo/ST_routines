# AIM ---------------------------------------------------------------------
# conver h5ad object to rds object

# python version ----------------------------------------------------------
Sys.setenv(RETICULATE_PYTHON = "/home/edo/micromamba/envs/env_sc/bin/python")
library(reticulate)
reticulate::use_condaenv(condaenv = "/home/edo/micromamba/envs/env_sc")

py_config()

# libraries ---------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(Seurat)
library(sceasy)

# convert the file --------------------------------------------------------
fileh5ad <- dir("../../data/Schirmer_ST", full.names = TRUE, recursive = TRUE) %>%
  str_subset("h5ad")

filerds <- fileh5ad %>%
  str_replace_all(".h5ad",".rds")

# generate all the rds files from the h5ad files
pmap(list(fileh5ad,filerds), function(h5ad,rds){
  sceasy::convertFormat(h5ad, from="anndata", to="seurat",outFile=rds)
})
