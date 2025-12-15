# AIM ---------------------------------------------------------------------
# conver h5ad object to rds object

# python version ----------------------------------------------------------
Sys.setenv(RETICULATE_PYTHON = "/home/edo/micromamba/envs/env_scvi/bin/python")
library(reticulate)
reticulate::use_condaenv(condaenv = "/home/edo/micromamba/envs/env_scvi")

py_config()

# libraries ---------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(Seurat)
library(sceasy)

# convert the file --------------------------------------------------------
fileh5ad <- dir("../data/Schirmer_ST", full.names = TRUE, recursive = TRUE) %>%
  str_subset("h5ad")

filerds <- fileh5ad %>%
  str_replace_all(".h5ad",".rds")

# test one file
# sceasy::convertFormat("../data/Schirmer_ST/GSM8563697/misc/GSM8563697_visium_CO37.h5ad",
#                       from="anndata",
#                       to="seurat",
#                       outFile="../data/Schirmer_ST/GSM8563697/misc/GSM8563697_visium_CO37_test.rds")

# generate all the rds files from the h5ad files
pmap(list(fileh5ad,filerds), function(h5ad,rds){
  sceasy::convertFormat(h5ad, from="anndata", to="seurat",outFile=rds)
})
