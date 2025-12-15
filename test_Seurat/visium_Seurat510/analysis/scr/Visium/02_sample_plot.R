# AIM ---------------------------------------------------------------------
# sample load the full raw dataset

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(hdf5r)
library(tidyverse)
library(DropletUtils)
library(patchwork)

# read in the file --------------------------------------------------------
# save the list of raw imported objects
list_sobj <- readRDS(file = "../out/object/list_spobj_filtered.rds")


SpatialFeaturePlot(list_sobj$`1142243F`, features = "nCount_Spatial") + theme(legend.position = "right")
