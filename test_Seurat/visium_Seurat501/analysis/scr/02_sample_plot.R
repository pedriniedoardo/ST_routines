# AIM ---------------------------------------------------------------------
# sample load the full raw dataset

# libraries ---------------------------------------------------------------
library(Seurat)

# read in the file --------------------------------------------------------
# save the list of raw imported objects
list_sobj <- readRDS(file = "../../visium/out/object/list_spobj_filtered.rds")


SpatialFeaturePlot(list_sobj$`1142243F`, features = "nCount_Spatial") + theme(legend.position = "right")
