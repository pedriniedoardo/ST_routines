# libraries ---------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(Seurat)

# read in the data --------------------------------------------------------
# test for the specific file
img <- Read10X_Image(image.dir = "../../data/GSE210616_RAW/",
                     image.name = "tissue_hires_image.png")

spobj <- Load10X_Spatial("../../data/GSE210616_RAW/",filename = "filtered_feature_bc_matrix.h5",image = img)

# change the scale factor since the image is highres
spobj@images$slice1@scale.factors$lowres <- spobj@images$slice1@scale.factors$hires

# plot the sample image
plot01 <- Seurat::SpatialPlot(spobj,pt.size.factor = 0,crop = F)+theme(legend.position = "none")
plot02 <- Seurat::SpatialPlot(spobj)
plot03 <- SpatialFeaturePlot(spobj, features = "nCount_Spatial",image.alpha = 0.2,pt.size.factor = 5) + theme(legend.position = "right")
plot04 <- VlnPlot(spobj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

wrap_plots(list(plot01,plot02,plot03,plot04))
ggsave("../../out/plot/01_test_GSE210616_RAW.pdf",width = 10,height = 10)

# SpatialFeaturePlot(spobj, features = "nCount_Spatial",interactive = F) + theme(legend.position = "right")
# SpatialFeaturePlot(spobj, features = "nFeature_Spatial",interactive = F) + theme(legend.position = "right")
# SpatialFeaturePlot(spobj, features = "nCount_Spatial",interactive = T) + theme(legend.position = "right")
# 
# SpatialFeaturePlot(spobj,interactive = T,features = "nCount_Spatial") + theme(legend.position = "right")
