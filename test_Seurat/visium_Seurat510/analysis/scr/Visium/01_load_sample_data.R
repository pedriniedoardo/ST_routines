# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(hdf5r)

# read in the data --------------------------------------------------------
# # test absinta files
# img_test <- Read10X_Image(image.dir = "data/Absinta/spatial/",
#                           image.name = "tissue_hires_image.png")
# 
# spobj_test <- Load10X_Spatial("data/Absinta/",filename = "filtered_feature_bc_matrix.h5",image = img_test)
# 
# # change the scale factor since the image is highres
# spobj_test@images$slice1@scale.factors$lowres <- spobj_test@images$slice1@scale.factors$hires
# 
# # plot the sample image
# SpatialFeaturePlot(spobj_test, features = "nCount_Spatial",interactive = T) + theme(legend.position = "right")

# test for the specific file
img <- Read10X_Image(image.dir = "../../data/GSE222322_visium_SpinalCord/GSM6919905/spatial/",
                     image.name = "tissue_hires_image.png")

spobj <- Load10X_Spatial("../../data/GSE222322_visium_SpinalCord/GSM6919905/",filename = "filtered_feature_bc_matrix.h5",image = img)

# change the scale factor since the image is highres
spobj@images$slice1@scale.factors$lowres <- spobj@images$slice1@scale.factors$hires

# plot the sample image
SpatialFeaturePlot(spobj, features = "nCount_Spatial",interactive = F) + theme(legend.position = "right")
# SpatialFeaturePlot(spobj, features = "nCount_Spatial",interactive = T) + theme(legend.position = "right")

SpatialFeaturePlot(spobj,interactive = T,features = "nCount_Spatial") + theme(legend.position = "right")


# # specify the high res image, which is the only one available in the project
# img <- Read10X_Image(image.dir = "../../data/GSE222322_visium_SpinalCord/GSM6919905/spatial",
#                      image.name = "tissue_hires_image.png")
# 
# img@scale.factors$lowres = img@scale.factors$hires # it is better to set the scale factors this way.
# spobj2 <- Load10X_Spatial("../../data/GSE222322_visium_SpinalCord/GSM6919905/",filename = "filtered_feature_bc_matrix.h5",image = img)
# SpatialFeaturePlot(spobj2, features = "nCount_Spatial",interactive = F) + theme(legend.position = "right")
