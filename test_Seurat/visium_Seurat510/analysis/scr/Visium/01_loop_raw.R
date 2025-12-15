# AIM ---------------------------------------------------------------------
# sample load the full raw dataset

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(hdf5r)
library(tidyverse)
library(DropletUtils)
library(patchwork)

# generate the .h5 file ---------------------------------------------------
# # define the sample id to programmatically read in all the folders
# sample_id <- dir("../data/misc/test_edo/raw/") %>%
#   str_remove(pattern = "sample_")
# 
# # loop the generation of the .h5 file regular folder
# # samp <- sample_id[1]
# lapply(sample_id,function(samp){
#   # track the progress
#   print(samp)
#   
#   # define the paths in input and output
#   path_in <- paste0("../data/misc/test_edo/raw/sample_",samp,"/",samp,"_raw_feature_bc_matrix/")
#   path_out <- paste0("../data/misc/test_edo/raw/sample_",samp,"/",samp,".h5")
#   
#   # read in the counts table
#   counts <- Read10X(data.dir = path_in,  gene.column=1)
#   
#   # generate the .h5 file
#   DropletUtils::write10xCounts(x = counts,path = path_out,type = "HDF5")
#   
# })
# generate the object -----------------------------------------------------

# define the sample id to programmatically read in all the folders
sample_id <- dir("../data/misc/test_edo/raw/") %>%
  str_remove(pattern = "sample_")

# loop the generation of the .h5 file regular folder
# samp <- sample_id[1]
list_sobj <- lapply(sample_id,function(samp){
  # track the progress
  print(samp)
  
  # define the paths in input and output
  path_in_image <- paste0("../data/misc/test_edo/raw/sample_",samp,"/",samp,"_spatial/")
  path_in_reads <- paste0("../data/misc/test_edo/raw/sample_",samp,"/")
  sample_h5 <- paste0(samp,".h5")
  path_in_meta <- paste0("../data/misc/test_edo/raw/sample_",samp,"/",samp,"_metadata.csv")
  
  # read in the image file. load the highres image
  img <- Read10X_Image(image.dir = path_in_image,
                       image.name = "tissue_hires_image.png")
  
  # load the full matrix of reads
  spobj <- Load10X_Spatial(path_in_reads,
                           filename = sample_h5,
                           image = img,
                           filter.matrix = F)
  
  # change the scale factor for the highres image
  spobj@images$slice1@scale.factors$lowres <- spobj@images$slice1@scale.factors$hires
  
  # read in the full meta
  meta_test <- read_csv(path_in_meta) %>%
    setNames(c("barcode","nCount_RNA","nFeature_RNA","patientid","subtype","Classification")) %>%
    # select(-c("nCount_RNA","nFeature_RNA")) %>%
    column_to_rownames("barcode")
  
  # add the metadata to the final object
  spobj_final <- AddMetaData(spobj,meta_test)
  
  return(spobj_final)
}) %>%
  setNames(sample_id)

# save the list of raw imported objects
saveRDS(list_sobj,file = "../out/object/list_spobj_raw.rds")

# explore the objects -----------------------------------------------------
#
sample_id <- names(list_sobj)

list_plot01 <- lapply(sample_id[1:2],function(samp){
  # track the progress 
  print(samp)
  
  # laod the file
  spobj <- list_sobj[[samp]]
  
  # generate the plots
  p1 <- SpatialFeaturePlot(spobj, features = "nCount_Spatial",interactive = F,image.alpha = 0,pt.size.factor = 3,max.cutoff = "q95") + theme(legend.position = "right")
  p2 <- SpatialFeaturePlot(spobj, features = "nCount_Spatial",interactive = F,image.alpha = 100,pt.size.factor = 0,max.cutoff = "q95") + theme(legend.position = "right")
  p3 <- SpatialDimPlot(spobj, group.by = "Classification",image.alpha = 0,pt.size.factor = 3) + theme(legend.position = "right")
  p4 <- SpatialFeaturePlot(spobj,crop = T,features = "nCount_Spatial",interactive = F,image.alpha = 50,pt.size.factor = 3,alpha = 0.5,max.cutoff = "q95") + theme(legend.position = "right")
  
  # assamble the plots
  (p1 + p2 + p3 + p4) +
    plot_annotation(samp,theme=theme(plot.title=element_text(hjust=0.5)))
})

list_plot02 <- lapply(sample_id[3:6],function(samp){
  # track the progress 
  print(samp)
  
  # laod the file
  spobj <- list_sobj[[samp]]
  
  # generate the plots
  p1 <- SpatialFeaturePlot(spobj, features = "nCount_Spatial",interactive = F,image.alpha = 0,pt.size.factor = 2,max.cutoff = "q95") + theme(legend.position = "right")
  p2 <- SpatialFeaturePlot(spobj, features = "nCount_Spatial",interactive = F,image.alpha = 100,pt.size.factor = 0,max.cutoff = "q95") + theme(legend.position = "right")
  p3 <- SpatialDimPlot(spobj, group.by = "Classification",image.alpha = 0,pt.size.factor = 2) + theme(legend.position = "right")
  p4 <- SpatialFeaturePlot(spobj,crop = T,features = "nCount_Spatial",interactive = F,image.alpha = 50,pt.size.factor = 2,alpha = 0.5,max.cutoff = "q95") + theme(legend.position = "right")
  
  # assamble the plots
  (p1 + p2 + p3 + p4) +
    plot_annotation(samp,theme=theme(plot.title=element_text(hjust=0.5)))
})



pdf("../out/plot/01_plot_raw.pdf",width = 22,height = 20)
c(list_plot01,list_plot02)
dev.off()


