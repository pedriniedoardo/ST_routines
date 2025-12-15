# AIM ---------------------------------------------------------------------
# sample load ST data
# this script relies on running the 02_ST_folder_to_h5_conversion.R and 02_ST_h5ad_to_rds_conversion.R

# libraries ---------------------------------------------------------------
library(Seurat)
library(patchwork)
library(tidyverse)


# define the folders ------------------------------------------------------
folder <- dir("../../data/Schirmer_ST", full.names = TRUE)
sample_name <- str_remove(folder, "../../data/Schirmer_ST/")

# load the raw seq data ---------------------------------------------------
# this data have all the raw data

# x <- "../../data/Schirmer_ST/GSM8563697"
list_obj <- lapply(folder, function(x){
  # follow the process
  print(x)
  
  # read in the coordinates of the spots
  LUT_barcode <- read_csv(file.path(x,"spatial/tissue_positions.csv"))
  
  # read in the image
  img <- Read10X_Image(image.dir = file.path(x,"spatial/"),
                       image.name = "tissue_hires_image.png")
  
  # load the spatial data
  spobj <- Load10X_Spatial(file.path(x),filename = "filtered_feature_bc_matrix.h5",image = img)
  
  # Add array_row and array_col to metadata
  spobj@meta.data <- spobj@meta.data %>%
    rownames_to_column("barcode") %>%
    left_join(LUT_barcode,by = "barcode") %>%
    column_to_rownames("barcode")
  
  # change the scale factor since the image is highres
  spobj@images$slice1@scale.factors$lowres <- spobj@images$slice1@scale.factors$hires
  
  return(spobj)
}) %>%
  setNames(sample_name)

# read in all the analyzed data to recover the annotation of the autors
sample_name_rds <- dir("../../data/Schirmer_ST", full.names = TRUE,recursive = T) %>%
  str_subset("rds")

# define the LUT for the file name to sample conversion
LUT_sample_rds <- data.frame(sample_name_rds) %>%
  mutate(sample_name = str_extract(sample_name_rds,"GSM[0-9]+"))

# rds <- LUT_sample_rds$sample_name_rds[1]
# name <- LUT_sample_rds$sample_name[1]
# add the anntation to the raw data
list_obj_annotated <- pmap(list(LUT_sample_rds$sample_name_rds,
                                LUT_sample_rds$sample_name), function(rds,name){
  # follow the process
  print(name)
  
  # sample read in a processed dataset
  data <- readRDS(rds)
  data@meta.data
  
  # add the metadata to the final object
  test<- list_obj[[name]]
  
  # update the metadata
  test@meta.data <- 
    test@meta.data %>%
    rownames_to_column("barcode") %>%
    left_join(data@meta.data %>% 
                rownames_to_column("barcode") %>%
                select(-c(array_row,array_col)),by = "barcode") %>%
    column_to_rownames("barcode")
  
  return(test)
}) %>%
  setNames(LUT_sample_rds$sample_name)

# sanity check ------------------------------------------------------------
dim(list_obj_annotated$GSM8563703)
dim(list_obj$GSM8563703)

dim(list_obj_annotated$GSM8563703@meta.data)
dim(list_obj$GSM8563703@meta.data)

# plot the sample image ---------------------------------------------------
# this will plot only the reads per spot
pdf("../../out/plot/panel_schirmer_all.pdf",width = 10,height = 5)
pmap(list(list_obj_annotated,names(list_obj_annotated)),function(x,name){
  # check the progress
  # print(name)
  # save the plot
  # plot <- SpatialFeaturePlot(x, features = "nCount_Spatial") + theme(legend.position = "top") + plot_annotation(name)
  # print(plot)
  p1 <- SpatialFeaturePlot(x, features = "nCount_Spatial",alpha=0,crop = F) + theme(legend.position = "none")
  p2 <- SpatialFeaturePlot(x, features = "nCount_Spatial",crop = F) + theme(legend.position = "right")
  
  wrap_plots(p1, p2) + plot_annotation(name)
})
dev.off()

# this will plot the annotation from the authors with relative count per area
pdf("../../out/plot/panel_schirmer_analyzed.pdf",width = 10,height = 5)
pmap(list(list_obj_annotated,names(list_obj_annotated)),function(x,name){
  
  # add it to the count dataset
  p1 <- Seurat::SpatialPlot(x,group.by = "areas") + theme(legend.position = "none")
  # violin plot reads per area
  p2 <- Seurat::VlnPlot(x,features = "nCount_Spatial",group.by = "areas") + theme(legend.position = "right")
  # SpatialFeaturePlot(x, features = "nCount_Spatial",alpha=0,crop = F) + theme(legend.position = "none")
  # SpatialFeaturePlot(x, features = "nCount_Spatial",crop = F) + theme(legend.position = "right")
  
  wrap_plots(p1, p2) + plot_annotation(name)
})
dev.off()

# focus on 703 ------------------------------------------------------------
test_rds <- list_obj_annotated$GSM8563703
p1 <- Seurat::SpatialPlot(test_rds,group.by = "areas") + theme(legend.position = "none")
p2 <- Seurat::VlnPlot(test_rds,features = "nCount_Spatial",group.by = "areas") + theme(legend.position = "right")
p12 <- SpatialFeaturePlot(test_rds, features = "nCount_Spatial",alpha=0,crop = T) + theme(legend.position = "none")
p22 <- SpatialFeaturePlot(test_rds, features = "nCount_Spatial",crop = T) + scale_fill_viridis_c(option = "turbo",breaks = c(0, 5000, 10000),limits = c(0, 10000)) + theme(legend.position = "right")

((p12+p1+p22)/(p2))+ plot_annotation("GSM8563703")

# try to fill in the NA classification
test_rds_fix <- test_rds

test_rds@meta.data %>%
  filter(is.na(areas) & array_row < 50 & array_col < 50) %>%
  ggplot(aes(x=array_col,y=array_row)) +geom_point()

test_rds@meta.data %>%
  filter(is.na(areas) & array_row < 50 & array_col < 50) %>%
  filter(is.na(areas) & array_row > 26 & array_col < 35) %>%
  ggplot(aes(x=array_col,y=array_row)) +geom_point()

# fis the meta
meta_test_fix <- test_rds@meta.data %>%
  # filter(is.na(areas)) %>%
  mutate(area_fix = case_when(is.na(areas) & array_row < 50 & array_col < 50 ~ "PPWM",
                              T ~ areas)) %>%
  mutate(area_fix = case_when(is.na(areas) & array_row >26 & array_col < 35 ~ NA,
                              T ~ area_fix))

test_rds_fix@meta.data <- meta_test_fix
# Seurat::SpatialPlot(test_rds_fix,group.by = "area_fix")
p1 <- Seurat::SpatialPlot(test_rds_fix,group.by = "area_fix") + theme(legend.position = "none")
p2 <- Seurat::VlnPlot(test_rds_fix,features = "nCount_Spatial",group.by = "area_fix") + theme(legend.position = "right")
wrap_plots(p1, p2) + plot_annotation("GSM8563703")

test_rds_fix_subset <- subset(test_rds_fix,subset = area_fix %in% c("PPWM","LR","LC","GM","Claustrum"))
p12 <- Seurat::SpatialPlot(test_rds_fix_subset,group.by = "area_fix") + theme(legend.position = "none")
p22 <- Seurat::VlnPlot(test_rds_fix_subset,features = "nCount_Spatial",group.by = "area_fix") + theme(legend.position = "right")
wrap_plots(p1, p2) + plot_annotation("GSM8563703")

# SpatialFeaturePlot(list_obj[[1]], features = "nCount_Spatial",alpha=0,crop = F) + theme(legend.position = "none")
# SpatialFeaturePlot(list_obj[[1]], features = "nCount_Spatial") + theme(legend.position = "right")
# SpatialFeaturePlot(list_obj[[1]], features = "nCount_Spatial",) + theme(legend.position = "right")
# SpatialFeaturePlot(spobj, features = "nCount_Spatial",interactive = F) + theme(legend.position = "right")
# SpatialFeaturePlot(spobj,interactive = T,features = "nCount_Spatial") + theme(legend.position = "right")
