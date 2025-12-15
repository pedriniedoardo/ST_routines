# AIM ---------------------------------------------------------------------
# manual draw the regions where to calculate the gradients

# libraries ---------------------------------------------------------------
library(Seurat)
library(SPATA2)
library(tidyverse)

# read in the object ------------------------------------------------------
list_obj <- readRDS("/media/edo/sandiskSSD/work/HSR/project_absinta/spatial_visium/220924_Visium_BrainMS_Martina/analysis/GitHub/Visium_brainMS_single_sample/out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")

# create the spata list of objects
list_SPATA2 <- lapply(list_obj, function(x) {
  # check the progress
  print(x)
  test <- SPATA2::asSPATA2(x,
                           sample_name = "SeuratProject",
                           platform = "VisiumSmall",
                           img_name = "slice1",
                           img_scale_fct = "lowres",
                           assay_name = "Spatial",
                           assay_modality = "gene")
  return(test)
})

# manual drawing ----------------------------------------------------------
# V01 3 mm
list_SPATA2$V01 <- createSpatialTrajectories(list_SPATA2$V01)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V01, 
  ids = "test", 
  color_by = "manual_segmentation")

# V02 3 mm
list_SPATA2$V02 <- createSpatialTrajectories(list_SPATA2$V02)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V02, 
  ids = "test", 
  color_by = "manual_segmentation")

# V03 3 mm
list_SPATA2$V03 <- createSpatialTrajectories(list_SPATA2$V03)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V03, 
  ids = "test", 
  color_by = "manual_segmentation")

# V04 3 mm
list_SPATA2$V04 <- createSpatialTrajectories(list_SPATA2$V04)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V04, 
  ids = "test", 
  color_by = "manual_segmentation")

# V05 3 mm
list_SPATA2$V05 <- createSpatialTrajectories(list_SPATA2$V05)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V05, 
  ids = "test", 
  color_by = "manual_segmentation")

# V06 3 mm
list_SPATA2$V06 <- createSpatialTrajectories(list_SPATA2$V06)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V06, 
  ids = "test", 
  color_by = "manual_segmentation")

# V07 3 mm
list_SPATA2$V07 <- createSpatialTrajectories(list_SPATA2$V07)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V07, 
  ids = "test", 
  color_by = "manual_segmentation")

# V08 3 mm
list_SPATA2$V08 <- createSpatialTrajectories(list_SPATA2$V08)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V08, 
  ids = "test", 
  color_by = "manual_segmentation")

# V09 3 mm
list_SPATA2$V09 <- createSpatialTrajectories(list_SPATA2$V09)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V09, 
  ids = "test", 
  color_by = "manual_segmentation")

# V10 3 mm
list_SPATA2$V10 <- createSpatialTrajectories(list_SPATA2$V10)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V10, 
  ids = "test", 
  color_by = "manual_segmentation")

# V11 3 mm
list_SPATA2$V11 <- createSpatialTrajectories(list_SPATA2$V11)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V11, 
  ids = "test", 
  color_by = "manual_segmentation")

# V11 3 mm
list_SPATA2$V12 <- createSpatialTrajectories(list_SPATA2$V12)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V12, 
  ids = "test", 
  color_by = "manual_segmentation")

# V11 3 mm
list_SPATA2$V14 <- createSpatialTrajectories(list_SPATA2$V14)
# check the manual segmentation
plotSpatialTrajectories(
  object = list_SPATA2$V14, 
  ids = "test", 
  color_by = "manual_segmentation")

# save the full list
saveRDS(list_SPATA2,file = "../../out/object/100_list_brain_SPATA2.rds")

# plot all the lists
list_plot <- lapply(list_SPATA2,function(x){
  p <- plotSpatialTrajectories(
    object = x, 
    ids = "test", 
    color_by = "manual_segmentation")
  
  return(p)
})
wrap_plots(list_plot)
