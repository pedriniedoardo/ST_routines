# AIM ---------------------------------------------------------------------
# manual draw the regions where to calculate the gradients for a specific slide

# libraries ---------------------------------------------------------------
library(Seurat)
library(SPATA2)
library(tidyverse)
library(cowplot)
library(patchwork)

# -------------------------------------------------------------------------
list_obj2 <- readRDS("/media/edo/sandiskSSD/work/HSR/project_absinta/spatial_visium/220924_Visium_BrainMS_Martina/analysis/GitHub/Visium_brainMS_single_sample/out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")

list_test2 <- lapply(list_obj2["V04"], function(x) {
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

# draw the line
test_alt <- createSpatialTrajectories(list_test2$V04)

# check the manual segmentation
plotSpatialTrajectories(
  object = test_alt, 
  ids = "test_2", 
  color_by = "manual_segmentation")

list_sig <- readRDS("../../data/senescence_pathways.rds")
genes <- list_sig$senmayo$Genes %>% unique()

# both, opt1 and opt2, have the same effect (addSignature() just allows to specify the
# assay of interest, which is fixed to 'gene' for addGeneSet()).
# opt1
test_alt <- addSignature(test_alt, name = "senmayo", class = "Custom", molecules = genes)
plotStsLineplot(test_alt, variables ="Custom_senmayo", id = "test_2")+ggtitle("V04")

ggsave("../../out/plot/100_SPATA2_Trajectory_Screening_senmayo_V04.pdf", width = 4, height = 4)

# save the object
saveRDS(test_alt, "../../out/object/100_SPATA2_Trajectory_Screening_seurat_plot_V04.rds")
