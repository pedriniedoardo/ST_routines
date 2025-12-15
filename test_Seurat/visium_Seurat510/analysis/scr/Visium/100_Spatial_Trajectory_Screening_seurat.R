# AIM ---------------------------------------------------------------------
# sample processing of a seurat object following the SPATA2 workflow

# libraries ---------------------------------------------------------------
library(Seurat)
library(SPATA2)
library(tidyverse)

# read in the object ------------------------------------------------------
list_obj <- readRDS("/media/edo/sandiskSSD/work/HSR/project_absinta/spatial_visium/220924_Visium_BrainMS_Martina/analysis/GitHub/Visium_brainMS_single_sample/out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")

# create a sample SPATA2 object to test the gradient
test <- SPATA2::asSPATA2(list_obj$V01,
                         sample_name = "SeuratProject",
                         platform = "VisiumSmall",
                         img_name = "slice1",
                         img_scale_fct = "lowres",
                         assay_name = "Spatial",
                         assay_modality = "gene")

# show image 
plotImage(test)

# manually draw the trajectory using createSpatialTrajectories(). Make sure to save the trajectory.
# test <- createSpatialTrajectories(test)

# created with code 
# plotSpatialTrajectories(
#   object = object_t269, 
#   ids = "horizontal_mid", 
#   color_by = "histology")

# plot the trajectory
plotSpatialTrajectories(
  object = test, 
  ids = "test_1", 
  color_by = "manual_segmentation")

# 2. Run the algorithm ====================================================
# The function to use is called spatialTrajectoryScreening(). The parameter variables takes the numeric variables that are supposed to be included in the screening process. Since all sorts of numeric variables can be included in the screening, the argument for the input is simply called variables. Here, we are using the genes that were already identified as spatially variable by SPARKX. The goal is to further analyze which of the genes are expressed in a non-random and biologically meaningful way along the trajectory.

# this is a wrapper around SPARK::sparkx(), to identify features with non-random spatial expression pattern.
test <- runSPARKX(object = test)

# keep genes with a sparkx pvalue of 0.01 or lower
spark_df <- getSparkxGeneDf(object = test, threshold_pval = 0.01)

# show results
spark_df

# name of the trajectory
id <- "test_1"

# pull the gene ids. `getSparkxGenes()` would work, too
input_genes <- spark_df[["genes"]]

# the function spatialTrajectoryScreening, screens the sample for numeric variables that follow specific expression changes along the course of the spatial trajectory.
# we can use the genes ids or any other numeric variable in the medatata
# note: the results are NOT stored in the SPATA2 object but in a separate object
# in this case focus on a variable in the metadata
sts_out_metadata <- 
  spatialTrajectoryScreening(
    object = test, 
    # ID of the spatial trajectory
    id = id,
    # the variables/genes to include in the screening
    # variables = input_genes,
    variables = c("nCount_Spatial"),
    # keep all the genes expressed genes to calculate the trend
    sign_var = "fdr",
    sign_threshold = 1,
    rm_zero_infl = F
  )

# use the gene filterered from sparkx
sts_out_gene <- 
  spatialTrajectoryScreening(
    object = test, 
    # ID of the spatial trajectory
    id = id,
    # the variables/genes to include in the screening
    variables = input_genes,
    # keep all the genes expressed genes to calculate the trend
    sign_var = "fdr",
    sign_threshold = 1,
    rm_zero_infl = F
  )

# Note: The output of spatialTrajectoryScreening() is not saved in the SPATA2 object but returned in a separate S4 object of class SpatialTrajectoryScreening.
# Do not overwrite the SPATA2 object
sts_out_metadata@results$significance
sts_out_metadata@results$model_fits
sts_out_metadata@results$model_fits %>%
  group_by(variables) %>% 
  slice_min(mae, n = 1)

# check if SPP1 is there
sts_out_gene@results$significance %>% 
  filter(variables == "SPP1")

# 3. Results =============================================================
# The first step of the screening identifies pattern that are unlikely due to random gene expression. The second step fits the non-random gene expression pattern to predefined models which guides in interpretation and screening for specific gene expression pattern.

# Slot @results$significance contains a data.frame with one row for each screened variable which provides information regarding the degree of randomness the inferred pattern contains as quantified by the total variation (tot_var). The p-value gives the probability to obtain such a total variation under complete randomness and indicates the degree of significance. Column fdr contains the adjusted p-value according to the False Discovery Rate.

sign_df <- 
  sts_out_gene@results$significance %>% 
  filter(fdr < 0.05)

# show significance data.frame
sign_df


sts_out_gene@results$model_fits %>%
  filter(variables %in% c("SPP1"))

# 3.1 Non-random gene expression gradients +++++++++++++++++++++++++++++++++
# The figures below show examples for gene expression gradient identified as non-random.

# extract variables names
# non_random <- getSgsResultsVec(sts_out) %>% head(4)
non_random <- getSgsResultsVec(sts_out_gene)

trajectory_add_on <- ggpLayerSpatialTrajectories(object = test, ids = "test_1")

plotSurfaceComparison(
  object = test,
  color_by = non_random,
  pt_clrsp = "Reds 3",
  outline = T,
  nrow = 6
) + 
  trajectory_add_on

plotStsLineplot(test, variables = non_random, id = "test_1", line_color = "red", nrow = 6) 

# 3.2 Random gene expression gradients +++++++++++++++++++++++++++++++++++++
# The figures below show examples for gene expression gradient identified as random.

# extract random variable names 
random <- 
  sts_out@results$significance %>% 
  filter(fdr > 0.05) %>% 
  slice_max(tot_var, n = 4) %>% 
  pull(variables) %>% 
  head(4)

plotSurfaceComparison(
  object = test,
  color_by = random,
  pt_clrsp = "BuPu",
  outline = T, 
  nrow = 1
) + 
  trajectory_add_on

plotStsLineplot(test, variables = random, id = "test_1", line_color = "blue", nrow = 1) 

# Random vs. non-random
# Looking at the gradient of variable LEPROT (random, blue) and SHISA5 (non-random, red) both feature a descending pattern along the trajectory. How come that one is identified as most likely random while the other is not? Spatial gradient screening decides which variables are most likely random or not random by computing the variance along the gradient, which is stored in the variable tot_var in sts_out@results$significance. Consider both gradient plots with the expression estimates plotted as black points, too. While SHISA5 features a comparatively smooth decline, LEPROT does not. The variance along its gradient is too high to be deemed significant.

# left plot
# plotStsLineplot(test, variables = "LEPROT", line_color = "blue", id = "test_1") + 
#   ggplot2::geom_point()

# right plot
# plotStsLineplot(test, variables = "SHISA5", line_color = "red", id = "test_1") + 
#   ggplot2::geom_point()


# 3.2 Model fits +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The second step uses predefined models and fits them to the inferred gradients of the pattern identified as non-random in the previous step. The figure below shows the default models used by SPATA2. They can be extended by the user for specific querries with the argument model_add.
showModels(nrow = 4) + 
  labs(x = "Distance along Trajectory [%]")


# Slot @results$model_fits contains the model fitting results. It is a data.frame where each row corresponds to a variable ~ model pair. The columns mae (mean absolute error) and rmse (root mean squared error) indicate the quality of the fit. The lower the value the better.

# only the ones with fdr < 0.05
best_fits <- sts_out@results$model_fits %>% 
  filter(variables %in% sign_df[["variables"]]) %>% 
  group_by(variables) %>% 
  slice_min(mae, n = 1)

# pull all the genes
best_fits_all <-
  sts_out@results$model_fits %>% 
  filter(variables %in% sts_out@results$significance[["variables"]]) %>% 
  group_by(variables) %>% 
  slice_min(mae, n = 1)

best_fits %>%
  filter(variables %in% c("SPP1"))

best_fits_all %>%
  filter(variables %in% c("SPP1","NUPR1","LUC7L3"))

# save the table for all the features with the assigned best model
best_fits_all %>%
  write_tsv("../out/table/best_fits_all.tsv")

# The following code chunk extracts the genes that followed each model best.
best_fits_by_model <- 
  group_by(best_fits, models) %>% 
  slice_min(mae, n = 1)
# filter(rmse < 0.2) # threshold suggestions for root mean squared error

group_by(best_fits, models)

best_fits_by_model

plotSurfaceComparison(
  object = test, 
  color_by = c(best_fits_by_model[["variables"]]), 
  outline = TRUE, 
  display_image = FALSE, 
  pt_clrsp = "Reds 3", 
  nrow = 2
) + 
  trajectory_add_on


plotStsLineplot(
  object = test, ,
  variables = c(best_fits_by_model[["variables"]]), 
  id = "test_1", 
  line_color = "red", 
  nrow = 2
)


plotStsLineplot(
  object = test, ,
  variables = c("SPP1","NUPR1","LUC7L3","CST3","C5orf63","NUPR1","MBP"), 
  id = "test_1", 
  line_color = "red", 
  nrow = 2
)

plot <- plotStsLineplot(
  object = test, ,
  variables = c("HMGB1","IGFBP5","IGFBP7","IL6ST","JUN","PECAM1","SERPINE1","IGFBP6","NUPR1","SPP1","CHIT1","IGKC","SYN1","SV2A"), 
  id = "test_1", 
  line_color = "red", 
  nrow = 4
)

plotSurfaceComparison(
  object = test, 
  color_by = c("HMGB1","IGFBP5","IGFBP7","IL6ST","JUN","PECAM1","SERPINE1","IGFBP6","NUPR1","SPP1","CHIT1","IGKC","SYN1","SV2A"), 
  outline = TRUE, 
  display_image = FALSE, 
  pt_clrsp = "Reds 3", 
  nrow = 4
) + 
  trajectory_add_on

plot$data %>%
  ggplot(aes(x = dist, y = values)) + geom_point()+facet_wrap(~variables,scales = "free")+geom_smooth(n=100)+theme_cowplot()+theme(strip.background = element_blank())

test %>%
  saveRDS("../out/object/SPATA2_object_V01.rds")

plotStsLineplot(test, variables ="nCount_Spatial", id = "test_1")

# test add geneset --------------------------------------------------------
# example genes of an hypothetical gene signature (gene set)
# read in the signature file
list_sig <- readRDS("../../data/senescence_pathways.rds")

genes <- list_sig$senmayo$Genes %>% unique()

# both, opt1 and opt2, have the same effect (addSignature() just allows to specify the
# assay of interest, which is fixed to 'gene' for addGeneSet()).
# opt1
test <- addSignature(test, name = "senmayo", class = "Custom", molecules = genes)

# opt2
# object <- addGeneSet(object, name = "EXCITATORY_SYNAPSE", class = "Custom", genes = genes)

gs <- "Custom_senmayo"

# visualize signature expression
plotSurface(test, color_by = gs)

# extract genes of signature
getGenes(test, signatures = gs)

# ?addGeneSet()
# ?concept_molecular_signatures
plotStsLineplot(test, variables ="Custom_senmayo", id = "test_1")
