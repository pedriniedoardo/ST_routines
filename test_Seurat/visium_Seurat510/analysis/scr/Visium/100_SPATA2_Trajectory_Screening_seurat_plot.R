# AIM ---------------------------------------------------------------------
# manual draw the regions where to calculate the gradients

# libraries ---------------------------------------------------------------
library(Seurat)
library(SPATA2)
library(tidyverse)
library(cowplot)
library(patchwork)

# read in the object ------------------------------------------------------
list_obj <- readRDS("../../out/object/100_list_brain_SPATA2.rds")

# read in the list of signatures ------------------------------------------
list_sig <- readRDS("../../data/senescence_pathways.rds")

# apply the same processing to all the slides
# sl <- list_obj[[1]]
list_SPATA2_out <- pmap(list(list_obj,names(list_obj)),function(sl,nm){
  
  # check the progress
  print(nm)
  
  # this is a wrapper around SPARK::sparkx(), to identify features with non-random spatial expression pattern.
  test <- runSPARKX(object = sl)
  
  # add gene set
  # example genes of an hypothetical gene signature (gene set)
  # read in the signature file
  genes <- list_sig$senmayo$Genes %>% unique()
  
  # both, opt1 and opt2, have the same effect (addSignature() just allows to specify the
  # assay of interest, which is fixed to 'gene' for addGeneSet()).
  # opt1
  test <- addSignature(test, name = "senmayo", class = "Custom", molecules = genes)
  
  # keep genes with a sparkx pvalue of 0.01 or lower
  spark_df <- getSparkxGeneDf(object = test, threshold_pval = 0.01)
  
  # show results
  # spark_df
  
  # name of the trajectory
  id <- "test"
  
  # pull the gene ids. `getSparkxGenes()` would work, too
  input_genes <- spark_df[["genes"]]
  
  # the function spatialTrajectoryScreening, screens the sample for numeric variables that follow specific expression changes along the course of the spatial trajectory.
  # we can use the genes ids or any other numeric variable in the medatata
  # note: the results are NOT stored in the SPATA2 object but in a separate object
  # in this case focus on a variable in the metadata
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
  # sts_out_gene@results$significance
  # sts_out_gene@results$model_fits
  # sts_out_gene@results$model_fits %>%
  #   group_by(variables) %>% 
  #   slice_min(mae, n = 1)
  
  # check if SPP1 is there
  sts_out_gene@results$significance %>% 
    filter(variables == "SPP1")
  
  sign_df <- 
    sts_out_gene@results$significance %>% 
    filter(fdr < 0.05)
  
  # show significance data.frame
  # sign_df
  
  # sts_out_gene@results$model_fits %>%
  #   filter(variables %in% c("SPP1"))
  
  # only the ones with fdr < 0.05
  best_fits <- sts_out_gene@results$model_fits %>% 
    filter(variables %in% sign_df[["variables"]]) %>% 
    group_by(variables) %>% 
    slice_min(mae, n = 1) %>%
    ungroup()
  
  # pull all the genes
  best_fits_all <-
    sts_out_gene@results$model_fits %>% 
    filter(variables %in% sts_out_gene@results$significance[["variables"]]) %>% 
    group_by(variables) %>% 
    slice_min(mae, n = 1) %>%
    ungroup()
  
  # # save the table for all the features with the assigned best model
  # best_fits_all %>%
  #   write_tsv("../out/table/best_fits_all.tsv")
  
  # The following code chunk extracts the genes that followed each model best.
  best_fits_by_model <- best_fits %>%
    group_by(models) %>% 
    slice_min(mae, n = 1)
  
  list_out <- list(slide = test, 
                   sts_out_gene = sts_out_gene, 
                   best_fits = best_fits, 
                   best_fits_all = best_fits_all, 
                   best_fits_by_model = best_fits_by_model)
  
  return(list_out)
})

saveRDS(list_SPATA2_out,"../../out/object/100_list_brain_SPATA2_processed.rds")
# list_SPATA2_out <- readRDS("../../out/object/100_list_brain_SPATA2_processed.rds")

# list_test <- list_SPATA2_out$V01
# plotStsLineplot(list_test$slide, variables ="Custom_senmayo", id = "test")+ggtitle("V01")

# plot all the senmayo traces
list_plot_senmayo <- pmap(list(list_SPATA2_out,names(list_SPATA2_out)),function(sl,nm){
  
  # check the progress
  print(nm)
  
  # plot the traces
  plot <- plotStsLineplot(sl$slide, variables ="Custom_senmayo", id = "test")+ggtitle(nm)
  
  return(plot)
  
})

wrap_plots(list_plot_senmayo)
ggsave("../../out/plot/100_senmayo_traces.pdf",width = 12, height = 12)

# plot a specific set of genes --------------------------------------------

# list_SPATA2_out$V01$sts_out_gene
# list_SPATA2_out$V01$slide
# list_SPATA2_out$V01$best_fits_all
# list_SPATA2_out$V01$sts_out_gene@results$significance[["variables"]]


# best_fits_all <-
#   sts_out@results$model_fits %>% 
#   filter(variables %in% sts_out@results$significance[["variables"]]) %>% 
#   group_by(variables) %>% 
#   slice_min(mae, n = 1)

list_SPATA2_out$V01$best_fits_all %>%
  filter(variables %in% c("NUPR1","HMGB1","IGFBP5","Custom_senmayo","SV2A","CHIT1","MBP","GFAP"))

plotStsLineplot(
  object = list_SPATA2_out$V01$slide,
  variables = c("NUPR1","HMGB1","IGFBP5","Custom_senmayo","SV2A","CHIT1","MBP","GFAP"), 
  id = "test", 
  line_color = "red", 
  nrow = 2
)
ggsave("../../out/plot/100_short_panel.pdf",width = 12, height = 6)
