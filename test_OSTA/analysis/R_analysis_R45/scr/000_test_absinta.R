# AIM ---------------------------------------------------------------------
# accessory function for the Snakemake pipeline.
# generate the spe object per sample to be processed by snakeamke

# libraries ---------------------------------------------------------------
library(VisiumIO)
library(tidyverse)
library(ggspavis)
library(SpatialExperiment)
library(patchwork)
library(scater)
library(scran)
library(pheatmap)
library(spatialLIBD)
library(markdown)
library(SpotSweeper)
library(BayesSpace)
library(ComplexHeatmap)
library(scuttle)
library(SpotSweeper)
library(Banksy)

# define the inputs -------------------------------------------------------
# define the list of samples
folder <- "../../data/test/"
sample <- dir(folder) %>%
  str_subset(negate = T,pattern = ".log")

sample_short <- str_replace(sample,pattern = "ABSINTA_spinalcord_",replacement = "S")

# loop processing ---------------------------------------------------------
# loop the standard processing
# samp <- "ABSINTA_spinalcord_00_045"
# samp_short <- "S00_045"

pmap(list(sample,sample_short),function(samp,samp_short){
  # track the progress of the processing
  print(samp_short)
  
  # read in the data --------------------------------------------------------
  # define the input files
  input <- paste0(folder,samp,"/outs")
  
  # read in data from spaceranger
  vis_test <- TENxVisium(
    spacerangerOut = input,
    processing = "filtered", 
    format="h5", 
    images="hires") %>%
    import()
  
  class(vis_test)
  dim(vis_test)
  
  # save the individual object
  base::saveRDS(vis_test,paste0("../../out/object/test/spe_output/000_",samp_short,".rds"))
})
