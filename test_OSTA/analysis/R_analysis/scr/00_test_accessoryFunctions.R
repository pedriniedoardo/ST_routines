# -------------------------------------------------------------------------
# test h5ad_to_seurat implemented in 00_accessory_functions.R
# stardist/Sample_HDP031_C1_expanded_nuclei.h5ad
# stardist/Sample_HDP031_C1_nuclei_grouped.h5ad
                                                                                                             
source("scr/00_accessoryFunctions.R")
test <- h5ad_to_seurat("../../data/stardist/Sample_HDP031_C1_expanded_nuclei.h5ad")
test

test@meta.data

test2 <- h5ad_to_seurat("../../data/stardist/Sample_HDP031_C1_nuclei_grouped.h5ad")
test2

test2@meta.data


# -------------------------------------------------------------------------


