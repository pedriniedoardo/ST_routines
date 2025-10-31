renv::install("tidyomics/tidySpatialWorkshop")

# In May 2025, the following packages should be installed from github repositories, to use the latest features. In case you have them pre installed, run the following command
renv::install(c("bioc::ggspavis",
                "stemangiola/tidySummarizedExperiment",
                "william-hutchison/tidySpatialExperiment",
                "stemangiola/tidybulk",
                "stemangiola/tidygate",
                "stemangiola/CuratedAtlasQueryR"))

renv::install("ggcorrplot")

# Then build the vignettes
renv::install("tidyomics/tidySpatialWorkshop", build_vignettes = TRUE, force=TRUE)

# To view vignette
library(tidySpatialWorkshop)
vignette("Introduction")

# -------------------------------------------------------------------------
renv::install("spatialLIBD")
renv::install("ExperimentHub")
renv::install("SpatialExperiment")
renv::install("MoleculeExperiment")
renv::install("SubcellularSpatialData")
renv::install("ExperimentHub")
renv::install("CuratedAtlasQueryR")
renv::install("tidyverse")
renv::install("scater")
renv::install("scran")
renv::install("scuttle")
renv::install("Seurat")
renv::install("SPOTlight")
renv::install("Banksy")
renv::install("hoodscanR")
renv::install("scRNAseq")

renv::install("scatterpie")
renv::install("zellkonverter")
renv::install("MangiolaLaboratory/cellNexus")
