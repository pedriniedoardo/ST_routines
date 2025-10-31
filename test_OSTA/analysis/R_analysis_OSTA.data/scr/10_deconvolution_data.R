# libraries ---------------------------------------------------------------
library(OSTA.data)

# read in the data --------------------------------------------------------
id <- "Visium_HumanBreast_Janesick"
pa <- OSTA.data_load(id)
unzip(pa, exdir="/media/edo/sandiskSSD/work/training/spatial_transcriptomic/test_OSTA/data/Visium_HumanBreast_Janesick")

id <- "Chromium_HumanBreast_Janesick"
pa <- OSTA.data_load(id)
unzip(pa, exdir="/media/edo/sandiskSSD/work/training/spatial_transcriptomic/test_OSTA/data/Chromium_HumanBreast_Janesick")
