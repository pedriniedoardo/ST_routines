# AIM ---------------------------------------------------------------------
# sample download raw files from geoquery

# libraries ---------------------------------------------------------------
library(GEOquery)
library(tidyverse)

# wrangling ---------------------------------------------------------------
gse <- getGEO('GSE222322',GSEMatrix = F,destdir = "../data/")

folder <- names(GSMList(gse))

# create the folders to save the data
lapply(folder, function(x){
  dir.create(paste0("../data/sample_GSE222322/",x))
})

# pull the link to the file
gse@gsms$GSM6919905@header$supplementary_file_2

# sample save one file
download.file(gse@gsms$GSM6919905@header$supplementary_file_1,destfile = paste0("../data/sample_GSE222322/GSM6919905/filtered_feature_bc_matrix.h5"))

# save all the files
pmap(list(gse@gsms,names(gse@gsms)),function(x,name){
  download.file(x@header$supplementary_file_1,destfile = paste0("../data/sample_GSE222322/",name,"/filtered_feature_bc_matrix.h5"))
  # gse@gsms$GSM6919905@header$supplementary_file_1
})

file_name <- basename(gse@gsms$GSM6919905@header$supplementary_file_1) %>% 
  str_remove_all("_filtered_feature_bc_matrix.h5")

write_tsv(x = data.frame("test"),file = paste0("test/",name,"/",file_name,".txt"))

