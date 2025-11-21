# libraries ---------------------------------------------------------------
library(OSTA.data)

# read in the data --------------------------------------------------------
# check the location fo the cache
tools::R_user_dir("BiocFileCache", which = "cache")

bfc <- BiocFileCache::BiocFileCache()
BiocFileCache::bfcinfo()
bfc@cache

# run the regular process
id <- "Visium_HumanBreast_Janesick"
pa <- OSTA.data_load(id)
unzip(pa, exdir="/media/edo/sandiskSSD/work/training/spatial_transcriptomic/test_OSTA/data/Visium_HumanBreast_Janesick")

id <- "Chromium_HumanBreast_Janesick"
pa <- OSTA.data_load(id)
unzip(pa, exdir="/media/edo/sandiskSSD/work/training/spatial_transcriptomic/test_OSTA/data/Chromium_HumanBreast_Janesick")


# -------------------------------------------------------------------------
# in case of issues with zip, use the following

# location fo the cache before changing it
tools::R_user_dir("BiocFileCache", which = "cache")
bfc <- BiocFileCache()
bfc@cache

# decouple the OSTA.data_load() call to avoid the zip issue
id <- "Visium_HumanBreast_Janesick"
bfc <- BiocFileCache()
bfcinfo()

url <- "https://osf.io/5n4q3"

# verify 'id' with informative error if not
id <- match.arg(id, OSTA.data_list(url))

# check if already cached
q <- bfcquery(bfc, id)

# retrieve, zip, cache & return path
no <- osfr::osf_retrieve_node(url)
df <- osfr::osf_ls_files(no, id)

# dir.create(td <- tempfile())
osfr::osf_download(df, path = "../../data/Visium_HumanBreast_Janesick", recurse=TRUE)

# zip(fnm <- paste0(id, ".zip"), dir())
zip_filename <- paste0("../../data/Visium_HumanBreast_Janesick.zip") # Define the name
# file_id <- dir("../../data/") %>% str_subset(paste0(df$name,collapse = "|"))
files_to_zip <- paste0("../../data/Visium_HumanBreast_Janesick/")              # Get the list of files

# this trigger the issue
# utils::zip(zipfile = zip_filename, files = files_to_zip)
zip::zip(zipfile = zip_filename, files = files_to_zip)

# add to BiocCace
bfcadd(bfc, zip_filename, fpath=file.path(".", zip_filename))
# -------------------------------------------------------------------------


