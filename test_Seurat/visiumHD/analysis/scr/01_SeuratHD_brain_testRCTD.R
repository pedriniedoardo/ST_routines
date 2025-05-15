# AIM ---------------------------------------------------------------------
# this is just to isolate the RCTD step from the main workflow to run some tests.

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(spacexr)

# read in the data --------------------------------------------------------
# read in the main object
object <- readRDS(file="../out/object/0506final_brain.rds")

# read in the subsetting file
cortex.coordinates <- as.data.frame(read.csv("../data/vignette_seuratHD/cortex-hippocampus_coordinates.csv"))

# load in the reference scRNA-seq dataset
ref <- readRDS("../data/vignette_seuratHD/allen_scRNAseq_ref.Rds")

# wrangling ---------------------------------------------------------------
# subset the main object
cortex <- CreateSegmentation(cortex.coordinates)
# library(sf)
object[["cortex"]] <- Overlay(object[["slice1.008um"]], cortex)
cortex <- subset(object, cells=Cells(object[['cortex']]))

# dimension of the full subset
dim(cortex@assays$Spatial.008um)

#sketch the cortical subset of the Visium HD dataset
DefaultAssay(cortex) <- "Spatial.008um"
cortex <- FindVariableFeatures(cortex)
cortex <- SketchData(
  object = cortex,
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch")

# confirm the dimension of the sketched slot
dim(cortex@assays$sketch)

# standard pre-process the sketched assay
DefaultAssay(cortex) <- "sketch"
cortex <- ScaleData(cortex) %>%
  RunPCA(assay="sketch", reduction.name = "pca.cortex.sketch", verbose = T) %>%
  FindNeighbors(reduction = "pca.cortex.sketch", dims = 1:50) %>%
  RunUMAP(reduction = "pca.cortex.sketch", reduction.name = "umap.cortex.sketch", return.model = T, dims = 1:50, verbose = T)

# save the object
saveRDS(cortex,"../out/object/01_cortex_testRCTD_small.rds")

# extract the count from the reference
Idents(ref) <- "subclass_label"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$subclass_label)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

# extract the count from the query
counts_hd <- cortex[["sketch"]]$counts
cortex_cells_hd <- colnames(cortex[["sketch"]])
coords <- GetTissueCoordinates(cortex)[cortex_cells_hd,1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 28)
saveRDS(RCTD,"../out/object/01_RCTD01_cortex_small.rds")

# this is quite an long process
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
saveRDS(RCTD,"../out/object/01_RCTD02_cortex_small.rds")

# notice that the result of the deconvolution is the size of the sketch
dim(RCTD@results$results_df)

# add results back to Seurat object
cortex <- AddMetaData(cortex, metadata = RCTD@results$results_df)

#
dim(cortex@meta.data)

# notice that at this point the annotation is sparse. we have added annotation only to the cells that have been tested in the sketched subset. all the others are NA
test_meta <- cortex@meta.data %>%
  group_by(first_type) %>%
  summarise(n = n()) %>%
  mutate(test = case_when(is.na(first_type)~NA,
                          T~"tested")) %>%
  group_by(test) %>%
  summarise(n = sum(n))

#
test_meta
dim(RCTD@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
cortex$first_type <- as.character(cortex$first_type)
cortex$first_type[is.na(cortex$first_type)] <- "Unknown"
cortex <- ProjectData(
  object = cortex,
  assay = "Spatial.008um",
  full.reduction = "pca.cortex",
  sketched.assay = "sketch",
  sketched.reduction = "pca.cortex.sketch",
  umap.model = "umap.cortex.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)

cortex@meta.data

# now after the Project data, we have access to a column full_first_type

# plot the data -----------------------------------------------------------
DefaultAssay(cortex) <- "sketch"
Idents(cortex) <- "first_type"
p1 <- DimPlot(object, reduction = "umap.sketch", label=F) + ggtitle("Sketched clustering (5000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(cortex) <- "Spatial.008um"
Idents(cortex) <- "first_type"
p2 <- DimPlot(cortex, reduction = "full.umap.sketch", label=F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2




