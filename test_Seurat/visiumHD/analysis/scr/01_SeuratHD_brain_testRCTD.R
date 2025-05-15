# AIM ---------------------------------------------------------------------
# this is just to isolate the RCTD step from the main workflow to run some tests.

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(spacexr)

# read in the data --------------------------------------------------------
# read in the main object
cortex <- readRDS(file="../out/object/00_cortex.rds")

# load in the reference scRNA-seq dataset
ref <- readRDS("../data/vignette_seuratHD/allen_scRNAseq_ref.Rds")

# wrangling ---------------------------------------------------------------
# dimension of the full subset
dim(cortex@assays$Spatial.008um)

# check the metadata. make sure there is no intermediate annotation
cortex@meta.data %>%
  head()

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

# save the object to avoid rerunning the above steps
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

# check the metadata
dim(cortex@meta.data)
head(cortex@meta.data) 

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


# -------------------------------------------------------------------------
# this is a test
cortex$spot_class <- as.character(cortex$spot_class)
cortex$spot_class[is.na(cortex$spot_class)] <- "miss"
# -------------------------------------------------------------------------



# this way we have coerced all the value in the first type column to a placeholder. for some reason the column cannot have NA for the ProjectData step.
head(cortex@meta.data)
cortex@meta.data %>%
  group_by(first_type) %>%
  summarise(n = n()) %>%
  mutate(test = case_when(is.na(first_type)~NA,
                          T~"tested")) %>%
  group_by(test) %>%
  summarise(n = sum(n))

# run the projection
cortex <- ProjectData(
  object = cortex,
  assay = "Spatial.008um",
  full.reduction = "pca.cortex",
  sketched.assay = "sketch",
  sketched.reduction = "pca.cortex.sketch",
  umap.model = "umap.cortex.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type",
                 # this is a test
                 full_spot_class = "spot_class")
)

# now there is a new column called full_first_type
cortex@meta.data

# if the Unknown label is passed automatically to to all the missing cells, I would expect that all the cells that were not part of the sketching process, would have been assigned to unknown.
# use the first_class variable as filter of defining the cells that were missing from the sketch process.
# define the LUT of the sketched barcodes
LUT_sketched <- cortex@meta.data %>%
  mutate(test = case_when(is.na(first_class)~"no_sketched",
                          T~"sketched")) %>%
  select(test) %>%
  rownames_to_column("barcodes")

# define the LUT of the sketched barcodes
LUT_sketched %>%
  group_by(test) %>%
  summarise(n = n())

# count the barcodes annotation per sketching condition
LUT_annotation <- cortex@meta.data %>%
  select(full_first_type) %>%
  rownames_to_column("barcodes")

summary_table <- LUT_annotation %>%
  left_join(LUT_sketched,by = "barcodes") %>%
  group_by(test,full_first_type) %>%
  summarise(n = n(),.groups = "drop") %>%
  group_by(test) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/tot)

summary_table %>%
  mutate(full_first_type = fct_reorder(full_first_type,-prop)) %>%
  ggplot(aes(x=full_first_type,y=prop,fill=test))+geom_col(position = "dodge")+
  scale_fill_manual(values = c("sketched" = "#E69F00", "no_sketched" = "#56B4E9"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# do the same for the spot_class
LUT_annotation2 <- cortex@meta.data %>%
  select(full_spot_class) %>%
  rownames_to_column("barcodes")

summary_table2 <- LUT_annotation2 %>%
  left_join(LUT_sketched,by = "barcodes") %>%
  group_by(test,full_spot_class) %>%
  summarise(n = n(),.groups = "drop") %>%
  group_by(test) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/tot)

summary_table2 %>%
  mutate(full_spot_class = fct_reorder(full_spot_class,-prop)) %>%
  ggplot(aes(x=full_spot_class,y=prop,fill=test))+geom_col(position = "dodge")+
  scale_fill_manual(values = c("sketched" = "#E69F00", "no_sketched" = "#56B4E9"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot the data -----------------------------------------------------------
DefaultAssay(cortex) <- "sketch"
Idents(cortex) <- "first_type"
p1 <- DimPlot(cortex, reduction = "umap.sketch", label=F) + ggtitle("Sketched clustering (5000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(cortex) <- "Spatial.008um"
Idents(cortex) <- "first_type"
p2 <- DimPlot(cortex, reduction = "full.umap.sketch", label=F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2





