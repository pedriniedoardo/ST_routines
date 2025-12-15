# AIM ---------------------------------------------------------------------
# sample processing of a single sample of visium using Seurat

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(tidyverse)

library(glmGamPoi)
library(Rfast2)
# library(presto)

# make sure to keep using v3 assays
options(Seurat.object.assay.version = "v3")

# Dataset -----------------------------------------------------------------
# Here, we will be using a recently released dataset of sagital mouse brain slices generated using the Visium v1 chemistry. There are two serial anterior sections, and two (matched) serial posterior sections.

# You can download the data here, and load it into Seurat using the Load10X_Spatial() function. This reads in the output of the spaceranger pipeline, and returns a Seurat object that contains both the spot-level expression data along with the associated image of the tissue slice. You can also use our SeuratData package for easy data access, as demonstrated below. After installing the dataset, you can type ?stxBrain to learn more.
# InstallData("stxBrain")
# brain <- LoadData("stxBrain", type = "anterior1")
# brain2 <- LoadData("stxBrain", type = "posterior1")
# # save the object for future uses
# saveRDS(brain,"../../data/misc/stxBrain.SeuratData_0.1.2_anterior1.rds")
# saveRDS(brain2,"../../data/misc/stxBrain.SeuratData_0.1.2_posterior1.rds")

# v3 objects
brain <- readRDS("../../data/misc/stxBrain.SeuratData_0.1.1_anterior1.rds")
brain2 <- readRDS("../../data/misc/stxBrain.SeuratData_0.1.1_posterior1.rds")

# v5 objects
# brain <- readRDS("../../data/misc/stxBrain.SeuratData_0.1.2_anterior1.rds")
# brain2 <- readRDS("../../data/misc/stxBrain.SeuratData_0.1.2_posterior1.rds")
allen_reference <- readRDS("../../data/misc/allen_cortex.rds")

# define the scale factor for the points. this is needed for older objects assay v3 running on Seurat 5.1.0
# ref: https://github.com/satijalab/seurat/issues/9049
spot_size <- brain@images$anterior1@scale.factors$fiducial/brain@images$anterior1@scale.factors$hires
SpatialFeaturePlot(brain, features = "nCount_Spatial",pt.size.factor = spot_size) + theme(legend.position = "right")

# check the scale factor of the images if using v3 with Seurat 5.1.0
# image <- brain[["anterior1"]]
# scale_factors <- ScaleFactors(image)
# 
# hires_scale_factor <- scale_factors[["hires"]]
# spot_radius <- scale_factors[["spot"]]

# if needed fix the scale factor
# if hires_scale_factor and spot_radius are equal - this should not be the case. Fortunately, there's an easy workaround:
# brain@images$anterior1@scale.factors$spot <- 89.4725

class(brain@assays$Spatial)
class(brain@assays$SCT)
class(allen_reference@assays$RNA)

# Data preprocessing ------------------------------------------------------
# The initial preprocessing steps that we perform on the spot by gene expression data are similar to a typical scRNA-seq experiment. We first need to normalize the data in order to account for variance in sequencing depth across data points. We note that the variance in molecular counts / spot can be substantial for spatial datasets, particularly if there are differences in cell density across the tissue. We see substantial heterogeneity here, which requires effective normalization.
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial",pt.size.factor = spot_size) + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# These plots demonstrate that the variance in molecular counts across spots is not just technical in nature, but also is dependent on the tissue anatomy. For example, regions of the tissue that are depleted for neurons (such as the cortical white matter), reproducibly exhibit lower molecular counts. As a result, standard approaches (such as the LogNormalize() function), which force each data point to have the same underlying ‘size’ after normalization, can be problematic.

# As an alternative, we recommend using sctransform (Hafemeister and Satija, Genome Biology 2019), which which builds regularized negative binomial models of gene expression in order to account for technical artifacts while preserving biological variance. For more details on sctransform, please see the paper here and the Seurat vignette here. sctransform normalizes the data, detects high-variance features, and stores the data in the SCT assay.
brain <- SCTransform(brain, assay = "Spatial", verbose = T)

# Dimensionality reduction, clustering, and visualization -----------------
# We can then proceed to run dimensionality reduction and clustering on the RNA expression data, using the same workflow as we use for scRNA-seq analysis.
brain <- brain %>%
  RunPCA(assay = "SCT", verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(verbose = T) %>%
  RunUMAP(reduction = "pca", dims = 1:30)

# We can then visualize the results of the clustering either in UMAP space (with DimPlot()) or overlaid on the image with SpatialDimPlot().
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3,pt.size.factor = spot_size)
p1 + p2

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,5, 8)), facet.highlight = TRUE, ncol = 3,pt.size.factor = spot_size)

# Identification of Spatially Variable Features ---------------------------
# Seurat offers two workflows to identify molecular features that correlate with spatial location within a tissue. The first is to perform differential expression based on pre-annotated anatomical regions within the tissue, which may be determined either from unsupervised clustering or prior knowledge. This strategy works will in this case, as the clusters above exhibit clear spatial restriction.
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3,pt.size.factor = spot_size)

# An alternative approach, implemented in FindSpatiallyVariables(), is to search for features exhibiting spatial patterning in the absence of pre-annotation. The default method (method = 'markvariogram), is inspired by the Trendsceek, which models spatial transcriptomics data as a mark point process and computes a ‘variogram’, which identifies genes whose expression level is dependent on their spatial location. More specifically, this process calculates gamma(r) values measuring the dependence between two spots a certain “r” distance apart. By default, we use an r-value of ‘5’ in these analyses, and only compute these values for variable genes (where variation is calculated independently of spatial location) to save time.

# We note that there are multiple methods in the literature to accomplish this task, including SpatialDE, and Splotch. We encourage interested users to explore these methods, and hope to add support for them in the near future.
brain <- FindSpatiallyVariableFeatures(brain,
                                       assay = "SCT",
                                       features = VariableFeatures(brain)[1:1000],
                                       selection.method = "moransi")

# this trigger an error
SpatiallyVariableFeatures(brain, method = "moransi")

# NA introduced in the meta table for the missing morinsi parameters on the missing features. on previous version of seurat it was a warning. on Seurat v5.1.0 it fails
brain@assays$SCT@meta.features %>%
  filter(!is.na(moransi.spatially.variable)) %>% dim()

# workaround extract the features with a moran parameter from the object
# brain_test <- brain
# 
# str(brain_test@assays$SCT)
# 
# brain_test@assays$SCT@SCTModel.list$counts@feature.attributes
# 
# brain_test@assays$SCT@meta.features <- brain_test@assays$SCT@meta.features %>%
#   mutate(moransi.spatially.variable.rank = case_when(is.na(moransi.spatially.variable) ~ nrow(brain),
#                                                      TRUE ~ moransi.spatially.variable))
# 
# brain_test@assays$SCT@meta.features <- brain_test@assays$SCT@meta.features %>%
#   mutate(moransi.spatially.variable = case_when(is.na(moransi.spatially.variable) ~ F,
#                                                      TRUE ~ moransi.spatially.variable))
# 
# brain_test@assays$SCT@meta.features <- brain_test@assays$SCT@meta.features %>%
#   mutate(MoransI_p.value = case_when(is.na(MoransI_p.value) ~ 1,
#                                      TRUE ~ MoransI_p.value))
# 
# brain_test@assays$SCT@meta.features <- brain_test@assays$SCT@meta.features %>%
#   mutate(MoransI_observed = case_when(is.na(MoransI_observed) ~ 0,
#                                      TRUE ~ MoransI_observed))
# 
# SpatiallyVariableFeatures(brain_test, method = "moransi")
top.features <- brain@assays$SCT@meta.features %>%
  filter(! is.na(moransi.spatially.variable.rank)) %>%
  top_n(6, desc(moransi.spatially.variable.rank)) %>%
  rownames()

# Now we visualize the expression of the top 6 features identified by this measure.
# top.features <- head(SpatiallyVariableFeatures(brain, method = "moransi"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1),pt.size.factor = spot_size)

# Subset out anatomical regions -------------------------------------------
# As with single-cell objects, you can subset the object to focus on a subset of data. Here, we approximately subset the frontal cortex. This process also facilitates the integration of these data with a cortical scRNA-seq dataset in the next section. First, we take a subset of clusters, and then further segment based on exact positions. After subsetting, we can visualize the cortical cells either on the full image, or a cropped image.

# this one is not working on a v3 object. with v5 there are some issue
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400 | image_imagecol < 150))
# make sure the coordinates are available
GetTissueCoordinates(cortex)

cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE, pt.size.factor = spot_size)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, label.size = 3,pt.size.factor = spot_size)
p1 + p2

# Integration with single-cell data ---------------------------------------
# At ~50um, spots from the visium assay will encompass the expression profiles of multiple cells. For the growing list of systems where scRNA-seq data is available, users may be interested to ‘deconvolute’ each of the spatial voxels to predict the underlying composition of cell types. In preparing this vignette, we tested a wide variety of decovonlution and integration methods, using a reference scRNA-seq dataset of ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, generated with the SMART-Seq2 protocol. We consistently found superior performance using integration methods (as opposed to deconvolution methods), likely because of substantially different noise models that characterize spatial and single-cell datasets, and integration methods are specifiically designed to be robust to these differences. We therefore apply the ‘anchor’-based integration workflow introduced in Seurat v3, that enables the probabilistic transfer of annotations from a reference to a query set. We therefore follow the label transfer workflow introduced here, taking advantage of sctransform normalization, but anticipate new methods to be developed to accomplish this task.

# We first load the data (download available here), pre-process the scRNA-seq reference, and then perform label transfer. The procedure outputs, for each spot, a probabilistic classification for each of the scRNA-seq derived classes. We add these predictions as a new assay in the Seurat object.

# allen_reference <- readRDS("/brahms/shared/vignette-data/allen_cortex.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
# plan("multisession", workers = 4)
# increase the RAM usage
options(future.globals.maxSize = 15 * 1000 * 1024^2)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = T) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = T) %>%
  RunPCA(verbose = FALSE)

# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

# Now we get prediction scores for each spot for each class. Of particular interest in the frontal cortex region are the laminar excitatory neurons. Here we can distinguish between distinct sequential layers of these neuronal subtypes, for example:

DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), ncol = 2, crop = TRUE,pt.size.factor = spot_size)

# Based on these prediction scores, we can also predict cell types whose location is spatially restricted. We use the same methods based on marked point processes to define spatially variable features, but use the cell type prediction scores as the “marks” rather than gene expression.
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "moransi",
                                        features = rownames(cortex), r.metric = 5, slot = "data")

# as before this will fail on Seurat v5.1.0
# top.clusters <- head(SpatiallyVariableFeatures(cortex, selection.method = "moransi"), 4)

top.clusters <- cortex@assays$predictions@meta.features %>%
  filter(! is.na(moransi.spatially.variable.rank)) %>%
  top_n(4, desc(moransi.spatially.variable.rank)) %>%
  rownames()

SpatialPlot(object = cortex, features = top.clusters, ncol = 2,pt.size.factor = spot_size)

# Finally, we show that our integrative procedure is capable of recovering the known spatial localization patterns of both neuronal and non-neuronal subsets, including laminar excitatory, layer-1 astrocytes, and the cortical grey matter.
SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT",
                                        "L6b", "Oligo"), ncol = 2, crop = FALSE, alpha = c(0.1, 1),pt.size.factor = spot_size)

# show the distribution of the score per annotation
meta_cortex <- cortex@meta.data %>%
  data.frame() %>%
  rownames_to_column("barcode")

LUT_cortex_score <- cortex@assays$predictions@data %>%
  data.frame() %>%
  rownames_to_column("cellID") %>%
  pivot_longer(names_to = "barcode",values_to = "score",-c(cellID)) %>%
  mutate(barcode = str_replace_all(barcode,pattern = "\\.",replacement = "-")) %>%
  filter(cellID != "max")

LUT_cortex_score_summary <- LUT_cortex_score %>%
  group_by(barcode) %>%
  nest() %>%
  mutate(top_cellID = map(data,function(x){
    order_cellID <- x %>%
      arrange(desc(score)) %>%
      pull(cellID)
    
    return(order_cellID[1])
  })) %>%
  unnest(top_cellID) %>%
  select(-data) %>%
  ungroup()

meta_cortex_full <- left_join(meta_cortex,
                              LUT_cortex_score %>%
                                pivot_wider(names_from = cellID,values_from = score) %>%
                                left_join(LUT_cortex_score_summary,by = "barcode"),
                              by = "barcode") %>%
  column_to_rownames("barcode")

# update the meta
cortex@meta.data <- meta_cortex_full

# plot the new annotation
Idents(cortex) <- "top_cellID"
SpatialDimPlot(cortex, crop = TRUE, label = TRUE,pt.size.factor = spot_size)

# Working with multiple slices in Seurat ----------------------------------
# This dataset of the mouse brain contains another slice corresponding to the other half of the brain. Here we read it in and perform the same initial normalization.

# brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)
# fix the scale factor
# brain2@images$posterior1@scale.factors$spot <- 89.4725
spot_size2 <- brain2@images$posterior1@scale.factors$fiducial/brain2@images$posterior1@scale.factors$hires

SpatialDimPlot(brain2, label = TRUE,pt.size.factor = spot_size2)

# In order to work with multiple slices in the same Seurat object, we provide the merge function.
brain.merge <- merge(brain, brain2)

# This then enables joint dimensional reduction and clustering on the underlying RNA expression data.
DefaultAssay(brain.merge) <- "SCT"

VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- brain.merge %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

# Finally, the data can be jointly visualized in a single UMAP plot. SpatialDimPlot() and SpatialFeaturePlot() will by default plot all slices as columns and groupings/features as rows.
DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))

SpatialDimPlot(brain.merge)

SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"),pt.size.factor = spot_size2)

# save the final objects --------------------------------------------------
saveRDS(brain,"../../out/object/Visium_mouse_brain_full_v3.rds")
saveRDS(cortex,"../../out/object/Visium_mouse_cortex_full_v3.rds")
saveRDS(brain.merge,"../../out/object/Visium_mouse_brainMerge_full_v3.rds")
