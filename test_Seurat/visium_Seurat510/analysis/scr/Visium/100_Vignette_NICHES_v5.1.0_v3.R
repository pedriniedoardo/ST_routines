# title: "Estimating Microenvironment from Spatial Data"
# author: "Junchen Yang & Micha Sam Brickman Raredon"
# date: "`r Sys.Date()`"

# NICHES is a toolset which transforms single-cell atlases into single-cell-signaling atlases. It is engineered to be computationally efficient and very easy to run. It interfaces directly with Seurat from Satija Lab. The cell-signaling outputs from NICHES may be analyzed with any single-cell toolset, including Seurat, Scanpy, Monocle, or others.

# Here, we show how NICHES may be used to estimate individual cellular microenvironment from spatial transcriptomic data.

# First, let's load dependencies.
library(Seurat)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(NICHES)
library(viridis)
library(tidyverse)

# Next, we load the data, perform basic pre-processing, and cluster the data so that we can visualize patterns of interest. For this vignette we will use basic Seurat clustering annotations to avoid the work of labeling celltypes, which are not necessary for this demonstration.

# InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
brain <- UpdateSeuratObject(brain) # JC: need to update seurat obj

# for Seurat v5.1.0 there is the need to update the size of the spots
spot_size <- brain@images$anterior1@scale.factors$fiducial/brain@images$anterior1@scale.factors$hires

# Normalization 
brain <- SCTransform(brain, assay = "Spatial", verbose = T)
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"),pt.size.factor = spot_size)

# Dimensional reduction with all cells
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
p1 <- DimPlot(brain, reduction = "umap",group.by = 'seurat_clusters', label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE,group.by = 'seurat_clusters', label.size = 3,pt.size.factor = spot_size)
p1 + p2

# These numeric annotations are satisfactory for the demonstration of NICHES. Any and all metadata can be carried over when NICHES is run, allowing multiple levels of clustering and sub-clustering to be leveraged downstream after micorenvironment calculation.

# Next, we will format the spatial coordinate metadata so that every cell has an explicitly labeled x and y coordinate.

## Format Spatial Coordinates and Normalize
brain@meta.data$x <- brain@images$anterior1@coordinates$row
brain@meta.data$y <- brain@images$anterior1@coordinates$col

DefaultAssay(brain) <- "Spatial"
brain <- NormalizeData(brain)

# NICHES can be run on imputed or non-imputed data. Here, we will use imputed data.

## Impute and Run NICHES
brain <- SeuratWrappers::RunALRA(brain)

NICHES_output <- RunNICHES(object = brain,
                           LR.database = "fantom5",
                           species = "mouse",
                           assay = "alra",
                           position.x = 'x',
                           position.y = 'y',
                           k = 4, 
                           cell_types = "seurat_clusters",
                           min.cells.per.ident = 0,
                           min.cells.per.gene = NULL,
                           meta.data.to.map = c('orig.ident','seurat_clusters'),
                           CellToCell = F,CellToSystem = F,SystemToCell = F,
                           CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = T)


# sample routine for a custom LR database ---------------------------------
# # is it possible to change the database?
# # LR.database	string. Default: "fantom5". Currently accepts "fantom5","omnipath", or "custom".
# ?LoadCustom
# LoadCustom()
# db_mouse <- NICHES::LoadFantom5(species = "mouse")
# str(db_mouse)
# 
# # custom_LR_database	
# # data.frame. Each row is a ligand-receptor mechanism where the first column corresponds to the source genes that express the ligands subunits (separated by '_') and the second column corresponds to the receptor genes that express the receptor subunits (separated by '_').
# db_mouse_df <- data.frame(source.subunits = db_mouse$source.subunits[,1],
#                           target.subunits = db_mouse$target.subunits[,1])
# 
# head(db_mouse_df)
# dim(db_mouse_df)
# 
# NICHES_output_test <- RunNICHES(object = brain,
#                                 LR.database = "custom",
#                                 custom_LR_database = db_mouse_df,
#                                 species = "mouse",
#                                 assay = "alra",
#                                 position.x = 'x',
#                                 position.y = 'y',
#                                 k = 4, 
#                                 cell_types = "seurat_clusters",
#                                 min.cells.per.ident = 0,
#                                 min.cells.per.gene = NULL,
#                                 meta.data.to.map = c('orig.ident','seurat_clusters'),
#                                 CellToCell = F,CellToSystem = F,SystemToCell = F,
#                                 CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = T)
# 
# niche_test <- NICHES_output_test[['NeighborhoodToCell']]
# Idents(niche_test) <- niche_test@meta.data$ReceivingType
# 
# # Scale and visualize
# niche_test <- niche_test %>%
#   ScaleData() %>%
#   FindVariableFeatures(selection.method = "disp") %>%
#   RunPCA()
# ElbowPlot(niche_test,ndims = 50)
# niche_test <- RunUMAP(niche_test,dims = 1:10)
# DimPlot(niche_test,reduction = 'umap',pt.size = 0.5,shuffle = T, label = T) +ggtitle('Cellular Microenvironment')+NoLegend()
# 
# niches_test.data <- GetAssayData(object =  niche_test[['NeighborhoodToCell']], slot = 'data')

# Which allows direct visualization of niche interactions of interest in a spatial context:

# Plot celltype specific niche signaling
# SpatialFeaturePlot(brain,
#                    features = c('Bmp2—Bmpr2','Efna1—Ephb6','Fgf1—Fgfr2'),
#                    slot = 'scale.data')
# -------------------------------------------------------------------------

# NICHES outputs a list of objects. Each object contains a certain style of cell-system signaling atlas. Above, we have only calculated a single one of interest, namely, individual cellular microenvironment. We next isolate this output and embed using UMAP to visualize the microenvironemnt of each cell.
niche <- NICHES_output[['NeighborhoodToCell']]
Idents(niche) <- niche@meta.data$ReceivingType

# Scale and visualize
niche <- ScaleData(niche)
niche <- FindVariableFeatures(niche,selection.method = "disp")
niche <- RunPCA(niche)
ElbowPlot(niche,ndims = 50)
niche <- RunUMAP(niche,dims = 1:10)
DimPlot(niche,reduction = 'umap',pt.size = 0.5,shuffle = T, label = T) +ggtitle('Cellular Microenvironment')+NoLegend()

# We can already see, from this plot, some notable overlap between the microenvironments of celltypes 1 & 7 and celltypes 6 & 3. Let's explore this more deeply by finding signaling mechanisms specific to each celltype niche, plotting some of the results in heatmap form:

# Find markers
mark <- FindAllMarkers(niche,min.pct = 0.25,only.pos = T,test.use = "roc")
GOI_niche <- mark %>% group_by(cluster) %>% top_n(5,myAUC)

DoHeatmap(niche,features = unique(GOI_niche$gene))+ 
  scale_fill_gradientn(colors = c("grey","white", "blue"))

# This confirms that celltypes 1 & 7 and 6 & 3 do indeed have some shared character.
# We can further confirm that identified celltype specific signaling mechanisms are indeed specific to tissue regions in which those cells are found, by plotting matched ligand and receptor pairs:

# Check that these make sense and print little plots
DefaultAssay(brain) <- 'alra'
p1 <- SpatialFeaturePlot(brain, crop = TRUE, features = "Fgf1",slot = "data",min.cutoff =  'q1',
                         max.cutoff = 'q99',pt.size.factor = spot_size)+ggtitle("Ligand")+theme(legend.position = "right")
p2 <- SpatialFeaturePlot(brain, crop = TRUE, features = "Fgfr2",slot = "data",min.cutoff =  'q1',
                         max.cutoff = 'q99',pt.size.factor = spot_size)+ggtitle("Receptor")+theme(legend.position = "right")

ggpubr::ggarrange(p1,p2)

# Further, and perhaps more usefully, we can map over the output from NICHES onto the original spatial object as follows:

# Add Niches output as an assay
niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
brain[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
DefaultAssay(brain) <- "NeighborhoodToCell"
brain <- ScaleData(brain)

# Which allows direct visualization of niche interactions of interest in a spatial context:

# Plot celltype specific niche signaling
SpatialFeaturePlot(brain,
                   features = c('Bmp2—Bmpr2','Efna1—Ephb6','Fgf1—Fgfr2'),
                   slot = 'scale.data',pt.size.factor = spot_size)

# add the spatial niche clustering to the main object
# dim(niche@meta.data)
# dim(brain)
# 
# brain@meta.data <- brain@meta.data %>% rownames_to_column("barcodes") %>%
#   left_join(niche@meta.data %>%
#               select(barcodes = ReceivingCell,niches_cluster = seurat_clusters,ReceivingType),by = "barcodes") %>%
#   column_to_rownames("barcodes")
# 
# identical(brain@meta.data$seurat_clusters,brain@meta.data$ReceivingType)
# identical(brain@meta.data$seurat_clusters,brain@meta.data$niches_cluster)
# 
# p3 <- SpatialDimPlot(brain, label = TRUE,group.by = 'niches_cluster', label.size = 3,pt.size.factor = spot_size)
# p4 <- SpatialDimPlot(brain, label = TRUE,group.by = 'seurat_clusters', label.size = 3,pt.size.factor = spot_size)
# 
# p3 + p4
# 
# Idents(brain) <- 'niches_cluster'
# p32 <- SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain), facet.highlight = TRUE,pt.size.factor = spot_size)+plot_annotation('Niche Clusters')
# 
# Idents(brain) <- 'seurat_clusters'
# p42 <- SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain), facet.highlight = TRUE,pt.size.factor = spot_size)+plot_annotation('seurat clusters')
# 
# wrap_plots(list(p32,p42))
