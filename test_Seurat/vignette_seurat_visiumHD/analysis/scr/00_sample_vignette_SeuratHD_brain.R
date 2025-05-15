# Visium HD support in Seurat ---------------------------------------------
# We have previously released support Seurat for sequencing-based spatial transcriptomic (ST) technologies, including 10x visium and SLIDE-seq. We have now updated Seurat to be compatible with the Visium HD technology, which performs profiling at substantially higher spatial resolution than previous versions.
# 
# Users can install the Visium HD-compatible release from Github. Existing Seurat workflows for clustering, visualization, and downstream analysis have been updated to support both Visium and Visium HD data.
# 
# We note that Visium HD data is generated from spatially patterned olignocleotides labeled in 2um x 2um bins. However, since the data from this resolution is sparse, adjacent bins are pooled together to create 8um and 16um resolutions. 10x recommends the use of 8um binned data for analysis, but Seurat supports in the simultaneous loading of multiple binnings - and stores them in a single object as multiple assays.
# 
# In this vignette, we provide an overview of some of the spatial workflows that Seurat supports for analyzing Visium HD data, in particular:
#   
#   Unsupervised clustering
# Identification of spatial tissue domains
# Subsetting spatial regions
# Integration with scRNA-seq data
# Comparing the spatial localization of different cell types
# Please note that Visium HD is a new data type, and we expect to update this vignette as we test additional methods for spatial data analysis. We strongly encourage users to explore how different parameter settings affect their results, to analyze data iteratively (and in collaboration with biological experts), and to orthogonally validate unexpected or surprising biological findings.
# 
# We focus our analysis on a Visium HD dataset from the mouse brain, available to download here but also run clustering workflow on a dataset from the mouse intestine.

# packages required for Visium HD
# install.packages("hdf5r")
# install.packages("arrow")

library(Banksy)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(gridExtra)
library(pals)
library(tidyverse)
library(presto)
library(sf)

# Load Visium HD data -----------------------------------------------------
# * Visium HD mouse brain dataset is available for download [here](https://support.10xgenomics.com/spatial-gene-expression/datasets)
# * The Seurat can store multiple binnings/resolutions in different assays
# * `bin.size` parameter specifies resolutions to load (8 and 16um are loaded by default)
# * Users can switch between resolutions by [changing the assay](https://satijalab.org/seurat/articles/multimodal_vignette)

localdir <- "../data/mouse_brain/" 
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))

# Setting default assay changes between 8um and 16um binning
Assays(object)

# specify one of the two
DefaultAssay(object) <- "Spatial.008um"

# 
vln.plot <- VlnPlot(object, features = 'nCount_Spatial.008um', pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(object, features = 'nCount_Spatial.008um') + theme(legend.position = "right")

# note that many spots have very few counts, in-part 
# due to low cellular density in certain tissue regions
vln.plot | count.plot

# Normalize datasets ------------------------------------------------------
# In this vignette we use standard log-normalization for spatial data. We note that the best normalization methods for spatial data are still being developed and evaluated, and encourage users to read manuscripts from the [Phipson/Davis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03241-7) and [Fan](https://www.biorxiv.org/content/10.1101/2023.08.30.555624v2) labs to learn more about potential caveats for spatial normalization.

# normalize both 8um and 16um bins
DefaultAssay(object) <- "Spatial.008um"
object <- NormalizeData(object)

DefaultAssay(object) <- "Spatial.016um"
object <- NormalizeData(object)

# Visualize gene expression -----------------------------------------------
# * Adjusting `pt.size.factor` (set to 1.2 by default) helps to visualize molecular and histological info in this HD dataset  
# * You can also adjust the `shape`, and `stroke` (outline) parameters for visualization

# switch spatial resolution to 16um from 8um
DefaultAssay(object) <- "Spatial.016um"
p1 <- SpatialFeaturePlot(object, features = "Rorb") + ggtitle("Rorb expression (16um)")

# switch back to 8um
DefaultAssay(object) <- "Spatial.008um"
p2 <- SpatialFeaturePlot(object, features = "Hpca") + ggtitle("Hpca expression (8um)")

p1 | p2

# Unsupervised clustering -------------------------------------------------
# While the standard scRNA-seq clustering workflow can also be applied to spatial datasets - we have observed that when working with Visium HD datasets, the [Seurat v5 sketch clustering workflow](https://satijalab.org/seurat/articles/seurat5_sketch_analysis) exhibits improved performance, especially for identifying rare and spatially restricted groups.
# 
# As described in [Hao et al, Nature Biotechnology 2023](https://www.nature.com/articles/s41587-023-01767-y) and [Hie et al](https://www.sciencedirect.com/science/article/pii/S2405471219301528), sketch-based analyses aim to 'subsample' large datasets in a way that preserves rare populations. Here, we sketch the Visium HD dataset, perform clustering on the subsampled cells, and then project the cluster labels back to the full dataset.
# 
# Details of the sketching procedure and workflow are described in [Hao et al, Nature Biotechnology 2023](https://www.nature.com/articles/s41587-023-01767-y) and the [Seurat v5 sketch clustering vignette](https://satijalab.org/seurat/articles/seurat5_sketch_analysis). Since the full Visium HD dataset fits in memory, we do not use any of the on-disk capabilities of Seurat v5 in this vignette.

# note that data is already normalized
DefaultAssay(object) <- "Spatial.008um"
object <- FindVariableFeatures(object)
object <- ScaleData(object)

# object before sketching
object
# An object of class Seurat 
# 38118 features across 492460 samples within 2 assays 
# Active assay: Spatial.008um (19059 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: Spatial.016um
# 2 spatial fields of view present: slice1.008um slice1.016um

# we select 50,0000 cells and create a new 'sketch' assay

object <- SketchData(
  object = object,
  # in case of Seurat 5.3.0: https://github.com/satijalab/seurat/issues/7171
  # features = VariableFeatures(object),
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# now there is a new assay in the object
object
# An object of class Seurat 
# 57177 features across 492460 samples within 3 assays 
# Active assay: sketch (19059 features, 2000 variable features)
# 2 layers present: counts, data
# 2 other assays present: Spatial.008um, Spatial.016um
# 2 spatial fields of view present: slice1.008um slice1.016um

# switch analysis to sketched cells
DefaultAssay(object) <- "sketch"

# perform clustering workflow
object <- FindVariableFeatures(object) %>%
  ScaleData() %>%
  RunPCA(assay="sketch", reduction.name = "pca.sketch") %>%
  FindNeighbors(assay="sketch", reduction = "pca.sketch", dims = 1:50) %>%
  FindClusters(cluster.name="seurat_cluster.sketched", resolution = 3) %>%
  RunUMAP(reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

# Now we can project the cluster labels, and dimensional reductions (PCA and UMAP) that we learned from the 50,000 sketched cells - to the entire dataset, using the `ProjectData` function. 

# In the resulting object, for all cells:
#   
# * cluster labels will be stored in `object$seurat_cluster.projected`
# * Projected PCA embeddings will be stored in `object[["pca.008um"]]`
# * Projected UMAP embeddings will be stored in `object[["umap.sketch"]]`

object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# We can visualize the clustering results for the sketched cells, as well as the projected clustering results for the full dataset:

DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label=F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label=F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2

# Of course, we can now also visualize the unsupervised clusters based on their spatial location. Note that running `SpatialDimPlot(object, interactive = TRUE)`, also enables interactive visualization and exploration.
SpatialDimPlot(object, label=T, repel=T, label.size = 4)

# When there are many different clusters (some of which are spatially restricted and others are mixed), plotting the spatial location of all clusters can be challenging to interpret. We find it helpful to plot the spatial location of different clusters individually. For example, we highlight the spatial localization of a few clusters below, which happen to correspond to different cortical layers:
Idents(object) <- "seurat_cluster.projected"
cells <- CellsByIdentities(object, idents=c(0,4,32,34,35))
p <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], 
                    cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T) + NoLegend()
p

# We can also find and visualize the top gene expression markers for each cluster:

# Crete downsampled object to make visualization either
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
object_subset <- subset(object, cells = Cells(object[['Spatial.008um']]), downsample=1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = 'Spatial.008um', only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

# Identifying spatially-defined tissue domains ----------------------------
# While the previous analyses consider each bin independently, spatial data enables cells to be defined not just by their neighborhood, but also by their broader spatial context.
# 
# In [Singhal et al.](https://www.nature.com/articles/s41588-024-01664-3), the authors introduce BANKSY, Building Aggregates with a Neighborhood Kernel and Spatial Yardstick (BANKSY). BANKSY performs multiple tasks, but we find it particularly valuable for identifying and segmenting tissue domains. When performing clustering, BANKSY augments a spot's expression pattern with both the mean and the gradient of gene expression levels in a spot's broader neighborhood.
# 
# We thank the authors for enabling BANKSY to be compatible with Seurat via the [`SeuratWrappers`](https://github.com/satijalab/seurat-wrappers) framework, which requires separate installation of the BANKSY package:
# library(SeuratWrappers)
# library(Banksy)

# Before running BANKSY, there are two important model parameters that users should consider:
# * `k_geom` : Local neighborhood size. Larger values will yield larger domains
# * `lambda` : Influence of the neighborhood. Larger values yield more spatially coherent domains
# 
# The `RunBanksy` function creates a new `BANKSY` assay, which can be used for dimensional reduction and clustering:
object <- RunBanksy(object,
                    lambda = 0.8,
                    verbose=TRUE, 
                    assay = 'Spatial.008um',
                    slot = 'data',
                    features = 'variable',
                    k_geom = 50)

# object after running BANSKY
object
# An object of class Seurat 
# 61177 features across 492460 samples within 4 assays 
# Active assay: BANKSY (4000 features, 0 variable features)
# 2 layers present: data, scale.data
# 3 other assays present: Spatial.008um, Spatial.016um, sketch
# 4 dimensional reductions calculated: pca.sketch, umap.sketch, full.pca.sketch, full.umap.sketch
# 2 spatial fields of view present: slice1.008um slice1.016um

DefaultAssay(object) <- "BANKSY"
object <- RunPCA(object, assay = 'BANKSY', reduction.name = "pca.banksy", features = rownames(object), npcs = 30) %>%
  FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)

Idents(object) <- "banksy_cluster"
p <- SpatialDimPlot(object, group.by="banksy_cluster",label=T, repel=T, label.size = 4,images = "slice1.008um") 
p

# As with unsupervised clustering, we can highlight the spatial location of each tissue domain individually:
banksy_cells <- CellsByIdentities(object)
p <- SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00","grey50"),facet.highlight = T, combine=T) + NoLegend()
p

# Subset out anatomical regions  ------------------------------------------
# Users may wish to segment or subset out a restricted region for further downstream analysis. For example, here we create a coordinate-defined segmentation mask marking cortical and hippocampal regions from the entire dataset using the `CreateSegmentation` function, and then identify cells that fall into this region with the `Overlay` function.
# 
# The list of coordinates is available for download [here](https://www.dropbox.com/scl/fi/qbs3j1alq33f0qz892ub3/cortex-hippocampus_coordinates.csv?rlkey=lsxglb15jhjdrircy9lb6n0rd&dl=0), and users can identify these boundaries when exploring their own datasets using the `interactive=TRUE` argument to `SpatialDimPlot`.
cortex.coordinates <- as.data.frame(read.csv("../data/vignette_seuratHD/cortex-hippocampus_coordinates.csv"))
cortex <- CreateSegmentation(cortex.coordinates)

# library(sf)
object[["cortex"]] <- Overlay(object[["slice1.008um"]], cortex)
cortex <- subset(object, cells=Cells(object[['cortex']]))

# Integration with scRNA-seq data (deconvolution) -------------------------
# Seurat v5 also includes support for [Robust Cell Type Decomposition](https://www.nature.com/articles/s41587-021-00830-w), a computational approach to deconvolve spot-level data from spatial datasets, when provided with an scRNA-seq reference. RCTD has been shown to accurately annotate spatial data from a variety of technologies, including SLIDE-seq, Visium, and the 10x Xenium in-situ spatial platform. We observe good performance with Visium HD as well.
# 
# To run RCTD, we first install the `spacexr` package from GitHub which implements RCTD. When running RCTD, we follow the instructions from the [RCTD vignette](https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html).

library(spacexr)

# RCTD takes an scRNA-seq dataset as a reference, and a spatial dataset as a query. For a reference, we use a mouse scRNA-seq dataset from the Allen Brain Atlas, available for download [here](https://www.dropbox.com/scl/fi/r1mixf4eof2cot891n215/allen_scRNAseq_ref.Rds?rlkey=ynr6s6wu1efqsjsu3h40vitt7&dl=0). The reference scRNAs-eq dataset has been reduced to 200,000 cells (and rare cell types <25 cells have been removed).

# We use the cortex Visium HD object as the spatial query. For computational efficiency, we sketch the spatial query dataset, apply RCTD to deconvolute the 'sketched' cortical cells and annotate them, and then project these annotations to the full cortical dataset.

#sketch the cortical subset of the Visium HD dataset
DefaultAssay(cortex) <- "Spatial.008um"
cortex <- FindVariableFeatures(cortex)
cortex <- SketchData(
  object = cortex,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch")

DefaultAssay(cortex) <- "sketch"
cortex <- ScaleData(cortex) %>%
  RunPCA(assay="sketch", reduction.name = "pca.cortex.sketch", verbose = T) %>%
  FindNeighbors(reduction = "pca.cortex.sketch", dims = 1:50) %>%
  RunUMAP(reduction = "pca.cortex.sketch", reduction.name = "umap.cortex.sketch", return.model = T, dims = 1:50, verbose = T)

# load in the reference scRNA-seq dataset
ref <- readRDS("../data/vignette_seuratHD/allen_scRNAseq_ref.Rds")

# check the ref
ref
# An object of class Seurat 
# 31053 features across 199993 samples within 1 assay 
# Active assay: RNA (31053 features, 2000 variable features)
# 1 layer present: counts

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
saveRDS(RCTD,"../out/object/RCTD01_vignette_SeuratHD_cortex.rds")

# this is quite an long process
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
saveRDS(RCTD,"../out/object/RCTD02_vignette_SeuratHD_cortex.rds")

# notice that the result of the deconvolution is the size of the sketch
dim(RCTD@results$results_df)

# add results back to Seurat object
cortex <- AddMetaData(cortex, metadata = RCTD@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
cortex$first_type <- as.character(cortex$first_type)
cortex$first_type[is.na(cortex$first_type)] <- 'Unknown'
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

DefaultAssay(object) <- "Spatial.008um"

# we only ran RCTD on the cortical cells
# set labels to all other cells as "Unknown"
object[[]][, "full_first_type"] <- "Unknown"
object$full_first_type[Cells(cortex)] <- cortex$full_first_type[Cells(cortex)]
Idents(object) <- 'full_first_type'

# now we can spatially map the location of any scRNA-seq cell type
# start with Layered (starts with L), excitatory neurons in the cortex
cells <- CellsByIdentities(object)
excitatory_names <- sort(grep("^L.* CTX",names(cells),value = TRUE))
p <- SpatialDimPlot(object, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T, ncol=4)
p

# We can now look for associations between the scRNA-seq labels of individual bins, and their tissue domain identity (as assigned by BANKSY). By asking which domains the excitatory neuron cells fall in, we can rename the BANKSY clusters as neuronal layers:

plot_cell_types <- function(data, label) {
  p <- ggplot(data, aes(x = get(label), y = n, fill = full_first_type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = ifelse(n >= min_count_to_show_label, full_first_type, "")), position = position_stack(vjust = 0.5), size = 2) +
    xlab(label) +
    ylab("# of Spots") +
    ggtitle(paste0("Distribution of Cell Types across ", label)) +
    theme_minimal()
}

cell_type_banksy_counts <- object[[]] %>%
  dplyr::filter(full_first_type %in% excitatory_names) %>%
  dplyr::count(full_first_type, banksy_cluster)

min_count_to_show_label <- 20

p <- plot_cell_types(cell_type_banksy_counts, "banksy_cluster")
p

# Based on this plot, we can now assign cells (even if they are not excitatory neurons) to individual neuronal layers. 

Idents(object) <- 'banksy_cluster'
object$layer_id <- 'Unknown'
object$layer_id[WhichCells(object,idents = c(7))] <- "Layer 2/3"
object$layer_id[WhichCells(object,idents = c(15))] <- "Layer 4"
object$layer_id[WhichCells(object,idents = c(5))] <- "Layer 5"
object$layer_id[WhichCells(object,idents = c(1))] <- "Layer 6"

# Finally, we can visualize the spatial distribution of other cell types, and ask which cortical layers they fall in. For example, in contrast to excitatory neurons, inhibitory (GABAergic) interneurons in the cortex are not spatially restricted to individual layers - but they do show biases. 
# 
# In our [previous analysis of STARmap data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687398/), and [consistent with previous work](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6456269/), we found that SST and PV interneuron classes tend to be restricted to layers 4-6, while VIP and Lamp5 interneurons tend to be located in layers 2/3. These results were based on an in-situ imaging technology which captures single-cell profiles - and here we ask whether we can find the same result in Visium HD spot-based data.

# set ID to RCTD label
Idents(object) <- 'full_first_type'

# Visualize distribution of 4 interneuron subtypes
inhibitory_names <- c("Sst","Pvalb","Vip","Lamp5")
cell_ids <- CellsByIdentities(object, idents = inhibitory_names)
p <- SpatialDimPlot(object, cells.highlight = cell_ids, cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T, ncol=4)
p

# create barplot to show proportions of cell types of interest
layer_table <- table(object$full_first_type, object$layer_id)[inhibitory_names,1:4]

neuron_props <- reshape2::melt(prop.table(layer_table), margin = 1)
ggplot(neuron_props, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cell type", y = "Proportion", fill = "Layer") +
  theme_classic()

saveRDS(object, file="../out/object/0506final_brain.rds")

# We recapitulate the same findings, previously identified in in-situ imaging data, in the Visium HD dataset. This highlights that the 8um binning of Visium HD, even though it does not represent true single cell resolution, is capable of accurately localizing scRNA-seq-defined cell types, although we strongly encourage users to orthogonally validate unexpected or surprising biological findings.
