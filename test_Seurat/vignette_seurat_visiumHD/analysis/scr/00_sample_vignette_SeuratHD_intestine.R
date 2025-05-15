# libraries ---------------------------------------------------------------
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

# Unsupervised clustering:  mouse intestine -------------------------------
# We briefly demonstrate  our sketch-clustering workflow on a second Visium HD dataset, from the Mouse Small Intestine (FFPE), available for download [here](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-intestine). We identify clusters, visualize their spatial locations, and report their top gene expression markers:  

localdir <- "../data/mouse_intestine/"
object <- Load10X_Spatial(data.dir = localdir, bin.size = 8)

DefaultAssay(object) <- "Spatial.008um"
object <- NormalizeData(object) %>%
  FindVariableFeatures() %>%
  ScaleData()

object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch")

DefaultAssay(object) <- "sketch"
object <- FindVariableFeatures(object) %>%
  ScaleData() %>%
  RunPCA(assay="sketch", reduction.name = "pca.sketch") %>%
  FindNeighbors(assay="sketch", reduction = "pca.sketch", dims = 1:50) %>%
  FindClusters(cluster.name="seurat_cluster.sketched", resolution = 3) %>% 
  RunUMAP(reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched"))

Idents(object) <- "seurat_cluster.projected"
DefaultAssay(object) <- "Spatial.008um"

p1 <- DimPlot(object, reduction = "umap.sketch", label=F) + theme(legend.position = "bottom")
p2 <- SpatialDimPlot(object, label=F) + theme(legend.position = "bottom")
p1 | p2

# We visualize the location of each cluster individually:

Idents(object) <- "seurat_cluster.projected"
cells <- CellsByIdentities(object, idents=c(1,5,18,26))
p <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T) + NoLegend()
p

DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
object_subset <- subset(object, cells = Cells(object[['Spatial.008um']]), downsample=1000)

DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features=top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

saveRDS(object, file="../out/object/0506final_intestine.rds")
