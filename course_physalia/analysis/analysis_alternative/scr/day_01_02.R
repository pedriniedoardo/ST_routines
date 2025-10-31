# 3. Downloading Example Dataset ------------------------------------------
# We'll work with data from the `spatialLIBD` package, which contains spatial transcriptomics data from human dorsolateral prefrontal cortex. This dataset provides an excellent example for learning spatial analysis techniques as it includes:
# - Multiple tissue sections
# - Known anatomical layers
# - Rich molecular information
# - High-quality imaging data

# libraries
library(SpatialExperiment)
library(ggspavis)
library(scuttle)
library(scater)
library(scran)
library(spatialLIBD)
library(ExperimentHub)

# In this section, we explore the handling and processing of spatial transcriptomics data using the `spatialLIBD` and `ExperimentHub` packages. The following R code block retrieves a specific dataset from the `ExperimentHub`, a Bioconductor project designed to manage and distribute large biological data sets. The code efficiently fetches the data, removes any existing dimensional reductions, and filters the dataset to include only selected samples. This approach is essential for analysing spatial patterns in gene expression across multiple samples, and the code below exemplifies how to manipulate these datasets in preparation for further analysis. This process is adapted from the `Banksy` package's vignette, which provides advanced methods for multi-sample spatial transcriptomics.
# Maynard and Torres et al., doi: [10.1038/s41593-020-00787-0](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8095368/), [tutorial](https://www.bioconductor.org/packages/release/data/experiment/vignettes/spatialLIBD/inst/doc/spatialLIBD.html) 

# To avoid error for SPE loading 
# https://support.bioconductor.org/p/9161859/#9161863
setClassUnion("ExpData", c("matrix", "SpatialExperiment"))

spatial_data <- 
  ExperimentHub::ExperimentHub() |> 
  spatialLIBD::fetch_data( eh = _, type = "spe")

names(libd_layer_colors) <- gsub("ayer", "", names(libd_layer_colors))

# Clear the reductions
reducedDims(spatial_data) <- NULL 

# check the total number of samples before filtering
colData(spatial_data)$sample_id %>% table()
# see there are 12 individual slices
# 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676 
# 4226   4384   4789   4634   3661   3498   4110   4015   3639   3673   3592   3460

# Select only 3 samples
spatial_data <- spatial_data[,spatial_data$sample_id %in% c("151673", "151675", "151676")]

# Display the object
spatial_data

# check the dimension now
colData(spatial_data)$sample_id %>% table()
# as expected only three samples are available now
# 151673 151675 151676 
# 3639   3592   3460

# We shows metadata for each cell, helping understand the dataset's structure.
col_data <- colData(spatial_data)
head(col_data)

# We access and display feature-related information from the dataset.
row_data <- rowData(spatial_data)
head(row_data)

# Here, we perform a preliminary examination of the assay data contained within the spatial dataset.
assay(spatial_data)[1:20, 75:100]

# 4. Data Visualisation and Manipulation ----------------------------------
# Spatial transcriptomics data requires specialized visualization approaches to understand both molecular and spatial aspects simultaneously. The `ggspavis` package provides powerful tools for this purpose:

# image data
imgData(spatial_data)

# Simple visualization of spatial data
ggspavis::plotSpots(spatial_data) + 
  facet_wrap(~sample_id)

# **Exercise 1.0** --------------------------------------------------------
# Calculate how many spots have been profiled, per sample. And calculate what is the % of the total available spots in Visium low resolution.
spatial_data@colData %>%
  data.frame() %>%
  group_by(sample_id) %>%
  summarise(n = n()) %>%
  mutate(total_spot = 4992) %>%
  mutate(prop_area = n/total_spot)

# # A tibble: 3 × 4
# sample_id     n total_spot prop_area
# <chr>     <int>      <dbl>     <dbl>
#   1 151673     3639       4992     0.729
# 2 151675     3592       4992     0.720
# 3 151676     3460       4992     0.693

# -------------------------------------------------------------------------

# We can enhance our understanding by adding layer annotations. In this dataset, layers L1-6 represent different cortical layers, while WM indicates white matter:
# Plot spots with anatomical annotations
ggspavis::plotSpots(
  spatial_data, 
  annotate = "spatialLIBD"
) + 
  facet_wrap(~sample_id)

# Explore additional visualisation features offered by the Visium platform, exposing the H&E (hematoxylin and eosin) image.
ggspavis::plotVisium(spatial_data, point_size = 0.5)

# This visualisation focuses on specific tissue features within the dataset, emphasising areas of interest.
# 1. **Mitochondrial Content**: High mitochondrial gene expression often indicates stressed or dying cells
# 2. **Library Size**: Total RNA content per spot
# 3. **Number of Detected Genes**: Diversity of gene expression per spot

ggspavis::plotVisium(
  spatial_data, 
  annotate = "spatialLIBD", 
  highlight = "in_tissue", 
  point_size =0.5
) + 
  facet_wrap(~sample_id)

# 5. Quality control and filtering ----------------------------------------
# We will use the `scater` package [McCarthy et al. 2017](https://academic.oup.com/bioinformatics/article/33/8/1179/2907823?login=true) to compute the three primary QC metrics we discussed earlWe'llUsing the scater Package for QC Metrics: We'll apply the `scater` package to compute three primary quality control metrics. We'll also use `ggspavis` for visualisation along with some custom plotting techniques.
# Previously, we visualised both on- and off-tissue spots. Moving forward, we focus on on-tissue spots for more relevant analyses. This block shows how to filter out off-tissue spots to refine the dataset.
# Source [OSTAWorkshopBioc2021](https://lmweber.org/OSTAWorkshopBioc2021/articles/Vignette03_Analysis_workflow.html)
# dim(spatial_data)


# Dataset dimensions before the filtering
dim(spatial_data)
# [1] 33538 10691

# Filtering Dataset to Retain Only On-Tissue Spots: We then refine our dataset by keeping only those spots that are on-tissue, aligning with our focus for subsequent analysis. The dimensions of the dataset after filtering are displayed to confirm the changes.

# Subset to keep only on-tissue spots
spatial_data <- spatial_data[, colData(spatial_data)$in_tissue == 1]
dim(spatial_data)

# Mitochondrial -----------------------------------------------------------
# Next, we identify mitochondrial genes, as they are indicative of cell health. Cells with high mitochondrial gene expression typically indicate poor health or dying cells, which we aim to exclude.

# Classify genes as "mitochondrial" (is_mito == TRUE) or not (is_mito == FALSE)
is_gene_mitochondrial <- grepl("(^MT-)|(^mt-)", rowData(spatial_data)$gene_name)
# plot the genes
rowData(spatial_data)$gene_name[is_gene_mitochondrial]

# After identifying mitochondrial genes, we apply quality control metrics to further clean the dataset. This involves adding per-cell QC measures and setting a threshold to exclude cells with high mitochondrial transcription.

# see the metadata before adding the metrics
spatial_data@colData

# Add per-cell QC metrics to spatial data using identified mitochondrial genes
spatial_data <- scater::addPerCellQC(
  spatial_data, 
  subsets = list(mito = is_gene_mitochondrial)
)

# Select expressed genes threshold
qc_mitochondrial_transcription <- colData(spatial_data)$subsets_mito_percent > 30

# after adding the metrics
spatial_data@colData


# Check how many spots are filtered out
table(qc_mitochondrial_transcription)

# After applying the QC metrics, it's crucial to visually assess their impact. This step involves checking the spatial pattern of the spots removed based on high mitochondrial transcription, helping us understand the distribution and quality of the remaining dataset.

# Add threshold in colData
colData(spatial_data)$qc_mitochondrial_transcription <- qc_mitochondrial_transcription

# Visualize spatial pattern of filtered spots
plotSpotQC(
  spatial_data, 
  plot_type = "spot",  
  annotate = "qc_mitochondrial_transcription"
) + 
  facet_wrap(~sample_id)

# Library Size Analysis ---------------------------------------------------
# This analysis focuses on examining the distribution of library sizes across different spots. It uses a histogram and density plot to visualise the range and commonality of library sizes in the dataset.

# see the per sample metric
data.frame(colData(spatial_data)) |> 
  ggplot(aes(x = sum)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60) +
  # facet_wrap(~sample_id)+
  geom_density() +
  scale_x_log10() +
  xlab("Library size") + 
  ylab("Density") + 
  theme_classic()

library(ggridges)
# see the metric split by sample to see if there is a batch
data.frame(colData(spatial_data)) |> 
  ggplot(aes(x = sum,y=sample_id)) +
  ggridges::geom_density_ridges(alpha = 0.5) +
  # geom_histogram(aes(y = after_stat(density)), bins = 60) +
  # facet_wrap(~sample_id)+
  # geom_density(alpha = 0.5) +
  scale_x_log10() +
  xlab("Library size") + 
  ylab("Density") + 
  theme_classic()

# EXERCISE ----------------------------------------------------------------
# Use the isOutlier() function on the mitochondrial QC, and we will ask whether you think our thresholding is more stringent than the outlier method for filtering spots.
library(scuttle)
df_QC_mito <- perCellQCMetrics(spatial_data,
                               subsets = list(mito = is_gene_mitochondrial))
# the output is a dataframe of QC fro all the spots
df_QC_mito

low_QC_mito <- isOutlier(df_QC_mito$sum, type="lower", log=TRUE)
table(low_QC_mito)

high_QC_mito <- isOutlier(df_QC_mito$sum, type="high", log=TRUE)
table(high_QC_mito)

high_QC_mito <- isOutlier(df_QC_mito$sum, type="high", log=TRUE)

## Add threshold in colData
colData(spatial_data)$qc_isOutlierMitoLow <- low_QC_mito
colData(spatial_data)$qc_isOutlierMitoHigh <- high_QC_mito

## Visualize spatial pattern of filtered spots
p1 <- plotSpotQC(
  spatial_data, 
  plot_type = "spot",  
  annotate = "qc_mitochondrial_transcription") + 
  facet_wrap(~sample_id)

p2 <- plotSpotQC(
  spatial_data, 
  plot_type = "spot",  
  annotate = "qc_isOutlierMitoLow") + 
  facet_wrap(~sample_id)

p3 <- plotSpotQC(
  spatial_data, 
  plot_type = "spot",  
  annotate = "qc_isOutlierMitoHigh") + 
  facet_wrap(~sample_id)

library(patchwork)
p1/p2/p3

# does it make any difference if we run the isOutlier function per sample rather than on the overall object

# -------------------------------------------------------------------------


#### Library Size Analysis

# This analysis focuses on examining the distribution of library sizes across different spots. It uses a histogram and density plot to visualise the range and commonality of library sizes in the dataset.

## Visualize library size distribution
data.frame(colData(spatial_data)) |> 
  ggplot(aes(x = sum)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60) +
  geom_density() +
  scale_x_log10() +
  xlab("Library size") + 
  ylab("Density") + 
  theme_classic()

# Setting Library Size Threshold: After examining the library sizes, a threshold is applied to identify spots with library sizes below 700, which are considered for potential exclusion from further analysis.

## Select library size threshold
qc_total_counts <- colData(spatial_data)$sum < 700

## Check how many spots are filtered out
table(qc_total_counts)

# Incorporating Library Size Threshold in Dataset: This step involves adding the library size threshold to the dataset's metadata and examining the spatial pattern of the spots that have been removed based on this criterion.

## Add threshold in colData
colData(spatial_data)$qc_total_counts <- qc_total_counts

## Check for putative spatial pattern of removed spots
plotSpotQC(
  spatial_data, 
  plot_type = "spot",  
  annotate = "qc_total_counts", 
) + 
  facet_wrap(~sample_id)

#### Detected genes
# This analysis examines how many genes are expressed per spot, using a histogram and density plot to visualise the distribution of gene counts across the dataset.
## Density and histogram of library sizes
data.frame(colData(spatial_data) ) |> 
  ggplot(aes(x = detected)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60) +
  geom_density() +
  scale_x_log10() +
  xlab("Number of genes with > 0 counts") + 
  ylab("Density") + 
  theme_classic()

# Setting Gene Expression Threshold: This block applies a threshold to identify spots with fewer than 500 detected genes, considering these for exclusion to ensure data quality.

## Select expressed genes threshold
qc_detected_genes <- colData(spatial_data)$detected < 500
## Check how many spots are filtered out
table(qc_detected_genes)

# Incorporating Gene Expression Threshold in Dataset: After setting the gene expression threshold, it is added to the dataset's metadata. The spatial pattern of spots removed based on this threshold is then examined.

## Add threshold in colData
colData(spatial_data)$qc_detected_genes <- qc_detected_genes

## Check for putative spatial pattern of removed spots
plotSpotQC(
  spatial_data, 
  plot_type = "spot",  
  annotate = "qc_detected_genes", 
) + 
  facet_wrap(~sample_id)

# Exploring the Relationship Between Library Size and Number of Genes Detected: This analysis explores the correlation between library size and the number of genes detected in each spot, providing insights into data quality and sequencing depth.
## Density and histogram of library sizes
data.frame(colData(spatial_data)) |> 
  ggplot(aes(sum, detected)) +
  geom_point(shape=".") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Library size") + 
  ylab("Number of genes with > 0 counts") + 
  theme_classic()

#### Combined filtering
# After applying all QC filters, this block combines them and stores the results in the dataset. The spatial patterns of all discarded spots are then reviewed to ensure comprehensive quality control.

## Store the set in the object
colData(spatial_data)$discard <- qc_total_counts | qc_detected_genes | qc_mitochondrial_transcription

## Check the spatial pattern of combined set of discarded spots
plotSpotQC(
  spatial_data, 
  plot_type = "spot",  
  annotate = "discard", 
) + 
  facet_wrap(~sample_id)

# The final step in data preprocessing involves removing all spots identified as low-quality based on the previously applied thresholds, refining the dataset for subsequent analyses.
spatial_data = spatial_data[,!colData(spatial_data)$discard ]

### 6. Dimensionality reduction
# Dimensionality reduction is essential in spatial transcriptomics due to the high-dimensional nature of the data, which includes vast gene expression profiles across various spatial locations. Techniques such as PCA (Principal Component Analysis) and UMAP (Uniform Manifold Approximation and Projection) are particularly valuable. PCA helps to reduce noise and highlight the most significant variance in the data, making it simpler to uncover underlying patterns and correlations. UMAP, ofen calculated from principal components (and not directly from features) preserves both global and local data structures, enabling more nuanced visualisations of complex cellular landscapes. Together, these methods facilitate a deeper understanding of spatial gene expression, helping to reveal biological insights such as cellular heterogeneity and tissue structure, which are crucial for both basic biological research and clinical applications.

#### Variable gene identification
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(spatial_data))
dec <- scran::modelGeneVar(spatial_data, subset.row = genes) 

# Visualisation
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)

# Get top variable genes 
dec <- scran::modelGeneVar(spatial_data, subset.row = genes, block = spatial_data$sample_id) 
hvg <- scran::getTopHVGs(dec, n = 1000)

rowData(spatial_data[head(hvg),])[,c("gene_id", "gene_name")]

# here the metadata for the genes
rowData(spatial_data)

#### PCA
# With the highly variable genes, we perform PCA to reduce dimensionality, followed by UMAP to visualise the data in a lower-dimensional space, enhancing our ability to observe clustering and patterns in the data.

spatial_data <- 
  spatial_data |> 
  scuttle::logNormCounts() |> 
  scater::runPCA(subset_row = hvg)

## Check correctness - names
reducedDimNames(spatial_data)

reducedDim(spatial_data, "PCA")[1:5, 1:5]

# As for single-cell data, we need to verify that there is not significant batch effect. If so we need to adjust for it (a.k.a. integration) before calculating principal component. Many adjustment methods to output adjusted principal components directly. 

#### UMAP
# You can appreciate that, in this case, selecting within-sample variable genes, we do not see major batch effects across samples. We see two major pixel clusters.
# We can appreciate that there are no major batch effects across samples, and we don't see grouping driven by sample_id.

set.seed(42)
spatial_data <- scater::runUMAP(spatial_data, dimred = "PCA")

scater::plotUMAP(spatial_data, colour_by = "sample_id", point_size = 0.2) 

# **Exercise 1.1** --------------------------------------------------------
# Visualise where the two macro clusters are located spatially. We will take a very pragmatic approach and get cluster label from splitting the UMAP coordinated in two (`colData()` and `reducedDim()` will help us, see above), and then visualise it with `ggspavis`.

# - Modify the `SpatialExperiment` object based on the UMAP1 dimension so to divide those 2 cluster
# - Visualise the UMAP colouring by the new labelling
# - Visualise the Visium slide colouring by the new labelling

UMAP_coord <- reducedDim(spatial_data,type = "UMAP") %>%
  data.frame() %>%
  mutate(manual_cluster = case_when(UMAP1>2~"clu_01",
                                    T~"clu_02"))


# add the metadata to the objetc
colData(spatial_data)$manual_cluster <- UMAP_coord$manual_cluster
scater::plotUMAP(spatial_data, colour_by = "manual_cluster", point_size = 0.2) 

plotVisium(spatial_data,annotate = "manual_cluster")

### 7. Clustering
# Clustering in spatial transcriptomics is crucial for understanding the intricate cellular composition of tissues. This technique groups cells/pixels based on similar gene expression profiles, enabling the identification of distinct cell types and states within a tissue's spatial context. Clustering reveals patterns of cellular organisation and differentiation, and interactions in the microenvironment.

#### Transcriptome based nearest neighbours
# First, we establish the number of nearest neighbors to use in the k-NN graph. This graph forms the basis for clustering, using the Walktrap algorithm to detect community structures that suggest natural groupings or clusters in the data. `buildSNNGraph` is from the `scan` package.

## Set number of Nearest-Neighbours (NNs)
k <- 10

## Build the k-NN graph
g_walk <- 
  spatial_data |> 
  scran::buildSNNGraph( 
    k = 10, 
    use.dimred = "PCA"
  ) |> 
  igraph::cluster_walktrap()

clus <- g_walk$membership
## Check how many
table(clus)

# Applying Clustering Labels and Visualising Results: After determining clusters, we apply these labels back to the spatial data and visualise the results using UMAP. This allows us to observe how the data clusters in a reduced dimension space, and further visualise how these clusters map onto the physical tissue sections for context.

# We can appreciate here that we get two main pixel clusters.

colLabels(spatial_data) <- factor(clus)

scater::plotUMAP(spatial_data, colour_by = "label") + scale_color_brewer(palette = "Paired")

# Those two clusters group the white matter from the rest of the layers.

## Plot in tissue map
ggspavis::plotSpots(spatial_data, annotate = "label") + 
  facet_wrap(~sample_id) +
  scale_color_brewer(palette = "Paired")

# As for comparison, we show the manually annotated regions. We can see that while the single cell style clustering catchers, the overall tissue, architecture, a lot of details are not retrieved. We clusters cannot faithfully recapitulate the tissue morphology. However, they might represent specific cell types within morphological regions.

## Plot ground truth in tissue map
ggspavis::plotSpots(spatial_data, annotate = "spatialLIBD") +
  facet_wrap(~sample_id) +
  scale_color_manual(values = libd_layer_colors)

#### Spatially-aware clustering
# To cluster spatial regions (i.e. tissue domain) rather than single-cell types, the clustering algorithms need to take spatial context into account. For example what is the transcriptional profile of the neighbouring pixels or neighbouring cells.

# BANKSY combines molecular and spatial information. BANKSY leverages the fact that a cell's state can be more fully represented by considering both its own transcriptome "nd that of its local microenvironment.This algorithm embeds cells within a combined space that incorporates their own transcriptome and that of their locell'svironment, representing both the cell state and the surrounding microenvironment.

# Overview of the algorithm 

# - \* Construct a neighborhood graph between cells in physical space (k-nearest neighbors or radius nearest neighbors). 
# - \* We use neighborhood graph to compute two matrices: 
# -- \*\* Average neighborhood expression matrix 
# -- \*\* "Azimuthal Gabor filter" matrix. It represents the transcriptomic microenvironment around each cell. It measures the gradient of gene expression in each cell's neighborhood. 
# 
# - \* These matrices are then scaled on the basis of a mixing parameter λ, which controls their relative weighting 
# - \* Concatenate these two matrices with the original gene–cell expression matrix 
# - \* Combine these three matrices by direct product
# [Singhal et al., 2025](https://www.nature.com/articles/s41588-024-01664-3)
# [Material source](https://bioconductor.org/packages/release/bioc/vignettes/Banksy/inst/doc/multi-sample.html)

library(Banksy)

# scale the counts, without log transformation
spatial_data <- spatial_data |> logNormCounts(log=FALSE, name = "normcounts")

# **Highly-variable genes**
# The Banksy documentation, suggest the use of `Seurat` for the detection of highly variable genes.

library(Seurat)

# Convert to list
spatial_data_list_for_seurat <- lapply(unique(spatial_data$sample_id), function(x) 
  spatial_data[,  spatial_data$sample_id == x  ]
)

seu_list <- lapply(spatial_data_list_for_seurat, function(x) {
    x <- as.Seurat(x, data = NULL)
    NormalizeData(x, scale.factor = 5000, normalization.method = 'RC')
})

# Compute HVGs
hvgs <- lapply(seu_list, function(x) {
    VariableFeatures(FindVariableFeatures(x, nfeatures = 2000))
})
hvgs <- Reduce(union, hvgs)
rm(seu_list, spatial_data_list_for_seurat)

length(hvgs)

rowData(spatial_data[head(hvgs),])[,c("gene_id", "gene_name")]

# We now split the data by sample, to compute the neighbourhood matrices.
# Convert to list
spatial_data_list <- lapply(unique(spatial_data$sample_id), function(x) 
  spatial_data[
    hvgs, 
    spatial_data$sample_id == x
    ]
)

spatial_data_list <- lapply(
  spatial_data_list, 
  computeBanksy, # Compute the component neighborhood matrices
  assay_name = "normcounts"
)

# Rejoin the datasets
spe_joint <- do.call(cbind, spatial_data_list)

# Here, we perform PCA using the BANKSY algorithm on the joint dataset. The group argument specifies how to treat different samples, ensuring that features are scaled separately per sample group to account for variations among them.

# Note: this step takes long time
spe_joint <- runBanksyPCA( # Run PCA on the Banskly matrix
  spe_joint, 
  lambda = 0.2, # spatial weighting parameter. Larger values (e.g. 0.8) incorporate more spatial neighborhood
  group = "sample_id", # Features belonging to the grouping variable will be z-scaled separately. 
  seed = 42
)

# Once the dimensional reduction is complete, we cluster the spots across all samples and use `connectClusters` to visually compare these new BANKSY clusters against manual annotations.

spe_joint <- clusterBanksy( # clustering on the principal components computed on the BANKSY matrix
  spe_joint, 
  lambda = 0.2, # spatial weighting parameter. Larger values (e.g. 0.8) incorporate more spatial neighborhood
  resolution = 0.7, # numeric vector specifying resolution used for clustering (louvain / leiden).
  seed = 42
)
colData(spe_joint)$clust_annotation  = colData(spe_joint)$Cluster

spe_joint <- connectClusters(spe_joint)

# As an optional step, we smooth the cluster labels for each sample independently, which can enhance the visual coherence of clusters, especially in heterogeneous biological samples.

scater::plotUMAP(spe_joint, colour_by = "clust_annotation") + scale_color_brewer(palette = "Paired")
ggspavis::plotSpots(spe_joint, annotate = "clust_annotation") + 
  facet_wrap(~sample_id) +
  scale_color_brewer(palette = "Paired")

# From SpiceMix paper [Chidester et al., 2023](https://www.nature.com/articles/s41588-022-01256-z)
spatial_data_list <- lapply(
  unique(spe_joint$sample_id), 
  function(x) 
    spe_joint[, spe_joint$sample_id == x]
)

spatial_data_list <- lapply(
  spatial_data_list, 
  smoothLabels, 
  cluster_names = "clust_M0_lam0.2_k50_res0.7",
  k = 6L, 
  verbose = FALSE
)
names(spatial_data_list) <- paste0("sample_", unique(spe_joint$sample_id))

# The raw and smoothed cluster labels are stored in the `colData` slot of each `SingleCellExperiment` or `SpatialExperiment` object.
cluster_metadata <- colData(spatial_data_list$sample_151673)[, c("clust_M0_lam0.2_k50_res0.7", "clust_M0_lam0.2_k50_res0.7_smooth")]

head(cluster_metadata)

# Using cluster comparison metrics like the adjusted Rand index (ARI) we evaluate the performance of our clustering approach. This statistical analysis helps validate the clustering results against known labels or pathologies.

# The Adjusted Rand Index (ARI) is a measure of the similarity between two data clusterings. Measures degree of overlapping between two partitions.
compareClusters(spatial_data_list$sample_151673, func = 'ARI')

# We calculate the ARI for each sample to assess the consistency and accuracy of our clustering across different samples.
ari <- sapply(spatial_data_list, function(x) as.numeric(tail(compareClusters(x, func = "ARI")[, 1], n = 1)))
ari

# Visualising Clusters and Annotations on Spatial Coordinates: We utilise the ggspavis package to visually map BANKSY clusters and manual annotations onto the spatial coordinates of the dataset, providing a comprehensive visual overview of spatial and clustering data relationships.
# Use scater:::.get_palette('tableau10medium')
library(cowplot)

pal <- c(
  "#1abc9c", "#3498db", "#9b59b6", "#e74c3c", "#34495e",
"#f39c12", "#d35400", "#7f8c8d", "#2ecc71", "#e67e22"
)

p1 <- ggspavis::plotSpots(
  do.call(cbind, spatial_data_list), 
  annotate = sprintf("%s_smooth", "clust_M0_lam0.2_k50_res0.7"), 
  pal = pal
) +
  facet_wrap(~sample_id) +
  theme(legend.position = "none") +
  labs(title = "BANKSY clusters")

p2 <- ggspavis::plotSpots(
  do.call(cbind, spatial_data_list), 
  annotate = sprintf("%s", "clust_M0_lam0.2_k50_res0.7"), 
  pal = pal
) +
  facet_wrap(~sample_id) +
  theme(legend.position = "none") +
  labs(title = "BANKSY clusters")

p3 <- ggspavis::plotSpots(spatial_data, annotate = "spatialLIBD") +
  facet_wrap(~sample_id) +
  scale_color_manual(values = libd_layer_colors) +
  theme(legend.position = "none") +
  labs(title = "spatialLIBD regions")

p1/p2/p3

# **Exercise 1.2** --------------------------------------------------------
# We have applied cluster smoothing using `smoothLabels`. How much do you think this operation has affected the cluster labels. To find out,
# -   Plot the non smoothed cluster
# -   identify the pixel that have been smoothed, and
# -   visualise them using `plotSpotQC` that we have used above.


# isolate the object
test <- do.call(cbind, spatial_data_list)

# wrangling the metadata
LUT_test <- colData(test) %>%
  data.frame() %>%
  mutate(BANSKY_concordant = clust_M0_lam0.2_k50_res0.7_smooth == clust_M0_lam0.2_k50_res0.7)
# add the meta to the object
colData(test)$BANSKY_concordant <- LUT_test$BANSKY_concordant

# plotting
p1 <- ggspavis::plotSpots(
  do.call(cbind, spatial_data_list), 
  annotate = sprintf("%s_smooth", "clust_M0_lam0.2_k50_res0.7"), 
  pal = pal
) +
  facet_wrap(~sample_id) +
  theme(legend.position = "none") +
  labs(title = "BANKSY clusters")

p2 <- ggspavis::plotSpots(
  do.call(cbind, spatial_data_list), 
  annotate = sprintf("%s", "clust_M0_lam0.2_k50_res0.7"), 
  pal = pal
) +
  facet_wrap(~sample_id) +
  theme(legend.position = "none") +
  labs(title = "BANKSY clusters")

p3 <- ggspavis::plotSpots(test,annotate = "BANSKY_concordant") +
  facet_wrap(~sample_id) +
  labs(title = "BANKSY clusters")

library(patchwork)
p1/p2/p3
  
### 8. Deconvolution of pixel-based spatial data
# One of the popular algorithms for spatial deconvolution is SPOTlight. [Elosua-Bayes et al., 2021](https://academic.oup.com/nar/article/49/9/e50/6129341).

# [Source](https://bioconductor.org/packages/devel/bioc/vignettes/SPOTlight/inst/doc/SPOTlight_kidney.html)
# 
# SPOTlight uses a seeded non-negative matrix factorization regression, initialized using cell-type marker genes and non-negative least squares.

#### Producing the reference for single-cell databases
# Here, we retrieve and prepare a single-cell RNA reference. The dataset in question, zhong-prefrontal-2018, originates from a study by Zhong et al. (2018), which offers a comprehensive single-cell transcriptomic survey of the human prefrontal cortex during development . Utilising the scRNAseq package, the dataset is fetched and subsequently processed to aggregate counts across cells sharing the same sample and cell type, thereby reducing data complexity and enhancing interpretability. Further filtering steps ensure the removal of empty columns and entries with missing cell type annotations. Finally, the logNormCounts function from the scuttle package is applied to perform log-normalisation, a crucial step for mitigating technical variability and preparing the data for accurate comparative analyses .

# Get reference
library(scRNAseq)
brain_reference <- fetchDataset("zhong-prefrontal-2018", "2023-12-22")

brain_reference <- 
  brain_reference |> 
  scuttle::aggregateAcrossCells(ids = paste(brain_reference$sample, brain_reference$cell_types, sep = "_")) 

brain_reference <- brain_reference[, brain_reference |> assay() |> colSums() > 0]
brain_reference <- brain_reference[, !brain_reference$cell_types |> is.na()]

brain_reference <- 
  brain_reference |> 
  logNormCounts()

head(colData(brain_reference))

# These are the cell types included in our reference, and the number of pseudobulk samples we have for each cell type.
table(brain_reference$cell_types)

# These are the number of samples we have for each of the three data sets.
table(brain_reference$sample)

# Now, we identify the variable genes within each dataset, to not capture technical effects, and identify the union of variable genes for further analysis.
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(brain_reference))

# Convert to list
brain_reference_list <- lapply(unique(brain_reference$dataset_id), function(x) brain_reference[, brain_reference$dataset_id == x])

dec = scran::modelGeneVar(brain_reference, subset.row = genes, block = brain_reference$sample_id)
hvg_CAQ = scran::getTopHVGs(dec, n = 1000)

hvg_CAQ = unique( unlist(hvg_CAQ))

head(hvg_CAQ)

# Initially, the code prepares the spatial data by setting gene names as row identifiers.
spatial_data_gene_name <- spatial_data
rownames(spatial_data_gene_name) <- rowData(spatial_data_gene_name)$gene_name
spatial_data_gene_name = logNormCounts(spatial_data_gene_name)

# We then identify and score relevant marker genes based on their expression and significance across different cell types. The result is a curated list of high-potential marker genes, organised and ready for deeper analysis and interpretation in the context of spatial gene expression patterns.

# This function provides a convenience wrapper for marker gene identification between groups of cells, based on running `pairwiseTTests.` 

# This function represents a simpler and more intuitive summary of the differences between the groups. We do this by realizing that the p-values for these types of comparisons are largely meaningless; individual cells are not meaningful units of experimental replication, while the groups themselves are defined from the data. Thus, by discarding the p-values, we can simplify our marker selection by focusing only on the effect sizes between groups.

mgs <- scran::scoreMarkers(
  brain_reference, 
  groups = brain_reference$cell_types,
  
  # Omit mitochondrial genes and keep all the genes in spatial
  subset.row = 
    grep("(^MT-)|(^mt-)|(\\.)|(-)", rownames(brain_reference), value=TRUE, invert=TRUE) |> 
    intersect(rownames(spatial_data_gene_name))
)

# Select the most informative markers
mgs_df <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.8, ]
  
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_df)

head(mgs_df)

# We now utilise `SPOTlight` to deconvolve tisslet'smposition from our independent mouse brain reference. The result is visualised through a scatter pie plot, which provides a graphical representation of the spatial distribution of cell types within a let'sfic sample. This visualisation aids in understanding the spatial heterogeneity.

library(SPOTlight)

res <- SPOTlight(
  x = brain_reference |> assay("logcounts"),
  y = spatial_data_gene_name |> assay("logcounts"),
  groups = brain_reference$cell_types,
  group_id = "cluster",
  mgs = mgs_df,
  hvg =  intersect(hvg_CAQ, rownames(spatial_data_gene_name)),
  weight_id = "mean.AUC",
  gene_id = "gene"
)

cell_first_sample = colData(spatial_data_gene_name)$sample_id=="151673"

plotSpatialScatterpie(
  x = spatial_data_gene_name[,cell_first_sample],
  y = res$mat[cell_first_sample,],
  cell_types = colnames(res$mat[cell_first_sample,]),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4
) 

# Let's visualise without pit'syte_cell and endothelial cells, which oclet'sa lot of the visual spectrum.
plotSpatialScatterpie(
  x = spatial_data_gene_name[,cell_first_sample],
  y = res$mat[cell_first_sample,-c(2,9)],
  cell_types = colnames(res$mat[cell_first_sample,-c(2,9)]),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4
) 

# Let's also exclude without muscle_cell, which occupy a lot of the visual spectrum.
plotSpatialScatterpie(
  x = spatial_data_gene_name[,cell_first_sample],
  y = res$mat[cell_first_sample,-c(2, 9, 5)],
  cell_types = colnames(res$mat[cell_first_sample,-c(2, 9, 5)]),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4
) 

# No, let's look at the correlation matrices to see which cell type are most often occurring rather than mutually exclusive within our data set.
plotCorrelationMatrix(res$mat)
mat_df = as.data.frame(res$mat)

# **Exercise 1.4** --------------------------------------------------------
# Rather than looking at the correlation matrix, overall, let's observe whether the correlation structure amongst cell types is consistent across samples. Do you think it's consistent or noticeably different?

# **Exercise 1.5** --------------------------------------------------------
## Exercise 1.5 (adapted to your current cell types)

# Some of the most positive correlations in the new matrix are seen between:

# - **Microglia** and **Neurons**  
# - **Astrocytes** and **Stem.cells**
# 
# > **Microglia** are the resident immune cells of the central nervous system, constantly surveying the parenchyma and clearing debris.  
# > **Neurons** are the electrically excitable cells that transmit and process information via synaptic connections.  
# > **Astrocytes** are star-shaped glia that support neuronal metabolism, regulate extracellular ions and neurotransmitter uptake.  
# > **Stem.cells** denote undifferentiated progenitors capable of self-renewal and differentiation into multiple neural lineages.
# 
# Let us now **visualise** where these pairs of cell types most co-occur in your spatial map. For **each** pair, carry out the following:
# 
# 1. **Label** any pixel where both cell types exceed 10 % abundance (i.e. > 0.1).  
# 2. **Label** any pixel where the _sum_ of their abundances exceeds 40 % (i.e. > 0.4).  
# 3. **Plot** the spatial coordinates of all pixels, **colouring** them by this new label (for example:  
#    - `0` = neither condition met  
#    - `1` = both abundances > 0.1  
#    - `2` = summed abundance > 0.4  
# 
# You should end up with two analogous visualisations:
# 
# - **Microglia + Neurons**  
# - **Astrocytes + Stem.cells**
# 
# Feel free to reuse your previous code, simply substituting the cell-type columns and updating the thresholds as above.

# -------------------------------------------------------------------------
# Bonus - Alternative reference from the Human Cell Atlas - using cellNexus
# [cellNexus](https://stemangiola.github.io/cellNexus/) is a query interface that allow the programmatic exploration and retrieval of the harmonised, curated and reannotated CELLxGENE single-cell human cell atlas. Data can be retrieved at cell, sample, or dataset levels based on filtering criteria.

# Harmonised data is stored in the ARDC Nectar Research Cloud, and most cellNexus functions interact with Nectar via web requests, so a network connection is required for most functionality.

# Mangiola et al., 2025 doi [doi.org/10.1101/2023.06.08.542671](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)

# Get reference
library(cellNexus)
library(HDF5Array)

tmp_file_path = tempfile()

brain_reference =
  
  # Query metadata across 30M cells
  get_metadata() |>
  
  # Filter your data of interest
  dplyr::filter(tissue_groups=="cerebral lobes and cortical areas", disease == "Normal") |> 
  
  # Collect pseudobulk as SummarizedExperiment
  get_pseudobulk() |> 
  
  # Normalise for Spotlight
  scuttle::logNormCounts() |> 
  
  # Save for fast reading
  HDF5Array::saveHDF5SummarizedExperiment(tmp_file_path, replace = TRUE)

library(HDF5Array)
brain_reference = HDF5Array::loadHDF5SummarizedExperiment(tmp_file_path)
my_metadata = colData(brain_reference)
head(my_metadata)

# These are the cell types included in our reference, and the number of pseudobulk samples we have for each cell type.

table(brain_reference$cell_type_harmonised)

# These are the number of samples we have for each of the three data sets.
table(brain_reference$dataset_id)

# The `collection_id` can be used to gather information on the cell database. e.g. <https://cellxgene.cziscience.com/collections/><collection_id>

table(brain_reference$collection_id)
