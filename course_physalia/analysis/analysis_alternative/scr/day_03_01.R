# Session 3 – Spatial analyses of imaging data
# In this session we will learn the basics of imaging-derived spatial transcriptomic data. We will learn how to visualise, manipulate and analyse single molecule data.
# We will maintain the use of `tidyomics` that we learned in `Session 2`. The programming style, in contrast of `Session 1` will make use of the `|>` (pipe) operator.

## Experimental technologies
# Here we briefly describe some of the major technologies. This section is contributed by Dr Luciano Martellotto.

### Xenium
# The Xenium platform from 10x Genomics translates the principles of padlock-probe chemistry and rolling-circle amplification into a fully automated, imaging-based workflow that delivers subcellular maps of RNA within tissue sections. Each sample is secured in a proprietary glass cassette that holds up to several serial sections in a rigid fiducial frame, ensuring precise registration between fluorescent images and subsequent histological stains. Once loaded, the instrument carries out probe hybridisation and ligation steps in situ, converting each target transcript into a circular DNA molecule that serves as the template for highly localised rolling-circle amplification.

# As amplification proceeds, each RNA-derived padlock probe generates a dense bundle of amplified DNA—known as an amplification product dot—at the original site of the transcript. These dots are then illuminated in successive rounds of fluorescent detection and cleavage, producing a unique spatial barcode for each gene target. By capturing high-resolution images after every cycle, Xenium reads out the identity and localisation of hundreds of transcripts at true subcellular resolution (often down to 200 nm), while preserving tissue morphology throughout the experiment.

# Following completion of the imaging sequence, sections can undergo immunofluorescence or standard haematoxylin & eosin staining without loss of register, enabling seamless integration of protein, RNA and anatomical data. The accompanying Xenium Explorer software performs automated cell segmentation—typically using DAPI-stained nuclei—and assigns each amplification dot to its host cell, yielding single-cell, spatially resolved gene expression matrices ready for downstream analysis.

### CoxMx

# The CosMx Single Molecular Imager, commercialised by NanoString, brings truly single-molecule mapping of RNA—and even proteins—into intact tissue sections with subcellular precision. In this approach, each target transcript is first recognised by a bespoke in situ hybridisation probe carrying a unique “readout domain” of oligonucleotides. These readout domains each host multiple photocleavable sites. During the experiment, a series of up to sixteen fluorescent reporter sets are sequentially hybridised to the readout domains, imaged in high-resolution three-dimensional fields, then cleaved away and washed off. The presence or absence of fluorescence across the sixteen cycles produces a binary barcode that identifies each individual molecule, while the precise x, y and z coordinates captured by the microscope pin down its location within or between cells.

# CosMx is compatible with both fresh-frozen and formalin-fixed paraffin-embedded samples and can scan as much as 300 mm² of tissue by stitching together up to 384 fields of view. Beyond RNA, the same cyclic detection chemistry can be applied to barcoded antibodies for multiplexed protein measurement—in some panels up to 68 different targets in the same section. After imaging, the AtoMx software platform integrates DAPI or membrane-marker images to segment individual cells and assign each molecular dot to its host cell, yielding spatially resolved count matrices at single-cell, subcellular resolution.

### MERFISH
# The Vizgen MERSCOPE system brings the power of MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridisation) into a user-friendly, end-to-end instrument for high-plex, subcellular mapping of RNA within intact tissue sections. In practice, a fixed tissue slice—whether fresh-frozen or formalin-fixed paraffin-embedded—is first permeabilised and then incubated with a library of “encoding” probes. Each encoding probe carries a short, target‐specific sequence that binds to its RNA of interest, followed by a readout sequence that serves as the scaffold for later detection.

# Once the probes are hybridised, the MERSCOPE instrument carries out a series of imaging rounds that reveal the unique binary barcode of each transcript. In each cycle, a set of fluorescent “readout” probes is flowed over the sample and allowed to bind to their complementary readout sequences. A high-resolution image is captured, the fluorophores are chemically cleaved and washed away, and the next readout cycle begins. By repeating this process—typically over a dozen or more rounds—each RNA species is assigned a distinct pattern of “on” and “off” signals across the images. Sophisticated error-robust encoding ensures that even if a spot is missed or a signal fluctuates, the correct gene identity can still be recovered.

# Following image acquisition, the built-in MERSCOPE Visualizer software aligns the hundreds of image stacks, decodes each fluorescent barcode into a specific RNA identity, and maps thousands of individual molecules back to their precise positions within cells. Researchers can then overlay the decoded RNA map with standard histology or immunofluorescence stains, revealing how gene expression patterns relate to tissue morphology.

### 1. Working with imaging-based data in Bioconductor with MoleculeExperiment

# Tidyverse library(tidyverse)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(purrr)
library(glue) # sprintf
library(stringr)
library(forcats)
library(tibble)

# Plotting
library(colorspace)
library(dittoSeq)
library(ggspavis)
library(RColorBrewer)
library(ggspavis)

# Analysis
library(scuttle)
library(scater)
library(scran)

# Data download
library(ExperimentHub)
library(SubcellularSpatialData)

# Tidyomics
library(tidySingleCellExperiment)
library(tidySummarizedExperiment)
library(tidySpatialExperiment)

# Niche analysis
library(hoodscanR)
library(scico)

#### SubcellularSpatialData
# This [data package](https://bioconductor.org/packages/release/data/experiment/html/SubcellularSpatialData.html) contains annotated datasets localized at the sub-cellular level from the STOmics, Xenium, and CosMx platforms, as analyzed in the publication by [Bhuva et al., 2025](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03241-7). It includes raw transcript detections and provides functions to convert these into `SpatialExperiment` objects.
# The data in this workshop we will be analyzing is the Xenium Mouse Brain dataset. The dataset has 3 serial sections of fresh frozen mouse brain. Raw transcript level data is provided with region annotations for each detection.
# The data is stored in the `ExperimentHub` package, and can be downloaded using the following queries:

# To avoid error for SPE loading 
# https://support.bioconductor.org/p/9161859/#9161863
setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))

eh = ExperimentHub()
query(eh, "SubcellularSpatialData")

# Brain Mouse data
tx = eh[["EH8230"]]
tx |> filter(sample_id=="Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs") |> nrow()
# 62,744,602

#### An overview of the data
# The data is however very large and thus we will work with a small subset of the data.
tx_small =  tx[sample(seq_len(nrow(tx)), size = nrow(tx)/500),]


# -------------------------------------------------------------------------
# However, since the data is very large, for the convenience of the workshop, we will directly download the small subset of the data.

# To avoid error for SPE loading 
# https://support.bioconductor.org/p/9161859/#9161863
setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))

options(timeout = max(300, getOption("timeout")))
tx_small_file = tempfile() 
utils:: download.file("https://zenodo.org/records/11213118/files/tx_small.rda?download=1", destfile = tx_small_file)
load(tx_small_file)
# -------------------------------------------------------------------------

tx_small = tx_small |> as_tibble()
tx_small

# Let's preview the object. The data is contained in a simple data frame.
tx_small |> 
  head()

# This dataset have been annotated for regions. Here we plot the regions in the sample. We can appreciate how, even subsampling the data 1 in 500, we still have a vast amount of data to visualise.
tx_small |>
    ggplot(aes(x, y, colour = region)) +
    geom_point(pch = ".") +
    facet_wrap(~sample_id, ncol = 2) +
    coord_fixed() +
    theme_minimal() +
    theme(legend.position = "none")

# Let's have a look how many regions have been annotated
tx_small |> 
  distinct(region)

# From this large dataset, we select a small region for illustrative purposes. 
tx_small_region =
  tx |>
  filter(x |> between(3700, 4200), y |> between(5000, 5500))

# how many genes are expressed per cell
tx_small_region %>%
  filter(cell %in% c(4426)) %>%
  group_by(cell,gene) %>%
  summarise(tot_gene = sum(counts))

tx_small_region %>%
  filter(cell %in% c(4426)) %>%
  group_by(cell,gene) %>%
  summarise()
  
tx_small_region %>%
  filter(cell %in% c(4426)) %>%
  group_by(cell) %>%
  summarise(tot_count = sum(counts))

tx_small_region %>%
  filter(cell %in% c(4426)) %>%
  filter(gene %in% c("Acsbg1"))

tx_small_region %>%
  group_by(cell) %>%
  summarise(n_count = n()) %>%
  ggplot(aes(x=n_count))+geom_histogram()+scale_x_continuous(trans = "log1p",breaks = c(0,1,5,10,20,50,100,500,1000,5000,10000))+theme_bw()

tx_small_region %>%
  group_by(cell,gene) %>%
  summarise(n_gene = n()) %>%
  ggplot(aes(x=n_gene))+geom_histogram()+scale_x_continuous(trans = "log1p")+theme_bw()

tx_small_region %>%
  group_by(cell,gene) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

tx_small_region %>%
  group_by(cell) %>%
  summarise(n_gene = n(),
            tot_count = sum(counts)) %>%
  ggplot(aes(x=tot_count))+geom_histogram()+scale_x_continuous(trans = "log1p",breaks = c(0,1,5,10,20,50,100,500,1000,5000,10000))

# -------------------------------------------------------------------------
# If you do not have tx loaded from before, load the pre-saved data:
tx_small_region_file = tempfile() 
utils::download.file("https://zenodo.org/records/11213155/files/tx_small_region.rda?download=1", destfile = tx_small_region_file)
load(tx_small_region_file)
# -------------------------------------------------------------------------

### 2. MoleculeExperiment
# The R package MoleculeExperiment includes functions to create and manipulate objects from the newly introduced MoleculeExperiment class, designed for analysing molecule-based spatial transcriptomics data from platforms such as Xenium by 10X, CosMx SMI by Nanostring, and Merscope by Vizgen, among others.

# `MoleculeExperiment` class uses cell boundary information instead of cell identifiers. And thus we won't use `MoleculeExperiment` directly. However, as it is an important part of bioconductor we briefly introduce this package.

# We show how we would import our table of probe location into a  `MoleculeExperiment`. For this section, we will go through the example code given in the [vignette material](https://www.bioconductor.org/packages/release/bioc/vignettes/MoleculeExperiment/inst/doc/MoleculeExperiment.html).
library(MoleculeExperiment)

repoDir = system.file("extdata", package = "MoleculeExperiment")
repoDir = paste0(repoDir, "/xenium_V1_FF_Mouse_Brain")

me = readXenium(repoDir, keepCols = "essential")
me

# In this object, besides the single molecule location, we have cell segmentation boundaries. We can use these boundaries to understand subcellular localisation of molecules and to aggregate molecules in cells.
ggplot_me() +
  geom_polygon_me(me, assayName = "cell", fill = "#F8DE7E", color="grey") +
  geom_point_me(me) +
  # zoom in to selected patch area
  coord_cartesian(
    xlim = c(4900, 4919.98),
    ylim = c(6400.02, 6420)
  )

# In this object we don't only have the cell segmentation but the nucleus segmentation as well. 
boundaries(me, "nucleus") = readBoundaries(
  dataDir = repoDir,
  pattern = "nucleus_boundaries.csv",
  segmentIDCol = "cell_id",
  xCol = "vertex_x",
  yCol = "vertex_y",
  keepCols = "essential",
  boundariesAssay = "nucleus",
  scaleFactorVector = 1
)

boundaries(me, "cell")
showMolecules(me)

bds_colours = setNames(
  c("#aa0000ff", "#ffaaffff"),
  c("Region 1", "Region 2")
)

ggplot_me() +
  # add cell segments and colour by cell id
  geom_polygon_me(me, byFill = "segment_id", colour = "black", alpha = 0.1) +
  # add molecule points and colour by feature name
  geom_point_me(me, byColour = "feature_id", size = 0.1) +
  # add nuclei segments and colour the border with red
  geom_polygon_me(me, assayName = "nucleus", fill = NA, colour = "red") +
  # zoom in to selected patch area
  coord_cartesian(xlim = c(4900, 4919.98), ylim = c(6400.02, 6420))

# remove the object 
rm(me)
gc()

# -------------------------------------------------------------------------
# `MoleculeExperiment` also has functions such as `dataframeToMEList()` and then `MoleculeExperiment()` where we can  organise our large data frame containing single molecules into a more efficient `MoleculeExperiment` object.
library(MoleculeExperiment)

tx_small_me = 
  tx_small |> 
  select(sample_id, gene, x, y) |> 
  dataframeToMEList(
    dfType = "molecules",
    assayName = "detected",
    sampleCol = "sample_id",
    factorCol = "gene",
    xCol = "x",
    yCol = "y"
  ) |> 
  MoleculeExperiment()

tx_small_me

# Here, we can appreciate the difference in size between the redundant data frame 
tx_small |> 
  object.size() |> 
  format(units = "auto")

# and the `MoleculeExperiment`.
tx_small_me |> 
  object.size() |> 
  format(units = "auto")

rm(tx_small)
rm(tx_small_me)
gc()

#### A preview of a zoomed in section of the tissue
# Now let's try to visualise just a small section, `tx_small_region`, that we downloaded earlier. You can appreciate, single molecules are coloured by cell. You can also appreciate the difference in density between regions. An aspect to note, is that not all probes are within cells. This depends on the segmentation process.
brewer.pal(7, "Set1")

tx_small_region |>
  filter(!is.na(cell)) |> 
  slice_sample(prop = 0.3) |> 
  ggplot(aes(x, y, colour = factor(cell))) +
  geom_point(shape=".") +
  
  facet_wrap(~sample_id, ncol = 2) +
  scale_color_manual(values = sample(colorRampPalette(brewer.pal(8, "Set2"))(1800))) +
  coord_fixed() +
  theme_minimal() +
  theme(legend.position = "none")

# Let's have a look to the probes that have not being unassigned to cells.
tx_small_region |>
  filter(is.na(cell)) |> 
  ggplot(aes(x, y, colour = factor(cell))) +
  geom_point(shape=".") +
  facet_wrap(~sample_id, ncol = 2) +
  scale_color_manual(values = sample(colorRampPalette(brewer.pal(8, "Set2"))(1800))) +
  coord_fixed() +
  theme_minimal() +
  theme(legend.position = "none")


# **Exercise 3.1** --------------------------------------------------------
# We want to understand how much data we are discarding, that does not have a cell identity.
# - Using base R grammar calculate what is the ratio of outside-cell vs within-cell, probes
# - Reproduce the same calculation with `tidyverse` 
# - Calculate the percentage of probes are within the cytoplasm but outside the nucleus

tx_small_region %>%
  mutate(in_cell = !is.na(cell)) %>%
  group_by(in_cell) %>%
  summarise(n = n(),.groups = "drop") %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

tx_small_region %>%
  mutate(in_cell = !is.na(cell)) %>%
  filter(in_cell==T) %>%
  group_by(overlaps_nucleus) %>%
  summarise(n = n(),.groups = "drop") %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)



### 3. Aggregation and analysis
# We will convert our cell by gene count to a `SpatialExperiment`. This object stores a cell by gene matrix with relative XY coordinates.
# `SubcellularSpatialData` package has a utility function that aggregated the single molecules in cells, where these cell ID have been identified with segmentation.
tx_spe = SubcellularSpatialData::tx2spe(tx)
tx_spe = tx_spe |> mutate(in_tissue = TRUE) 

# If you do not have tx loaded from before, load the pre-saved data converted to `SpatialExperiment`:
tx_spe_file = tempfile() 
utils::download.file("https://zenodo.org/records/11213166/files/tx_spe.rda?download=1", destfile = tx_spe_file)
# load("~/Downloads/tx_spe.rda")
load(tx_spe_file)

# Keep just the annotated regions.
tx_spe = tx_spe |> filter(!is.na(region))

# Let have a look to the `SpatialExperiment`.
tx_spe

# Here we introduce the `ggspavis` package to visualize spatial transcriptomics data. This package requires a column called `in_tissue` to be present in the `SpatialExperiment` object. Here we edit our data include this column.
tx_spe = tx_spe |> mutate(in_tissue = TRUE) 

# Let's have a look to our `SpatialExperiment`.
tx_spe

# Let's have a look at how many cells have been detected for each region
tx_spe |> 
  add_count(region) |> 
  ggplot(aes(fct_reorder(region, n, .desc = TRUE))) +
  geom_bar() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, hjust=1, size = 2))

# We normalise the `SpatialExperiment` using `scater`.
tx_spe = 
  tx_spe |> 
  # Scaling and tranformation
  scater::logNormCounts() 

# We then visualise what is the relationship between variance and total expression across cells.
tx_spe |> 
  
  # Gene variance
  scran::modelGeneVar(block = tx_spe$sample_id) |> 
  
  # Reformat for plotting
  as_tibble(rownames  = "feature") |> 
  
  # Plot
  ggplot(aes(mean, total)) +
  geom_point() +
  geom_smooth(color="red")+
  xlab("Mean log-expression") + 
  ylab("Variance") +
  theme_bw()

# For further analysis, we subset the dataset to allow quicker calculations.
tx_spe_sample_1 = 
  tx_spe |>
  filter(sample_id=="1") |> 
  slice_sample(prop = 0.2)

# As we have done previously, we calculate variable informative genes, for further analyses.
genes <- !grepl(pattern = "NegControl.+|BLANK.+", x = rownames(tx_spe_sample_1))

# Get the top 2000 genes.
top.hvgs = 
  tx_spe_sample_1 |>
  
  scran::modelGeneVar(subset.row = genes) |> 
  
  # Model gene variance and select variable genes per sample
  getTopHVGs(n=200) 

top.hvgs

# The selected subset of genes can then be passed to the subset.row argument (or equivalent) in downstream steps. 
tx_spe_sample_1 =  
  tx_spe_sample_1 |> 
  fixedPCA(subset.row=top.hvgs)

# We then use the gene expression to cluster sales based on their similarity and represent these clusters in a two dimensional embeddings (UMAP)
# Louvain clustering is a popular method used in single-cell RNA sequencing (scRNA-seq) data analysis to identify groups of cells with similar gene expression profiles. This method is based on the Louvain algorithm, which is widely used for detecting community structures in large networks.

# The Louvain algorithm is designed to maximize a metric known as modularity, which measures the density of edges inside communities compared to edges outside communities.

# It operates in two phases: 
# - first, it looks for small communities by optimizing modularity locally, and
# - second it aggregates nodes belonging to the same community and repeats the process.
cluster_labels = 
  tx_spe_sample_1 |> 
  scran::clusterCells(
    use.dimred="PCA", 
    BLUSPARAM=bluster::NNGraphParam(k=20, cluster.fun="louvain")
  ) |> 
  as.character()

cluster_labels |> 
  head()

# Now we add this cluster column to our `SpatialExperiment`
tx_spe_sample_1 = 
  tx_spe_sample_1 |> 
  mutate(clusters = cluster_labels)

tx_spe_sample_1 |> select(.cell, clusters)

# As we have done before, we caclculate UMAPs for visualisation purposes.
# This step takes long time.
## Check how many
tx_spe_sample_1 = 
  tx_spe_sample_1 |>
  runUMAP() 

# Now, let's visualise the clusters in UMAP space.
tx_spe_sample_1 |> 
  plotUMAP(colour_by = "clusters") +
  scale_color_discrete(
    colorRampPalette(brewer.pal(9, "Set1"))(30)
  )

# **Exercise 3.2**  -------------------------------------------------------
# Let's try to understand the identity of these clusters performing gene marker detection.
# In the previous sections we have seen how to do gene marker selection for sequencing-based spatial data. We just have to adapt it to our current scenario.
# - Score the markers (scran::scoreMarkers or tx_spe_sample_1)
# - Filter top markers (filter mean.AUC > 0.8)
# - Focus on Cluster 1 and try to guess the cell type (subset first element in the list, copy and paste the first 5 genes, and quickly look in public resources about what cell type those gene are markers of)
# - Plot the umap colouring by the top marker of cluster 1 (plotReducedDim())

mgs2 <- scran::scoreMarkers(x = tx_spe_sample_1,groups = tx_spe_sample_1$clusters)
mgs_df2 <- lapply(names(mgs2), function(i) {
  x <- mgs2[[i]]
  
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.8, ]
  
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})

# bind the column
mgs_df2 <- purrr::reduce(mgs_df2,bind_rows)

# filter only genes from cluster 1
mgs_df2 %>%
  filter(cluster == 1)

# glutamatergic neurons
list_plot <- map(c("Hpcal1","Bhlhe22","Cpne6","Cdh4","Slc17a7"),function(x){
  plotReducedDim(tx_spe_sample_1, "UMAP", colour_by=x)
})

library(patchwork)
wrap_plots(list_plot)

# To understand whether the cell clusters explain morphology as opposed to merely cell identity, we can color cells according to annotated region. As we can see we have a lot of regions. We have more regions that cell clusters.
tx_spe_sample_1 |> 
  plotUMAP(colour_by = "region") +
  scale_color_discrete(
    brewer.pal(n = 30, name = "Set1")
  ) +
  guides(color="none")

# Let's try to understand the morphological distribution of cell clusters in space. 
# Plot ground truth in tissue map. 
tx_spe_sample_1 |> 
    ggspavis::plotSpots(annotate = "clusters") + 
    guides(color = "none")

# For comparison the annotated regions
tx_spe_sample_1 |> 
  ggspavis::plotSpots(annotate = "region") + 
      scale_color_manual(values = colorRampPalette( brewer.pal(9,"Set1") )(150) ) +
  guides(color = "none")

# **Exercise 3.3** --------------------------------------------------------
# **Spatial-aware clustering:** Apply the spatial aware clustering method BANKSY. Taking as example the code run for session 2.
#
library(Banksy)

# scale the counts, without log transformation
spatial_data <- tx_spe_sample_1
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
hvgs <- hvgs[!grepl(pattern = "NegControl.+|BLANK.+", x = hvgs)]
rm(seu_list, spatial_data_list_for_seurat)

length(hvgs)

rowData(spatial_data[head(hvgs),])[,c("gene", "genetype")]

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
colData(spe_joint)$clust_annotation  = colData(spe_joint)$clusters

spe_joint <- connectClusters(spe_joint)

scater::plotUMAP(spe_joint, colour_by = "clust_annotation") + scale_color_brewer(palette = "Paired")
ggspavis::plotSpots(spe_joint, annotate = "clust_annotation") + 
  facet_wrap(~sample_id) +
  scale_color_brewer(palette = "Paired")


# -------------------------------------------------------------------------
### 4. Neighborhood analyses
# hoodscanR [Liu et al., 2025](https://www.bioconductor.org/packages/release/bioc/vignettes/hoodscanR/inst/doc/Quick_start.html)
# [Source](https://divingintogeneticsandgenomics.com/post/neighborhood-cellular-niches-analysis-with-spatial-transcriptome-data-in-seurat-and-bioconductor/)
# 
# Algorithm:
# - Nearest cells detection by Approximate Nearest Neighbor (ANN) search algorithm
# - Calculating euclidean distance matrix between cells and their k-nearest neighbors
# - Cell-level annotations provided by users are used to construct a cell annotation matrix
# - Identify cellular neighborhoods uses the SoftMax function, enhanced by a "shape" parameter that governs the "influence radius". This measures probability of a cell type to be found in a neighbour.
# - The K-means clustering algorithm finds recurring neighbours
# In order to perform neighborhood scanning, we need to firstly identify k (in this example, k = 100) nearest cells for each cells. The searching algorithm is based on Approximate Near Neighbor (ANN) C++ library from the RANN package.
tx_spe_neighbours <- 
  tx_spe_sample_1 |> 
  readHoodData(anno_col = "clusters") |> 
  findNearCells(k = 100)

# The output of findNearCells function includes two matrix, an annotation matrix and a distance matrix.
# the $cells slot contains the the cell id of the neighbours
tx_spe_neighbours$cells[1:10, 1:5]
dim(tx_spe_neighbours$cells)
tx_spe_neighbours$distance[1:10, 1:5]

# We can then perform neighborhood analysis using the function scanHoods. This function incldue the modified softmax algorithm, aimming to genereate a matrix with the probability of each cell associating with their 100 nearest cells.
# Calculate neighbours
pm <- scanHoods(tx_spe_neighbours$distance)

# We can then merge the probabilities by the cell types of the 100 nearest cells. We get the probability distribution of each cell all each neighborhood. 
hoods <- mergeByGroup(pm, tx_spe_neighbours$cells)
hoods[1:2, 1:10]

# We plot randomly plot 50 cells to see the output of neighborhood scanning using plotHoodMat. In this plot, each value represent the probability of the each cell (each row) located in each cell type neighborhood. The rowSums of the probability maxtrix will always be 1.
hoods |> 
  as.data.frame() |>
  rownames_to_column(var = "cell") |>
  mutate(
      cell = str_replace(cell, ".*outs_(\\d+)$", "Xenium_\\1")
    ) |>
  column_to_rownames(var = "cell") |> 
  as.matrix() |>
  plotHoodMat(n = 50)

# We can then merge the neighborhood results with the `SpatialExperiment` object using `mergeHoodSpe` so that we can conduct more neighborhood-related analysis.
tx_spe_sample_1 =  tx_spe_sample_1 |> mergeHoodSpe(hoods)
tx_spe_sample_1

# We can see what are the neighborhood distributions look like in each cluster using `plotProbDist`. Here we only plot 10 clusters
tx_spe_sample_1 |> 
  plotProbDist(
    pm_cols = colnames(hoods),
    by_cluster = TRUE, 
    plot_all = TRUE, 
    show_clusters = as.character(seq(10))
    )

# The clusters can then be plot on the tissue using `plotissue`
tx_spe_sample_1 |> plotTissue(color = clusters)
