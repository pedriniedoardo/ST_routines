# Intermediate processing {#sec-seq-intermediate-processing}

## Preamble

# Introduction ------------------------------------------------------------
# This chapter demonstrates methods for several intermediate processing steps -- normalization, feature selection, and dimensionality reduction -- that are required before downstream analysis methods can be applied.
# We use methods from the `r BiocStyle::Biocpkg("scater")` [@McCarthy2017-scater] and `r BiocStyle::Biocpkg("scran")` [@Lun2016-scran] Bioconductor packages, which were originally developed for scRNA-seq data, with the (simplified) assumption that these methods can be applied to spot-based ST data by treating spots as equivalent to cells. In addition, we discuss alternative spatially-aware methods, and provide links to later parts of the book where these methods are described in more detail.
# Following the processing steps, this chapter also includes short demonstrations for subsequent steps, including clustering and identification of marker genes. For more details on these steps, see the later analysis chapters or workflow chapters.


# Dependencies
library(scran)
library(scater)
library(ggspavis)
library(pheatmap)
library(patchwork)
library(BayesSpace)
library(SpatialExperiment)

# set seed for reproducibility
set.seed(100)

# Load data ---------------------------------------------------------------
# Here, we load 10x Genomics Visium data from one postmortem human brain tissue section from the dorsolateral prefrontal cortex (DLPFC) brain region [@Maynard2021-DLPFC], which has previously undergone quality control and filtering (see also @sec-seq-quality-control).
# <!--
#   Here, we load a 10x Genomics Visium dataset that will also be used in several of the following chapters. This dataset has previously been preprocessed using tools outside R and saved in `SpatialExperiment` format, and is available for download from the `r BiocStyle::Biocpkg("STexampleData")` package.
# 
# This dataset consists of one sample (Visium capture area) of postmortem human brain tissue from the dorsolateral prefrontal cortex (DLPFC) brain region, measured with the 10x Genomics Visium platform. The dataset is described in @Maynard2021-DLPFC.
# -->

# load data
spe <- readRDS("../../out/object/seq-spe_qc.rds")

# Normalization {#sec-seq-intermediate-processing-norm} -------------------
### Library size normalization
# A simple and fast approach for normalization is "library size normalization", which consists of calculating log-transformed normalized counts ("logcounts") using library size factors, treating each spot as equivalent to a single cell. This method is simple, fast, and generally provides a good baseline. It can be calculated using the `r BiocStyle::Biocpkg("scater")` [@McCarthy2017-scater] and `r BiocStyle::Biocpkg("scran")` [@Lun2016-scran] packages.
# [See also @sec-ind-normalization.]{.aside}

# However, library size normalization does not make use of any spatial information. In some datasets, the simplified assumption that each spot can be treated as equivalent to a single cell is not appropriate, and can create problems during analysis. (*"Library size confounds biology in spatial transcriptomics data"* [@Bhuva2024-library-size].)

# Some alternative non-spatial methods from scRNA-seq workflows (e.g. normalization by deconvolution) are also less appropriate for spot-based ST data, since spots can contain multiple cells from one or more cell types, and datasets can contain multiple samples (e.g. multiple tissue sections), which may result in sample-specific clustering.

# calculate library size factors
spe <- computeLibraryFactors(spe)

summary(sizeFactors(spe))
hist(sizeFactors(spe), breaks = 50, main = "Histogram of size factors")

# calculate logcounts and store in new assay
spe <- logNormCounts(spe)

assayNames(spe)

# Spatially-aware normalization
# `r BiocStyle::Biocpkg("SpaNorm")` [@Salim2025-SpaNorm] was developed for ST data, using a gene-wise model (e.g. negative binomial) in which variation is decomposed into a library size-dependent (technical) and -independent (biological) component.

# Feature selection -------------------------------------------------------
# Highly variable genes (HVGs)
# Identifying a set of top "highly variable genes" (HVGs) is a standard step for feature selection in many scRNA-seq workflows, which can also be used as a simple and fast baseline approach in spot-based ST data. This makes the simplified assumption that spots can be treated as equivalent to cells.

# Here, we take a standard approach for selecting HVGs; see @sec-ind-feature-selection-testing-hvgs for more details. We first remove mitochondrial genes, since these tend to be very highly expressed and are not of main biological interest.

# identify mitochondrial genes
nm <- rowData(spe)$gene_name
mt <- grepl("^MT-", nm, ignore.case = TRUE)
table(mt)

# remove them
spe <- spe[!mt, ]

# Next, apply methods from `r BiocStyle::Biocpkg("scran")`. This gives a list of top HVGs, which can be used as the input for subsequent steps. The parameter `prop` defines the proportion of HVGs to select -- for example, `prop = 0.1` returns the top 10%. (Alternatively, argument `n` can be used to selected a fixed number of top HVGs, say 1000 or 2000.)

# fit mean-variance relationship
dec <- modelGeneVar(spe)

# plot
plot(dec$mean, dec$total, xlab = "mean (logexpr)", ylab = "variance (logexpr)")
curve(metadata(dec)$trend(x), add = TRUE, col = "dodgerblue")

# select top HVGs
sel <- getTopHVGs(dec, prop = 0.1)

# number of HVGs selected
length(sel)

# Spatially variable genes (SVGs)
# Alternatively, spatially-aware methods can be used to identify a set of "spatially variable genes" (SVGs). These methods take the spatial coordinates of the measurements into account, which can generate a more biologically informative ranking of genes associated with biological structure within the tissue area. The set of SVGs may then be used either instead of or complementary to HVGs in subsequent steps.
# A number of methods have been developed to identify SVGs. For more details and examples, see @sec-ind-feature-selection-testing-svgs.

# Dimensionality reduction ------------------------------------------------
# A standard next step is to perform dimensionality reduction using principal component analysis (PCA), applied to the set of top HVGs (or SVGs). The set of top principal components (PCs) can then be used as the input for subsequent steps. See @sec-ind-dimensionality-reduction for more details.
# Note that we have set a random seed at the start of the chapter for reproducibility, since the fast implementation in `scater::runPCA()` uses a random initialization.
spe <- runPCA(spe, subset_row = sel)

# In addition, we can perform nonlinear dimensionality reduction using the UMAP algorithm, applied to the set of top PCs (default 50). The first two UMAP dimensions can be used for visualization purposes by plotting them on the x and y axes.
spe <- runUMAP(spe, dimred = "PCA")

# Clustering {#sec-seq-intermediate-processing-clustering}
# <!-- To do: consider moving sections from here onwards into later chapters (either workflows or downstream steps) - or new workflow chapter (before DLPFC) -->
# After completing the intermediate processing steps above, we next demonstrate some short examples of downstream analysis steps.
# Here, we run a spatially-aware clustering algorithm, `r BiocStyle::Biocpkg("BayesSpace")` [@Zhao2021-BayesSpace], to identify spatial domains. BayesSpace was developed for sequencing-based ST data, and takes the spatial coordinates of the measurements into account.
# For more details on clustering, see @sec-ind-clustering.

# run BayesSpace clustering
.spe <- spatialPreprocess(spe, skip.PCA = TRUE)
.spe <- spatialCluster(.spe, nrep = 1000, burn.in = 100, q = 10, d = 20)

# cluster labels
table(spe$BayesSpace <- factor(.spe$spatial.cluster))

# Visualizations ----------------------------------------------------------
# We can visualize the cluster labels in x-y space as spatial domains, alongside the manually annotated reference labels (`ground_truth`) available for this dataset.

# using plotting functions from ggspavis package
# and formatting using patchwork package
plotCoords(spe, annotate = "ground_truth") + geom_point(size = 2)
plotCoords(spe, annotate = "BayesSpace") + geom_point(size = 2) +
  plot_layout() & 
  theme(legend.key.size = unit(0, "lines")) & 
  scale_color_manual(values = unname(pals::trubetskoy()))
table(spe$BayesSpace, spe$ground_truth, useNA = "ifany")

# Inspecting the spatial distribution of the top PCs, we can observe that the main driver of expression variability is the distinction between white matter (WM) and non-WM cortical layers.
pcs <- reducedDim(spe, "PCA")
pcs <- pcs[, seq_len(4)]
lapply(colnames(pcs), \(.) {
  spe[[.]] <- pcs[, .]
  plotCoords(spe, annotate = .)
}) |> 
  wrap_plots(nrow = 1) & coord_equal() & 
  geom_point(shape = 16, stroke = 0, size = 2) & 
  scale_color_gradientn(colors = pals::jet(), n.breaks = 3) & 
  theme_void() & theme(
    plot.title = element_text(hjust = 0.5), 
    legend.key.width = unit(0.2, "lines"), 
    legend.key.height = unit(0.8, "lines"))

# Differential expression -------------------------------------------------
# Having clustered our spots to identify spatial domains, we can test for genes that are differentially expressed (DE) between spatial domains, which can be interpreted as marker genes for the spatial domains; see also @sec-ind-feature-selection-testing-degs.
# Here, we will use pairwise t-tests and specifically test for upregulation (as opposed to downregulation), i.e. expression should be higher in the cluster for which a gene is reported to be a marker. See @sec-ind-clustering for additional details.
# [Alternatively, we could also test for spatially variable genes (SVGs) within domains to identify genes with spatial expression patterns that vary independently of the spatial arrangement of domains; see @sec-ind-feature-selection-testing-svgs for details.]{.aside}

# using scran package
mgs <- findMarkers(spe, groups = spe$BayesSpace, direction = "up")
top <- lapply(mgs, \(df) rownames(df)[df$Top <= 2])
length(top <- unique(unlist(top)))

# We can visualize selected markers as a heatmap (of cluster-wise expression means):
# compute cluster-wise averages
pbs <- aggregateAcrossCells(spe, 
                            ids = spe$BayesSpace, subset.row = top, 
                            use.assay.type = "logcounts", statistics = "mean")
# use gene symbols as feature names
mtx <- t(assay(pbs))
colnames(mtx) <- rowData(pbs)$gene_name
# using pheatmap package
pheatmap(mat = mtx, scale = "column")

# Or spatial plots (of spot-level expression values):
# gene-wise spatial plots
gs <- c("MBP", "PLP1", "NRGN", "SNAP25", "NEFL", "HPCAL1")
ps <- lapply(gs, \(.) {
  plotCoords(spe, 
             annotate = ., 
             feature_names = "gene_name", 
             assay_name = "logcounts") })
# figure arrangement
wrap_plots(ps, nrow = 2) & 
  theme(legend.key.width = unit(0.4, "lines"), 
        legend.key.height = unit(0.8, "lines")) & 
  scale_color_gradientn(colors = rev(hcl.colors(9, "Rocket")))

# Save data {.unnumbered} -------------------------------------------------
# Save data object for re-use within later chapters.

colLabels(spe) <- spe$BayesSpace
saveRDS(spe, "../../out/object/seq-spe_cl.rds")
