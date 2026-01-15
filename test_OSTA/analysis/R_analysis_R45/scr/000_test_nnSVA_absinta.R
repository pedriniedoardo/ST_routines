# AIM ---------------------------------------------------------------------
# test the nnSVN package for the selection of Spatially Variable Genes

# renv integration --------------------------------------------------------

# to load the packages
# source(".Rprofile")

# in the config specify the following
# renv_library_path: "renv/library/linux-rocky-9.5/R-4.5/x86_64-conda-linux-gnu"

# libraries ---------------------------------------------------------------
library(VisiumIO)
library(tidyverse)
library(ggspavis)
library(SpatialExperiment)
library(patchwork)
library(scater)
library(scran)
library(pheatmap)
library(spatialLIBD)
library(markdown)
library(SpotSweeper)
library(BayesSpace)
library(ComplexHeatmap)
library(scuttle)
library(SpotSweeper)
library(Banksy)
library(SpaNorm)
library(Voyager)
library(nnSVG)

# Snakemake integation ----------------------------------------------------
# define the input
vis_input_id <- "../../out/object/test/spe_output/000_S2014_072.rds"
# vis_input_id <- snakemake@input$spe_input
message("input object: ", vis_input_id)

# define the output
# save the filtered dataset
vis_output_filtered_id <- "../../out/object/test/spe_output/000_S2014_072_SpaNorm_filtered.rds"
# vis_output_filtered_id <- snakemake@output$spe_output_filtered
message("output object: ", vis_output_filtered_id)

# save the unfiltered dataset
vis_output_unfiltered_id <- "../../out/object/test/spe_output/000_S2014_072_SpaNorm_unfiltered.rds"
# vis_output_unfiltered_id <- snakemake@output$spe_output_unfiltered
message("output object: ", vis_output_unfiltered_id)

# read in the data --------------------------------------------------------
# read in the raw spe object generated in the step before
vis_test <- readRDS(vis_input_id)

# QC ----------------------------------------------------------------------

# calculate QC metrics ----------------------------------------------------
# We calculate quality control (QC) metrics using the `r BiocStyle::Biocpkg("scater")` package [@McCarthy2017], and apply simple global thresholding-based QC methods to identify any low-quality spots, as described in @sec-seq-qc-global.

# For more details on QC methods, including more advanced QC approaches, see @sec-seq-quality-control.
# subset to keep only spots over tissue
vis_test <- vis_test[, colData(vis_test)$in_tissue == 1]

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(vis_test)$Symbol)

# Calculate QC metrics using `r BiocStyle::Biocpkg("scater")` [@McCarthy2017].
# calculate per-spot QC metrics and store in colData
vis_test <- addPerCellQC(vis_test, subsets = list(mito = is_mito))

# add the log sum metrics
vis_test$log_sum <- log(vis_test$sum)

# manual filters ----------------------------------------------------------
# define some thresholds fro removal of low QC spots
# make the call and add it to the metadata
vis_test$qc_ManualGeneLow <- vis_test$detected < 500
vis_test$qc_ManualCountLow <- vis_test$sum < 600
vis_test$qc_ManulaMitoHigh <- vis_test$subsets_mito_percent > 30

# combined set of identified spots for manual approach
vis_test$discard_manual <- vis_test$qc_ManualCountLow | vis_test$qc_ManualGeneLow | vis_test$qc_ManulaMitoHigh

# automatic filters -------------------------------------------------------
# test automatic outliers detection
# Use the isOutlier() function on the mitochondrial QC, and we will ask whether you think our thresholding is more stringent than the outlier method for filtering spots.
df_QC_mito <- perCellQCMetrics(vis_test,
                               subsets = list(mito = is_mito))

# make the call
high_QC_mito <- isOutlier(df_QC_mito$subsets_mito_percent, type="high", log=TRUE)
low_QC_count <- isOutlier(vis_test$sum, type="low", log=TRUE)
low_QC_gene <- isOutlier(vis_test$detected, type="low", log=TRUE)

# Add threshold in colData
vis_test$qc_isOutlierCountLow <- low_QC_count
vis_test$qc_isOutlierGeneLow <- low_QC_gene
vis_test$qc_isOutlierMtHigh <- high_QC_mito

# combined set of identified spots for isOutlier approach
vis_test$discard_isOutlier <- vis_test$qc_isOutlierCountLow | vis_test$qc_isOutlierGeneLow | vis_test$qc_isOutlierMtHigh

# Local outlier detection -------------------------------------------------
# make the call and add it to the metadata
vis_test <- localOutliers(vis_test, metric = "sum", direction = "lower", log = TRUE)
vis_test <- localOutliers(vis_test, metric = "detected", direction = "lower", log = TRUE)
vis_test <- localOutliers(vis_test, metric = "subsets_mito_percent",direction = "higher", log = FALSE)

# combined set of identified spots for local outlier approach
vis_test$discard_LocalOutlier <- vis_test$sum_outliers | vis_test$detected_outliers | vis_test$subsets_mito_percent_outliers

# artifact detection ------------------------------------------------------
# fill the routine for the artifact detection


# filter the slide --------------------------------------------------------
# pick one/none filtering

# remove combined set of low-quality spots
# for simplicity, in this automatic processing
vis_test_filter <- vis_test[, !vis_test$discard_LocalOutlier]

# remove features with all 0 counts
vis_test_filter <- vis_test_filter[rowSums(counts(vis_test_filter))>0, ]

# Normalization -----------------------------------------------------------
# calculate library size factors
vis_test_filter <- computeLibraryFactors(vis_test_filter)

# calculate logcounts
vis_test_filter <- logNormCounts(vis_test_filter)

# Feature selection (HVGs) ------------------------------------------------
# potentially remove uninteresting features here

# fit mean-variance relationship, decomposing variance into 
# technical and biological components
dec <- modelGeneVar(vis_test_filter)

# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)

# filter low-expressed using default filtering parameters, in this case keep mito genes we know are important to define a cell type
# it is fine to create another object as we are going to use just the ranking for the genes
spe_nnSVG <- filter_genes(vis_test_filter,filter_mito = F)

# removing genes only does not affect the logNormCounst calculation provided we don't rerun the computeLibraryFactors funciton.
# removing genes will change the total number of reads per sample, therefore the size factor per sample
# changing the size factor per sample will change the normalization process
# I think there is no reason to rerun the computeLibraryFactors after removal of genes.

# calculate logcounts
# spe_nnSVG <- logNormCounts(spe_nnSVG)

# in this case, instead of the default HVG use the selection of the top SVG for the downstream analysis
# try to use the same range of genes
set.seed(2144)
test_SVG <- nnSVG(spe_nnSVG,
                  assay_name = "logcounts",
                  n_threads = 20,
                  verbose = F)

# extract gene-level results
res_nnSVG <- rowData(test_SVG)
saveRDS(res_nnSVG,"../../out/object/test/000_res_nnSVG_S2014_072.rds")

# -------------------------------------------------------------------------
# all.equal(sizeFactors(vis_test_filter),sizeFactors(spe_nnSVG))
# 
# id_gene <- rownames(assay(spe_nnSVG,i = "logcounts"))
# 
# all.equal(assay(vis_test_filter,i = "logcounts")[id_gene,],
#           assay(spe_nnSVG,i = "logcounts")[id_gene,])
