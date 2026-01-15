# AIM ---------------------------------------------------------------------
# quick testing of nnSVG for spatially variable genes

# libraries ---------------------------------------------------------------
library(SpatialExperiment)
library(STexampleData)
library(scran)
library(nnSVG)
library(ggplot2)

# read in the data --------------------------------------------------------
# load example dataset from STexampleData package
spe <- Visium_humanDLPFC()
dim(spe)

# preprocessing steps -----------------------------------------------------
# keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)

# filter low-expressed and mitochondrial genes
# using default filtering parameters
spe <- filter_genes(spe)
dim(spe)

# for a faster runtime reduce the size of the dataset 
set.seed(123)
n <- 100
sub <- spe[, sample(ncol(spe), n)]

# using library size factors
sub <- computeLibraryFactors(sub)
sub <- logNormCounts(sub)
assayNames(sub)

# select small set of random genes and several known SVGs for faster runtime in this example
set.seed(123)
ix_random <- sample(seq_len(nrow(sub)), 10)

known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
ix_known <- which(rowData(sub)$gene_name %in% known_genes)

ix <- c(ix_known, ix_random)

sub <- sub[ix, ]
dim(sub)

# run nnSVG
# set seed for reproducibility
set.seed(123)

# using a single thread in this example
start_time <- Sys.time()

sub <- nnSVG(sub,
             assay_name = "logcounts",
             n_threads = 1,
             verbose = F)

end_time <- Sys.time()

# Calculate difference
duration <- end_time - start_time
print(duration)


# using a single thread in this example
start_time <- Sys.time()

sub <- nnSVG(sub,
             assay_name = "logcounts",
             n_threads = 8,
             verbose = F)

end_time <- Sys.time()

# Calculate difference
duration <- end_time - start_time
print(duration)