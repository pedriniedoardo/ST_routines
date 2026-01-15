# read in a sample file
library(spatialDE)

set.seed(123)
vis_test_filter <- mockSVG(size = 30, tot_genes = 200, de_genes = 1, return_SPE = TRUE)
vis_test_filter <- computeLibraryFactors(vis_test_filter)
vis_test_filter <- logNormCounts(vis_test_filter)

vis_test_filter

# test 01 -----------------------------------------------------------------
# rerun lognormcounts on the object
test <- vis_test_filter

dim(vis_test_filter)
dim(test)

all.equal(assay(test,i = "logcounts") %>% matrix(),assay(vis_test_filter,i = "logcounts") %>% matrix())
all.equal(sizeFactors(test),sizeFactors(vis_test_filter))

# does rerunning and logNormCounts affect the results?
test <- logNormCounts(test)

all.equal(assay(test,i = "logcounts") %>% matrix(),assay(vis_test_filter,i = "logcounts") %>% matrix())
all.equal(sizeFactors(test),sizeFactors(vis_test_filter))

# rerunning on the same object do not change the result

# what if I also rerun the size factor calculation
test <- computeLibraryFactors(test)
test <- logNormCounts(test)

all.equal(assay(test,i = "logcounts") %>% matrix(),assay(vis_test_filter,i = "logcounts") %>% matrix())
all.equal(sizeFactors(test),sizeFactors(vis_test_filter))

# rerunning on the same object do not change the result

# test 02 -----------------------------------------------------------------
# what if I remove samples
set.seed(123)
n <- 100
sub <- vis_test_filter[, sample(ncol(vis_test_filter), n)]
test2 <- sub

dim(vis_test_filter)
dim(test2)

all.equal(assay(test2,i = "logcounts") %>% matrix(),assay(sub,i = "logcounts") %>% matrix())
all.equal(sizeFactors(test2),sizeFactors(sub))

# does rerunning and logNormCounts affect the results?
test2 <- logNormCounts(test2)

all.equal(assay(test2,i = "logcounts") %>% matrix(),assay(sub,i = "logcounts") %>% matrix())
all.equal(sizeFactors(test2),sizeFactors(sub))

# change the result

# removing samples will change the total number of reads per sample, therefore the size factor per sample
# changing the size factor per sample will change the normalization process

# what if I also rerun the size factor calculation?
test22 <- computeLibraryFactors(test2)
test22 <- logNormCounts(test22)

all.equal(assay(test22,i = "logcounts") %>% matrix(),assay(sub,i = "logcounts") %>% matrix())
all.equal(sizeFactors(test22),sizeFactors(sub))

all.equal(assay(test22,i = "logcounts") %>% matrix(),assay(test2,i = "logcounts") %>% matrix())
all.equal(sizeFactors(test22),sizeFactors(test2))

# change the result
# it looks like it has recognized samples have been removed and therefore has recalcualted the LibraryFactors

# test 03 -----------------------------------------------------------------
# what if I remove genes
set.seed(123)
gene_id <- sample(rownames(vis_test_filter),100,replace = F)
sub2 <- vis_test_filter[gene_id,]
test3 <- sub2

dim(vis_test_filter)
dim(test3)

all.equal(assay(test3,i = "logcounts") %>% matrix(),assay(sub2,i = "logcounts") %>% matrix())
all.equal(sizeFactors(test3),sizeFactors(sub2))

# removign genes only does not affect the logNormCounst calculation provided we don't rerun the computeLibraryFactors funciton.
# removing genes will change the total number of reads per sample, therefore the size factor per sample
# changing the size factor per sample will change the normalization process

# does rerunning and logNormCounts affect the results?
test3 <- logNormCounts(test3)

all.equal(assay(test3,i = "logcounts") %>% matrix(),assay(sub2,i = "logcounts") %>% matrix())
all.equal(sizeFactors(test3),sizeFactors(sub2))

# rerunning on the same object do not change the result

# removing genes will change the total number of reads per sample, therefore the size factor per sample
# changing the size factor per sample will change the normalization process

# what if I also rerun the size factor calculation?
test33 <- computeLibraryFactors(test3)
test33 <- logNormCounts(test33)

all.equal(assay(test33,i = "logcounts") %>% matrix(),assay(sub2,i = "logcounts") %>% matrix())
all.equal(sizeFactors(test33),sizeFactors(sub2))

all.equal(assay(test33,i = "logcounts") %>% matrix(),assay(test3,i = "logcounts") %>% matrix())
all.equal(sizeFactors(test33),sizeFactors(test3))

# change the result


# -------------------------------------------------------------------------
# This result occurs because of how size factors are handled and standardized within the logNormCounts function in the scater/Bioconductor ecosystem.
# Size factors are relative values that must be centered around 1 for the specific set of cells present in the object.

# By default, the function logNormCounts() performs a sanity check on your size factors. It ensures that the mean of the size factors across all columns (cells/spots) in the object is equal to 1.If the mean is not 1, logNormCounts internally rescales them before calculating the log-counts.

# Why removing SAMPLES changes results (Test 02)
# When you subset samples (columns), you are taking a slice of the original size factors.
# Before Subsetting: The size factors for the full 900 cells had a mean of exactly 1.
# After Subsetting (sub): You selected a random 100 cells. The size factors for these specific 100 cells will likely not have a mean of exactly 1 (it might be 0.98 or 1.02 purely by chance).
# Running logNormCounts (test2): The function detects that the mean of these 100 size factors is not 1. It re-centers them (e.g., divides them all by their new mean).Because the size factors changed (they were shifted/scaled), the resulting normalized expression values (logcounts) also changed.
# Key Takeaway: Size factors are relative to the cohort of cells. If you change the cohort, the relative scaling changes.

# Why removing GENES does not change results (Test 03)
# Size factors are a property of the columns (cells), not the rows (genes). They represent how much "total signal" a cell has compared to the average.
# After Subsetting (sub2): You removed rows, but you kept the exact same 900 cells.
# The Size Factors: Since the cells are identical to the original object, the size factors are untouched. Their mean is still exactly 1.
# Running logNormCounts (test3): The function sees that the size factors are already properly centered (mean = 1). It does not need to change them. It simply applies the normalization formula using the existing values.
# Note: In your code, you eventually ran computeLibraryFactors on the gene-subsetted data (test33). This did change the results. This is because calculating size factors from scratch relies on the total counts of genes. If you delete genes, the total counts change, and new size factors are born. But simply running logNormCounts on existing factors (your specific question) does not change anything.

# should I rerun the computeLibraryFactors and logNormCounts In the case of selecting genes only?
# No, you should NOT rerun computeLibraryFactors and logNormCounts in this case. You should skip it and use the size factors calculated on the full dataset.