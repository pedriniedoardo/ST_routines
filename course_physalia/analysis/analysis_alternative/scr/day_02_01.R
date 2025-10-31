# Session 2: Tidying spatial data

## Introduction to tidyomics

# `tidyomics` represents a significant advancement in bioinformatics analysis by bridging the gap between Bioconductor and the tidyverse ecosystem. This integration provides several key benefits:
  
# 1. **Unified Analysis Framework**: Combines the power of Bioconductor's specialized biological data structures with tidyverse's intuitive data manipulation
# 2. **Maintained Compatibility**: Preserves original data containers and methods, ensuring long-term support
# 3. **Enhanced Workflow Efficiency**: Enables streamlined analysis pipelines using familiar tidyverse syntax
# 
# The ecosystem includes several specialized packages:
# - `tidySummarizedExperiment`: For bulk RNA-seq analysis
# - `tidySingleCellExperiment`: For single-cell data
# - `tidySpatialExperiment`: For spatial transcriptomics
# - Additional tools: `plyranges`, `nullranges`, `tidyseurat`, `tidybulk`, `tidytof`
# 
# 
# [tidySpatialWorkshop](https://github.com/tidyomics/tidySpatialWorkshop) 
# [tidy transcriptomic manifesto](https://tidyomics.github.io/tidyomicsBlog/post/2021-07-07-tidy-transcriptomics-manifesto/)
# 
# `tidyomics` is an interoperable software ecosystem that bridges Bioconductor and the tidyverse. `tidyomics` is installable with a single homonymous meta-package. This ecosystem includes three new packages: tidySummarizedExperiment, tidySingleCellExperiment, and tidySpatialExperiment, and five publicly available R packages: `plyranges`, `nullranges`, `tidyseurat`, `tidybulk`, `tidytof`. Importantly, `tidyomics` leaves the original data containers and methods unaltered, ensuring compatibility with existing software, maintainability and long-term Bioconductor support. 
# 
# `tidyomics` is presented in "The tidyomics ecosystem: Enhancing omic data analyses" [Hutchison and Keyes et al., 2025](https://www.biorxiv.org/content/10.1101/2023.09.10.557072v1)

library(SpatialExperiment)

# Tidyverse library(tidyverse)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(purrr)
library(glue) # sprintf
library(stringr)

# Plotting
library(colorspace)
library(dittoSeq)
library(ggspavis)

# Analysis
library(scuttle)
library(scater)
library(scran)

# Similarly to **Section 2**, this section uses `spatialLIBD` and `ExperimentHub` packages to gather spatial transcriptomics data.

library(spatialLIBD)
library(ExperimentHub)

# To avoid error for SPE loading 
# https://support.bioconductor.org/p/9161859/#9161863
setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))

spatial_data <- 
  ExperimentHub::ExperimentHub() |> 
  spatialLIBD::fetch_data( eh = _, type = "spe")

# Clear the reductions
reducedDims(spatial_data) = NULL 

# Make cell ID unique
colnames(spatial_data)  = paste0(colnames(spatial_data), colData(spatial_data)$sample_id)
rownames(spatialCoords(spatial_data)) = colnames(spatial_data) # Bug?

# Display the object
spatial_data

# If `ExperimentHub` should not work. The `spatial_data` object from the previous code block can be downloaded from [Zenodo - 10.5281/zenodo.11233385](https://zenodo.org/records/11233385/files/tidySpatialWorkshop_spatial_data.rds?download=1)

# Working with tidySpatialExperiment --------------------------------------
# The `tidySpatialExperiment` package creates a bridge between `SpatialExperiment` objects and the tidyverse ecosystem. It provides:
# 
# 1. A tidy data view of `SpatialExperiment` objects
# 2. Compatible dplyr, tidyr, ggplot and plotly functions
# 3. Seamless integration with existing SpatialExperiment functionality


# 1. tidySpatialExperiment package ----------------------------------------
# `tidySpatialExperiment` provides a bridge between the `SpatialExperiment` single-cell package and the tidyverse [@wickham2019welcome]. It creates an invisible layer that enables viewing the `SpatialExperiment` object as a tidyverse tibble, and provides `SpatialExperiment`-compatible `dplyr`, `tidyr`, `ggplot` and `plotly` functions.

# If we load the `tidySpatialExperiment` package and then view the single cell data, it now displays as a tibble. 
library(tidySpatialExperiment)
spatial_data

# Data interface, display -------------------------------------------------
# If we want to revert to the standard `SpatialExperiment` view we can do that.
options("restore_SpatialExperiment_show" = TRUE)
spatial_data

# If we want to revert back to tidy SpatialExperiment view we can.
options("restore_SpatialExperiment_show" = FALSE)
spatial_data
# Note that **rows** in this context refers to rows of the abstraction, not **rows** of the SpatialExperiment which correspond to genes **tidySpatialExperiment** prioritizes cells as the units of observation in the abstraction, while the full dataset, including measurements of expression of all genes, is still available "in the background".

# Original behaviour is preserved -----------------------------------------
# The tidy representation behaves exactly as a native `SpatialExperiment`. It can be interacted with using [SpatialExperiment commands](https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html) such as `assays`.

assays(spatial_data)

# 2. Tidyverse commands ---------------------------------------------------
# We can also interact with our object as we do with any tidyverse tibble. We can use `tidyverse` commands, such as `filter`, `select` and `mutate` to explore the `tidySpatialExperiment` object. Some examples are shown below and more can be seen at the `tidySpatialExperiment` [website](https://stemangiola.github.io/tidySpatialExperiment/articles/introduction.html#tidyverse-commands-1).

# Select ------------------------------------------------------------------
# We can use `select` to view columns, for example, to see the filename, total cellular RNA abundance and cell phase. 
# If we use `select` we will also get any view-only columns returned, such as the UMAP columns generated during the preprocessing.
spatial_data |> select(.cell, sample_id, in_tissue, spatialLIBD)

# Note that some columns are always displayed no matter whet. These column include special slots in the objects such as reduced dimensions, spatial coordinates (mandatory for `SpatialExperiment`), and sample identifier (mandatory for `SpatialExperiment`). 

# Although the select operation can be used as a display tool, to explore our object, it updates the `SpatialExperiment` metadata, subsetting the desired columns.
colData(spatial_data) <- 
    spatial_data |> 
    select(.cell, sample_id, in_tissue, spatialLIBD)

# to modify the object
spatial_data <-
  spatial_data |>
  select(.cell, sample_id, in_tissue, spatialLIBD)

# To select columns of interest, we can use `tidyverse` powerful pattern-matching tools. For example, using the method `contains` to select 
spatial_data |> 
  select(.cell, contains("sum")) 

# Filter ------------------------------------------------------------------
# We can use `filter` to subset rows, for example, to keep our three samples we are going to work with.
# We just display the dimensions of the dataset before filtering
ncol(spatial_data)

spatial_data <- 
  spatial_data |> 
  filter(sample_id %in% c("151673", "151675", "151676"))

spatial_data

# Here we confirm that the tidy R manipulation has changed the underlining object.
ncol(spatial_data)

# In comparison the base-R method recalls the variable multiple times
spatial_data <- spatial_data[,spatial_data$sample_id %in% c("151673", "151675", "151676")]

# Or for example, to see just the rows for the cells in spatialLIBD region L1.
spatial_data |> dplyr::filter(sample_id == "151673", spatialLIBD == "L1")


# Flexible, more powerful filters with `stringr`
spatial_data |> 
  dplyr::filter(
    subject |> str_detect("Br[0-9]1"), 
    spatialLIBD == "L1"
  )

# Summarise ---------------------------------------------------------------
# The integration of all spot/pixel/cell-related information in one table abstraction is very powerful to speed-up data exploration ana analysis.
spatial_data |> 
  filter(sum_umi < 1000) |> 
  count(sample_id)

# Mutate ------------------------------------------------------------------
# We can use `mutate` to create a column. For example, we could create a new `Phase_l` column that contains a lower-case version of `Phase`. 

# Note that the special columns `sample_id`, `pxl_col_in_fullres`, `pxl_row_in_fullres`, `PC*` are view only and cannot be mutated.
spatial_data |>
  mutate(spatialLIBD_lower = tolower(spatialLIBD)) |>
  select(.cell, spatialLIBD, spatialLIBD_lower)

# We can update the underlying `SpatialExperiment` object, for future analyses. And confirm that the `SpatialExperiment` metadata has been mutated.
spatial_data <- 
  spatial_data |>
  mutate(spatialLIBD_lower = tolower(spatialLIBD))

spatial_data |> 
  colData() |>
  _[,c("spatialLIBD", "spatialLIBD_lower")]

# We can mutate columns for on-the-fly analyses and exploration. Let's suppose one column has capitalisation inconsistencies, and we want to apply a unique filter.
spatial_data |>
  mutate(spatialLIBD = tolower(spatialLIBD)) |>
  filter(spatialLIBD == "wm")

# Extract -----------------------------------------------------------------
# We can use tidyverse commands to polish an annotation column. We will extract the sample, and group information from the file name column into separate columns. 
# Simulate file path
spatial_data <- spatial_data  |> mutate(file_path = glue("../data/single_cell/{sample_id}/outs/raw_feature_bc_matrix/"))

# First take a look at the file column
spatial_data |> select(.cell, file_path)

# Extract specific identifiers from complex data paths, simplifying the dataset by isolating crucial metadata. This process allows for clearer identification of samples based on their file paths, improving data organization.

# Create column for sample
spatial_data <- spatial_data |>
  # Extract sample ID from file path and display the updated data
  tidyr::extract(file_path, "sample_id_from_file_path", "\\.\\./data/single_cell/([0-9]+)/outs/raw_feature_bc_matrix/", remove = FALSE)

# Take a look
spatial_data |> select(.cell, sample_id_from_file_path, everything())

# Unite -------------------------------------------------------------------
# We could use tidyverse `unite` to combine columns, for example to create a new column for sample id combining the sample and subject id (BCB) columns.
spatial_data <- spatial_data |> unite("sample_subject", sample_id, subject, remove = FALSE)

# Take a look
spatial_data |> select(.cell, sample_id, sample_subject, subject)
spatial_data |> separate("sample_subject", into = c("t1","t2"))

meta_obj <- colData(spatial_data) %>%
  data.frame() %>%
  separate(sample_subject,into = c("f1","f2"),sep = "_")

colData(spatial_data)$f1 <- meta_obj$f1
colData(spatial_data)$f2 <- meta_obj$f2

colData(spatial_data)

spatial_data <- spatial_data %>%
  separate("sample_subject", into = c("t1","t2"))

spatial_data

meta_obj %>%
  separate

# 3. Advanced filtering/gating and pseudobulk -----------------------------
# `tidySpatialExperiment` provide a interactive advanced tool for gating region of interest for streamlined exploratory analyses.

# This capability is powered by `tidygate`. We show how you can visualise your data and manually drawing gates to select one or more regions of interest using an intuitive tidy grammar. From https://bioconductor.org/packages/devel/bioc/vignettes/tidySpatialExperiment/inst/doc/overview.html

# Let's draw an arbitrary gate interactively
spatial_data
spatial_data |> 
  
  # Filter one sample
  filter(in_tissue, sample_id=="151673") |> 
  
  # Gate based on tissue morphology
  tidySpatialExperiment::gate(alpha = 0.1, colour = "spatialLIBD") 

spatial_data_gated <- tidygate_env$gates

# You can reload a pre-made gate for reproducibility
data(spatial_data_gated)

spatial_data <- 
  spatial_data |> 
  # Filter one sample
  filter(in_tissue, sample_id=="151673") |> 
  # Gate based on tissue morphology
  tidySpatialExperiment::gate(alpha = 0.1, colour = "spatialLIBD", programmatic_gates = tidySpatialWorkshop::spatial_data_gated) 

# `tidySpatialExperiment` added a `.gated` column to the `SpatialExperiment` object. We can see this column in its tibble abstraction.
spatial_data |> select(.cell, .gated)

# We can count how many pixels we selected with simple `tidyverse` grammar
spatial_data |> count(.gated)

# To have a visual feedback of our selection we can plot the slide annotating by our newly created column.
spatial_data |> 
  ggspavis::plotVisium(annotate = ".gated")

# We can also filter, for further analyses
spatial_data |> 
  filter(.gated == 1)


# **Exercise 2.1** --------------------------------------------------------
# Gate roughly the white matter layer of the tissue (bottom-left) and visualise in UMAP reduced dimensions where this manual gate is distributed.
# - Calculate PCA, UMAPs as we did for Session 1
# - Gate the area of white matter
# - Plot UMAP dimensions according to the gating

genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(spatial_data))
dec <- scran::modelGeneVar(spatial_data, subset.row = genes) 

# Get top variable genes 
dec <- scran::modelGeneVar(spatial_data, subset.row = genes, block = spatial_data$sample_id) 
hvg <- scran::getTopHVGs(dec, n = 1000)

rowData(spatial_data[head(hvg),])[,c("gene_id", "gene_name")]

# here the metadata for the genes
rowData(spatial_data)

spatial_data <- 
  spatial_data |> 
  scuttle::logNormCounts() |> 
  scater::runPCA(subset_row = hvg)

set.seed(42)
spatial_data <- scater::runUMAP(spatial_data, dimred = "PCA")

# do the gating
spatial_data <- 
  spatial_data |> 
  # Filter one sample
  filter(in_tissue, sample_id=="151673") |> 
  # Gate based on tissue morphology
  tidySpatialExperiment::gate(alpha = 0.1, colour = "spatialLIBD") 

# plot UMAP
scater::plotUMAP(spatial_data, colour_by = ".gated", point_size = 0.2) 

# test more object
spatial_data <- 
  spatial_data |> 
  # Gate based on tissue morphology
  tidySpatialExperiment::gate(alpha = 0.1, colour = "spatialLIBD") 

table(spatial_data$sample_id,spatial_data$.gated)
table(spatial_data$sample_id)

# 4. Work with features ---------------------------------------------------
# By default `tidySpatialExperiment` (as well as `tidySingleCellExperiment`) focus their tidy abstraction on pixels and cells, as this is the key analysis and visualisation unit in spatial and single-cell data. This has proven to be a practical solution to achieve elegant `tidy` analyses and visualisation.

# In contrast, bulk data focuses to features/genes for analysis. In this case its tidy representation with `tidySummarizedExperiment` prioritise features, exposing them to the user.
# If you want to interact with features, the method `join_features` will be helpful. For example, we can add one or more features of interest to our abstraction.

# Let's add the astrocyte marker GFAP 

# Find out ENSEMBL ID
rowData(spatial_data) |> 
  as_tibble() |> 
  filter( gene_name == "GFAP")

# Join the feature to the metadata
# add the gene expression to the metadata
# this can be a problem if run multiple times.
spatial_data <- 
  spatial_data |> 
  join_features("ENSG00000131095", shape="wide")

spatial_data |> 
  select(.cell, ENSG00000131095)

# **Exercise 2.2**
#   Join the endothelial marker PECAM1 (CD31, look for ENSEMBL ID), and plot in space the pixel that are in the 0.75 percentile of EPCAM1 expression. Are the PECAM1-positive pixels (endothelial?) spatially clustered?
#   
# - Get the ENSEMBL ID
# - Join the feature to the tidy data abstraction
# - Calculate the 0.75 quantile across all pixels `mutate()`
# - Label the cells with high PECAM1
# - Plot the slide colouring for the new label 

rowData(spatial_data) |> 
  as_tibble() |> 
  filter(gene_name == "PECAM1")

# Join the feature to the metadata
# add the gene expression to the metadata
spatial_data


test <- spatial_data |> 
  join_features("ENSG00000261371", shape="wide") %>%
  mutate(quantile_PECAM1 = pnorm(q = ENSG00000261371,sd = sd(ENSG00000261371),mean = mean(ENSG00000261371))) %>%
  mutate(high_PECAM1 = quantile_PECAM1>0.75)

spatial_data |> 
  join_features("ENSG00000261371", shape="wide") %>%
  mutate(quantile_PECAM1 = pnorm(q = ENSG00000261371,sd = sd(ENSG00000261371),mean = mean(ENSG00000261371))) %>% pull(quantile_PECAM1)

spatial_data |> 
  join_features("ENSG00000261371", shape="wide") %>%
  mutate(quantile_PECAM1_02 = quantile(ENSG00000261371,0.75)) %>% table(.$quantile_PECAM1_02)
  mutate(high_PECAM1 = quantile_PECAM1>0.75)

test %>%
  ggspavis::plotSpots(annotate = "high_PECAM1")+facet_wrap(~sample_id)

# 5. Summarisation/aggregation --------------------------------------------
# Distinct
# We can quickly explore the elements of a variable with distinct

spatial_data |> 
  distinct(sample_id)

# We can `distinct` across multiple variables
spatial_data |> 
  distinct(sample_id, Cluster)

# Count
# We can gather more information counting the instances of a variable
spatial_data |> 
  count(Cluster) |> 
  arrange(desc(n))

# We calculate summary statistics of a subset of data
spatial_data |> 
  filter(Cluster==1) |> 
  count(sample_id) |> 
  arrange(desc(n))

# Aggregate
# For summarised analyses, we can aggregate pixels/cells as pseudobulk with the function `aggregate_cells`. This also works for `SingleCellExeriment`.We obtain a `SummarizedExperiment`. 
spe_regions_aggregated <-
  spatial_data |>
  aggregate_cells(c(sample_id, spatialLIBD))

spe_regions_aggregated

spe_regions_aggregated@assays@data$counts[1:5,1:5]

# `tidyomics` allows to cross spatial, single-cell (Bioconductor and seurat), and bulk keeping a consistent interface.
library(tidySummarizedExperiment)

spe_regions_aggregated

# You will be able to apply the familiar `tidyverse` operations
spe_regions_aggregated |> 
  filter(sample_id == "151673")

# 6. tidyfying your workflow ----------------------------------------------
# We will take workflow used in **Session 2**, performed using mostly base R syntax and convert it to tidy R syntax. We will show you how the readability and modularity of your workflow will improve. 

# Subset to keep only on-tissue spots.

# **Base R approach:**
spatial_data_base <- spatial_data[, colData(spatial_data)$in_tissue == 1]

# **Tidyverse Approach:**
spatial_data_tidy <- 
  spatial_data |> 
  filter(in_tissue == 1) 

# **Specific Differences and Advantages:**
# The `tidyverse` `filter()` function clearly states the intent to filter the dataset, whereas the Base R approach uses subsetting which might not be immediately clear to someone unfamiliar with the syntax.

# The `tidyverse` approach inherently supports chaining further operations without manually checking dimensions, assuming that users trust the operation to behave as expected.

# Manipulating feature information
# For `SingleCellExperiment` there is no tidy API for manipulating feature wise data yet, on the contrary for `SummarizedExperiment`, because gene-centric the abstraction  allow for direct gene information manipulation. Currently, `tidySingleCellExperiment` and `tidySpatialExperiment` do not prioritize the manipulation of features (genes). 
# 
# While these functions can employ genes for cell manipulation and visualisation, as demonstrated in `join_features()`, they lack tools for altering feature-related information. Instead, their primary focus is on cell information, which serves as the main observational unit in single-cell data. This contrasts with bulk RNA sequencing data, where features are more central.
# The tidy API for `SingleCellExperiment` has feature-manipulation API among our plans. See [tidyomics challenges](https://github.com/orgs/tidyomics/projects/1)

# **Base R approach:**
is_gene_mitochondrial <- grepl("(^MT-)|(^mt-)", rowData(spatial_data)$gene_name)
rowData(spatial_data)$gene_name[is_gene_mitochondrial]

# Quality Control: --------------------------------------------------------
# Apply quality control measures to exclude cells based on mitochondrial content and read/gene count, a common indicator of cell health and viability.

# **Base R approach:**
spatial_data <- addPerCellQC(spatial_data, subsets = list(mito = is_gene_mitochondrial))

# Select expressed genes threshold
qc_mitochondrial_transcription <- colData(spatial_data)$subsets_mito_percent > 30
colData(spatial_data)$qc_mitochondrial_transcription <- qc_mitochondrial_transcription

# **Tidyverse Approach:**
spatial_data <- 
  spatial_data |> 
  # Add QC
  addPerCellQC(subsets = list(mito = is_gene_mitochondrial)) |> 
  # Add threshold in colData
  mutate(qc_mitochondrial_transcription = subsets_mito_percent > 30)

spatial_data |> select(.cell, qc_mitochondrial_transcription)


# **Specific Differences and Advantages:**
# `tidyverse` pipelines these operations without storing intermediate results, directly updating the dataset. Base R separates these steps, requiring manual tracking of variables and updating the dataset in multiple steps, increasing complexity and potential for errors.

# Direct Data Mutation: Tidyverse directly mutates the dataset within the pipeline, whereas Base R extracts, computes, and then reassigns values, which can be more verbose and less efficient in terms of workflow clarity and execution.

#### Group-specific analyses

# **Base R approach:**
# get gene for subset
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(spatial_data))

# Convert to list
spatial_data_list <- lapply(unique(spatial_data$sample_id), function(x) spatial_data[, spatial_data$sample_id == x])

# Detect sample-specific hughly-variable genes
marker_genes = 
  lapply( spatial_data_list,
          function(x){
            dec = scran::modelGeneVar(x, subset.row = genes)
            scran::getTopHVGs(dec, n = 1000)
          }
  ) 

head(unique(unlist(marker_genes)))

# **Tidyverse Approach: group_split**
# get gene for subset
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(spatial_data))

marker_genes = 
  spatial_data |> 
  
  # Grouping
  group_split(sample_id) |> 
  
  # Loop across the list elements
  map(~ .x |> 
        scran::modelGeneVar(subset.row = genes) |> 
        scran::getTopHVGs(n = 1000)
  ) |> 
  reduce(union)

marker_genes |> head()

# **Tidyverse Approach: nest**
spatial_data |> 
  nest(sample_data = -sample_id) |> 
  mutate(marker_genes = map(sample_data, ~ 
                              .x |> 
                              scran::modelGeneVar(subset.row = genes) |> 
                              scran::getTopHVGs(n = 1000)
  )) 


# **Specific Differences and Advantages:**
# `tidyverse` neatly handles grouping and plotting within a single chain, using `nest()` or `group_split()` and `map()` for compartmentalized operations, which organizes the workflow into a coherent sequence. 
# tidyverse's `map()` is a powerful functional language tool, which can return arbitrary types, such as `map_int`, `map_char`, `map_lgl`.It is integrated into the data manipulation workflow, making it part of the data pipeline.

# Multi-parameter filtering

# **Base R approach:** 
## # Mitochondrial transcription
qc_mitochondrial_transcription <- colData(spatial_data)$subsets_mito_percent > 30
colData(spatial_data)$qc_mitochondrial_transcription <- qc_mitochondrial_transcription

# ## Select library size threshold
qc_total_counts <- colData(spatial_data)$sum < 700
colData(spatial_data)$qc_total_counts <- qc_total_counts

# ## Select expressed genes threshold
qc_detected_genes <- colData(spatial_data)$detected < 500
colData(spatial_data)$qc_detected_genes <- qc_detected_genes

# ## Find combination to filter
colData(spatial_data)$discard <- qc_total_counts | qc_detected_genes | qc_mitochondrial_transcription

# # Filter
spatial_data = spatial_data[,!colData(spatial_data)$discard ]

# **Tidyverse Approach:**
spatial_data_filtered = 
  spatial_data |> 

  mutate(
    discard = 
      subsets_mito_percent > 30 |
      sum < 700 |
      detected < 500
  ) |> 
  filter(!discard)

# **Specific Differences and Advantages:**
# 
# **Tidyverse:** The code directly applies multiple filtering conditions within a single filter() function, making it highly readable and concise. The conditions are clearly laid out, and the operation directly modifies the spatial_data dataframe. This approach is more intuitive for managing complex filters as it condenses them into a singular functional expression.
# 
# **Base R:** The approach first calculates each condition and stores them within the colData of the dataset. These conditions are then combined to create a logical vector that flags rows to discard. Finally, it subsets the data by removing rows that meet any of the discard conditions. This method is more verbose and requires manually handling intermediate logical vectors, which can introduce errors and complexity in tracking multiple data transformations.
# 
# **Why tidyverse might be better in this context:**
# 
# **Coding efficiency:** `tidyverse` chains operations, reducing the need for intermediate variables and making the code cleaner and less error-prone.
# 
# **Readability:** The filter conditions are all in one place, which simplifies understanding what the code does at a glance, especially for users familiar with the tidyverse syntax.
# 
# **Maintainability:** Fewer and self-explanatory lines of code and no need for intermediate steps make the code easier to maintain and modify, especially when conditions change or additional filters are needed.

# -------------------------------------------------------------------------

library(scuttle)

# scater::isOutlier(colData(spatial_data)$subsets_mito_percent, type="higher") %>% table()

is_gene_mitochondrial <- grepl("(^MT-)|(^mt-)", rowData(spatial_data)$gene_name)

spatial_data |> 
  nest(sample_data = -sample_id) |> 
  mutate(prop_mito_isOutlier = map(sample_data, function(x){
    
    scater::isOutlier(x$subsets_mito_percent, type="higher") %>% table()
  } 
  )) %>%
  # tidy(prop_mito_isOutlier)
  unnest(prop_mito_isOutlier)

# 7. Visualisation --------------------------------------------------------
# Here, we will show how to use ad-hoc spatial visualisation, as well as `ggplot` to explore spatial data we will show how `tidySpatialExperiment` allowed to alternate between tidyverse visualisation, and any visualisation compatible with `SpatialExperiment`. 

# Ad-hoc visualisation: Plotting the regions
# Let's visualise the regions that spatialLIBD labelled across three Visium 10X samples.
spatial_data_filtered |> 
  ggspavis::plotSpots(annotate = "spatialLIBD") +
  facet_wrap(~sample_id) +
  scale_color_manual(values = libd_layer_colors |> str_remove("ayer")) +
  theme(legend.position = "none") +
  labs(title = "spatialLIBD regions")

# Custom visualisation: Plotting the regions
spatial_data_filtered |> 
  ggplot(aes(array_row, array_col)) +
  geom_point(aes(color = spatialLIBD)) +
  facet_wrap(~sample_id) +
  coord_fixed() +
  theme(legend.position = "none") +
  labs(title = "spatialLIBD regions")

# Custom visualisation: Plotting RNA output
# Now, let's observe what is the difference in total transcriptional cell output across regions. We can appreciate that different regions of these Visium slide is characterised by significantly different total RNA output. For example, the region one has a low R&D output, on the contrary regions to an L3, characterised by a high RNA output.
# We could conclude that when we use thresholding to filter "low-quality" pixels we have to be careful about possible biological and spatial effects.
spatial_data_filtered |> 
  ggplot(aes(sum_umi, color = spatialLIBD)) +
  geom_density() + 
  facet_wrap(~sample_id) +
  scale_color_manual(values = libd_layer_colors |> str_remove("ayer")) +
  scale_x_log10() +
  theme_bw()

# We provide another example of how the use of tidy. Spatial experiment makes custom visualisation, very easy and intuitive, leveraging `ggplot` functionalities. We will observe the relationship between mitochondrial transcription percentage, and total gene counts. We expect this relationship to be inverse as cells with higher mitochondrial transcription percentage tent to have a more limited transcriptional gene pool (e.g. for dieying or damaged cells).
spatial_data_filtered |> 
  ggplot(aes(subsets_mito_percent, sum_gene)) + 
  geom_point(aes(color = spatialLIBD), size=0.2) +  
  stat_ellipse(aes(group = spatialLIBD), alpha = 0.3) +
  scale_color_manual(values = libd_layer_colors |>
  str_remove("ayer")) +
  scale_y_log10() +
  theme_bw()

# Interestingly, if we plot the correlation between these two quantities we observe heterogeneity among regions, with L1 showing a very small association.
spatial_data_filtered |> 
  ggplot(aes(subsets_mito_percent, sum_gene)) + 
  geom_point(aes(color = spatialLIBD), size=0.2) +  
  scale_color_manual(values = libd_layer_colors |>    str_remove("ayer")) +
  geom_smooth(method="lm") + 
  facet_wrap(~spatialLIBD) + 
  scale_y_log10() +
  theme_bw()

# Let's take a step further and group the correlations according to samples, to see whether different samples show different correlations.
spatial_data |> 
  ggplot(aes(subsets_mito_percent, sum_gene)) + 
  geom_point(aes(color = spatialLIBD), size=0.2) +  
  scale_color_manual(values = libd_layer_colors |> str_remove("ayer")) +
  geom_smooth(aes(group = sample_id), method="lm") + 
  facet_wrap(~spatialLIBD) + 
  scale_y_log10() +
  theme_bw()

# As you can appreciate, the relationship between the number of genes, probed Purcell and their mitochondrial prescription abundance it's quite  consistent.
# **Excercise 2.3 (assisted)**
# To to practice the use of `tidyomics` on spatial data, we propose a few exercises that connect manipulation, calculations and visualisation. These exercises are just meant to be simple use cases that exploit tidy R streamlined language.
# We assume that the cells we filtered as non-alive or damaged, characterised by being enriched uniquely for mitochondrial, genes, and genes, linked to up apoptosis. It is good practice to check these assumption. This exercise aims to estimate what genes are differentially expressed between filtered and unfiltered cells. Then visualise the results.
# Use `tidyomic`/`tidyverse` tools to label dead cells and perform differential expression within each region. Some of the comments you can use are: `mutate`, `nest`, `map`, `aggregate_cells`, `tidybulk:::test_differential_abundance`, 

# A hist:
# 
# - spatial_data |> 
#   mutate(
#     dead = ...
#     
# - Aggregate by sample, dead status, ad annotated region
# 
# - `nest` by annotated region
# 
# - use `map` to test DE
# :::
# 
# ::: {.note}
**Excercise 2.4**

Inspired by our audience, let's try to use `tidyomics` to identify potential Amyloid Plaques.

Amyloid plaques are extracellular deposits primarily composed of aggregated amyloid-beta (Aβ) peptides. They are a hallmark of Alzheimer's disease (AD) and are also found in certain other neurodegenerative conditions.

Amyloid plaques can be found in the brains of mice, particularly in transgenic mouse models that are engineered to develop Alzheimer's disease-like pathology. 
Although amyloid plaques themselves are extracellular, the presence and formation of these plaques are associated with specific gene expression changes in the surrounding and involved cells. These gene markers are indicative of the processes that contribute to amyloid plaque formation, as well as the cellular response to these plaques ([Ranman et al., 2021](https://molecularneurodegeneration.biomedcentral.com/articles/10.1186/s13024-021-00465-0).)

- join_features()
- mutate() 
- ggspavis::plotSpots()
:::
  
  ```{r}
marker_genes_of_amyloid_plaques = c("APP", "PSEN1", "PSEN2", "CLU", "APOE", "CD68", "ITGAM", "AIF1")

rownames(spatial_data) = rowData(spatial_data)$gene_name

```

The excercise includes

- Join the features

- Rescaling

- Summarising signature (sum), `mutate()`

- Plotting colousing by the signature
:::
  
  
  **Session Information**
  
  ```{r}
sessionInfo()
```

**References**
  
  ```{css echo=FALSE}
.note {
  margin: 30px;
  padding: 1em;
  background: #FFF8F0;
    border: 1px solid #EFE8E0;
  border-radius: 10px;
}
```