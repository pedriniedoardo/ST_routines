# 29  Feature selection & testing -----------------------------------------

# 29.1.1 Introduction -----------------------------------------------------
# Unsupervised clustering and label-transfer approaches yield (discrete, non-overlapping) grouping of cells by transcriptional -- and presumably functional -- similarity. 
# Identifying **differentially expressed genes (DEGs)**, i.e., genes that are up-/down-regulated in one or few subpopulation(s), helps characterize clusters (and, in an unsupervised setting, find biologically meaningful cluster labels).

# By contrast, methods to identify **spatially variable genes (SVGs)** aim to find genes with spatially correlated patterns of expression.
# Feature selection in terms of SVGs can be used as a spatially-aware alternative to mean-variance relationship-based HVGs (see @sec-seq-intermediate-processing), or to identify biologically informative genes as candidates for experimental follow-up.

# Here, we demonstrate selected methods to identify SVGs *de novo* (using `nnSVG`), and based on pre-computed spatial clusters (using `DESpace`).
# We further compare these to DEGs as well as HVGs, and discuss the conceptual similarities and differences between genes identified through these different approaches.

# 29.1.2 Dependencies -----------------------------------------------------
library(DESpace)
library(tidyverse)
library(ggspavis)
library(ggrepel)
library(nnSVG)
library(patchwork)
library(pheatmap)
library(scater)
library(scran)
library(SpatialExperiment)
library(ComplexHeatmap)
library(nnSVG)
library(spatialDE)

# set seed for random number generation
# in order to make results reproducible
set.seed(123)
# load data from previous chapters
# (post quality control & clustering)
spe <- readRDS("../../out/object/seq-spe_cl.rds")

# Non-spatial -------------------------------------------------------------
# We can also use methods originally developed for bulk and single-cell transcriptomics, such as those for detecting highly variable genes (HVGs) and differentially expressed genes (DEGs), on spatial omics data. 
# Note that HVGs are defined only based on molecular features (i.e., gene expression), and do not incorporate spatial information. 
# In the context of single-cell and ST, differentially expressed (DE) may refer to differences between groups of samples, or to differences between clusters (e.g., cell types, or spatial domains); below, we consider changes across spatial structures.

# Highly variable genes (HVGs) --------------------------------------------
# Feature selection is used as a preprocessing step to identify a subset of biologically informative features (e.g. genes), in order to reduce noise (due to both technical and biological factors) and improve computational performance.

# [Note that HVGs are defined based only on molecular features (i.e. gene expression), and do not take any spatial information into account. If the spatial patterns in gene expression in a dataset mainly reflect spatial distributions of cell types (defined by gene expression), then relying on HVGs for downstream analyses may be sufficient. However, if there is further biologically meaningful spatial structure that is not captured in this way, spatially-aware methods may be needed instead.]{.aside}

# Identifying a set of top "highly variable genes" (HVGs) is a standard step for feature selection in many scRNA-seq workflows, which can also be used as a simple and fast baseline approach in spot-based ST data, or high-plex imaging based ST data; for the former, this makes the simplified assumption that spots can be treated as equivalent to cells.

# The set of top HVGs can then be used as the input for subsequent steps, such as dimensionality reduction and clustering. For a comprehensive discussion on HVG selection, we refer readers to [OSCA](https://bioconductor.org/books/release/OSCA.basic/feature-selection.html#hvg-selection).

# In this example, we use the `r BiocStyle::Biocpkg("scran")` package [@Lun2016-scran] to identify HVGs in a two-step procedure. We first model the mean-variance relationship, which decomposes variance into a technical component (smooth fit) and biological component (deviation thereof).
# Secondly, HVGs may be selected based on their rank, using a fixed number `n` or proportion `p` of genes:

# [A more stringent selection may be obtained by also adding an `FDR.threshold` on the significance of large variances *relative* to other genes.]{.aside}
                                                                                                                                                                                                     
# fit mean-variance relationship, decomposing 
# variance into technical & biological components
dec <- modelGeneVar(spe)

# select top 2,000 ranked genes (in terms
# of their biological variance component)
hvg <- getTopHVGs(dec, n=2e3)

# visualize mean-variance relationship
# par(mar = c(4, 4, 0, 0))
# fit <- metadata(dec)
# plot(fit$mean, fit$var, cex = 0.5, xlab = "mean expression", ylab = "variance")
# points(dec[hvg, "mean"], dec[hvg, "total"], cex = 0.5, col = "dodgerblue")
# curve(fit$trend(x), add = TRUE, lwd = 2, col = "tomato")

# ggplot implementation
# model the variance by the expression
fit <- metadata(dec)

# save all the value of mean expression and variance per gene
df_all <- data.frame(
  mean = fit$mean,
  var  = fit$var) %>%
  rownames_to_column("gene")

# Highlighted the hvg 
df_hvg <- df_all %>%
  filter(gene %in% hvg)

# generate the plot of the varinace by expression
ggplot(df_all, aes(x = mean, y = var)) +
  # Base layer: all points (black, equivalent to first plot() call)
  geom_point(size = 0.5) + 
  
  # Highlight layer: HVGs (blue, equivalent to points() call)
  geom_point(data = df_hvg, aes(x = mean, y = var), 
             color = "dodgerblue", size = 0.5) +
  
  # Trend line (equivalent to curve() call)
  stat_function(fun = fit$trend, color = "tomato", linewidth = 1) +
  
  # Labels and theme
  labs(x = "mean expression", y = "variance") +
  theme_bw()

# Differentially expressed genes (DEGs) -----------------------------------
# Having clustered our spots (or cells), we can test for differentially expressed genes (DEGs) between clusters. These can be interpreted as marker genes and used to interpret clusters in terms of their function, and to help annotate them (i.e. assign biologically meaningful labels).

# For details on identifying DE genes between groups of cells from scRNA-seq data, we refer readers to [OSCA](https://bioconductor.org/books/release/OSCA.basic/marker-detection.html).

# Here, we will use pairwise t-tests and specifically test for upregulation (as opposed to downregulation), i.e. expression should be higher in the cluster for which a gene is reported to be a marker; see @sec-ind-clustering for details.

# differential gene expression analysis
mgs <- findMarkers(spe,
                   groups=spe$BayesSpace,
                   direction="up")
# select for a few markers per cluster
deg <- lapply(mgs, function(df){
  rownames(df)[df$Top <= 3]
})
length(deg <- unique(unlist(deg)))

# We can visualize selected marker genes as a heatmap where bins represent the average expression (here, log-transformed library size-normalized counts) of a given gene (= columns) in a given cluster (= rows).
# [To visually amplify differences, we use `scale = "column"`. This will, for every gene, subtract the mean and divide by the standard deviation across per-cluster means, bringing genes of potentially very different expression levels to a comparable scale.]

# compute cluster-wise averages
pbs <- aggregateAcrossCells(spe, 
                            ids = spe$BayesSpace, subset.row = deg, 
                            use.assay.type = "logcounts",
                            statistics = "mean")

# use gene symbols as feature names
assay(pbs)[1:5,1:5]

mtx <- t(assay(pbs))
colnames(mtx) <- rowData(pbs)$gene_name

# # using pheatmap package
# pheatmap(mat = mtx, scale = "column")

# use the complexHeatmap implementaiton
mtx_scaled <- scale(mtx)

# confirm the correct scaling by gene (in this case the columns)
colSums(mtx_scaled)
colSds(mtx_scaled)

# plot the heatmap
Heatmap(mtx_scaled)

# In addition, we can plot in x-y space, i.e. coloring spots by their expression of a given marker gene:
# select top-3 markers for each cluster & get gene symbols
gs <- unique(unlist(lapply(mgs, function(df){
  head(rownames(df), 3)
  })))

gs <- rowData(spe)$gene_name[match(gs, rownames(spe))]

# gene-wise spatial plots
ps <- lapply(gs, function(gene) {
  plotCoords(spe, 
             annotate = gene, 
             feature_names = "gene_name", 
             assay_name = "logcounts") })

# figure arrangement
wrap_plots(ps, nrow = 4) & 
  theme(legend.key.width = unit(0.4, "lines"), 
        legend.key.height = unit(0.8, "lines")) & 
  scale_color_gradientn(colors = rev(hcl.colors(9, "Rocket")))

# Spatially-aware ---------------------------------------------------------
# Here, we demonstrate brief examples of how to identify a set of top SVGs using (`r BiocStyle::Biocpkg("nnSVG")` [@Weber2023-nnSVG] and `r BiocStyle::Biocpkg("DESpace")` [@Cai2024-DESpace]). These methods are available through Bioconductor and can be easily integrated into Bioconductor-based workflows. 

# Spatially-variable genes (SVGs) -----------------------------------------
# SVGs are usually identified integrating gene expression measurements with spatial coordinates, either with or without predefined spatial domains.
# SVGs approaches can be broadly categorized into three types: overall SVGs, spatial domain-specific SVGs, and cell type-specific SVGs [@Yan2024-categorization-SVGs]. 

# The detection of **overall SVGs** is sometimes used as a feature selection step for further downstream analyses such as spatially-aware clustering (see @sec-ind-clustering).
# Several methods have been proposed for detecting overall SVGs; below, we report some few notable examples:
  
# - `r BiocStyle::Biocpkg("spatialDE")` [@Svensson2018-SpatialDE] and `r BiocStyle::Biocpkg("nnSVG")` [@Weber2023-nnSVG], which are based on a Gaussian process model;
# - [Moran's I](https://en.wikipedia.org/wiki/Moran%27s_I) and [Geary's C](https://en.wikipedia.org/wiki/Geary%27s_C), which ranks genes according to their observed spatial autocorrelation (see @sec-ind-spatial-statistics); and,
# - `r BiocStyle::Githubpkg("xzhoulab/SPARK")` [@Sun2020-SPARK; @Zhu2021-SPARK-X], which uses a non-parametric test of the covariance matrices of the spatial expression data.

# **Spatial domain-specific** SVGs target changes in gene expression between spatial domains, which are used to summarize the whole spatial information.
# Genes displaying expression changes across spatial clusters indicate SVGs. 
# Spatial domains can be predefined based on morphology knowledge, or identified through spatially-aware clustering approaches such as `r BiocStyle::Biocpkg("BayesSpace")` [@Zhao2021-BayesSpace] and `r BiocStyle::Biocpkg("Banksy")` [@Singhal2024-BANKSY].
# Methods that belong to this category include `r BiocStyle::Biocpkg("DESpace")` [@Cai2024-DESpace] in R, and `r BiocStyle::Githubpkg("jianhuupenn/SpaGCN")` [@Hu2021-SpaGCN] in Python.

# **Cell type-specific** SVG methods leverage external cell type annotations to identify SVGs within cell types. 
# They analyze interaction effects between cell types and spatial coordinates [@Yan2024-categorization-SVGs].
# Examplary methods in R include `r BiocStyle::Biocpkg("CTSV")` [@Yu2022-CTSV], *C-SIDE* (implemented in `r BiocStyle::Biocpkg("spacexr")`; [@Cable2022-C-SIDE]), and `r BiocStyle::Githubpkg("shanyu-stat/spVC")` [@Yu2024-spVC].

# nnSVG -------------------------------------------------------------------
# In this example, we use a small subset of the dataset for faster runtime. We select a subset of the data, by subsampling the set of spots and including stringent filtering for lowly expressed genes. 
# [A full analysis with `nnSVG` using all spots for this dataset and default filtering parameters for an individual Visium sample from human brain tissue (available from `r BiocStyle::Biocpkg("spatialLIBD")`) takes around 45 minutes on a standard laptop.]{.aside}

# set random seed for number generation
# in order to make results reproducible
set.seed(123)
# sample 100 spots to decrease runtime in this demo
# (note: skip this step in full analysis)
n <- 100
sub <- spe[, sample(ncol(spe), n)]

# filter lowly expressed genes using stringent criteria to decrease runtime in this demo
# (note: use default criteria in full analysis)

sub <- filter_genes(sub,     # filter for genes with...
                    filter_genes_ncounts=10, # at least 10 counts in
                    filter_genes_pcspots=3)  # at least 3% of spots

# re-normalize counts post-filtering
sub <- logNormCounts(sub)

# run nnSVG
set.seed(123)
sub <- nnSVG(sub,
             assay_name = "logcounts",
             n_threads = 20,
             verbose = T)

# extract gene-level results
res_nnSVG <- rowData(sub)

# show results
head(res_nnSVG, 3)

# count significant SVGs (at 5% FDR significance level)
table(res_nnSVG$padj <= 0.05)

# identify the top-ranked SVGs
res_nnSVG$gene_name[res_nnSVG$rank == 1]

# store the top 6 SVGs
res_nnSVG %>%
  data.frame() %>%
  arrange(rank) %>%
  top_n(n = 6,wt = -rank)

# pull the genes id in order
top_nnSVG <- res_nnSVG %>%
  data.frame() %>%
  arrange(rank) %>%
  top_n(n = 6,wt = -rank) %>%
  pull(gene_id)
  
# DESpace -----------------------------------------------------------------
# `DESpace` relies on pre-computed spatial domains to summarize the primary spatial structures of the data; see @sec-ind-clustering on clustering.
plotCoords(spe, 
           annotate="BayesSpace") +
  theme(legend.key.size=unit(0, "lines")) +
  scale_color_manual(values=unname(pals::trubetskoy()))

# run DESpace
res <- svg_test(spe, cluster_col="BayesSpace")
head(res_DESpace <- res$gene_results)

# count significant SVGs (at 5% FDR significance level)
table(res_DESpace$FDR <= 0.05)

# Downstream analyses  ----------------------------------------------------
# The set of top SVGs may be further investigated, e.g. by plotting the spatial expression of several top genes and via gene pathway analyses (i.e., comparing significant SVGs with known gene sets associated with specific biological functions). 
top_HVGs <- getTopHVGs(dec, n=(n <- 6))
top_DEGs <- lapply(mgs, function(df){
  rownames(df)[df$Top == 1]
  })
top_DEGs <- unique(unlist(top_DEGs))[seq_len(6)]
top_DESpace <- res_DESpace$gene_id[seq_len(6)]

# make a unique list of genes and pick all the ranks
top_genes <- c(top_HVGs,top_DESpace,top_nnSVG) %>% unique()

# pull all the rankings
# for DEGs it is not clear if the ranking makes sense
df_ranking <- list(
  # data.frame(gene_id = top_HVGs,rank = seq_along(top_HVGs),method = "HVG"),
  # # data.frame(gene = top_DEGs,rank = seq_along(top_DEGs),method = "DEG"),
  # data.frame(gene_id = top_DESpace,rank = seq_along(top_DESpace),method = "DESpace"),
  # data.frame(gene_id = top_nnSVG,rank = seq_along(top_nnSVG),method = "SVG")
  data.frame(gene_id = getTopHVGs(dec)) %>%
    mutate(rank = row_number()) %>%
    filter(gene_id %in% top_genes) %>%
    mutate(method = "HVG") %>%
    select(gene_id,rank,method),
  
  res_DESpace %>%
    mutate(rank = row_number()) %>%
    filter(gene_id %in% top_genes) %>%
    mutate(method = "DESpace") %>%
    select(gene_id,rank,method) %>%
    remove_rownames(),
  
  res_nnSVG %>%
    data.frame() %>%
    filter(gene_id %in% top_genes) %>%
    mutate(method = "SVG") %>%
    select(gene_id,rank,method) %>%
    remove_rownames()
  
) %>%
  bind_rows() %>%
  # add the gene name from the raw annotation
  left_join(y = rowData(sub) %>%
              data.frame() %>%
              dplyr::select(gene_id,gene_name),by = "gene_id") %>%
  # if missing use the gene_id
  mutate(gene = case_when(is.na(gene_name)~gene_id,
                          T~gene_name))

# plot the rankings
ggplot(df_ranking, aes(x = method, y = rank, group = gene)) +
  geom_line(aes(color = gene), size = 1.2, alpha = 0.8) +
  geom_point(aes(color = gene), size = 3) +
  geom_text_repel(aes(label = gene),
                  size = 3,
                  direction = "y",
                  nudge_x = 0.1,
                  bg.color = "white", # Adds a white outline around letters
                  bg.r = 0.15        # Radius of the outline (thickness)
                  ) +
  scale_y_reverse() +
  theme_minimal() +
  coord_cartesian(ylim = c(1,20)) +
  labs(
    title = "Comparison of Gene Rankings",
    subtitle = "Tracking gene positions across HVG, DESpace, and SVG methods",
    y = "Rank (1 = Top)",
    x = "Method"
  ) +
  theme(legend.position = "none") # Hide legend if labels are used on plot

# build also some metrics to compare the rankings
# Convert to Wide Format (Genes x Methods)
df_wide <- df_ranking %>%
  select(-c(gene_id,gene_name)) %>%
  pivot_wider(names_from = method, values_from = rank)

df_matrix <- df_wide %>%
  column_to_rownames("gene")

# Calculate W fro globa correlation
library(DescTools)
w_results <- KendallW(df_matrix,test = T)
w_results

# use spearman and kendall rank correlation for pairwise comparison
df_ranking

df_crossing <- crossing(id_ref = unique(df_ranking$method),
                        id_query = unique(df_ranking$method)) %>%
  # filter the unique combinations
  filter(id_ref < id_query)

# in case of missing values what is the behavipous or cor.test
getOption("na.action")

# id_ref <- "HVG"
# id_query <- "SVG"

df_cor <- pmap(df_crossing, function(id_ref,id_query){
  # keep track of the processing
  print(c(id_ref,id_query))
  
  # run the correlation analysis
  df_id_ref <- df_ranking %>%
    filter(method == id_ref) %>%
    select(gene,rank)
  
  df_id_query <- df_ranking %>%
    filter(method == id_query) %>%
    select(gene,rank)
  
  df_full <- full_join(df_id_ref,df_id_query,by = "gene",suffix = c(".ref",".query"))
  

  df_spearman <- cor.test(df_full$rank.ref, df_full$rank.query, method = "spearman") %>%
    broom::tidy() %>%
    mutate(ref = id_ref,
           query = id_query)
  
  df_kendall <- cor.test(df_full$rank.ref, df_full$rank.query, method = "kendall") %>%
    broom::tidy() %>%
    mutate(ref = id_ref,
           query = id_query)
  
  bind_rows(df_spearman,df_kendall) %>%
    mutate(comparison = paste0(id_ref,"_",id_query))
}) %>%
  bind_rows()

df_cor

# in case I am not interested in the ranking itself, but in the presence/absence of the gene in the list, I could use the jaccard score.
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

# test
a <- c('potato', 'tomotto', 'chips', 'baloon')
b <- c('car', 'chips', 'bird', 'salt')

jaccard(a, b)

# df_crossing is the same
# pull the top 100 genes for all methods. in this case the ranking is not interesting

df_ranking2 <- list(
  # data.frame(gene_id = top_HVGs,rank = seq_along(top_HVGs),method = "HVG"),
  # # data.frame(gene = top_DEGs,rank = seq_along(top_DEGs),method = "DEG"),
  # data.frame(gene_id = top_DESpace,rank = seq_along(top_DESpace),method = "DESpace"),
  # data.frame(gene_id = top_nnSVG,rank = seq_along(top_nnSVG),method = "SVG")
  data.frame(gene_id = getTopHVGs(dec)) %>%
    mutate(rank = row_number()) %>%
    filter(rank %in% 1:100) %>%
    mutate(method = "HVG") %>%
    select(gene_id,rank,method),
  
  res_DESpace %>%
    mutate(rank = row_number()) %>%
    filter(rank %in% 1:100) %>%
    mutate(method = "DESpace") %>%
    select(gene_id,rank,method) %>%
    remove_rownames(),
  
  res_nnSVG %>%
    data.frame() %>%
    filter(rank %in% 1:100) %>%
    mutate(method = "SVG") %>%
    select(gene_id,rank,method) %>%
    remove_rownames()
  
) %>%
  bind_rows() %>%
  # add the gene name from the raw annotation
  left_join(y = rowData(sub) %>%
              data.frame() %>%
              dplyr::select(gene_id,gene_name),by = "gene_id") %>%
  # if missing use the gene_id
  mutate(gene = case_when(is.na(gene_name)~gene_id,
                          T~gene_name))


df_jaccard_score <- pmap(list(id_ref = df_crossing$id_ref,
                              id_query = df_crossing$id_query), function(id_ref,id_query){
                                
                                # calculate the jaccard score
                                a <- df_ranking2 %>%
                                  filter(method == id_ref) %>% pull(gene_id)
                                
                                b <- df_ranking2 %>%
                                  filter(method == id_query) %>% pull(gene_id)
                                
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame("id_ref" = id_ref,
                                                 "id_query" = id_query,
                                                 "jaccard_score" = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

# check the table
df_jaccard_score

# Visualization -----------------------------------------------------------
# To visualize the expression levels of selected genes in spatial coordinates on the tissue slide, we can use plotting functions from the `r BiocStyle::Biocpkg("ggspavis")` package.

# get top SVGs from each method
gs <- list(
  HVGs=top_HVGs, 
  DEGs=top_DEGs, 
  DESpace=top_DESpace,
  nnSVG=top_nnSVG)

# get gene symbols from ensembl identifiers
idx <- match(unlist(gs), rowData(spe)$gene_id)
.gs <- rowData(spe)$gene_name[idx]

# expression plots for each top gene
ps <- lapply(seq_along(.gs), function(.) {
  plotCoords(spe, 
             point_size=0,
             annotate=.gs[.], 
             assay_name="logcounts", 
             feature_names="gene_name") + 
    if (. %% 6 == 1) list(
      ylab(names(gs)[ceiling(./6)]), 
      theme(axis.title.y=element_text()))
}) 
wrap_plots(ps, nrow=4) & theme(
  legend.key.width=unit(0.4, "lines"),
  legend.key.height=unit(0.8, "lines")) &
  scale_color_gradientn(colors = pals::parula())

# Comparison --------------------------------------------------------------
# We can compare the ranks of the genes detected by each method and compute pairwise correlations, highlighting two known cortical layer-associated SVGs: *MOBP* and *SNAP25*.
# While all four methods - HVGs, DEGs, SVGs identified through `nnSVG` and `DESpace` - yield similar gene ranks, each method captures distinct aspects of gene expression, with varying pairwise correlations between them.

# subset data for the 113 genes nnSVG is based on
sub_gene <- rowData(sub)$gene_id
sub_DESpace <- res_DESpace[sub_gene, ] 
sub_HVG <- as.data.frame(dec[sub_gene, ])

# aggregate DEGs for each spatial domain
subset_mgs <- lapply(mgs, function(x){
  x[sub_gene, 1:3] %>%
    as.data.frame()
})
# sub_DEG <- bind_cols(subset_mgs) 
sub_DEG <- do.call(cbind, subset_mgs) 
top_cols <- sub_DEG %>% select(ends_with(".Top"))
sub_DEG$Top <- do.call(pmin, c(top_cols, na.rm = TRUE))

# compute ranks 
sub_DESpace$rank <- rank(sub_DESpace$FDR, ties.method = "first")
sub_HVG$rank <- rank(-1 * sub_HVG$bio, ties.method = "first")
sub_DEG$rank <- rank(sub_DEG$Top, ties.method = "first")

# combine 'rank' column for each method
res_all <- list(
  nnSVG = as.data.frame(res_nnSVG),
  DESpace = sub_DESpace, HVGs = sub_HVG, DEGs = sub_DEG
)
rank_all <- do.call(cbind, lapply(res_all, `[[`, "rank"))
rank_all <- data.frame(gene_id = sub_gene, rank_all)

# known SVGs for this data
known_genes <- c("MOBP", "SNAP25")

# method names
method_names <- c("nnSVG", "DESpace", "HVGs", "DEGs")

# convert data structure
rank_long <- rank_all %>%
  left_join(data.frame(rowData(sub)[, c("gene_id", "gene_name")]), by = "gene_id") %>%
  select(all_of(c("gene_name", method_names))) %>%
  pivot_longer(cols = c(nnSVG, DESpace, HVGs, DEGs), 
               names_to = "method", 
               values_to = "rank")
# all pairwise comparisons
# df_pairs <- t(combn(method_names, 2)) |> as.data.frame()
df_pairs <- crossing(V1 = method_names, V2 = method_names) %>%
  filter(V1>V2)

# function to plot each pairwise comparison
plot_pairwise_comparison <- function(m1, m2, df) {
  
  # filter the data for the two methods being compared
  df <- df %>%
    filter(method %in% c(m1, m2)) %>%
    pivot_wider(names_from = method, values_from = rank)
  
  # compute pearson correlation between methods
  cor_val <- cor(df[[m1]], df[[m2]], method = "spearman")
  
  # plot
  ggplot(df, aes(x = .data[[m1]], y = .data[[m2]])) +
    geom_point() + 
    geom_text_repel(data = df %>% filter(gene_name %in% known_genes), 
                    aes(label = gene_name), color = "red", size = 3.25, 
                    nudge_x = 50, nudge_y = 10, box.padding = 0.5) +
    labs(x = paste(m1, "rank"), y = paste(m2, "rank"),
         title = paste(m2, "vs.", m1, ": Cor = ", round(cor_val, 2))) +
    #scale_color_manual(values = c("darkorange", "firebrick3", "deepskyblue2")) +
    theme_bw() + coord_fixed() + 
    xlim(c(0, 120)) + ylim(c(0, 120))
}

# generate and display all pairwise comparison plots
# plots <- lapply(seq_len(nrow(df_pairs)), function(i) {
#   plot_pairwise_comparison(df_pairs[i, 2], df_pairs[i, 1],rank_long)
# })
plots <- pmap(df_pairs,function(V1,V2){
  plot_pairwise_comparison(V1, V2,rank_long)
})

wrap_plots(plots, ncol = 3)
