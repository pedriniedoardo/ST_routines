# Deconvolution {#sec-seq-deconvolution} ----------------------------------

# Introduction ------------------------------------------------------------
# Sequencing-based ST data can contain zero to multiple cells per spot, which might be fully or only partially covered by cells, depending on the spatial resolution of the platform and the tissue cell density (see also @sec-seq-introduction and the schematic Figure below). This aspect of the data implies that there may be a mixture of cell types in a spot and thus a mixture of transcriptional programs.

# To help understand these mixtures, at least 20 deconvolution techniques have been proposed for spot-level ST data. Some methods require borrowing insights from a scRNA-seq reference dataset, while others can be reference-free. Based on their underlying algorithms, @Li2023-benchmark-deconvolution have grouped methods into five categories:

#   - **probabilistic-based**: use Bayesian inference, likelihood estimation, or probabilistic modeling to estimate cell type compositions while incorporating uncertainty. Available tools include `r BiocStyle::Biocpkg("spacexr")` (i.e. `RCTD`), `r BiocStyle::Githubpkg("YMa-lab/CARD")`, `r BiocStyle::Biocpkg("SpatialDecon")`, and `r BiocStyle::Githubpkg("JEFworks-Lab/STdeconvolve")` in R, and `r BiocStyle::Githubpkg("SpatialTranscriptomicsResearch/std-poisson")`, `r BiocStyle::Githubpkg("almaan/stereoscope")`, `r BiocStyle::Githubpkg("scverse/scvi-tools")` (i.e. `DestVI`), `r BiocStyle::Githubpkg("DongqingSun96/STRIDE")`, and `r BiocStyle::Githubpkg("BayraktarLab/cell2location")` in Python.
# - **non-negative matrix factorization (NMF)-based**: decompose gene expression data into latent components representing different cell types. Available tools include `r BiocStyle::Githubpkg("RubD/Giotto")` (i.e. `SpatialDWLS`) and `r BiocStyle::Biocpkg("SPOTlight")` in R, and [NMFReg](https://github.com/tudaga/NMFreg_tutorial) in Python.
# - **graph-based**: use graph neural networks or graph-based optimization to model spatial relationships. Available tools include `r BiocStyle::Githubpkg("leihouyeung/SD2")` in R, and `r BiocStyle::Githubpkg("Su-informatics-lab/DSTG")` and `r BiocStyle::Githubpkg("ma-compbio/SpiceMIx")` in Python.
# - **optimal transport (OT)-based**: infer spatial gene expression distributions by mapping scRNA-seq and ST data. Available tools include `r BiocStyle::Githubpkg("zcang/SpaOTsc")` and `r BiocStyle::Githubpkg("rajewsky-lab/novosparc")` in Python.
# - **deep learning-based**: align and integrate single-cell and spatial transcriptomics data with neural networks. For example, `r BiocStyle::Githubpkg("broadinstitute/Tangram")` in Python.
# 
# Among these, `std-poisson`, `STdeconvolve`, and `SpiceMix` are reference-free methods. Methods that incorporate spatial location information are `CARD`, `DSTG`, `SD2`, `Tangram`, `cell2location`, `DestVI`, `std-poisson`, and `SpiceMix`.
# 
# In this section, we will demonstrate deconvolution of cell types per spot, using `RCTD` on Visium and VisiumHD datasets.

# Dependencies {#sec-seq-deconvolution-load-data} -------------------------
library(BiocParallel)
library(CARDspa)
library(DropletUtils)
library(ggplot2)
library(ggspavis)
library(OSTA.data)
library(patchwork)
library(pheatmap)
library(scran)
library(scater)
library(spacexr)
library(SpatialExperiment)
library(VisiumIO)
library(tidyverse)
library(zip)
library(ComplexHeatmap)

# In this example of Visium breast cancer data [@Janesick2023-high-res], we perform cell type deconvolution without a single-cell (Chromium) reference and compare the concordance with the provided Visium annotation provided by 10x Genomics.

# -------------------------------------------------------------------------
# library(BiocFileCache)
# 
# # location fo the cache before changing it
# tools::R_user_dir("BiocFileCache", which = "cache")
# bfc <- BiocFileCache()
# bfc@cache
# 
# # decouple the OSTA.data_load() call to avoid the zip issue
# id <- "Visium_HumanBreast_Janesick"
# bfc <- BiocFileCache()
# bfcinfo()
# 
# url <- "https://osf.io/5n4q3"
# 
# # verify 'id' with informative error if not
# id <- match.arg(id, OSTA.data_list(url))
# 
# # check if already cached
# q <- bfcquery(bfc, id)
# 
# # retrieve, zip, cache & return path
# no <- osfr::osf_retrieve_node(url)
# df <- osfr::osf_ls_files(no, id)
# 
# # dir.create(td <- tempfile())
# osfr::osf_download(df, path = "../../data/Visium_HumanBreast_Janesick", recurse=TRUE)
# 
# # zip(fnm <- paste0(id, ".zip"), dir())
# zip_filename <- paste0("../../data/Visium_HumanBreast_Janesick.zip") # Define the name
# # file_id <- dir("../../data/") %>% str_subset(paste0(df$name,collapse = "|"))
# files_to_zip <- paste0("../../data/Visium_HumanBreast_Janesick/")              # Get the list of files
# 
# # utils::zip(zipfile = zip_filename, files = files_to_zip)
# zip::zip(zipfile = zip_filename, files = files_to_zip)
# 
# # add to BiocCace
# bfcadd(bfc, zip_filename, fpath=file.path(".", zip_filename))

# id <- "Chromium_HumanBreast_Janesick"
# # check if already cached
# q <- bfcquery(bfc, id)
# 
# # retrieve, zip, cache & return path
# no <- osfr::osf_retrieve_node(url)
# df <- osfr::osf_ls_files(no, id)
# 
# # dir.create(td <- tempfile())
# osfr::osf_download(df, path = "../../data/Chromium_HumanBreast_Janesick", recurse=TRUE)
# 
# # zip(fnm <- paste0(id, ".zip"), dir())
# zip_filename <- paste0("../../data/Chromium_HumanBreast_Janesick.zip") # Define the name
# # file_id <- dir("../../data/") %>% str_subset(paste0(df$name,collapse = "|"))
# files_to_zip <- paste0("../../data/Chromium_HumanBreast_Janesick/")              # Get the list of files
# 
# # utils::zip(zipfile = zip_filename, files = files_to_zip)
# zip::zip(zipfile = zip_filename, files = files_to_zip)
# 
# # add to BiocCace
# bfcadd(bfc, zip_filename, fpath=file.path(".", zip_filename))
# -------------------------------------------------------------------------

# retrieve dataset from OSF repository
# id <- "Visium_HumanBreast_Janesick"
# pa <- OSTA.data_load(id)
# OSTA.data_list()
# dir.create(td <- tempfile())
# unzip(pa, exdir=td)

# read into 'SpatialExperiment'
vis <- TENxVisium(
  spacerangerOut="../../data/Visium_HumanBreast_Janesick/outs", 
  processing="filtered", 
  format="h5", 
  images="lowres") |> 
  import()

# retrieve spot annotations & add as metadata
df <- read.csv("../../data/Visium_HumanBreast_Janesick/annotation.csv")
# match the barcodes
cs <- match(colnames(vis), df$Barcode)
identical(colnames(vis), df$Barcode)

vis$anno <- factor(df$Annotation[cs])

# set gene symbols as feature names
rownames(vis) <- make.unique(rowData(vis)$Symbol)
vis

# check the metadata
colData(vis)

# plot the grid on the tissue
xy <- spatialCoords(vis) * scaleFactors(vis)
ys <- nrow(imgRaster(vis)) - range(xy[, 2])
xs <- range(xy[, 1])
box <- geom_rect(
  xmin=xs[1], xmax=xs[2], ymin=ys[1], ymax=ys[2], 
  col="black", fill=NA, linetype=2, linewidth=2/3)

# make the plot individula
plotVisium(vis, spots=FALSE, point_size=1) + box
plotVisium(vis, point_size=1, zoom=TRUE) + 
  plot_layout(nrow=1) & facet_null()

# plot the panel
plotVisium(vis, spots=FALSE, point_size=1) + box + 
  plotVisium(vis, point_size=1, zoom=TRUE) + 
  plot_layout(nrow=1) & facet_null()

# Deconvolution is performed after quality control, as detailed in @sec-seq-quality-control, and is usually performed on unnormalized and untransformed (i.e. raw) counts. Here, we quickly check some typically spot-level metrics.
sub <- list(mt=grep("^MT-", rownames(vis)))
vis <- addPerCellQCMetrics(vis, subsets=sub)

# plot the QC metrics
vis$log_sum <- log1p(vis$sum)

plotCoords(vis, 
           annotate="log_sum") + 
  ggtitle("log library size") + 
  theme(
    legend.key.width=unit(0.5, "lines"), 
    legend.key.height=unit(1, "lines")) +
  scale_color_gradientn(colors=pals::jet()) +
  geom_point(size = 1.5) 

plotCoords(vis, 
           annotate="subsets_mt_percent") + 
  ggtitle("% mitochondrial") +
  theme(
    legend.key.width=unit(0.5, "lines"), 
    legend.key.height=unit(1, "lines")) +
  scale_color_gradientn(colors=pals::jet()) +
  geom_point(size = 1.5) 

ggplot(
  data.frame(colData(vis)), 
  aes(x=sum, y=subsets_mt_percent)) + 
  geom_point() + geom_density_2d() +
  scale_x_log10() + scale_y_sqrt() +
  theme(aspect.ratio=2/3) +
  plot_layout(nrow=1) & theme(
    legend.key.width=unit(0.5, "lines"), 
    legend.key.height=unit(1, "lines")) & 
  scale_color_gradientn(colors=pals::jet())

# make a unique plot
plotCoords(vis, 
           annotate="log_sum") + 
  ggtitle("log library size") + 
  plotCoords(vis, 
             annotate="subsets_mt_percent") + 
  ggtitle("% mitochondrial") + 
  ggplot(
    data.frame(colData(vis)), 
    aes(x=sum, y=subsets_mt_percent)) + 
  geom_point() + geom_density_2d() +
  scale_x_log10() + scale_y_sqrt() +
  theme(aspect.ratio=2/3) +
  plot_layout(nrow=1) & theme(
    legend.key.width=unit(0.5, "lines"), 
    legend.key.height=unit(1, "lines")) & 
  scale_color_gradientn(colors=pals::jet())

# A few spots have low library sizes, and can be removed.
table(vis$sum > 1000)
vis <- vis[, vis$sum > 1000]

# We first visualize the spot-level cell type annotation provided by 10x Genomics.
plotCoords(vis, 
           annotate="anno", point_size=1, 
           pal=unname(pals::trubetskoy())) + 
  theme(legend.key.size=unit(0, "lines")) +
  geom_point(size = 1.5)

# Now, we load the single-cell (Chromium) reference data for the Visium dataset. To streamline the demonstration, we consolidate some of the cell type annotations provided by 10x Genomics (i.e. `Annotation`) into more generalized categories (i.e. `Annogrp`).

# retrieve dataset from OSF repo
# id <- "Chromium_HumanBreast_Janesick"
# pa <- OSTA.data_load(id)
# dir.create(td <- tempfile())
# unzip(pa, exdir=td)

# read into 'SingleCellExperiment'
# fs <- list.files(td, full.names=TRUE)
# h5 <- grep("h5$", fs, value=TRUE)
sce <- read10xCounts("../../data/Chromium_HumanBreast_Janesick/filtered_feature_bc_matrix.h5", col.names=TRUE)

# use gene symbols as feature names
rownames(sce) <- make.unique(rowData(sce)$Symbol)

# retrieve cell type labels
# csv <- grep("csv$", fs, value=TRUE)
cd <- read.csv("../../data/Chromium_HumanBreast_Janesick/cell_metadata.csv", row.names=1)

# ignore mixtures
lab <- cd$Annotation
lab[grepl("Hyb", lab)] <- NA 

# simplify annotations
pat <- c("B Cell"="B",
         "T Cell"="T",
         "Mac"="macro",
         "Mast"="mast",
         "DCs"="dendritic",
         "Peri"="perivas",
         "End"="endo",
         "Str"="stromal",
         "Inv"="tumor",
         "Myo"="myoepi",
         "Hyb" = "hybrid")

# -------------------------------------------------------------------------
# my implementation
df_meta <- cd %>%
  rownames_to_column("barcode")

df_pat <- data.frame(Annogrp = pat) %>%
  rownames_to_column("pattern")

# pat <- "B Cell"
# ann <- "B"
df_meta_full <- pmap(list(df_pat$pattern,df_pat$Annogrp), function(pat,ann){
  df_meta %>%
    filter(str_detect(Annotation,pattern = pat)) %>%
    mutate(Annotation2 = ann)
}) %>%
  bind_rows() %>%
  # collapse the hybrid annotations
  group_by(barcode) %>%  
  summarize(size = n(), 
            Annogrp = paste(Annotation2,collapse="|")) %>%
  # column_to_rownames("barcode") %>%
  # if hybrid make it as NA
  mutate(Annogrp = case_when(str_detect(Annogrp, "hybrid")~"hybrid",
                             T~Annogrp)) %>%
  # join the with the full metadata
  left_join(df_meta,y = .,by = c("barcode")) %>%
  column_to_rownames("barcode") %>%
  mutate(Annogrp_final = case_when(is.na(Annogrp)~Annotation,
                                   Annogrp == "hybrid" ~ NA,
                                   T~Annogrp))

# count the annotations
table(df_meta_full$Annogrp_final)

# check meta before update
colData(sce)

sce$Annogrp <- df_meta_full[colnames(sce),"Annogrp_final"]

# -------------------------------------------------------------------------
# # former implementation
# pat <- c(
#   "B Cell"="B", "T Cell"="T", "Mac"="macro", "Mast"="mast",
#   "DCs"="dendritic", "Peri"="perivas", "End"="endo",
#   "Str"="stromal", "Inv"="tumor", "Myo"="myoepi")
# for (. in names(pat))
#   lab[grep(., lab)] <- pat[.]
# lab <- gsub("\\s", "", lab)
# 
# meta_new_02 <- cd %>%
#   mutate(Annogrp = lab)
# 
# test_compare <- inner_join(meta_new_02 %>%
#              rownames_to_column(),
#            df_meta_full %>%
#              rownames_to_column(),by = "rowname")
# 
# table(test_compare$Annogrp.x,test_compare$Annogrp_final)
# -------------------------------------------------------------------------


# We only keep the Chromium data with an annotation and are not labeled as "Hybrid", as these correspond to mixed subpopulations.
sce <- sce[, !is.na(sce$Annogrp)]
dim(sce)

# RCTD --------------------------------------------------------------------
# Next, we perform deconvolution with `r BiocStyle::Biocpkg("spacexr")` (also known as RCTD)[@Cable2022-RCTD].
# [Note that RCTD can also be adapted to Visium HD data with `rctd_mode = "doublet"`, as demonstrated by [@deOliveira2025-high-def] and @sec-seq-workflow-visium-hd.]{.aside}
# By default, `runRctd()`'s `rctd_mode = "doublet"` specifies at most two subpopulations coexist in a data unit (i.e. within a spot); here, we set `rctd_mode = "full"` in order to allow for an arbitrary number of subpopulations to be fit instead.

rctd_data <- createRctd(vis,
                        sce,
                        cell_type_col="Annogrp")

res <- runRctd(rctd_data,
               max_cores=20,
               rctd_mode="full")

# Weights inferred by `RCTD` should be normalized such that proportions of cell types sum to 1 for each spot:
# scale weights such that they sum to 1
str(assay(res))
ws <- assay(res)
# ws <- sweep(ws, 2, colSums(ws), `/`)
# scale the weights
ws2 <- t(t(ws)/(colSums(ws)))
# confirm the scaling
colSums(ws2)

ws_rctd <- data.frame(t(as.matrix(ws2)))
round(ws_rctd[1:5, 1:5], 2)

# add proportion estimates to colData
id_col <- paste0(names(ws_rctd),".RCTD")
colData(vis)[id_col] <- ws_rctd[colnames(vis), ]

# confirm the addition of the metadata to the onject
colData(vis)

# CARD --------------------------------------------------------------------
# Another method that can be used is `CARD`. First, we rename the columns of spatial coordinates for `CARD`.
# <!-- First, we realized delayed count matrices and rename spatial coordinate columns to make the input data compatible with CARD [Note: RCTD handles the conversion internally with spacexr::check_counts()]{.aside} -->

# realize delayed matrices, as CARD does 
# not yet support delayed matrix handling
counts(sce) <- as(counts(sce), "sparseMatrix")
counts(vis) <- as(counts(vis), "sparseMatrix")
colnames(spatialCoords(vis)) <- c("x", "y")

# Next, we perform the `CARD` deconvolution. Here, we demonstrate `CARD`'s interoperability with `SingleCellExperiment` and `SpatialExperiment`. [Note: `CARD` can also take a reference matrix, a reference cell type annotation column, a spatial count matrix, and a spatial coordinates data frame as separate items in `sc_count`, `sc_meta`, `spatial_count`, and `spatial_location`, respectively. However, we encourage simplifying the process by using existing Bioconductor classes.]{.aside} The deconvolution result matrix is already normalized such that the sum of cell type proportions for each spot is equal to 1. 

set.seed(2025)
CARD_obj <- CARD_deconvolution(
  spe=vis,
  sce=sce,
  sc_count=NULL,
  sc_meta=NULL,
  spatial_count=NULL,
  spatial_location=NULL,
  ct_varname="Annogrp",
  ct_select=NULL,      # use all 'sce$Annogrp' cell types
  sample_varname=NULL, # use all 'sce' as one 'ref' sample 
  mincountgene=100,
  mincountspot=5)
ws_card <- CARD_obj$Proportion_CARD

# order cell type names alphabetically, as for RCTD
id_order <- colnames(ws_rctd) %>% str_replace_all(pattern = "\\.",replacement = " ")
ws_card <- data.frame(ws_card[,id_order])
round(ws_card[1:5, 1:5], 2)

# Visualization -----------------------------------------------------------
# First, we define a couple accessory functions.
.plt_xy <- function(ws, vis, col, point_size) {
  xy <- spatialCoords(vis)[rownames(ws), ]
  colnames(xy) <- c("x", "y")
  df <- cbind(ws, xy)
  ggplot(df, aes(x, y, col=.data[[col]])) + 
    coord_equal() + theme_void() + 
    geom_point(size=point_size)
}

.plt_decon <- function(ws, vis) {
  ps <- lapply(names(ws), \(.) .plt_xy(ws, vis, col=., point_size=0.3))
  ps |> wrap_plots(nrow=3) & theme(
    legend.key.width=unit(0.5, "lines"),
    legend.key.height=unit(1, "lines")) &
    scale_color_gradientn(colors=pals::jet())
}

# We can visualize deconvolution weights in x-y space, i.e., coloring by  the proportion of a given cell type estimated to fall within a given spot:
### RCDT
.plt_decon(ws=ws_rctd, vis)

### CARD
.plt_decon(ws=ws_card, vis)

# The deconvolution results can also be viewed as a heatmap, where rows = cells and columns = clusters:
plot_heat_ws <- function(ws, string){
  p <- pheatmap(ws, 
                show_rownames=FALSE, show_colnames=TRUE, main=string,
                cellwidth=12, treeheight_row=5, treeheight_col=5)
  return(p)
}
plot_heat_ws(ws_rctd, string="RCTD") 
plot_heat_ws(ws_card, string="CARD")

# compare the correlation between the scores
left_join(
  ws_rctd %>%
    rownames_to_column("barcode") %>%
    pivot_longer(names_to = "annotation",values_to = "score.rctd",-barcode),
  ws_card %>%
    rownames_to_column("barcode") %>%
    pivot_longer(names_to = "annotation",values_to = "score.card",-barcode)) %>%
  ggplot(aes(x=score.rctd,y=score.card)) + geom_point(shape = 1,alpha=0.5) +
  facet_wrap(~annotation,scale = "free") +
  theme_bw() +
  theme(strip.background = element_blank())

# In both methods, we see that more than half of the spots are estimated to have a stromal proportion of more than 50%. Few spots have an intense and distinct signal for cancerous subpopulations, DCIS1 and DCIS2. For the following analysis, we focus on `RCTD` as an example.

# For comparison with spot annotations provided by 10x Genomics, we include majority voted cell type from deconvolution by `RCTD`.
# Note that, because stromal cells show broad signals across the entire tissue, to better investigate immune cell signals, we remove stromal from the majority vote calculation for an alternative label: `RCTD_no_stroma`.
ws <- ws_rctd
# derive majority vote label
ids <- names(ws)[apply(ws, 1, which.max)]
names(ids) <- rownames(ws)
vis$RCTD.major <- factor(ids[colnames(vis)])
colData(vis)

# derive majority vote excluding stromal cells
ws_no_stroma <- ws[, colnames(ws) != "stromal"]
ids_no_stroma <- names(ws_no_stroma)[apply(ws_no_stroma, 1, which.max)]
names(ids_no_stroma) <- rownames(ws)
vis$RCTD.major.no_stroma <- factor(ids_no_stroma[colnames(vis)])

# We can visualize these three annotations spatially:
lapply(c("anno", "RCTD.major", "RCTD.major.no_stroma"), function(.){
  plotCoords(vis, annotate=.) + geom_point(size = 1.5)
}) %>%
  wrap_plots(nrow=1) &
  theme(legend.key.size=unit(0, "lines")) &
  scale_color_manual(values=unname(pals::trubetskoy())) 

# Note the strong stromal signals and macrophages being the second most common cell type for stromal cells.
# To help characterize subpopulations from deconvolution, we can view their distribution against the provided annotation:
cd <- data.frame(colData(vis))
df <- as.data.frame(with(cd, table(RCTD.major, anno)))
fd <- as.data.frame(with(cd, table(RCTD.major.no_stroma, anno)))

ggplot(df, aes(Freq, RCTD.major, fill=anno)) + 
  ggtitle("RCTD") +
  ggplot(fd, aes(Freq, RCTD.major.no_stroma, fill=anno)) + 
  ggtitle("RCTD_no_stroma") +
  plot_layout(nrow=1, guides="collect") &
  labs(x="Proportion", y=NULL) &
  coord_cartesian(expand=FALSE) &
  geom_col(width=1, col="white", position="fill") &
  scale_fill_manual(values=unname(pals::trubetskoy())) &
  theme_minimal() & theme(
    aspect.ratio=1,
    legend.key.size=unit(2/3, "lines"),
    plot.title=element_text(hjust=0.5))

# Next, we can investigate the agreement between the provided annotation against the two deconvolution majority vote labels.
hm <- function(mat, string){
  pheatmap(
    mat, show_rownames=TRUE, show_colnames=TRUE, main=string,
    cellwidth=10, cellheight=10, treeheight_row=5, treeheight_col=5)
}

hm(prop.table(table(vis$anno, vis$RCTD.major), 2), string="RCTD")
hm(prop.table(table(vis$anno, vis$RCTD.major.no_stroma), 2), string="RCTD_no_stroma") 

# -------------------------------------------------------------------------
# try the same approach but using the jaccard metrics
# define the jaccard score function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

# test
a <- c('potato', 'tomotto', 'chips', 'baloon')
b <- c('car', 'chips', 'bird', 'salt')

jaccard(a, b)

df_meta_full2 <- colData(vis) %>% as.data.frame() %>% rownames_to_column("barcode")

df_crossing <- crossing(ref_anno = unique(df_meta_full2$anno),
                        ref_RCTD.major = unique(df_meta_full2$RCTD.major))
df_crossing

# build the scatter plot
# ref_anno <- df_crossing$ref_anno[1]
# ref_RCTD.major <- df_crossing$ref_RCTD.major[1]
df_jaccard_score <- pmap(list(ref_anno = df_crossing$ref_anno,
                              ref_RCTD.major = df_crossing$ref_RCTD.major), function(ref_anno,ref_RCTD.major){
                                
                                # calculate the jaccard score
                                a <- df_meta_full2 %>%
                                  filter(anno == ref_anno) %>% pull(barcode)
                                b <- df_meta_full2 %>%
                                  filter(RCTD.major == ref_RCTD.major) %>% pull(barcode)
                                
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame(ref_anno= ref_anno,
                                                 ref_RCTD.major = ref_RCTD.major,
                                                 jaccard_score = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

head(df_jaccard_score)

# shape it as a matrix
mat_jaccard_score <- df_jaccard_score %>%
  pivot_wider(names_from = ref_RCTD.major,values_from = jaccard_score) %>%
  column_to_rownames("ref_anno")

mat_jaccard_score

# plot the matrix
ht_02 <- Heatmap(mat_jaccard_score,
                 name = "Jaccard score",
                 # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                 col = viridis::viridis(option = "turbo",n = 20),
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize = 8),
                 column_names_side = "bottom",
                 column_names_gp = gpar(fontsize = 8),
                 row_dend_reorder = FALSE,
                 column_dend_reorder = FALSE,
                 row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 show_column_names = T,
                 show_row_names = T)

ht_02
# -------------------------------------------------------------------------

# Overall, we observe agreement between the provided spot labels and the `RCTD` deconvolution derived annotations.
# Before cleaning up stromal, some immune cell types, such as dendritic and mast, never had a chance to have the highest cell type proportion. 
# On the left panel, among among all the spots annotated by `RCTD` as T cells, nearly all of them are from the "immune" type in the provided annotation. Strong agreements are also observed for spots with cell type of "DCIS1", "DCIS2", and "Invasive tumor".

# Next, we prepare the principal components (PCs) needed to perform PC regression:
# log-library size normalization
vis <- logNormCounts(vis)
# feature selection 
dec <- modelGeneVar(vis)
hvg <- getTopHVGs(dec, prop=0.1)
# dimension reduction 
set.seed(1234)
vis <- runPCA(vis, ncomponents=20, subset_row=hvg)

# We fit the deconvolution result of each cell type against the first 10 PCs to obtain 10 regressions.
idx <- rownames(ws)
ids <- colnames(ws)
pcs <- reducedDim(vis, "PCA")
pcs <- pcs[idx, seq_len(10)]
pcr <- lapply(ids, \(id) {
  fit <- summary(lm(pcs ~ ws[[id]]))
  r2 <- sapply(fit, \(.) .$adj.r.squared)
  data.frame(id, pc=seq_along(r2), r2)
}) |> do.call(what=rbind)

# Here we plot the coefficient of determination of the first 10 PCs for each cell type.
pcr$id <- factor(pcr$id, ids)
pal <- pals::trubetskoy()
ggplot(pcr, aes(pc, r2, col=id)) +
  geom_line(show.legend=FALSE) + geom_point() +
  scale_color_manual("predictor", values=unname(pal)) +
  scale_x_continuous(breaks=c(1, seq(5, 20, 5))) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(x="principal component", y="coeff. of determination") +
  guides(col=guide_legend(override.aes=list(size=2))) +
  coord_cartesian(xlim=c(1, 10)) +
  theme_minimal() + theme(
    panel.grid.minor=element_blank(),
    legend.key.size=unit(0, "lines"))

# Let's inspect the key drivers of (expression) variability in terms of PCs. 
# Considering deconvolution results from above, we can see that, e.g.:
# -   PC1 distinguishes stromal, tumor, macrophage from the rest of the tissue
# -   PC3, PC4 and PC5 separate DCIS1, T and endothelial cells, respectively
pcs <- reducedDim(vis, "PCA")
pcs <- pcs[rownames(ws), seq_len(10)]

# specify subpopulations & PCs to visualize
var <- c("DCIS.1", "T", "endo")
var <- c(var, colnames(pcs)[3:5])

# visualize deconvolution weights alongside PCs
lapply(var, \(.) {
  .plt_xy(
    cbind(ws, pcs), vis, col=., point_size=0.3) +
    scale_color_gradientn(., colors=pals::jet())
}) |>
  wrap_plots(nrow=2) & theme(
    plot.title=element_blank(),
    legend.key.width=unit(0.5, "lines"),
    legend.key.height=unit(1, "lines"))

# Note that the direction of each PC is irrelevant from how much variation it explains.

# In conclusion, deconvolution-based cell type proportion estimates  are able to recapitulate PCs and, in turn, expression variability. [Apart from being a tool for spot deconvolution, `RCTD` can be used as a label transfer tool to annotate imaging-based ST data, such as for Xenium and MERSCOPE. For this, the default `doublet_mode = "doublet"` should be used, and a certainty score would be returned to indicate doublets with two predicted cell types.]{.aside}

# Appendix ----------------------------------------------------------------
# Benchmarks {.unnumbered}
# Benchmarking studies of deconvolution methods often require generating synthetic spots to establish a ground truth for cell type proportions. This process involves either simulating artificial tissue patterns or aggregating counts from scRNA-seq or imaging-based spatial transcriptomics data into spots. However, it is equally important to evaluate how these methods perform in tissues with highly spatially-heterogeneous cell type compositions, such as real cancer samples. Below are three comprehensive benchmarking studies, two of which that incorporate both artificial and real datasets.
# - @Sang-aram2023-Spotless developed a pipeline to benchmark 11 deconvolution methods, including `RCTD`, across 63 synthetic, 3 binned, and 2 real datasets. `RCTD` and `cell2location` were the most recommended methods. Figure 2 gives an overview of the benchmarking results.
# - @Li2023-benchmark-deconvolution benchmarked 18 deconvolution methods, including `RCTD`, across 50 simulated and real datasets. Among these methods, `CARD`, `cell2location`, and `Tangram` are highly recommended. Figure 1 summarizes the method performance, and Figure 4 gives a flowchart of how to decide on which method to use.
# - @Gaspard-Boulinc2025-deconvolution review and compare available cell-type deconvolution methods, and provide a continuously updated web-based summary table.
