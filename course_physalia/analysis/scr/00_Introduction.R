# This R script was converted from Introduction.Rmd

# Original R Markdown setup chunk
# knitr::opts_chunk$set(echo = TRUE, cache = FALSE)

## Instructors

# Dr. Stefano Mangiola is leading the Computational Cancer immunology group at the South Australian immunoGENomics Cancer Institute (SAiGENCI). He uses single-cell and spatial technologies to investigate the tumor microenvironment and the immune system. Beyong data production, his focus in on the integration and modelling of large-scale single-cell data resources. He is the author of `tidytranscriptiomics` and co-leads the `tidyomics` endevour.
#
# - BLUESKY: https://bsky.app/profile/stemang.bsky.social South Australian immunoGENomics Cancer Institute (SAiGENCI).
#
# - TWITTER/X: https://x.com/steman_research
#
#
# Malvica Kharbanda expert of spatial analyses at the .
#
# - TWITTER/X: https://x.com/Malvikharbanda

## Workshop partner: Physalia

# Original R Markdown chunk for displaying an image
# library(here)
# knitr::include_graphics(here("inst/images/physalia-min.png"))
# Note: knitr::include_graphics is for R Markdown. In a plain R script, you might load and display images differently if needed.

## Workshop goals and objectives

### What you will learn

# -   The basics of spatial profiling technologies
# -   Analysis and manipulation of sequencing-based spatial data.
# -   The basics of tidy R analyses of biological data with `tidyomics`
# -   How to interface `SpatialExperiment` with tidy R manipulation and visualisation
# -   Analysis and manipulation of imaging-based spatial data.

## Getting started

### Local

# You can view the material at the workshop webpage
# [here](https://tidyomics.github.io/tidySpatialWorkshop/index.html).

## Workshop package installation

# If you want to install the packages and material post-workshop, the
# instructions are below. The workshop is designed for R `4.4` and
# Bioconductor 3.19.

# Original R Markdown chunk, eval=FALSE
# # Install workshop package
# #install.packages('BiocManager')
# BiocManager::install("tidyomics/tidySpatialWorkshop", dependencies = TRUE)
#     
# # Then build the vignettes
# BiocManager::install("tidyomics/tidySpatialWorkshop", build_vignettes = TRUE, force=TRUE)
# 
# # To view vignette
# library(tidySpatialWorkshop)
# vignette("Introduction")

## Interactive execution of the vignettes

# From command line, and enter the tidySpatialWorkshop directory.

# Original R Markdown chunk for shell commands
# # Open the command line
# git clone git@github.com:tidyomics/tidySpatialWorkshop.git
# Note: This is a shell command, not R code.

# Alternatively download the [git zipped package](https://github.com/tidyomics/tidySpatialWorkshop/archive/refs/heads/devel.zip). Uncompress it. And enter the directory.

# Announcements

# Tidyomics is now published in (Nature Methods)[https://www.nature.com/articles/s41592-024-02299-2]. And availabel for (free) here[https://www.biorxiv.org/content/10.1101/2023.09.10.557072v3].

# Introduction to Spatial Omics

### Objective

# Provide a foundational understanding of spatial omics, covering different technologies and the distinctions between imaging and
# sequencing in experimental and analytical contexts.

### Workshop Structure

#### Day 1

##### 1. Welcome and Introduction

# -   Introduction of the instructor
# -   Introduction of the crowd
# -   Overview and goals of the workshop.

##### 2. What is Spatial Omics?

# -   Definition and significance in modern biology.
# -   Key applications and impact.
# -   Overview of different spatial omics technologies.
# -   Comparison of imaging-based vs sequencing-based approaches.

##### 3. Sequencing Spatial Omics

# -   Detailed comparison of methodologies.
# -   Experimental design considerations.
# -   Data analysis challenges and solutions.

##### 5. Analysis of sequencing based spatial data

# -   Getting Started with SpatialExperiment.
# -   Data Visualisation and Manipulation.
# -   Quality control and filtering.
# -   Dimensionality reduction.
# -   Spatial Clustering.
# -   Deconvolution of pixel-based spatial data.

#### Day 2

##### 1. Introduction to tidyomics

# -   Use tidyverse on spatial, single-cell, pseudobulk and bulk genomic data   

##### 2. Working with tidySpatialExperiment

# -   tidySpatialExperiment package
# -   Tidyverse commands
# -   Advanced filtering/gating and pseudobulk
# -   Work with features
# -   Summarisation/aggregation
# -   tidyfying your workflow
# -   Visualisation

#### Day 3

##### 1. Imaging Spatial Omics

# -   Detailed comparison of methodologies.
# -   Experimental design considerations.
# -   Data analysis challenges and solutions.

##### 2. Spatial analyses of imaging data

# -   Working with imaging-based data in Bioconductor with MoleculeExperiment
# -   Aggregation and analysis
# -   Clustering
# -   Neighborhood analyses

# End of converted Introduction.Rmd