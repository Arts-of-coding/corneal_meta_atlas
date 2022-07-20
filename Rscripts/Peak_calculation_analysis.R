# Note: if a library does not load, then use this below
# install.packages('PACKAGE_NAME)
# And if that does not work, then:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("PACKAGE_NAME")
#if (!require("pacman")) install.packages("pacman")
pacman::p_load(here,
               tidyverse,
               janitor, # Cleaning column names
               scales, # Transform axis scales
               ggrepel)

# loading in all important libraries
require("devtools")
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvoutlier)
library(limma)
library(knitr)
library(SingleCellExperiment)
library(scater)
library(RColorBrewer)
library(plot3D)
library(stringr)
library(SAVER)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(ggpubr)
library(circlize)
library(cowplot)
library(clustree)
library(grid)
library(gridExtra)
library(ape)
library(ggplot2)
library(DO.db)
library(clusterProfiler)
library(BiocParallel)
library(scANANSESeurat)
library(viridis)
library(DESeq2)
#install.packages("ggalluvial")
library(ggalluvial)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

###########################################################################
# Storing results in directories of your choice
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_peaks/"

# setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

# load in the region file where the peaks were found for the motifs
counts <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE_11032022/ATAC_motif/motif.count.txt', sep= '\t', header = T, row.names=1,comment.char = "")
counts <- counts[sum(rowSums(counts)==0),]
# test with one motif p53
