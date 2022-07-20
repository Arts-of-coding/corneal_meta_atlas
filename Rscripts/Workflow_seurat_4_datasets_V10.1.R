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
library(fgsea)

###########################################################################
# Storing results in directories of your choice
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/"

# setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

# Loading in all functions
source("modify_vlnplot_v1.R")
source("stacked_vlnplot_v1.R")
source("QC_seurobj_v1.R")
source("PCA_seurobj_v1.R")
source("QC_seurobj_clustering_v1.R")
source("clustres_seurobj_v1.R")
source("clustres_seurobj_int_v1.R")

###########################################################################
# Importing the separate datasets and removing the doublets
reference <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_Collin_2021/20220113 CollinRNAannotated.rds")
query <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_Catala2021/20220114 CatalaRNAannotated.rds")
query2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_Li2021/20220114 LiRNAannotated.rds")
query3 <-readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_Gautam2022/...")

reference <- AddMetaData(reference, metadata="Co", col.name="Condition")

query <-AddMetaData(query, metadata="Ca", col.name="Condition")

query2 <- AddMetaData(query2, metadata="Li", col.name="Condition")

query3 <- AddMetaData(query3, metadata="Ga", col.name="Condition")

###########################################################################
# Remove doublets for all datasets

# Collin
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
seu_Collin <- reference
sweep.res.list_Collin <- paramSweep_v3(seu_Collin, PCs = 1:30, sct = FALSE)
sweep.stats_Collin <- summarizeSweep(sweep.res.list_Collin, GT = FALSE)
bcmvn_Collin <- find.pK(sweep.stats_Collin)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu_Collin@meta.data$costum_clustering)           ## ex: annotations <- seu_Collin@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_Collin@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

# Choose pK 0.06 based on the highest BC metrix
seu_Collin <- doubletFinder_v3(seu_Collin, PCs = 1:30, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_Collin <- doubletFinder_v3(seu_Collin, PCs = 1:30, pN = 0.25, pK = 0.06, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

DimPlot(seu_Collin, label=F,label.size = 8, group.by = 'DF.classifications_0.25_0.06_1806') + ggtitle("Doublet estimation")

reference <- subset(x=seu_Collin, subset = (DF.classifications_0.25_0.06_1806 == "Singlet"))
saveRDS(reference, file = paste0(resultsdir,"/Collin_singlets.rds"))

# Catala
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
seu_Catala <- query
sweep.res.list_Catala <- paramSweep_v3(seu_Catala, PCs = 1:18, sct = FALSE)
sweep.stats_Catala <- summarizeSweep(sweep.res.list_Catala, GT = FALSE)
bcmvn_Catala <- find.pK(sweep.stats_Catala)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu_Catala@meta.data$costum_clustering)           ## ex: annotations <- seu_Catala@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_Catala@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

# Choose pK 0.005 based on the highest BC metrix
seu_Catala <- doubletFinder_v3(seu_Catala, PCs = 1:18, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_Catala <- doubletFinder_v3(seu_Catala, PCs = 1:18, pN = 0.25, pK = 0.06, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

DimPlot(seu_Catala, label=F,label.size = 8, group.by = 'DF.classifications_0.25_0.005_1289') + ggtitle("Doublet estimation")

query <- subset(x=seu_Catala, subset = (DF.classifications_0.25_0.005_1289 == "Singlet"))
saveRDS(query, file = paste0(resultsdir,"/Catala_singlets.rds"))

# Li
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
seu_Li <- query2
sweep.res.list_Li <- paramSweep_v3(seu_Li, PCs = 1:16, sct = FALSE)
sweep.stats_Li <- summarizeSweep(sweep.res.list_Li, GT = FALSE)
bcmvn_Li <- find.pK(sweep.stats_Li)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu_Li@meta.data$costum_clustering)           ## ex: annotations <- seu_Li@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_Li@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

# Choose pK 0.18 based on the highest BC metrix
seu_Li <- doubletFinder_v3(seu_Li, PCs = 1:16, pN = 0.25, pK = 0.18, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_Li <- doubletFinder_v3(seu_Li, PCs = 1:16, pN = 0.25, pK = 0.18, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

DimPlot(seu_Li, label=F,label.size = 8, group.by = 'DF.classifications_0.25_0.18_1347') + ggtitle("Doublet estimation")

query2 <- subset(x=seu_Li, subset = (DF.classifications_0.25_0.18_1347 == "Singlet"))
saveRDS(query2, file = paste0(resultsdir,"/Li_singlets.rds"))

# Gautam
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
seu_Gautam <- query3
sweep.res.list_Gautam <- paramSweep_v3(seu_Gautam, PCs = 1:24, sct = FALSE)
sweep.stats_Gautam <- summarizeSweep(sweep.res.list_Gautam, GT = FALSE)
bcmvn_Gautam <- find.pK(sweep.stats_Gautam)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu_Gautam@meta.data$costum_clustering)           ## ex: annotations <- seu_Gautam@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_Gautam@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

# Choose pK 0.18 based on the highest BC metrix
seu_Gautam <- doubletFinder_v3(seu_Gautam, PCs = 1:24, pN = 0.25, pK = 0.14, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_Gautam <- doubletFinder_v3(seu_Gautam, PCs = 1:24, pN = 0.25, pK = 0.14, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

DimPlot(seu_Gautam, label=F,label.size = 8, group.by = 'DF.classifications_0.25_0.14_831') + ggtitle("Doublet estimation")

query3 <- subset(x=seu_Gautam, subset = (DF.classifications_0.25_0.14_831 == "Singlet"))
saveRDS(query3, file = paste0(resultsdir,"/Gautam_singlets.rds"))

############################################################
# Reload datasets with only Singlets
reference <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220207/Collin_singlets.rds")
query <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220207/Catala_singlets.rds")
query2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220207/Li_singlets.rds")
query3 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220211/Gautam_singlets.rds")

fullcornea_epi <- merge(reference, y = c(query,query2,query3), add.cell.ids = c("Collin", "Catala","Li","Gautam"), project = "fullcornea_epi")
unname(fullcornea_epi$Condition)

# Setting the correct parameters for the anndata conversion
fullcornea_epi$batch <-unname(fullcornea_epi$Condition)
fullcornea_epi@meta.data$batch <- unname(fullcornea_epi$Condition)
# generate training data for the machine learning model for every group

# Save it as a separate object raw integration of the four datasets
saveRDS(fullcornea_epi, file = paste0(resultsdir,"/fullcornea_epi.rds"))

############################################################
# generate a contingency table for 3 and 4 dataset integration
# import the datasets
scVI_4 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220214/fullcornea_epi.rds")
scVI_labels_4 <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/scVI_leiden_clusters_fc_epi_0.35_15022022.tsv', sep= '\t', header = F, row.names=1,comment.char = "")

# check if the rownames are the same after processing
sum(names(scVI_4@active.ident) == rownames(scVI_labels_4))

# add the IDs if correct
scVI_4$scVI_label <- as.factor(scVI_labels_4$V2)

# Generate the probaility table
cont_table_4 <- prop.table(table(scVI_4$scVI_label,scVI_4$costum_clustering),margin=2)*100

# Order the rows and columns to make it more clear
pdf(paste(resultsdir,'heatmap_cont_table_4_dataset.pdf',sep="/") ,width=10,height=5,paper='special')
f1 = colorRamp2(c(0, 100), c("white", "darkred"), space = "RGB")

# Order the rows and columns to make it more clear
cont_df <- as.data.frame.matrix(cont_table_4)
vec <- c(6,3,2,8,0,9,1,12,18,4,7,10,11,5,13,15,16,14,17,20,19)
vec_num<- vec
cont_df<-cont_df[order(match(rownames(cont_df), vec)), , drop = FALSE]
vec2 <- c("Ca_C10_CLC","Co_C6_LPC","Li_C6_PC","Ga_C5_Cj","Ga_C7_Cj","Li_C5_DC","Ga_C3_THCEC","Ga_C9_THCEC","Ca_C8_LSC","Co_C5_LNCP","Li_C1_DC","Li_C3_DC","Ca_C7_WSEL","Li_C8_TALSC","Li_C4_PC","Li_C2_Cj","Ga_C2_Cj",
          "Co_C3_CjB","Co_C2_CjS","Ga_C6_Cj","Ca_C5_BCE","Co_C1_SCE","Ca_C6_LCE","Ga_C4_EHCEC","Co_C10_BCE",
          "Co_C4_CSSC","Co_C7_LSt","Ca_C1_ASK","Ca_C4_TSK","Ca_C3_GSK","Ca_C2_GSK","Ga_C1_CF","Co_C8_LF","Ga_C10_CF","Ca_C11_CES",
          "Ca_C9_CEM","Li_C9_LSC","Co_C14_LV","Co_C9_BV","Ga_C8_Mel","Li_C7_Mel",
          "Ca_C12_UN","Co_C11_MEC","Ca_C13_ESD","Li_C10_UN","Ga_C12_CTC","Co_C13_IC","Ga_C11_Mon","Co_C12_FCEC")

cont_df<-cont_df[vec2]
matz <- as.matrix(cont_df)

vec3 <- as.character(vec)

ht1 = Heatmap(matz, col = f1, cluster_columns = F,cluster_rows = F, name = "Percentage_found"
              ,right_annotation  = rowAnnotation(foo = anno_text(vec3, location = 0.5, just = "center")))
htlist1 <- ht1

v = c("Limbal stem cell","Limbal stem cell","Limbal epithelial","Limbal epithelial","Conjunctival","Conjunctival"
      ,"Corneal epithelial","Corneal epithelial","Corneal epithelial","Stromal","Stromal","Stromal",
      "Stromal","Stromal","Stromal","Endothelial","Vessels","Melanocytes","Immune cells","Unknown","Unknown")
f2 <-viridis(10,option = "H",alpha = 0.8)
z = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

# Calculating the number of cells
num_cells <- c()
for (i in rownames(matz)){
  print(i)
  num_cells <- c(num_cells,sum(scVI_4@meta.data$scVI_label == i))
}

# Determining the colors of the subsets
show_col(viridis_pal(option = "H",alpha = 0.8)(10))
vec0<- viridis(10,option = "H",alpha = 0.8)
show_col(viridis_pal(option = "H",alpha = 0.8)(40))
vec <- viridis(40,option = "H",alpha = 0.8)
print(vec)

newcol <- c(vec[19:20],vec[15:16],vec[2:3],vec[5:7],vec[26:31],vec0[3],vec0[10],vec0[7],vec0[4],vec[38:39])


ha = rowAnnotation(foo = anno_text(as.character(num_cells), location = 0.5, just = "center",
                                   gp = gpar(fill = "black", col = "white", border = "black"),
                                   width = max_text_width(as.character(num_cells))*1.2))
#newcol <- c(vec[19:20],vec[15:16],vec[2:3],vec[5:7],vec[28:31],vec0[3],vec0[10],vec0[7],vec0[4],vec[38:39])

f3 <- newcol
small_df <- data.frame(num_cells,row.names = vec_num)
small_mat <- as.matrix(small_df)

ht_list = htlist1 + Heatmap(v, col = f2,name = "Cell origin", width = unit(0.3, "cm"), heatmap_legend_param = list(
  at = unique(v),
  labels = unique(v),
  title = "Cell origin",
  col = f2,
  legend_height = unit(4, "cm"))) +Heatmap(z, name = "# of cells", col = f3,cluster_rows = F, width = unit(1.5, "cm"), show_heatmap_legend = FALSE,
                                                                                              cell_fun = function(j, i, x, y, width, height, fill) {
                                                                                                grid.text(sprintf("%.1f", small_mat[i, j]), x, y, gp = gpar(fontsize = 11))
                                                                                              })
print(draw(ht_list))
dev.off()

z2 <- z-1
colmat <- matrix(c(z2, newcol), ncol = 2)
# Write the colors as a table
write.table(colmat, file = paste0(resultsdir,'/unannotated_meta_colors.tsv'), sep = '\t')

marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/LSC_Marker_Genes.csv"
marker_genes_df <- read.table(marker_gene_file, header = TRUE, sep = ',',comment.char = "#", stringsAsFactors = F)
seur_obj <- scVI_4
DefaultAssay(seur_obj) <- "RNA"

seur_obj@active.ident <- seur_obj$scVI_label
seur_obj@active.ident <- factor(seur_obj@active.ident, 
                                levels=vec_num)

fullmarkers2<-NULL

# Make dotplots for all datasets with their cell populations and marker genes
for (paper in unique(marker_genes_df$abreviation)){
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes_scVI_4_int/'),paper)
  print(marker_dir)
  system(paste("mkdir -p ", marker_dir))
  sub_marker_df <- marker_genes_df[marker_genes_df$abreviation == paper,]
  
  fullmarkers <- NULL
  for (cell_type in unique(sub_marker_df$name_plot)){
    marker_genes <- sub_marker_df[sub_marker_df$name_plot == cell_type,]$marker_genes
    marker_genes <- unlist(strsplit(marker_genes, ','))
    fullmarkers <- c(fullmarkers,marker_genes)
    m <- length(fullmarkers)
    n <- length(unique(sub_marker_df$name_plot))
    if (m<30){
      w <- 10
    }else{
      w <- m/5
    }
  }
  fullmarkers <- unique(fullmarkers)
  fullmarkers2<-c(fullmarkers2,fullmarkers)
  pdf(paste0(paste(paste(marker_dir, paper, sep = '/'), sep = '_'),'full_markers.pdf'), width =w, height = 5)
  print(DotPlot(object = seur_obj, features = fullmarkers,cluster.idents = T,)+ RotatedAxis())
  dev.off()
}

fullmarkers2 <- unique(fullmarkers2)
m <- length(fullmarkers2)

#cluster markergenes based on expression by complex heatmap
# Generate the pseudobulk-table for correlation between datasets
sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))
pseudobulk_df <- NULL
pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))
for (i in unique(seur_obj@active.ident)){
pseudobulk_df[[i]] <- as.vector(rowSums(counts(sce_qc)[,seur_obj@active.ident == i]))
}
pseudobulk_df <- pseudobulk_df[pseudobulk_df$`row.names(counts(sce_qc))`%in%fullmarkers2,]
rownames(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df<-pseudobulk_df[,-1]

matz <- as.matrix(pseudobulk_df)
ht1 = Heatmap(matz, col = f1, cluster_columns = T,cluster_rows = F)
colord <- row_order(ht1)
matz <- matz[colord,]
new_order <- rownames(matz)

pdf(paste0(resultsdir,'/full_markers_scVI4.pdf'), width =m/5, height = 5)
print(DotPlot(object = seur_obj, features = new_order,cluster.idents = F)+ RotatedAxis())
dev.off()

neworder2<- c("ACTA2","CDH19","CCL3","MITF","KIT","ACKR1","COL4A3","COL1A1","COL3A1","COL4A1","COL5A1","FBLN1","KERA","LUM","CD34","MMP2","MMP3","MMP1","PAX6","KRT3","S100A9","KRT7","MUC1","TP63","KRT6A","KRT14","CPVL","GPHA2")

seur_obj@active.ident <- seur_obj$scVI_label
vec4 <- rev(vec_num)
seur_obj@active.ident <- factor(seur_obj@active.ident, 
                                levels=vec4)
neworder3<-rev(neworder2)

pdf(paste0(resultsdir,'/sub_markers_scVI4.pdf'), width =7, height = 5)
print(DotPlot(object = seur_obj, features = neworder3,cluster.idents = F)+ RotatedAxis())
dev.off()

# mesenchymal stromal cells (cd73 = NT5E; thy1 = CD90; CD105 = END)

neworder4<- c("KRT19","PROM1","CD34","ALDH1A1","ALDH3A1","ACTA2","NT5E","THY1","ENG","MMP2",
              "MMP3","MMP1","MMP12","MMP","MMP14","COL1A1","COL3A1","COL4A1","COL5A1","B3GNT7",
              "CHST6","FN1","EDA","ABCG2","BMI1","ALCAM","NOTCH1","SIX2","CXADR","PTGDS","PDK4","LY6D","LY6E","LY6H","LY6K","LYPD2")
pdf(paste0(resultsdir,'/new_markers_scVI4.pdf'), width =8, height = 5)
print(DotPlot(object = seur_obj, features = neworder4,cluster.idents = F)+ RotatedAxis())
dev.off()

neworder4<- c("MUC1","KRT7","KRT8","MMP2","MMP9","MMP14")
pdf(paste0(resultsdir,'/new_markers_Cj_scVI4.pdf'), width =8, height = 5)
print(DotPlot(object = seur_obj, features = neworder4,cluster.idents = F)+ RotatedAxis())
dev.off()


neworder4<-c("KRT15","SOX9","ACTN1","FZD7","KRT17","ATF3","IFITM3","CD63","MT1A","SOCS3")
pdf(paste0(resultsdir,'/new_markers_LESC_scVI4.pdf'), width =8, height = 5)
print(DotPlot(object = seur_obj, features = neworder4,cluster.idents = F)+ RotatedAxis())
dev.off()


neworder4<-c("MMP10","KRT24","KRT15","GPHA2","KRT12","KRT3","KRT4","CEACAM7","PHLDA1","LAMA5","LAMA2","LAMA4","LAMB2","LAMB3",
             "LAMC2","FLG","IVL","GJB2","GJB6","GJA1","HSPG2","ITGA3","ITGB4","CAV1","CXCL14","CKS2","MOXD1","LAMA3","TNFRSF21")
pdf(paste0(resultsdir,'/new_markers_staining_scVI4.pdf'), width =8, height = 5)
print(DotPlot(object = seur_obj, features = neworder4,cluster.idents = F)+ RotatedAxis())
dev.off()

neworder2<- c("ACTA2","CDH19","CCL3","MITF","KIT","ACKR1","COL4A3","COL1A1","COL3A1","COL4A1","COL5A1","FBLN1","KERA","LUM","CD34","MMP2","MMP3","MMP1","PAX6","KRT24","KRT3","KRT12","S100A8","KRT7","MUC1","TXNIP","KRT6A","KRT14","CAV1","CPVL","GPHA2","KRT15","TP63")

neworder3<-rev(neworder2)

pdf(paste0(resultsdir,'/final_markers_scVI4.pdf'), width =11, height = 5)
print(DotPlot(object = seur_obj, features = neworder3,cluster.idents = F)+ RotatedAxis())
dev.off()

################################################################################
# Correlation analysis of cells determine if the clusters can be joined or not
# Generate the pseudobulk-table for correlation between datasets
sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc$sample <- seur_obj@active.ident
sce_qc$sample

pseudobulk_df <- NULL

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))
pseudobulk_df[["C2"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 2]))
pseudobulk_df[["C8"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 8]))

pseudobulk_df[["C0"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 0]))
pseudobulk_df[["C9"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 9]))

pseudobulk_df[["C1"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 1]))
pseudobulk_df[["C12"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample ==12]))
pseudobulk_df[["C18"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 18]))

pseudobulk_df[["C4"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 4]))
pseudobulk_df[["C7"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 7]))

pseudobulk_df[["C5"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 5]))
pseudobulk_df[["C13"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$sample == 13]))

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

# Generate a function to perform spearman correlation
var <- FindVariableFeatures(seur_obj)

pseudobulk_df <- pseudobulk_df[rownames(pseudobulk_df) %in%var@assays$RNA@var.features,]


# library(pspearman)
# 
# s <- unname(cor.test(pseudobulk_df$C2, pseudobulk_df$C13,  method = "spearman",exact = F, alternative = "greater")$statistic)
# n <- nrow(pseudobulk_df)
# 
# pspearman(s=s,n = 2, lower.tail = T,
#           approximation = c("exact", "AS89", "t-distribution"))


library(Kendall)

noms_x <- c()
noms_y  <- c()
vals <- c()
p.val <- c()
for (i in unique(colnames(pseudobulk_df))){
  for (j in unique(colnames(pseudobulk_df))){
    K <- Kendall(pseudobulk_df[[i]], pseudobulk_df[[j]])
    p.val<-c(p.val,K$sl[1])
    vals<-c(vals,K$tau[1])
    noms_x <- c(noms_x,i)
    noms_y <- c(noms_y,j)
    
  }
}

#colSums(matrix(vals, nrow=4))

df <- as.data.frame(x=vals)
df$cond1 <- noms_x
df$cond2 <- noms_y
df$cond3 <- noms_x==noms_y
# split the columns
df2 <-  df[df$cond3==T,]
df3 <-  df[df$cond3==F,]
df3 <- df3[match(unique(df3$vals), df3$vals),]
df3 <- df3[order(df3$cond1, df3$cond2,decreasing = T), ]
#df <- rbind(df2,df3)

M <- matrix(1, length(unique(noms_x)), length(unique(noms_x)))

transform(df3, Freq = ave(df3$f, df3$cond1, FUN = length))
df3$Freq = unname(table(df3$cond1)[df3$cond1])


# filling the matrix based on number of conditions
df3 <- df3[order(-ave(df3$Freq, df3$cond1, FUN = max), -df3$Freq), ]

df3$Freq2 = unname(table(df3$cond2)[df3$cond2])

df3 <- df3 %>% 
  arrange(desc(Freq),Freq2)

M[lower.tri(M)] <- df3[[1]] 


noms<-c(df3[1,2],unique(df3$cond2))

rownames(M)<-noms
colnames(M)<-noms

df3 <- df3 %>% 
  arrange(Freq,desc(Freq2))

M[upper.tri(M)]<- df3[[1]] 




pdf(paste(resultsdir,'Correlation_datasets_2_8.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C2, ylab = "2")
ggqqplot(pseudobulk_df$C8, ylab = "8")

#library(ConsRank)

#X<- as.matrix(pseudobulk_df)
#Tau_X(X, Y=NULL)

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C2, pseudobulk_df$C8,  method = "spearman",exact = F, alternative = "greater")['p.value']
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)
spearmanpval <- round(as.numeric(as.character(unname(spearmancortest['p.value'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=C2, y=C8)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C2,na.rm=T)/10*1,y = max(pseudobulk_df$C8,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()



pdf(paste(resultsdir,'Correlation_datasets_0_9.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C0, ylab = "0")
ggqqplot(pseudobulk_df$C9, ylab = "9")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C0, pseudobulk_df$C9,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=C0, y=C9)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C0,na.rm=T)/10*1,y = max(pseudobulk_df$C9,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()

pdf(paste(resultsdir,'Correlation_datasets_1_12.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C1, ylab = "1")
ggqqplot(pseudobulk_df$C12, ylab = "12")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C1, pseudobulk_df$C12,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=C1, y=C12)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C1,na.rm=T)/10*1,y = max(pseudobulk_df$C12,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()

pdf(paste(resultsdir,'Correlation_datasets_1_18.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C1, ylab = "1")
ggqqplot(pseudobulk_df$C18, ylab = "18")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C1, pseudobulk_df$C18,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=C1, y=C18)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C1,na.rm=T)/10*1,y = max(pseudobulk_df$C18,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()

pdf(paste(resultsdir,'Correlation_datasets_18_12.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C18, ylab = "18")
ggqqplot(pseudobulk_df$C12, ylab = "12")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C18, pseudobulk_df$C12,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=C18, y=C12)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C18,na.rm=T)/10*1,y = max(pseudobulk_df$C12,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()

pdf(paste(resultsdir,'Correlation_datasets_4_7.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C4, ylab = "4")
ggqqplot(pseudobulk_df$C7, ylab = "7")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C4, pseudobulk_df$C7,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=C4, y=C7)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C4,na.rm=T)/10*1,y = max(pseudobulk_df$C7,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()

pdf(paste(resultsdir,'Correlation_datasets_5_13.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$C5, ylab = "5")
ggqqplot(pseudobulk_df$C13, ylab = "13")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$C5, pseudobulk_df$C13,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)
#spearmancorstat <- round(as.numeric(as.character(unname(spearmancortest['statistic'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=C5, y=C13)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$C5,na.rm=T)/10*1,y = max(pseudobulk_df$C13,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()
################################################################################
# Continue with 4 datasets integration
# Re-name the clustering of the UMAP
cluster_order <- c('CE','CB','LE','LESC','CSSC','CF','LSC',
                   'SK','LE','Cj','SK','TSK','CE','CF',
                   'Mel','EC','Ves','IC','CE','MF','CDH19+')

# Generate new cluster labeling:
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)

new.cluster.ids <- c('CE','CB','LE','LESC','CSSC','CF','LSC',
                     'SK','LE','Cj','SK','TSK','CE','CF',
                     'Mel','EC','Ves','IC','CE','MF','CDH19+')

seur_obj$scVI_label <- plyr::mapvalues(x = as.factor(seur_obj$scVI_label), from = current.cluster.ids, to = new.cluster.ids)
seur_obj$scVI_label <- as.factor(seur_obj$scVI_label)

# Export labeling for visualization on the UMAP in Python
write.csv(seur_obj$scVI_label,paste0(resultsdir,"/labels_4datasets.csv"),row.names=F)

# Check clusters on  the first two PCs
seur_obj <- FindVariableFeatures(seur_obj)
seur_obj <- ScaleData(seur_obj)
seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = 2)

mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
pca <- seur_obj[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

# Checking the first PC
Stdev(object = seur_obj[["pca"]])[1]

seur_obj@active.ident <- seur_obj$scVI_label
pdf(paste0(resultsdir,'/scVI_4datasets.pdf') ,width=6,height=5,paper='special')
pc1_viz <-2
pc2_viz <-1
y_label = "PC2"#paste0(paste0(paste0("PC",pc2_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
x_label = "PC1"#paste0(paste0(paste0("PC",pc1_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
PC_dimred <- DimPlot(seur_obj, reduction = "pca",label = F, dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
print(PC_dimred)
dev.off()

#saveRDS(seur_obj,paste0(resultsdir,"/4datasets_annotated_joined.rds"))

seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")

################################################################################
# Plot literature marker genes on the names populations
seur_obj@active.ident <- seur_obj$scVI_label
vec_num <- c('LSC','LESC','LE','Cj','CE','CSSC',
                'SK','TSK','CF','MF',
                'EC','Ves','Mel','IC','CDH19+')
vec4 <- rev(vec_num)
seur_obj@active.ident <- factor(seur_obj@active.ident, 
                                levels=vec4)

neworder2<- c("ACTA2","CDH19","CCL3","MITF","ACKR1","COL4A3","COL4A1","COL5A1","FBLN1","KERA","LUM","CD34","MMP2","MMP3","MMP1","PAX6","AREG","KRT24","KRT3","KRT12","S100A8","MUC1","TXNIP","KRT6A","KRT14","CAV1","CPVL","GPHA2","MMP10","TP63")

neworder3<-rev(neworder2)

pdf(paste0(resultsdir,'/final_markers_meta.pdf'), width =9, height = 5)
print(DotPlot(object = seur_obj, features = neworder3,cluster.idents = F)+ RotatedAxis())
dev.off()

#Lets find marker genes for each cluster
cluster.markers <- FindAllMarkers(seur_obj, only.pos = TRUE)

# making a nice heatmap of expression in each cluster not significant yet
n_genes <- 100
heatmap.markers <- cluster.markers %>% group_by(cluster) %>% top_n(n_genes, avg_log2FC)
DoHeatmap(seur_obj, features = heatmap.markers$gene) + NoLegend()

#filter on only sigi
cluster.markers$gene_name_shorter <- str_sub(cluster.markers$gene,end = -1)
cluster.markers_fc <- cluster.markers[cluster.markers$p_val_adj < 0.01,]

cluster.markers_fc <- cluster.markers_fc[(cluster.markers_fc$avg_log2FC > 0.58) ,]
# avg_logFC...
table(cluster.markers_fc$cluster) # number of significantly expressed genes with a high log2FC (>0.58)

# GO terms
plot_list <- list()
cluster_GOs <- c()
#expressed_genes <- row.names(counts(sce_qc)[rowSums(counts(sce_qc))>20,])
expressed_genes <- rownames(cluster.markers_fc)

unique(expressed_genes)

# Write the markers as a table
write.table(cluster.markers_fc, file = paste0(resultsdir,'/cluster_markers_all.csv'), sep = ',')
write.table(heatmap.markers, file = paste0(resultsdir,'/heatmap_markers_all.csv'), sep = ',')

#saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated2.rds") # 10 PC
#saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated_4clus.rds")
#seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated2.rds") # 10 PC
#seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated_4clus.rds")

###########################################################################
# Appending the GO terms of the unannotated clusters
plot_list <- list()
cluster_GOs <- c()
#expressed_genes <- rownames(cluster.markers_fc)

cluster.markers <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/cluster_markers_all.csv', sep = ',', header = TRUE, row.names = 1)

df <- NULL
df_go <- NULL
for (cluster in unique(cluster.markers$cluster)){
  print(cluster)
  df_subset <- cluster.markers[cluster.markers$cluster == cluster,]
  cluster_genes <- unique(str_sub(df_subset$gene,end = -1))
  PC1_ego <- enrichGO(gene = cluster_genes,
                      #universe = expressed_genes,
                      OrgDb         = 'org.Hs.eg.db',
                      keyType       = "SYMBOL",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      #pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
  if(is.null(PC1_ego)){
    go_plot1 <- rectGrob(gp=gpar(col=NA))
  }else{
    go_plot1 <- barplot(PC1_ego, showCategory=10)
  }
  PC1_ego <- simplify(PC1_ego, cutoff = 0.8, by = "p.adjust", select_fun = min)
  df<-PC1_ego@result
  df$condition  <-rep((cluster),times=nrow(df))
  df_go <- rbind(df_go,df)
  plot_list <- c(plot_list, list(go_plot1))
  cluster_GOs <- append(cluster_GOs, (paste('cluster of cell type ', cluster, sep = '')))
}
dev.off()

# plotting the GO-terms
grob2 <- ggarrange(plotlist= plot_list, nrow = 21, ncol = 1, labels = cluster_GOs) #, align = 'v')
pdf(paste(resultsdir,'/all_DEGS_heatmap.pdf',sep="") ,width=10,height=40,paper='special')
print(ggarrange(grob2, ncol =1 , nrow = 18, widths= c(6, 4)))
print(grob2)
dev.off()

# Heatmap cluster marking genes:
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

heatmap.markers <- cluster.markers_fc
cluster.markers_fc$log2FC <- cluster.markers_fc$avg_log2FC
top10 <- cluster.markers_fc %>%  group_by(cluster) %>%  top_n(n = 5)
top5_genes <- DoHeatmap(slot = "data",seur_obj, features = top10$gene) + NoLegend()
#top10 <- cluster.markers_fc %>% group_by(cluster) %>% top_n(n = 10, wt = 'avg_log2FC')
#marker_top5 <- DoHeatmap(slot= "data",seur_obj, features = top10$gene) + NoLegend()

n_genes <- 100
heatmap.markers <- cluster.markers_fc %>% group_by(cluster) #%>% top_n(n_genes, avg_logFC)

cell_type_plot <- DimPlot(seur_obj, label=TRUE,pt.size = 2) + ggtitle("Timepoint")
cluster_plot <- DimPlot(seur_obj, label=TRUE, group.by = 'costum_clustering', pt.size = 2) + ggtitle("Louvain Clustering")
grob_clustering = ggarrange(plotlist = list(cell_type_plot, cluster_plot) , ncol =1) #labels = cluster_GOs) #align = 'v')
nclust<- length(unique(as.numeric(cluster.markers_fc$cluster)))
RGB_colours_ggplot <- as.list(as.character(gg_color_hue(nclust)))
names <- as.numeric(unique(heatmap.markers$cluster))
col_fun = colorRamp2(names, RGB_colours_ggplot)#make sure the rows correspond to ggplot colour mapping

seur_obj2 <- seur_obj
seur_obj@assays$SCT <- NULL
seur_obj@assays$integrated <- NULL

mat <- as.matrix(seur_obj@assays$RNA[heatmap.markers$gene,]) # removed scale data
mat <- rbind(mat,seur_obj@active.ident)
mat <- mat[,order(mat[nrow(mat),])]

column_ha <- HeatmapAnnotation(cluster = mat[nrow(mat),], col =list(cluster = col_fun))
breaks <- mat[nrow(mat),]
mat <- mat[-nrow(mat),]
row_ha = rowAnnotation(adj_p_vallue = heatmap.markers$p_val_adj)
clust_heatmap <- grid.grabExpr(draw(Heatmap(mat, column_split = breaks,row_split = as.numeric(heatmap.markers$cluster), cluster_columns = F, cluster_rows = F,show_row_names = F,show_column_names = F,row_names_gp = gpar(fontsize = 6), top_annotation = column_ha, left_annotation = row_ha, row_names_rot = -40)))

###########################################################################
# Plotting the GO terms of the un-annotated clusters and the heatmap

#expressed_genes <- rownames(cluster.markers_fc)

pdf(paste0(resultsdir,'/GO_terms_clusters_epi_sub.pdf') ,width=12,height=120,paper='special')
#cowplot::plot_grid(plotlist =grob2, ncol = 1, nrow = 1)#, rel_heights = c(1, 4,2), rel_widths = c(3,2,2))
print(grob2)
dev.off()

pdf(paste0(resultsdir,'/complex_heatmap_epi_sub.pdf') ,width=40,height=24,paper='special')
cowplot::plot_grid(top5_genes,grob_clustering, ncol = 3, nrow = 1, rel_heights = c(1, 4,2), rel_widths = c(3,2,2))
dev.off()

# Write the markers as a table if not done already
#write.table(cluster.markers, file = paste0(resultsdir,'/cluster_markers.csv'), sep = ',')


cluster.markers <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/cluster_markers_all.csv', sep = ',', header = TRUE, row.names = 1)
df <- NULL
df_go <- NULL
plot_list <- list()
cluster_GOs <- c()
for (cluster in unique(cluster.markers$cluster)){
  print(cluster)
  df_subset <- cluster.markers[cluster.markers$cluster == cluster,]
  cluster_genes <- unique(str_sub(df_subset$gene,end = -1))
  PC1_ego <- enrichGO(gene = cluster_genes,
                      #universe = expressed_genes,
                      OrgDb         = 'org.Hs.eg.db',
                      keyType       = "SYMBOL",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      #pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
  if(is.null(PC1_ego)){
    go_plot1 <- rectGrob(gp=gpar(col=NA))
  }else{
    go_plot1 <- barplot(PC1_ego, showCategory=10)
  }
  PC1_ego <- simplify(PC1_ego, cutoff = 0.8, by = "p.adjust", select_fun = min)
  df<-PC1_ego@result
  df$condition  <-rep((cluster),times=nrow(df))
  df_go <- rbind(df_go,df)
  plot_list <- c(plot_list, list(go_plot1))
  cluster_GOs <- append(cluster_GOs, (paste('cluster of cell type ', cluster, sep = '')))
}
dev.off()

# Go term dotplot
write.table(df_go, file = paste0(resultsdir,'/go_cluster_markers.csv'), sep = ',',row.names = F)


#df_go$gene.ratio <- foo$num$new
df_go2 <- df_go[df_go$Count>10 &df_go$p.adjust<0.05,]
df_go2 <- df_go2 %>% group_by(condition)  %>% slice_head(n=3)#%>% top_n(n = 10, wt = 'Count')

pdf(paste0(resultsdir,'/go_dot_epi_sub_sub.pdf') ,width=14,height=7,paper='special')

vec_cond<-c("LSC","LESC","LE","CE","Cj","Mel","IC","CSSC","CF","MF","SK","TSK","EC","Ves","CDH19+")
vec_cond <- rev(vec_cond)

df_go2<-df_go2[order(match(as.factor(df_go2$condition), vec_cond)), , drop = FALSE]

df_go2 <- df_go2[as.factor(str_order(df_go2$Description,decreasing = T)),]

ggplot(data = df_go2, aes(x = factor(condition,level=vec_cond), y = Description, 
                          color = `p.adjust`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis") + coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

pdf(paste0(resultsdir,'/go_dot_epi.pdf') ,width=7,height=5,paper='special')

vec_cond<-c("LSC","LESC","LE","Cj","CE")
df_go2 <- df_go[df_go$Count>10 &df_go$p.adjust<0.05,]
df_go2 <- df_go2[df_go2$condition %in% vec_cond,]

df_go2 <- df_go2 %>% group_by(condition) %>% arrange(desc(Count),.by_group = T)

df_go2 <- df_go2 %>% group_by(condition)  %>% slice_head(n=3)#%>% top_n(n = 10, wt = 'Count')
# df_go2<-df_go2[order(match(as.factor(df_go2$condition), vec_cond)), , drop = FALSE]
# 
# df_go2 <- df_go2[as.factor(str_order(df_go2$Description,decreasing = T)),]
#mid<-1.0e−06
p1 <-ggplot(data = df_go2, aes(x = factor(condition,level=vec_cond), y = Description, 
                          color = -log10(`p.adjust`), size = Count, show.legend = FALSE)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_size_area()
p1
dev.off()

#df_go$gene.ratio <- foo$num$new
df_go3 <- df_go[df_go$Count>10 &df_go$p.adjust<0.05,]
#df_go3 <- df_go3 %>% group_by(condition)  %>% slice_head(n=3)#%>% top_n(n = 10, wt = 'Count')

pdf(paste0(resultsdir,'/go_dot_non_epi.pdf') ,width=7,height=5,paper='special')

vec_cond2<-c("CSSC","CF","MF","SK","TSK","Mel","IC","EC","Ves","CDH19+")
#df_go2 <- df_go[df_go$Count>10 &df_go$p.adjust<0.0005,]
df_go3 <- df_go3[df_go3$condition %in% vec_cond2,]

df_go3 <- df_go3 %>% group_by(condition) %>% arrange(desc(Count),.by_group = T)

df_go3 <- df_go3 %>% group_by(condition)  %>% slice_head(n=3)#%>% top_n(n = 10, wt = 'Count')


p2<-ggplot(data = df_go3, aes(x = factor(condition,level=vec_cond2), y = Description, 
                          color = -log10(`p.adjust`), size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_size_area()
p2
dev.off()


# Joined go plots
library(patchwork)
design_go <- "AB"

pdf(paste0(resultsdir,'/go_wrap.pdf') ,width=14,height=5,paper='special')
wrap_plots(list(A=p1,B=p2),design=design_go)
dev.off()
write.table(df_go, file = paste0(resultsdir,'/go_cluster_markers.csv'), sep = ',',row.names = F)

###############################################################################
#Pahway analysis with progeny
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("progeny")

library(progeny)
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")

seur_obj <- progeny(seur_obj, scale=TRUE, organism="Human", top=500, perm=1, 
                return_assay = TRUE)

## We create a data frame with the specification of the cells that belong to 
## each cluster to match with the Progeny scores.
CellsClusters <- data.frame(Cell = names(Idents(seur_obj)), 
                            CellType = as.character(Idents(seur_obj)),
                            stringsAsFactors = FALSE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
seur_obj <- Seurat::ScaleData(seur_obj, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(seur_obj, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

mydata_full <- progeny_scores_df[,c("Pathway","Activity","CellType")]
for (i in unique(seur_obj$scVI_label)){
  
  pdf(paste0(resultsdir,"/",i,'Progeny.pdf') ,width=8,height=6,paper='special')
  mydata <- mydata_full[mydata_full$CellType==i,]
  mydata <- mydata[,c("Pathway","Activity")]
  names(mydata) <- c("group", "value")

# function for computing mean, DS, max and min values
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

p1 <- ggplot(aes(y = value, x = factor(group)), data = mydata)
p1 <- p1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") #+ geom_jitter(position=position_jitter(width=.2), size=3) + ggtitle("Boxplot con media, 95%CI, valore min. e max.") + xlab("Gruppi") + ylab("Valori")
print(p1)
dev.off()
}

pdf(paste0(resultsdir,"/heatmap_progeny_raw.pdf") ,width=5,height=4,paper='special')
mat <- t(summarized_progeny_scores_df)
mat <- mat[,c("LSC","LESC","LE","Cj","CE","CSSC","CF","MF","SK","TSK","Mel","IC","EC","Ves","CDH19+")]

f1 = colorRamp2(c(-1.5, 0, 1.5), c("blue", "#EEEEEE", "red"), space = "RGB")

ht1 = Heatmap(mat, col = f1, cluster_columns = F,cluster_rows = T, name = "progeny_raw_score")
print(ht1)
dev.off()

###############################################################################
#Gene set enrichment analysis with GSEA and epifactors database
epig_factors <- read.table(file = "genes.txt",sep = '\t',header = T)
library(plyr)

epig_factors <- epig_factors[epig_factors$Function!="-",]
table <- NULL
#i<-"PRC2"
for(i in unique(epig_factors$Function)){
  print(i)
  epi_sel <- epig_factors[epig_factors$Function==i,]
  col <- epi_sel$HGNC.approved.symbol
  print(col)
  table[[i]] <- col
}

max.length <- max(sapply(table, length))

l <- table
## Add NA values to list elements
l <- lapply(l, function(v) { c(v, rep(NA, max.length-length(v)))})
## Rbind
table2 <-do.call(rbind, l)

pseudobulk_df<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220309/pseudobulk_reps_DE_datasets_markers_split2.tsv',header = T,row.names = 1)

# For the deseq2 matrix vs ESC
countfile <- pseudobulk_df

#rm(pseudobulk_df)

lakovst <- countfile
#rm(countfile)
# Generate coldata dataframe
coldata <- NULL

# conditions
j <- unlist(strsplit(colnames(pseudobulk_df), split = "_"))
j <- j[seq(1,length(j),2)]

# Change _1 & _2 to condition 1 and _3 & _4 to condition 2
c <- c("1","2")
d <- paste(j,c,sep="_")

# reps
cols <- colnames(pseudobulk_df)

# type
xx <- rep("paired-end",length(cols))

coldata <- data.frame(cols=cols,condition=d,condition2=j,type=xx)

#coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/col2.tsv', sep = '\t', header = TRUE, row.names = 1)
#coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210802/col.tsv', sep = '\t', header = TRUE, row.names = 1)

# for complex heatmap
#coldata <- coldata[1:22,]
rownames(coldata)<- coldata$cols
coldata <- coldata[,c("condition","condition2","type")]

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# checking the data
head(lakovst,2)

coldata

# setting the rowdata in coldata similar to the coldata in the count matrix
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(lakovst))

all(rownames(coldata) == colnames(lakovst))

lakovst <- lakovst[, rownames(coldata)]
all(rownames(coldata) == colnames(lakovst))


pathways.hallmark <- gmtPathways("msigdb.v7.2.symbols.gmt")

# show a few lines from the pathways file
head(pathways.hallmark)

pathways.hallmark <- pathways.hallmark[1:1000]

pathways.epig <- table

pathways <- pathways.epig

df_gsea <- NULL
LSC_gsea <- NULL
LE_gsea <- NULL
CSSC_gsea <- NULL
res_up_LSC <- NULL
res_down_LSC <- NULL
ranks_LSC <- NULL
res_up_LE <- NULL
res_down_LE <- NULL
ranks_LE <- NULL
res_up_CSSC <- NULL
res_down_CSSC <- NULL
ranks_CSSC

for(i in unique(coldata$condition2)){
print(i)
coldata2 <- coldata
coldata2$condition2

coldata2$condition2[which(coldata2$condition2!=i)] <- "others"

dds <- DESeqDataSetFromMatrix(countData = round(lakovst),
                              colData = coldata2,
                              design = ~ condition2)
dds

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
rm(keep)
# data tranfromation
#vsd <- vst(dds, blind=FALSE)

#________________DE_analysis_____________#
dds <- DESeq(dds) #This would take some time
res <- results(dds, alpha=0.05)
summary(res)

#_________________GSEA___________________#
# Steps toward doing gene set enrichment analysis (GSEA):

# 1- obtaining stats for ranking genes in your experiment,
# 2- creating a named vector out of the DESeq2 result
# 3- Obtaining a gene set from mysigbd
# 4- doing analysis


# already we performed DESeq2 analysis and have statistics for working on it
res$SYMBOL <- rownames(res)

# # important notice: if you have not such stats in your result (say comming from edgeR),
# # you may need to create a rank metric for your genes. To do this:
# # metric = -log10(pvalue)/sign(log2FC)
# 
# 
# 
# # Map Ensembl gene IDs to the symbol. First, create a mapping table.
# ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
#                                     key=res$row, 
#                                     columns="SYMBOL")#,
#                                     #keytype="ENSEMBL")
# names(ens2symbol)[1] <- "row"
# 
# ens2symbol <- as_tibble(ens2symbol)
# ens2symbol
# # joining
# res <- merge(data.frame(res), ens2symbol, by=c("row"))

# This is nessesary for genes with GRCh38.p13
# remove the NAs, averaging statitics for a multi-hit symbol
# res <- as.data.frame(res)
# res2 <- res %>% 
#   select(stat,SYMBOL) %>% 
#   na.omit() %>% 
#   group_by(SYMBOL) %>% 
#   summarize(stat=mean(stat))
#res2

# WIP
# creating  a named vector [ranked genes]
ranks <- res$stat
names(ranks) <- res$SYMBOL

library(snow)
SnowParam(workers = 1)

#Running fgsea algorithm:
fgseaRes <- fgseaMultilevel(pathways=pathways, stats=ranks)

# Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)

# To see what genes are in each of these pathways:
# gene.in.pathway <- pathways.hallmark %>% 
#   enframe("pathway", "SYMBOL") %>% 
#   unnest(cols = c(SYMBOL)) %>% 
#   inner_join(res, by="SYMBOL")

# add gene_ratio
vec <- lengths(pathways)

df <- data.frame(lapply(vec, type.convert), stringsAsFactors=FALSE)
df2 <- as.data.frame(t(df))

names_pathway <- as.character(names(pathways))
names(names_pathway)<-NULL

df2$pathway <- names_pathway
names(df2) <- c("total","pathway")

fgseaResTidy_ratio <- merge(fgseaResTidy, df2, by = "pathway")
fgseaResTidy_ratio$ratio <- fgseaResTidy_ratio$size/fgseaResTidy_ratio$total

pathways.hallmark[2]

fgseaResTidy_ratio$condition  <-rep((i),times=nrow(fgseaResTidy_ratio))
res_up <- res[res$log2FoldChange>0,]
res_down <- res[res$log2FoldChange<0,]
if (i=="LSC"){
  LSC_gsea <- fgseaResTidy
  res_up_LSC <- res_up$SYMBOL
  res_down_LSC <- res_down$SYMBOL
  ranks_LSC <- ranks} else if ((i=="LE")){
    LE_gsea <- fgseaResTidy
    res_up_LE <- res_up$SYMBOL
    res_down_LE <- res_down$SYMBOL
    ranks_LE <- ranks} else if (i=="CSSC"){
      CSSC_gsea <- fgseaResTidy
      res_up_CSSC <- res_up$SYMBOL
      res_down_CSSC <- res_down$SYMBOL
      ranks_CSSC <- ranks}


df_gsea <- rbind(df_gsea,fgseaResTidy_ratio)
}

# Go term dotplot
str(df_gsea)
df_gsea$leadingEdge <- as.character(df_gsea$leadingEdge)
write.table(df_gsea, file = paste0(resultsdir,'/gsea_DEG_one_v_all.csv'), sep = ',',row.names = T)

df_gsea2 <- df_gsea[df_gsea$size>2 &df_gsea$padj<0.05,]
df_gsea2 <- df_gsea2 %>% group_by(condition)  %>% slice_head(n=10)#%>% top_n(n = 10, wt = 'Count')

pdf(paste0(resultsdir,'/gsea_dot.pdf') ,width=6,height=3,paper='special')
vec_cond<-c("LSC","LE","CSSC")
vec_cond <- rev(vec_cond)

df_gsea2<-df_gsea2[order(match(as.factor(df_gsea2$condition), vec_cond)), , drop = FALSE]
df_gsea2 <- df_gsea2[!is.na(df_gsea2$pathway),]
#df_gsea2 <- df_gsea2[as.factor(str_order(df_gsea2$pathway,decreasing = T)),]

ggplot(data = df_gsea2, aes(x = factor(condition,level=rev(vec_cond)), y = pathway, 
                          color = `padj`, size = size)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("fgsea enrichment analysis") #+ coord_flip()+
  #theme(axis.text.x = element_text(angle = 270, hjust=0))
dev.off()

df_gsea3 <- df_gsea2[,c("pathway","leadingEdge","condition")]

#########################################################################
# Look at the up and down regulation of associated factors

###
# LSC 
# polycomb
ranks <- ranks_LSC
pdf(paste0(resultsdir,'/gsea_PcG_LSC_es.pdf') ,width=7,height=4,paper='special')
plotGseaTable(pathways["Polycomb group (PcG) protein"], ranks, LSC_gsea, 
              gseaParam=0.5,)
dev.off()

LSC_up <- res_up_LSC[which(res_up_LSC %in% unlist(unname(pathways["Polycomb group (PcG) protein"])))]
LSC_down <-res_down_LSC[which(res_down_LSC %in% unlist(unname(pathways["Polycomb group (PcG) protein"])))]

PcG_up <- epig_factors[epig_factors$HGNC.approved.symbol%in%LSC_up,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]
PcG_down <- epig_factors[epig_factors$HGNC.approved.symbol%in%LSC_down,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]

# Histone mod erase co-factor
term <- "Histone modification erase cofactor, TF"
pdf(paste0(resultsdir,'/gsea_H_ER_TF_LSC_es.pdf') ,width=7,height=4,paper='special')
plotGseaTable(pathways[term], ranks, LSC_gsea, 
              gseaParam=0.5,)
dev.off()

LSC_up <-res_up_LSC[which(res_up_LSC %in% unlist(unname(pathways[term])))]
LSC_down <-res_down_LSC[which(res_down_LSC %in% unlist(unname(pathways[term])))]

HERTF_LSC_up <- epig_factors[epig_factors$HGNC.approved.symbol%in%LSC_up,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]
HERTF_LSC_down <- epig_factors[epig_factors$HGNC.approved.symbol%in%LSC_down,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]

# Histone mod erase co-factor
term <- "Histone chaperone"
pdf(paste0(resultsdir,'/gsea_H_CHAP_LSC_es.pdf') ,width=7,height=4,paper='special')
plotGseaTable(pathways[term], ranks, LSC_gsea, 
              gseaParam=0.5,)
dev.off()

LSC_up <-res_up_LSC[which(res_up_LSC %in% unlist(unname(pathways[term])))]
LSC_down <-res_down_LSC[which(res_down_LSC %in% unlist(unname(pathways[term])))]

chap_LSC_up <- epig_factors[epig_factors$HGNC.approved.symbol%in%LSC_up,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]
chap_LSC_down <- epig_factors[epig_factors$HGNC.approved.symbol%in%LSC_down,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]
###


###
# LE
term <- "Histone modification read"
ranks <- ranks_LE
pdf(paste0(resultsdir,'/gsea_H_READ_LE_es.pdf') ,width=7,height=4,paper='special')
plotGseaTable(pathways[term], ranks, LE_gsea, 
              gseaParam=0.5,)
dev.off()

LE_up <-res_up_LE[which(res_up_LE %in% unlist(unname(pathways[term])))]
LE_down <-res_down_LE[which(res_down_LE %in% unlist(unname(pathways[term])))]

READ_LE_up <- epig_factors[epig_factors$HGNC.approved.symbol%in%LE_up,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]
READ_LE_down <- epig_factors[epig_factors$HGNC.approved.symbol%in%LE_down,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]

term <- "Histone modification erase cofactor"
pdf(paste0(resultsdir,'/gsea_H_ERCO_LE_es.pdf') ,width=7,height=4,paper='special')
plotGseaTable(pathways[term], ranks, LE_gsea, 
              gseaParam=0.5,)
dev.off()

LE_up <-res_up_LE[which(res_up_LE %in% unlist(unname(pathways[term])))]
LE_down <-res_down_LE[which(res_down_LE %in% unlist(unname(pathways[term])))]

ERCO_LE_up <- epig_factors[epig_factors$HGNC.approved.symbol%in%LE_up,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]
ERCO_LE_down <- epig_factors[epig_factors$HGNC.approved.symbol%in%LE_down,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]

###
# CSSC
ranks <- ranks_CSSC
term <- "DNA modification, RNA modification"
pdf(paste0(resultsdir,'/gsea_DNA_RNA_CSSC_es.pdf') ,width=7,height=4,paper='special')
plotGseaTable(pathways[term], ranks, CSSC_gsea, 
              gseaParam=0.5,)
dev.off()

CSSC_up <-res_up_CSSC[which(res_up_CSSC %in% unlist(unname(pathways[term])))]
CSSC_down <-res_down_CSSC[which(res_down_CSSC %in% unlist(unname(pathways[term])))]

DR_CSSC_up <- epig_factors[epig_factors$HGNC.approved.symbol%in%CSSC_up,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]
DR_CSSC_down <- epig_factors[epig_factors$HGNC.approved.symbol%in%CSSC_down,c("HGNC.approved.symbol","Protein.complex","Target.entity","Product")]


# delete below?

# library(enrichplot)
# 
# gseaplot2(fgseaResTidy, geneSetID = 1, title = fgseaResTidy$ID[1],color="red",pvalue_table=T)
# 
# install.packages("remotes")
# remotes::install_github("AliSajid/BioPathNet")
# library(BioPathNet)
# 
# 
# 
# res_up <- res[res$log2FoldChange>0,]
# res_down <- res[res$log2FoldChange<0,]
# 
# res_up <- res_up$SYMBOL
# res_down <- res_down$SYMBOL
#   
# em <- GSEAResult(results=fgseaRes, pathways=pathways, lower=-40, upper=40, alpha=0.05, upreg=tbl_df(res_up), downreg=tbl_df(res_down))
# em@geneList <- ranks
# 
# gseMe3 <- gseGO(geneList     = geneList,
#                 OrgDb        = org.Mm.eg.db,
#                 keyType = 'SYMBOL',
#                 ont          = "BP",
#                 scoreType = "std",
#                 pvalueCutoff = 0.05,
#                 verbose      = FALSE)
# 
# ## feature 1: numeric vector
# geneList = res$log2FoldChange
# 
# ## feature 2: named vector
# names(geneList) = as.character(rownames(res))
# 
# ## feature 3: decreasing orde
# geneList = sort(geneList, decreasing = TRUE)
# head(geneList, 5)
# 
# gseMe3 <- gseGO(geneList     = geneList,
#       OrgDb        = pathways,
#       keyType = 'SYMBOL',
#       #ont          = "BP",
#       scoreType = "std",
#       pvalueCutoff = 0.05,
#       verbose      = FALSE)
# 
# test2 <- GSEA(geneList, exponent = 1, nPerm = 1000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE=pathways, TERM2NAME = NA, verbose = TRUE, seed = FALSE)
# test3<-do.call(rbind.data.frame, pathways)
# bind_rows(pathways)

df <- tibble::enframe(pathways) %>%
  dplyr::mutate(value = purrr::map_chr(value, toString))
df <- df %>% 
  mutate(value = strsplit(as.character(value), ",")) %>% 
  unnest(value)

names(df) <- c("gs_name","entrez_gene")
head(df)

res_new <- res
res_new$log <- -log10(res$pvalue)

#gl <- res_new$log
gl <- res_new$log2FoldChange

names(gl) <- res_new$SYMBOL
gl<- sort(gl, decreasing = TRUE)
library(enrichplot)

head(gl)
em <- GSEA(gl ,TERM2GENE = df,pvalueCutoff = 1,verbose = T)
em[1:5,1:5]
em <- GSEA(geneList=gl, exponent = 1, nPerm = 1000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE=df)
em2<- gseDO(gl)

vec_noms <- df[df$gs_name=="Polycomb group (PcG) protein",]
#vec_noms <- df[df$gs_name=="Histone chaperone",]

vec_noms <- vec_noms$entrez_gene

install.packages("liger")
library(liger)
gsea(values=gl,
     geneset=vec_noms,
     power = 1,
     rank = T,
     weight = rep(1, length(values)),
     n.rand = 10000,
     plot = TRUE,
     main = "",
     return.details = FALSE,
     quantile.threshold = min(100/n.rand, 0.1),
     #random.seed = 1,
     mc.cores = 1
)

data("org.Hs.GO2Symbol.list")  
universe <- unique(unlist(org.Hs.GO2Symbol.list))  # get universe
gs <- org.Hs.GO2Symbol.list[[1]]  # get a gene set
vals <- rnorm(length(universe), 0, 10)  
names(vals) <- universe
vals[gs] <- rnorm(length(gs), 100, 10)

gsea(values=gl, geneset=vec_noms, mc.cores=1, n.rand=100, main="Epi:Polycomb group (PcG) protein") 


head(vals)
str(gs)
gseDO(
  sort(ranks,decreasing = T),
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea")
em["results"] <- NULL

str(em["results"])
enricher(gl,
         TERM2GENE=df,
         pvalueCutoff = 1)

em

pdf(paste0(resultsdir,'/gsea_PcG_LSC.pdf') ,width=4,height=4.5,paper='special')
gseaplot2(em, geneSetID = names(em@geneSets[43]), title = names(em@geneSets[43]),color="green",pvalue_table=F)
dev.off()

library(DOSE)
data(geneList)
edo2 <- gseDO(geneList)
gseaplot2(edo2, geneSetID = 44, title = edo2$Description[44],color="green",pvalue_table=F)

#______________________VISUALIZATION______________________________#

#__________bar plot _______________#
# Plot the normalized enrichment scores. 
#Color the bar indicating whether or not the pathway was significant:
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")
###############################################################################
# Check the two first PCs of the corneal epithelial populations
seur_obj <-subset(x=seur_obj, subset = (scVI_label == "LSC"|scVI_label == "LESC"|scVI_label =="LE"| scVI_label =="CE"))

seur_obj <- FindVariableFeatures(seur_obj)
seur_obj <- ScaleData(seur_obj)
seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = 2)

mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
pca <- seur_obj[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

# Checking the first PC
Stdev(object = seur_obj[["pca"]])[1]
seur_obj$scVI_label <- droplevels(seur_obj$scVI_label)
seur_obj@active.ident <- seur_obj$scVI_label
pdf(paste0(resultsdir,'/scVI_4datasets_full.pdf') ,width=6,height=5,paper='special')
pc1_viz <-2
pc2_viz <-1
y_label = "PC2"#paste0(paste0(paste0("PC",pc2_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
x_label = "PC1"#paste0(paste0(paste0("PC",pc1_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
PC_dimred <- DimPlot(seur_obj, reduction = "pca",label = F, dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
print(PC_dimred)
dev.off()

# Save the subsetted object
seur_obj@active.ident <- droplevels(seur_obj@active.ident)
seur_obj$scVI_label<-droplevels(seur_obj$scVI_label)
saveRDS(seur_obj,paste0(resultsdir,"/4datasets_epi_subset_celnum.rds"))

################################################################################
# Import dataset for VIA pseudotime trajectory inference
# seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220308/4datasets_epi_subset_celnum.rds")
# 
# # 
# labels_cell<-as.character(names(seur_obj$orig.ident))
# labels_cluster<-as.character(seur_obj$scVI_label)
# 
# df <- data.frame(labels_cell,labels_cluster)
# names(df) <- c("cell_id","group_id")
# 
# write.csv(df,"/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220308/cornea_ids.csv", row.names = T)
# 
# counts.df <- seur_obj@assays$RNA@counts %>% as.matrix %>% t %>% as.data.frame
# write.csv(counts.df,"/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220308/cornea_counts.csv", row.names = T)


###############################################################################
# Generate barplot overview first (added to the matrix later)
newdf <- NULL

for (xx in unique(seur_obj$scVI_label)) {
  print(xx)
  for (j in unique(seur_obj$Condition)){
    val<-sum(seur_obj$scVI_label==xx & seur_obj$Condition==j)
    #col<- dfcol[rownames(dfcol)%in% xx,]
    df <- data.frame(values=val,Dataset=j,cell=xx)
    newdf <- rbind(newdf, df)
  }
}

pdf(paste(resultsdir,'Cell_dist_datasets_epi_full.pdf',sep="/") ,width=10,height=8,paper='special')
# Stacked barchart showing the composition
print(ggplot(newdf, aes(fill=Dataset, y=values, x=cell)) +
        geom_bar(position="stack", stat="identity") + labs(x="Cell type",y="Number of cells"))
dev.off()

# Re-order that the colors make sense for the bar_plot
newdf <- newdf[order(newdf$cell),]
#dfcol2 <- dfcol[order(rownames(dfcol)),]

pdf(paste(resultsdir,'Cell_dist_datasets_stacked_epi_full.pdf',sep="/") ,width=5,height=5,paper='special')
print(ggplot(newdf, aes(fill=cell, y=values, x=Dataset)) +
        geom_bar(position="fill", stat="identity", width = 0.3) + labs(x="Dataset",y="Cell proportions"))#+scale_fill_manual(values = dfcol2)
dev.off()
################################################################################
# Generate a z-score table to determine relative expression in-between cell 
# populations
# Splitting the pseudobulk datasets unbiased into two replicates for Z-score calculation and over pseudotime
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")
sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc$sample <- seur_obj$scVI_label

sce_qc$sample <- factor(sce_qc$sample)
sce_qc$sample

sce_qc$rep <- seur_obj@meta.data$orig.ident
sce_qc$rep

unique(seur_obj$scVI_label)

for (sample in unique(sce_qc$sample)){
  pseudobulk_df[ , paste0(sample)] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample])
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulkdf_full.tsv'), sep = '\t',quote = F, row.names = F)

pseudobulk_df <- NULL
pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

#sce_qc$reps <- "rep1"

# use rep2 for Catala and one subject from Li
#sce_qc$reps[sce_qc$rep == "GautamSRR11470710"| sce_qc$rep =="GautamSRR14742510"| sce_qc$rep =="GautamSRR14742511"] <- "rep2"
#sce_qc$reps[sce_qc$rep == "CatalaGSM5651509" | sce_qc$rep == "CatalaGSM5651511" | sce_qc$rep == "CatalaGSM5651513" | sce_qc$rep == "CatalaGSM56511515"| sce_qc$rep == "CatalaGSM5651117"| sce_qc$rep == "CatalaGSM5651119"] <- "rep2"
#sce_qc$reps[sce_qc$rep == "CatalaGSM5651510" | sce_qc$rep == "CatalaGSM5651512" | sce_qc$rep == "CatalaGSM5651514" | sce_qc$rep == "CatalaGSM56511516"| sce_qc$rep == "CatalaGSM5651118"|sce_qc$rep == "CatalaGSM5651120"] <- "rep2"

# putting values of artificial rep2 and rep1 into the pseudobulkdataframecolumns
for (sample in unique(sce_qc$sample)){
    #if(sample == "EC"|sample == "SK"|sample == "TSK"){
      #sce_qc$reps <- "rep2"
      
      # use rep2 for Catala and one subject from Li
      #sce_qc$reps[sce_qc$rep == "GautamSRR11470710"| sce_qc$rep =="GautamSRR14742510"| sce_qc$rep =="GautamSRR14742511"] <- "rep2"
      #sce_qc$reps[sce_qc$rep == "CatalaGSM5651510" | sce_qc$rep == "CatalaGSM5651512" | sce_qc$rep == "CatalaGSM5651513"| sce_qc$rep == "CatalaGSM5651514"] <- "rep2"
      #sce_qc$reps[sce_qc$rep == "CatalaGSM5651509" | sce_qc$rep == "CatalaGSM5651511"  | sce_qc$rep == "CatalaGSM56511515"| sce_qc$rep == "CatalaGSM5651117"| sce_qc$rep == "CatalaGSM5651119" | sce_qc$rep == "CatalaGSM56511516"| sce_qc$rep == "CatalaGSM5651118"|sce_qc$rep == "CatalaGSM5651120"] <- "rep1"  
      #sce_qc$reps <- transform(sce_qc$reps, group = sample(c("rep1", "rep2"), nrow(df),
           #                              replace = TRUE, prob = c(0.5, 0.5)))
      
      # sce_qc$reps <- sce_qc$reps %>%
      #   mutate(group = sample(c('rep1', 'rep2'), n(), 
      #                         replace = TRUE, prob = c(0.5, 0.5)))
     # pseudobulk_df[ , paste0(sample,"_1")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep1"])
     # pseudobulk_df[ , paste0(sample,"_2")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$rep == "rep2"])
    #}else{
      sce_qc$reps <- "rep1"
      
      # use rep2 for Catala and one subject from Li
      sce_qc$reps[sce_qc$rep == "GautamSRR11470710"| sce_qc$rep =="GautamSRR14742510"| sce_qc$rep =="GautamSRR14742511"] <- "rep2"
      sce_qc$reps[sce_qc$rep == "CatalaGSM5651509" | sce_qc$rep == "CatalaGSM5651511" | sce_qc$rep == "CatalaGSM5651513" | sce_qc$rep == "CatalaGSM56511515"| sce_qc$rep == "CatalaGSM5651117"| sce_qc$rep == "CatalaGSM5651119"] <- "rep2"
      sce_qc$reps[sce_qc$rep == "CatalaGSM5651510" | sce_qc$rep == "CatalaGSM5651512" | sce_qc$rep == "CatalaGSM5651514" | sce_qc$rep == "CatalaGSM56511516"| sce_qc$rep == "CatalaGSM5651118"|sce_qc$rep == "CatalaGSM5651120"] <- "rep2"
      
      pseudobulk_df[ , paste0(sample,"_1")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep1"])
      pseudobulk_df[ , paste0(sample,"_2")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep2"])
   # }
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulk_reps_DE_datasets_markers_split2.tsv'), sep = '\t',quote = F, row.names = F)

pseudobulk_df<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220309/pseudobulk_reps_DE_datasets_markers_split2.tsv',header = T,row.names = 1)

# For the deseq2 matrix vs ESC
lakocountfile <- pseudobulk_df
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score_datasets"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

lakovst <- lakocountfile

# Generate coldata dataframe
coldata <- NULL

# conditions
j <- unlist(strsplit(colnames(pseudobulk_df), split = "_"))
j <- j[seq(1,length(j),2)]

# Change _1 & _2 to condition 1 and _3 & _4 to condition 2
c <- c("1","2")
d <- paste(j,c,sep="_")

# reps
cols <- colnames(pseudobulk_df)

# type
xx <- rep("paired-end",length(cols))

coldata <- data.frame(cols=cols,condition=d,condition2=j,type=xx)

#coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/col2.tsv', sep = '\t', header = TRUE, row.names = 1)
#coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210802/col.tsv', sep = '\t', header = TRUE, row.names = 1)

# for complex heatmap
#coldata <- coldata[1:22,]
rownames(coldata)<- coldata$cols
coldata <- coldata[,c("condition","condition2","type")]

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# checking the data
head(lakovst,2)

coldata

# setting the rowdata in coldata similar to the coldata in the count matrix
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(lakovst))

all(rownames(coldata) == colnames(lakovst))

lakovst <- lakovst[, rownames(coldata)]
all(rownames(coldata) == colnames(lakovst))

dds <- DESeqDataSetFromMatrix(countData = round(lakovst),
                              colData = coldata,
                              design = ~ condition)
dds
###################################################################
# for complex heatmap
#dds2 <- DESeq(dds)

# transpose and scale the matrix per row Z-score per gene per sample is scaled upon
vsd <- assay(vst(dds,blind = T))
Z <- t(scale(t(vsd)))
Z
Z_score <- as.data.frame(Z,row.names = rownames(lakovst))

# joining the columns on the means sequential
n <- 2
Z_joined <- t(rowMeans(t(Z_score), as.integer(gl(ncol(Z_score), n, ncol(Z_score)))) / n)

df %>% mutate(mean_all = rowMeans(.),
              mean_sel = rowMeans(select(., select_vars)))

#  generate list of factors
vec <- rownames(coldata)
x <- split(vec, ceiling(seq_along(vec)/2))

scoretable <- as.data.frame(do.call(cbind, lapply(x, function(i) rowMeans(Z_score[, i]))), row.names = rownames(lakovst))
colnames(scoretable) <- unique(coldata$condition2)
write.table(scoretable, file = paste0(resultsdir, '/Zscoretable_markers_split2.tsv'),sep = "\t", quote = F,row.names = T,col.names = T)
scoretable1 <- scoretable

################################################################################
markers<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/cluster_markers_all.csv', sep = ',', header = TRUE, row.names = 1)
data_new2 <- markers %>%                                      # Top N highest values by group
  arrange(desc(avg_log2FC)) %>% 
  group_by(cluster) %>%
  slice(1:3) # now 3 first 5
data_new2    

#data_new2 <- data_new2[data_new2$cluster %in% unique(j),]
markers_vec <- unique(data_new2$gene_name_shorter)

scoretable1 <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score_datasets20220323/Zscoretable_markers_split2.tsv", sep="\t",row.names=1)
################################################################################
# Generate the complex heatmap of the z-score
pdf(paste(resultsdir,'heatmap_z-score_table_4_dataset.pdf',sep="/") ,width=5,height=11,paper='special')
f1 = colorRamp2(c(-2,0, 2), c("blue","white", "red"), space = "RGB")

newdf <- NULL

for (xx in unique(seur_obj$scVI_label)) {
  print(xx)
  for (j in unique(seur_obj$Condition)){
    val<-sum(seur_obj$scVI_label==xx & seur_obj$Condition==j)
    #col<- dfcol[rownames(dfcol)%in% xx,]
    df <- data.frame(values=val,Dataset=j,cell=xx)
    newdf <- rbind(newdf, df)
  }
}

# Order the rows and columns to make it more clear
cont_df <- as.data.frame.matrix(scoretable1)
# Re-order
cont_df <- cont_df[c("LSC","LESC","LE","CE","Cj","Mel","CSSC","CF","IC","Ves")]

cont_df<-na.omit(cont_df)

#cont_df <- cont_df[1:200,]
matz <- as.matrix(cont_df)
rownames(matz)<-rownames(cont_df)
vec3 <- markers_vec

# defining multiple vectors
newdf <- newdf[newdf$cell%in%colnames(cont_df),]

m3 <- newdf[newdf$Dataset == "Co",]$values
co <- m3
m3 <- newdf[newdf$Dataset == "Ca",]$values
ca <- m3
m3 <- newdf[newdf$Dataset == "Ga",]$values
ga <- m3
m3 <- newdf[newdf$Dataset == "Li",]$values
li <- m3

# Creating matrix
m <- cbind(co, ca, ga, li)
rownames(m) <- unique(newdf$cell)

m2 <- m/rowSums(m)

m2 <- m2[c("LSC","LESC","LE","CE","Cj","Mel","CSSC","CF","IC","Ves"),,drop=FALSE]

anno = anno_barplot(m2, gp = gpar(fill = 2:5), bar_width = 0.5, height = unit(1, "cm"))

ha = HeatmapAnnotation(Composition=anno,
                       show_legend = c("dataset" = TRUE))

# Barplots don't generate legends so need to do it manually
lgd_list = list(
  Legend(labels = c("Collin", "Catala","Gautam","Li"), title = "Dataset annotation", 
         legend_gp = gpar(fill = 2:5)))

# simplify the dendrogram
hc = hclust(dist(matz))
group = cutree(hc, k = 10)

ht1 = Heatmap(matz, col = f1, cluster_columns = F, name = "Relative expression",
              top_annotation  = ha, cluster_rows = cluster_within_group(t(matz), group), 
              row_split = 10, border = TRUE)

# annotate rownames
vec1 <- which(rownames(matz) %in% vec3, arr.ind = T)
vec3 <- rownames(matz[rownames(matz) %in% vec3 ,])

htlist1 <- ht1 + rowAnnotation(link = anno_mark(at =  vec1,labels = vec3))
draw(htlist1, 
     heatmap_legend_list = lgd_list)
dev.off()

################################################################################
# GSEA with epigenetic pathways

# Extract the gmx file
epig_factors <- read.table(file = "genes.txt",sep = '\t',header = T)

list <- NULL
extract_enrichment <-function(j,database) {
  list <- database[,database$j]$HGNC.approved.symbol
  return(list)
}

vec <-unique(epig_factors$Protein.complex)
sapply(X = vec,FUN = extract_enrichment(database=epig_factors))




################################################################################
# Perform single cell footprinting analysis

# Load in ATAC object from snapATAC analysis
library(Signac)
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("biovizBase")

library(biovizBase)

library(GenomicRanges)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("exomeCopy")
library(exomeCopy)

meta_all.atac <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/meta_all_atac_seurat.rds")

singlecellfootprint <- function(cell_types,motif,frag_path,cell_id_dir,hint_motif_path) {
  print(cell_types)
  
  # Read in the regions from footprint analysis
  df1 <- read.table(file = paste0(hint_motif_path,cell_types,"_mpbs.bed"))
  motif_sel <- motif
  df1 <- df1[df1$V4==motif_sel,]
  
  res <- sapply(split(df1, df1$V4), function(j){
    GRanges(seqnames = j$V1,
            ranges = IRanges(start = j$V2,
                             end = j$V3,
                             names = j$V4))
  })
  res <- res[[1]]
  
  # Add the ranges of the motif around the center -50 bp and up 50 bp
  start(res) <- start(res) - 42 
  end(res) <- end(res) + 42 
  res2<-subdivideGRanges(res,subsize=1)
  
  # Load in the fragment_data
  frag.path <- frag_path
  
  # Load in the corresponding cell IDs 
  cells <- read.csv(paste0(cell_id_dir,cell_types,"_atac_metadata_sample.csv"),header = F)
  cells <- cells$V1
  
  # Calculating the fragments for each position and fragment dataset
  fragments <- CreateFragmentObject(
    path = frag.path,
    cells = cells,
    validate.fragments = F
  )
  peak_matrix <- FeatureMatrix(
    fragments = fragments,
    features = res2,cells = cells
  )
  peak_mat2 <- as.matrix(peak_matrix)
  peak_mat2

  # Join positions from multiple regions around motifs to get reads per footprint
  range <- 1:101 #101
  peak_mat5 <- NULL

  calc_joined <- function(j) {
    peak_mat3<- peak_mat2[seq(j, nrow(peak_mat2), by=101),] # 101
    peak_mat4 <-colSums(peak_mat3)
    print(peak_mat4[130])
    peak_mat5 <- rbind(peak_mat5,peak_mat4)
    return(peak_mat)
  }
  peak_mat5 <-sapply(range, calc_joined)
  peak_mat5 <- t(peak_mat5)
  colnames(peak_mat5)<-colnames(peak_mat2)
  return(peak_mat5)
}

peak_mat2 <- d2_mat
hist(peak_mat2)



cell_type <- c("LSC","LESC")#,"LE","Cj","CE","CSSC","CF"
motif <-"GM.5.0.p53.0001"
cell_id_dir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE_11032022/ATAC_peaks/"
hint_motif_path <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching/"
frag_path1 <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/cellranger-atac-2.0.0/ataclako2/outs/fragments.tsv.gz"
frag_path2 <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/cellranger-atac-2.0.0/ataclako3/outs/fragments.tsv.gz"
frag_path3 <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/cellranger-atac-2.0.0/ataclako4/outs/fragments.tsv.gz"

d1_mat <- NULL
d2_mat <- NULL
d3_mat <- NULL

i<-"LSC"

for (i in cell_type){
  d1_mat <- NULL
  d2_mat <- NULL
  d3_mat <- NULL
  d1_mat <- singlecellfootprint(cell_types=i,motif=motif,cell_id_dir=cell_id_dir,frag_path=frag_path1,hint_motif_path=hint_motif_path)
  d2_mat <- singlecellfootprint(cell_types=i,motif=motif,cell_id_dir=cell_id_dir,frag_path=frag_path2,hint_motif_path=hint_motif_path)
  d3_mat <- singlecellfootprint(cell_types=i,motif=motif,cell_id_dir=cell_id_dir,frag_path=frag_path3,hint_motif_path=hint_motif_path)
  
  # Join matrices from frag_paths
  df0 <- as.data.frame(d1_mat)
  df2 <- as.data.frame(d2_mat)
  df3 <- as.data.frame(d3_mat)
  df4 <- cbind(df2 + df3)
  df4 <- df4[,colnames(df3)]
  df5 <- cbind(df4+df0)
  
  # Writing the colsums for all fragments for each cell population
  write.table(data.frame("pos"=rownames(df5),df5), file = paste0(resultsdir,'/',motif,cell_type,'.tsv'), sep = '\t',quote = F, row.names = F)
}

# Normalize the matrices by number of motifs
calc_motif <- function(i) {
  number <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching/",i,"_mpbs.bed"))
  sum(number$V4==motif)
}

# Print heatmaps with normalized counts
cells <- c("LSC","LESC")#,"LE","Cj","CE","CSSC","CF"
df2 <- sapply(cells, calc_motif)
vircols <- viridis(n=100)
fvir2 = colorRamp2(c(0,0.025,0.05,0.075,0.1), c(vircols[1],vircols[25],vircols[50],vircols[75],vircols[100]), space = "RGB")

pdf(paste0(paste(paste(resultsdir, sep = '/'), sep = '_'),motif,'.pdf'), width =6, height = 2)
for (i in cells){
  print(i)
  val <- unname(df2[i])
  mat<-read.table(paste0('/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220530','/',motif,i,'.tsv'),row.names = 1,header = T)
  mat <- mat/val
  print(Heatmap(t(as.matrix(mat[,])),cluster_columns = F,col = fvir2,show_row_names = F,name = paste0("counts/motif_num_",i)))
}
dev.off()

# 
# 
# i <- "LSC"
# 
# for (i in c("LSC","LESC","LE","Cj","CE","CSSC","CF")){
#   print(i)
#   df1 <- read.table(file = paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/hint_atac/MotifMatching/",i,"_mpbs.bed"))
#   motif_sel <- motif
#   df1 <- df1[df1$V4==motif_sel,]
#   
#  
#   res <- sapply(split(df1, df1$V4), function(j){
#   GRanges(seqnames = j$V1,
#           ranges = IRanges(start = j$V2,
#                            end = j$V3,
#                            names = j$V4))
#   })
# 
#   res <- res[[1]]
#   
#   # add the ranges of the motif around the center motif
#   start(res) <- start(res) - 42
#   end(res) <- end(res) + 42
# 
#   res2<-subdivideGRanges(res,subsize=1)
#   
#   frag.path <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/cellranger-atac-2.0.0/ataclako2/outs/fragments.tsv.gz"
#   
#   #meta_all.atac <- meta_all.atac[meta_all.atac@meta.data[1:3]]
#   cells <- read.csv(paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE_11032022/ATAC_peaks/",i,"_atac_metadata_sample.csv"),header = F)
#   cells <- cells$V1
#   
#   fragments <- CreateFragmentObject(
#   path = frag.path,
#   cells = cells,
#   validate.fragments = F
#   )
#   
#   peak_matrix <- FeatureMatrix(
#   fragments = fragments,
#   features = res2,cells = cells
#   
#   )
#   
#   peak_mat2 <- as.matrix(peak_matrix)
#   
#   #peak_mat2 <- peak_mat2[,colSums(peak_mat2)!=0] #rowSums(peak_mat2)!=0
#   
#   #First_pos
#   
#   range <- 1:101
#   peak_mat5 <- NULL
#   
#   calc_joined <- function(j) {
#   peak_mat3<- peak_mat2[seq(j, nrow(peak_mat2), by=101),]
#   peak_mat4 <-colSums(peak_mat3)
#   #rownames(peak_mat4) <- i
#   print(peak_mat4)
#   peak_mat5 <- rbind(peak_mat5,peak_mat4)
#   }
#   
#   peak_mat5 <-t(sapply(range, calc_joined))
#   colnames(peak_mat5)<-colnames(peak_mat2)
#   
#   #Heatmap(t(as.matrix(peak_mat5)),cluster_columns = F,col = viridis(n=100))
#   
#   d1_mat <- peak_mat5
#   
#   
#   frag.path <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/cellranger-atac-2.0.0/ataclako3/outs/fragments.tsv.gz"
#   
#   #meta_all.atac <- meta_all.atac[meta_all.atac@meta.data[1:3]]
#   # cells <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE_11032022/ATAC_peaks/LSC_atac_metadata_sample.csv")
#   # cells <- cells$ACTACGATCCATATCT.1
#   
#   fragments <- CreateFragmentObject(
#   path = frag.path,
#   cells = cells,
#   validate.fragments = F
#   )
#   
#   peak_matrix <- FeatureMatrix(
#   fragments = fragments,
#   features = res2,cells = cells
#   
#   )
#   
#   peak_mat2 <- as.matrix(peak_matrix)
#   
#   #peak_mat2 <- peak_mat2[,colSums(peak_mat2)!=0] #rowSums(peak_mat2)!=0
#   
#   #First_pos
#   
#   range <- 1:101
#   peak_mat5 <- NULL
#   calc_joined <- function(j) {
#   peak_mat3<- peak_mat2[seq(j, nrow(peak_mat2), by=101),]
#   peak_mat4 <-colSums(peak_mat3)
#   #rownames(peak_mat4) <- i
#   print(peak_mat4)
#   peak_mat5 <- rbind(peak_mat5,peak_mat4)
#   }
#   
#   peak_mat5 <-t(sapply(range, calc_joined))
#   colnames(peak_mat5)<-colnames(peak_mat2)
#   
#   #Heatmap(t(as.matrix(peak_mat5)),cluster_columns = F,col = viridis(n=100))
#   
#   d2_mat <- peak_mat5
#   
#   
#   
#   
#   frag.path <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/cellranger-atac-2.0.0/ataclako4/outs/fragments.tsv.gz"
#   
#   #meta_all.atac <- meta_all.atac[meta_all.atac@meta.data[1:3]]
#   # cells <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE_11032022/ATAC_peaks/LSC_atac_metadata_sample.csv")
#   # cells <- cells$ACTACGATCCATATCT.1
#   
#   fragments <- CreateFragmentObject(
#   path = frag.path,
#   cells = cells,
#   validate.fragments = F
#   )
#   
#   peak_matrix <- FeatureMatrix(
#   fragments = fragments,
#   features = res2,cells = cells
#   
#   )
#   
#   peak_mat2 <- as.matrix(peak_matrix)
#   
#   #peak_mat2 <- peak_mat2[,colSums(peak_mat2)!=0] #rowSums(peak_mat2)!=0
#   
#   #First_pos
#   
#   range <- 1:101
#   peak_mat5 <- NULL
#   calc_joined <- function(j) {
#   peak_mat3<- peak_mat2[seq(j, nrow(peak_mat2), by=101),]
#   peak_mat4 <-colSums(peak_mat3)
#   #rownames(peak_mat4) <- i
#   print(peak_mat4)
#   peak_mat5 <- rbind(peak_mat5,peak_mat4)
#   }
#   
#   peak_mat5 <-t(sapply(range, calc_joined))
#   colnames(peak_mat5)<-colnames(peak_mat2)
#   
#   #Heatmap(t(as.matrix(peak_mat5)),cluster_columns = F,col = viridis(n=100))
#   
#   d3_mat <- peak_mat5
#   
#   
#   
#   #$cbind(d2_mat[1], d2_mat[-1] + d3_mat[-1])
#   df0 <- as.data.frame(d1_mat)
#   df2 <- as.data.frame(d2_mat)
#   df3 <- as.data.frame(d3_mat)
#   
#   df4 <- cbind(df2 + df3)
#   
#   
#   df4 <- df4[,colnames(df3)]
#   
#   df5 <- cbind(df4+df0)
#   
#   Heatmap(t(as.matrix(df5)),cluster_columns = F,col = viridis(n=100),show_row_names = F)
#   
#   write.table(data.frame("pos"=rownames(df5),df5), file = paste0(resultsdir,'/p53_single',i,'.tsv'), sep = '\t',quote = F, row.names = F)
# 
# }



mat/100
############# DELETE BELOW?
library("magrittr")
df3 %>% select(-type) %>%
  add(df2 %>% select(-type)) %>%
  mutate(type = df3$type)




peak_mat3<- peak_mat2[seq(i, 15272, by=17),]
peak_mat4 <-t(as.data.frame(colSums(peak_mat3)))
rownames(peak_mat4) <- i


peak_mat3 <- peak_mat2[]
# not run
total_fragments <- CountFragments(frag.path)
fragments_test <- total_fragments[sum(total_fragments$CB==c("ACTACGATCCATATCT-1")), "frequency_count"]

pbmc <- FRiP(
  object = meta_all.atac,
  assay = peak_mat2,
  total.fragments = fragments_test
)





res <- res$GM.5.0.p53.0001
DefaultAssay(meta_all.atac) <- 'peaks'





# ATAC analysis add gene annotation information
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(meta_all.atac) <- annotations

# We exclude the first dimension as this is typically correlated with sequencing depth
meta_all.atac <- RunTFIDF(meta_all.atac)
meta_all.atac <- FindTopFeatures(meta_all.atac, min.cutoff = "q0")
meta_all.atac <- RunSVD(meta_all.atac)
meta_all.atac <- RunUMAP(meta_all.atac, reduction = "lsi", dims = 2:10, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#meta_all.atac < FindVariableFeatures(meta_all.atac)
#meta_all.atac <- RunPCA(meta_all.atac)
#meta_all.atac <- RunUMAP(meta_all.atac, dims = 1:10)
DimPlot(meta_all.atac, reduction = 'umap.atac', group.by = 'predict.id', label=F,label.size = 8,cols = colmat$V2)



# Generate new cluster labeling:
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)

cluster_order <- c('Cj','CE','LE','LESC','CSSC','CF','LSC',
                   'SK','LE','Cj','SK','TSK','CE','CF',
                   'Mel','EC','Ves','IC','CE','MF','CDH19+')
colmat <- matrix(c(z2, newcol), ncol = 2)
colmat <- as.data.frame(colmat)

new.cluster.ids <- c('Cj','CE','LE','LESC','CSSC','CF','LSC',
                     'SK','LE','Cj','SK','TSK','CE','CF',
                     'Mel','EC','Ves','IC','CE','MF','CDH19+')

colmat$V1 <- plyr::mapvalues(x = as.factor(colmat$V1), from = current.cluster.ids, to = new.cluster.ids)
colmat$V1 <- as.factor(colmat$V1)
colmat <- colmat[!duplicated(colmat$V1),]

meta_all.atac$predict.id <- as.factor(meta_all.atac$predict.id)

#colVector <- cols$V2
colmat<-colmat[order(match(colmat$V1,levels(meta_all.atac$predict.id) )), , drop = FALSE]
#cols$V1 <- as.vector(cols$V1)
#cols <- cols[!duplicated(cols$V1),]

#,cols = cols$V2




# Downstream analysis from the two datasets?
# ...

# Include all corneal epithelial, limbal and conjunctival clusters from the full cornea
# Add the other two datasets to the seurat object 

















































# query2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220207/Li_singlets.rds")
# #query3 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220207/....rds")
# 
# # generate training data for the machine learning model for every group
# 
# # Extract labels
# labels<-as.character(reference$costum_clustering)
# write.csv(labels,paste0(resultsdir,"/training_labels_collins.csv"),row.names=F)
# 
# # make counts file readable in python, convert to h5ad, repeat for testing & training data
# test_dat<-CreateSeuratObject(counts=reference@assays$RNA@counts)
# saveRDS(test_dat,paste0(resultsdir,"/test_dat_collins.rds"))
# 
# # Extract labels
# labels<-as.character(query$costum_clustering)
# write.csv(labels,paste0(resultsdir,"/training_labels_catala.csv"),row.names=F)
# 
# # make counts file readable in python, convert to h5ad, repeat for testing & training data
# test_dat<-CreateSeuratObject(counts=query@assays$RNA@counts)
# saveRDS(test_dat,paste0(resultsdir,"/test_dat_catala.rds"))
# 
# # Extract labels
# labels<-as.character(query2$costum_clustering)
# write.csv(labels,paste0(resultsdir,"/training_labels_li.csv"),row.names=F)
# 
# # make counts file readable in python, convert to h5ad, repeat for testing & training data
# test_dat<-CreateSeuratObject(counts=query2@assays$RNA@counts)
# saveRDS(test_dat,paste0(resultsdir,"/test_dat_li.rds"))

# Extract labels
#labels<-as.character(query3$costum_clustering)
#write.csv(labels,paste0(resultsdir,"/training_labels_gautam.csv"),row.names=F)

# make counts file readable in python, convert to h5ad, repeat for testing & training data
#test_dat<-CreateSeuratObject(counts=query3@assays$RNA@counts)
#saveRDS(test_dat,paste0(resultsdir,"test_dat_gautam.rds"))

#### Create 5- fold cross-validation data ####
#install.packages("rBayesianOptimization") # only run on first install
# library(rBayesianOptimization)
# 
# source("Cross_Validation.R")
# Cross_Validation("training_labels_collins.csv", OutputDir = paste0(resultsdir,"/collin_cv/"))
# Cross_Validation("training_labels_catala.csv", OutputDir = paste0(resultsdir,"/catala_cv/"))
# Cross_Validation("training_labels_li.csv", OutputDir = paste0(resultsdir,"/li_cv/"))
# #Cross_Validation("training_labels_li.csv", OutputDir = paste0(resultsdir,"gautam_cv/"))
# 
# DefaultAssay(reference) <- "RNA"
# DefaultAssay(query) <- "RNA"
# DefaultAssay(query2) <- "RNA"
# pdf(paste0(paste(paste(resultsdir, sep = '/'), sep = '_'),'test_antiparallel_elf3_TGFBi.pdf'), width =5, height = 5)
# 
# print(FeaturePlot(reference, features = "ELF3"))
# print(FeaturePlot(reference, features = "TGFBI"))
# print(FeaturePlot(query, features = "ELF3"))
# print(FeaturePlot(query, features = "TGFBI"))
# print(FeaturePlot(query2, features = "ELF3"))
# print(FeaturePlot(query2, features = "TGFBI"))
# dev.off()
# 
# ##############################################################################
# # predict cell types based on the machine learning algorithm
# catala_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/collins_catala_SVM/SVM_Pred_Labels.csv")
# query$ml_labels_collin <- catala_labels$X0
# catala_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/collins_catala_SVMrej/SVMrej_Pred_Labels.csv")
# query$ml_labels_rej_collin <- catala_labels$X0
# 
# pdf(paste(resultsdir, "catala_labels_collin_ML.pdf", sep = '/'), width = 6, height = 5)
# cols1 <- col_mat <- viridis(14,alpha = 0.5,option = "H",begin = 0.1)
# p1 <- DimPlot(query, label=F,label.size = 8, group.by = 'ml_labels_rej_collin',cols=cols1) + ggtitle("Louvain Clustering")
# p2 <- DimPlot(query, label=F,label.size = 8, group.by = 'ml_labels_collin') + ggtitle("Louvain Clustering")
# p1
# p2
# dev.off()
# 
# catala_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/li_catala_SVM/SVM_Pred_Labels.csv")
# query$ml_labels_li <- catala_labels$X0
# catala_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/li_catala_SVMrej/SVMrej_Pred_Labels.csv")
# query$ml_labels_rej_li <- catala_labels$X0
# 
# pdf(paste(resultsdir, "catala_labels_li_ML.pdf", sep = '/'), width = 6, height = 5)
# cols1 <- col_mat <- viridis(14,alpha = 0.5,option = "H")
# p1 <- DimPlot(query, label=F,label.size = 8, group.by = 'ml_labels_rej_li',cols=cols1) + ggtitle("Louvain Clustering")
# p2 <- DimPlot(query, label=F,label.size = 8, group.by = 'ml_labels_li') + ggtitle("Louvain Clustering")
# p1
# p2
# dev.off()
# 
# # predict cell types based on the machine learning algorithm
# li_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/collins_li_SVM/SVM_Pred_Labels.csv")
# query2$ml_labels_collin <- li_labels$X0
# li_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/collins_li_SVMrej/SVMrej_Pred_Labels.csv")
# query2$ml_labels_rej_collin <- li_labels$X0
# 
# pdf(paste(resultsdir, "li_labels_collins_ML.pdf", sep = '/'), width = 6, height = 5)
# cols1 <- col_mat <- viridis(14,alpha = 0.5,option = "H",begin = 0.1)
# p1 <- DimPlot(query2, label=F,label.size = 8, group.by = 'ml_labels_rej_collin',cols=cols1) + ggtitle("Louvain Clustering")
# p2 <- DimPlot(query2, label=F,label.size = 8, group.by = 'ml_labels_collin') + ggtitle("Louvain Clustering")
# p1
# p2
# dev.off()
# 
# # predict cell types based on the machine learning algorithm
# li_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/catala_li_SVM/SVM_Pred_Labels.csv")
# query2$ml_labels_catala <- li_labels$X0
# li_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/catala_li_SVMrej/SVMrej_Pred_Labels.csv")
# query2$ml_labels_rej_catala <- li_labels$X0
# 
# pdf(paste(resultsdir, "li_labels_catala_ML.pdf", sep = '/'), width = 6, height = 5)
# cols1 <- col_mat <- viridis(14,alpha = 0.5,option = "H",begin = 0.1)
# p1 <- DimPlot(query2, label=F,label.size = 8, group.by = 'ml_labels_rej_catala',cols=cols1) + ggtitle("Louvain Clustering")
# p2 <- DimPlot(query2, label=F,label.size = 8, group.by = 'ml_labels_catala') + ggtitle("Louvain Clustering")
# p1
# p2
# dev.off()
# 
# # predict cell types based on the machine learning algorithm
# collin_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/catala_collins_SVM/SVM_Pred_Labels.csv")
# reference$ml_labels_catala <- collin_labels$X0
# collin_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/catala_collins_SVMrej/SVMrej_Pred_Labels.csv")
# reference$ml_labels_rej_catala <- collin_labels$X0
# 
# pdf(paste(resultsdir, "collins_labels_catala_ML.pdf", sep = '/'), width = 6, height = 5)
# cols1 <- col_mat <- viridis(14,alpha = 0.5,option = "H",begin = 0.1)
# p1 <- DimPlot(reference, label=F,label.size = 8, group.by = 'ml_labels_rej_catala',cols=cols1) + ggtitle("Louvain Clustering")
# p2 <- DimPlot(reference, label=F,label.size = 8, group.by = 'ml_labels_catala') + ggtitle("Louvain Clustering")
# p1
# p2
# dev.off()
# 
# # predict cell types based on the machine learning algorithm
# collin_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/li_collins_SVM/SVM_Pred_Labels.csv")
# reference$ml_labels_li <- collin_labels$X0
# collin_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/li_collins_SVMrej/SVMrej_Pred_Labels.csv")
# reference$ml_labels_rej_li <- collin_labels$X0
# 
# pdf(paste(resultsdir, "collins_labels_li_ML.pdf", sep = '/'), width = 6, height = 5)
# cols1 <- col_mat <- viridis(14,alpha = 0.5,option = "H",begin = 0.1)
# p1 <- DimPlot(reference, label=F,label.size = 8, group.by = 'ml_labels_rej_li',cols=cols1) + ggtitle("Louvain Clustering")
# p2 <- DimPlot(reference, label=F,label.size = 8, group.by = 'ml_labels_li') + ggtitle("Louvain Clustering")
# p1
# p2
# dev.off()
# 
# # set ml labels from own clustering
# reference$ml_labels_collin <- as.character(droplevels(reference$costum_clustering))
# reference$ml_labels_rej_collin <- as.character(droplevels(reference$costum_clustering))
# 
# query$ml_labels_catala <- as.character(droplevels(query$costum_clustering))
# query$ml_labels_rej_catala <- as.character(droplevels(query$costum_clustering))
# 
# query2$ml_labels_li <- as.character(droplevels(query2$costum_clustering))
# query2$ml_labels_rej_li <- as.character(droplevels(query2$costum_clustering))
# 
# pdf(paste(resultsdir, "Orig_ID.pdf", sep = '/'), width = 6, height = 5)
# seur_obj <- reference
# p2 <- DimPlot(seur_obj, label=F,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p2
# seur_obj <- query
# p2 <- DimPlot(seur_obj, label=F,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p2
# seur_obj <- query2
# p2 <- DimPlot(seur_obj, label=F,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p2
# dev.off()

# Predict cell types of other cell populations based on Collin with seurat
# # perform standard preprocessing on each object
# reference <- NormalizeData(reference)
# reference <- FindVariableFeatures(reference)
# reference <- ScaleData(reference)
# 
# query <- NormalizeData(query)
# query <- FindVariableFeatures(query)
# query <- ScaleData(query)
# 
# # find anchors
# anchors <- FindTransferAnchors(reference = reference, query = query)
# 
# # transfer labels
# predictions <- TransferData(anchorset = anchors, refdata = reference$costum_clustering)
# query <- AddMetaData(object = query, metadata = predictions)
# 
# pdf(paste(resultsdir, "Query_anno_reference.pdf", sep = '/'), width = 6, height = 5)
# seur_obj <- query
# p2 <- DimPlot(seur_obj, label=F,label.size = 8, group.by = 'predicted.id') + ggtitle("Louvain Clustering")
# p2
# dev.off()
# 
# query2 <- NormalizeData(query2)
# query2 <- FindVariableFeatures(query2)
# query2 <- ScaleData(query2)
# 
# # find anchors
# anchors <- FindTransferAnchors(reference = reference, query = query2)
# 
# # transfer labels
# predictions <- TransferData(anchorset = anchors, refdata = reference$costum_clustering)
# query2 <- AddMetaData(object = query2, metadata = predictions)
# 
# pdf(paste(resultsdir, "Query2_anno_reference.pdf", sep = '/'), width = 8, height = 8)
# seur_obj <- query2
# p2 <- DimPlot(seur_obj, label=F,label.size = 8, group.by = 'predicted.id') + ggtitle("Louvain Clustering")
# p2
# dev.off()

# ##############################################################################
# # Integration
# your_list <- list(reference, query, query2)
# soref <- SCTransform(reference, verbose = FALSE)
# soquery <- SCTransform(query, verbose = FALSE)
# soquery2 <- SCTransform(query2, verbose = FALSE)
# 
# ########################
# ### Integrate
# ########################
# dtlist <- list(reference = soref, query = soquery, query2 = soquery2)
# intfts <- SelectIntegrationFeatures(object.list = dtlist, nfeatures = 2000) # maxes out at 4080 (why?)
# dtlist <- PrepSCTIntegration(object.list = dtlist,
#                              anchor.features = intfts)
# anchors <- FindIntegrationAnchors(object.list = dtlist, normalization.method = "SCT",
#                                   anchor.features = intfts,reduction = "cca")
# integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# 
# # test it how it works
# integrated <- RunPCA(integrated)
# 
# integrated <- RunUMAP(integrated, dims = 1:30)
# 
# cornea.combined <- integrated
# 
# # add the ml labels again
# vec1 <- c(reference$ml_labels_collin,query$ml_labels_collin,query2$ml_labels_collin)
# vec2 <- c(reference$ml_labels_rej_collin,query$ml_labels_rej_collin,query2$ml_labels_rej_collin)
# cornea.combined$ml_full_collin <- as.factor(vec1)
# cornea.combined$ml_full_rej_collin <- as.factor(vec2)
# 
# vec1 <- c(reference$ml_labels_catala,query$ml_labels_catala,query2$ml_labels_catala)
# vec2 <- c(reference$ml_labels_rej_catala,query$ml_labels_rej_catala,query2$ml_labels_rej_catala)
# cornea.combined$ml_full_catala <- as.factor(vec1)
# cornea.combined$ml_full_rej_catala <- as.factor(vec2)
# 
# vec1 <- c(reference$ml_labels_li,query$ml_labels_li,query2$ml_labels_li)
# vec2 <- c(reference$ml_labels_rej_li,query$ml_labels_rej_li,query2$ml_labels_rej_li)
# cornea.combined$ml_full_li <- as.factor(vec1)
# cornea.combined$ml_full_rej_li <- as.factor(vec2)
# 
# # Plot the tree and the final clustering
# pdf(paste(resultsdir, "clustering_costum_transfered_cca_SCT_30dims.pdf", sep = '/'), width = 6, height = 5)
# #p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
# p0 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_collin', label=F,label.size = 8)
# p1 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_collin', label=F,label.size = 8)
# p01 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_catala', label=F,label.size = 8)
# p11 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_li', label=F,label.size = 8)
# p02 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_li', label=F,label.size = 8)
# p12 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_catala', label=F,label.size = 8)
# p2 <- DimPlot(cornea.combined, label=F,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p3 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
# p4 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
# #ape::plot.phylo(x = data.tree, direction = "downwards",)
# #
# p0
# p1
# p01
# p12
# p11
# p02
# p2
# p3
# p4
# #print(multiplot(p1, p2, p3, cols = 3))
# dev.off()
# 
# saveRDS(cornea.combined, file = paste0(resultsdir,"/cornea_integrated_SCT_unannotated_cca_30dims_ML.rds"))
# 
# ################################################################################
# # SCT normalization with rpca integration
# dtlist2 <- lapply(X = dtlist, FUN = RunPCA, features = intfts)
# anchors <- FindIntegrationAnchors(object.list = dtlist2, normalization.method = "SCT",
#                                   anchor.features = intfts,reduction = "rpca")
# integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# integrated <- RunPCA(integrated)
# integrated <- RunUMAP(integrated, dims = 1:30)
# 
# DimPlot(integrated, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
# 
# cornea.combined <- integrated
# 
# # add the ml labels again
# vec1 <- c(reference$ml_labels_collin,query$ml_labels_collin,query2$ml_labels_collin)
# vec2 <- c(reference$ml_labels_rej_collin,query$ml_labels_rej_collin,query2$ml_labels_rej_collin)
# cornea.combined$ml_full_collin <- as.factor(vec1)
# cornea.combined$ml_full_rej_collin <- as.factor(vec2)
# 
# vec1 <- c(reference$ml_labels_catala,query$ml_labels_catala,query2$ml_labels_catala)
# vec2 <- c(reference$ml_labels_rej_catala,query$ml_labels_rej_catala,query2$ml_labels_rej_catala)
# cornea.combined$ml_full_catala <- as.factor(vec1)
# cornea.combined$ml_full_rej_catala <- as.factor(vec2)
# 
# vec1 <- c(reference$ml_labels_li,query$ml_labels_li,query2$ml_labels_li)
# vec2 <- c(reference$ml_labels_rej_li,query$ml_labels_rej_li,query2$ml_labels_rej_li)
# cornea.combined$ml_full_li <- as.factor(vec1)
# cornea.combined$ml_full_rej_li <- as.factor(vec2)
# 
# # Plot the tree and the final clustering
# pdf(paste(resultsdir, "clustering_costum_transfered_rPCA_SCT_30dims.pdf", sep = '/'), width = 6, height = 5)
# #p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
# p0 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_collin', label=F,label.size = 8)
# p1 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_collin', label=F,label.size = 8)
# p01 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_catala', label=F,label.size = 8)
# p11 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_li', label=F,label.size = 8)
# p02 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_li', label=F,label.size = 8)
# p12 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_catala', label=F,label.size = 8)
# p2 <- DimPlot(cornea.combined, label=F,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p3 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
# p4 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
# #ape::plot.phylo(x = data.tree, direction = "downwards",)
# #
# p0
# p1
# p01
# p12
# p11
# p02
# p2
# p3
# p4
# #print(multiplot(p1, p2, p3, cols = 3))
# dev.off()
# 
# saveRDS(cornea.combined, file = paste0(resultsdir,"/cornea_integrated_SCT_unannotated_rPCA_30dims_ML.rds"))
# 
# ###############################################################################
# # non-SCT integration (lognormalized)
# 
# your_list <- list(reference, query, query2)
# 
# SnowParam(workers = 2)
# your_list <- lapply(X = your_list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method ="vst", nfeatures = 2000)
# })
# 
# SnowParam(workers = 2)
# cornea.anchors <- FindIntegrationAnchors(object.list = your_list, dims = 1:30,reduction = 'cca')
# 
# SnowParam(workers = 2)
# cornea.combined <- IntegrateData(anchorset = cornea.anchors, dims = 1:30)
# 
# # original unmodified data still resides in the 'RNA' assay
# DefaultAssay(cornea.combined) <- "RNA"
# 
# # Run the standard workflow for visualization and clustering
# cornea.combined <- ScaleData(cornea.combined, verbose = FALSE)
# cornea.combined <- FindVariableFeatures(cornea.combined)
# cornea.combined <- RunPCA(cornea.combined, verbose = FALSE)
# cornea.combined <- RunUMAP(cornea.combined, reduction = "pca", dims = 1:30)
# 
# # add the ml labels again
# vec1 <- c(reference$ml_labels_collin,query$ml_labels_collin,query2$ml_labels_collin)
# vec2 <- c(reference$ml_labels_rej_collin,query$ml_labels_rej_collin,query2$ml_labels_rej_collin)
# cornea.combined$ml_full_collin <- as.factor(vec1)
# cornea.combined$ml_full_rej_collin <- as.factor(vec2)
# 
# vec1 <- c(reference$ml_labels_catala,query$ml_labels_catala,query2$ml_labels_catala)
# vec2 <- c(reference$ml_labels_rej_catala,query$ml_labels_rej_catala,query2$ml_labels_rej_catala)
# cornea.combined$ml_full_catala <- as.factor(vec1)
# cornea.combined$ml_full_rej_catala <- as.factor(vec2)
# 
# vec1 <- c(reference$ml_labels_li,query$ml_labels_li,query2$ml_labels_li)
# vec2 <- c(reference$ml_labels_rej_li,query$ml_labels_rej_li,query2$ml_labels_rej_li)
# cornea.combined$ml_full_li <- as.factor(vec1)
# cornea.combined$ml_full_rej_li <- as.factor(vec2)
# 
# # Plot the tree and the final clustering
# pdf(paste(resultsdir, "clustering_costum_transfered_cca_normalized_30dims.pdf", sep = '/'), width = 6, height = 5)
# #p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
# p0 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_collin', label=F,label.size = 8)
# p1 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_collin', label=F,label.size = 8)
# p01 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_catala', label=F,label.size = 8)
# p11 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_li', label=F,label.size = 8)
# p02 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_li', label=F,label.size = 8)
# p12 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_catala', label=F,label.size = 8)
# p2 <- DimPlot(cornea.combined, label=F,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p3 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
# p4 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
# #ape::plot.phylo(x = data.tree, direction = "downwards",)
# #
# p0
# p1
# p01
# p12
# p11
# p02
# p2
# p3
# p4
# #print(multiplot(p1, p2, p3, cols = 3))
# dev.off()
# 
# saveRDS(cornea.combined, file = paste0(resultsdir,"/cornea_integrated_normalized_unannotated_cca_30dims_ML.rds"))
# 
# 
# ###########################################################################
# # Log normalized rpca integration
# SnowParam(workers = 2)
# cornea.anchors <- FindIntegrationAnchors(object.list = your_list, dims = 1:30,reduction = 'rpca')
# 
# SnowParam(workers = 2)
# cornea.combined <- IntegrateData(anchorset = cornea.anchors, dims = 1:30)
# 
# # original unmodified data still resides in the 'RNA' assay
# DefaultAssay(cornea.combined) <- "integrated"
# 
# # Run the standard workflow for visualization and clustering
# cornea.combined <- ScaleData(cornea.combined, verbose = FALSE)
# cornea.combined <- RunPCA(cornea.combined, verbose = FALSE)
# cornea.combined <- RunUMAP(cornea.combined, reduction = "pca", dims = 1:30)
# 
# # add the ml labels again
# vec1 <- c(reference$ml_labels_collin,query$ml_labels_collin,query2$ml_labels_collin)
# vec2 <- c(reference$ml_labels_rej_collin,query$ml_labels_rej_collin,query2$ml_labels_rej_collin)
# cornea.combined$ml_full_collin <- as.factor(vec1)
# cornea.combined$ml_full_rej_collin <- as.factor(vec2)
# 
# vec1 <- c(reference$ml_labels_catala,query$ml_labels_catala,query2$ml_labels_catala)
# vec2 <- c(reference$ml_labels_rej_catala,query$ml_labels_rej_catala,query2$ml_labels_rej_catala)
# cornea.combined$ml_full_catala <- as.factor(vec1)
# cornea.combined$ml_full_rej_catala <- as.factor(vec2)
# 
# vec1 <- c(reference$ml_labels_li,query$ml_labels_li,query2$ml_labels_li)
# vec2 <- c(reference$ml_labels_rej_li,query$ml_labels_rej_li,query2$ml_labels_rej_li)
# cornea.combined$ml_full_li <- as.factor(vec1)
# cornea.combined$ml_full_rej_li <- as.factor(vec2)
# 
# # Plot the tree and the final clustering
# pdf(paste(resultsdir, "clustering_costum_transfered_rPCA_normalized_30dims.pdf", sep = '/'), width = 6, height = 5)
# #p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
# p0 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_collin', label=F,label.size = 8)
# p1 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_collin', label=F,label.size = 8)
# p01 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_catala', label=F,label.size = 8)
# p11 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_li', label=F,label.size = 8)
# p02 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_rej_li', label=F,label.size = 8)
# p12 <- DimPlot(cornea.combined, reduction = 'umap', group.by = 'ml_full_catala', label=F,label.size = 8)
# p2 <- DimPlot(cornea.combined, label=F,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p3 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
# p4 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
# #ape::plot.phylo(x = data.tree, direction = "downwards",)
# #
# p0
# p1
# p01
# p12
# p11
# p02
# p2
# p3
# p4
# #print(multiplot(p1, p2, p3, cols = 3))
# dev.off()
# 
# saveRDS(cornea.combined, file = paste0(resultsdir,"/cornea_integrated_normalized_unannotated_rPCA_30dims_ML.rds"))
# 
# ################################################################################
# # Check the first two PCs subsetted Collins
# cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20220126/cornea_integrated_normalized_unannotated_rPCA_30dims_ML.rds")
# 
# cells.to.subset <- c('Co_C1_SCE','Co_C10_BCE','Co_C2_CjS','Co_C3_CjB','Co_C5_LNCP','Co_C6_LPC')
# seur_obj <- cornea.combined
# seur_obj <-subset(x=seur_obj, subset = (ml_full_collin == "Co_C1_SCE"|ml_full_collin == "Co_C10_BCE"|ml_full_collin == "Co_C2_CjS"|ml_full_collin == "Co_C3_CjB"|
#                                           ml_full_collin == "Co_C5_LNCP"|ml_full_collin == "Co_C6_LPC"))
# 
# # Check clusters on  the first two PCs
# seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = 2)
# mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
# pca <- seur_obj[["pca"]]
# 
# # Get the total variance:
# total_variance <- sum(matrixStats::rowVars(mat))
# eigValues = (pca@stdev)^2  ## EigenValues
# varExplained = eigValues / total_variance
# 
# # Checking the first PC
# Stdev(object = seur_obj[["pca"]])[1]
# 
# Idents(seur_obj) <- "ml_full_collin"
# pdf(paste0(resultsdir,'/PCplot_ML_non-omitted.pdf') ,width=6,height=5,paper='special')
# pc1_viz <-2
# pc2_viz <- 1
# y_label = "PC2"#paste0(paste0(paste0("PC",pc2_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
# x_label = "PC1"#paste0(paste0(paste0("PC",pc1_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
# PC_dimred <- DimPlot(seur_obj, reduction = "pca",label = F, dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
# print(PC_dimred)
# Idents(seur_obj) <- "ml_full_rej_collin"
# PC_dimred <- DimPlot(seur_obj, reduction = "pca",label = F, dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
# print(PC_dimred)
# dev.off()
# 
# saveRDS(seur_obj,paste0(resultsdir,"/unfiltered_collins_ML_subset_26012022.rds"))
# 
# # Checking the PCs after removing unknown
# seur_obj <-subset(x=seur_obj, subset = (ml_full_rej_collin == "Co_C1_SCE"|ml_full_rej_collin == "Co_C10_BCE"|ml_full_rej_collin == "Co_C2_CjS"|ml_full_rej_collin == "Co_C3_CjB"|
#                                           ml_full_rej_collin == "Co_C5_LNCP"|ml_full_rej_collin == "Co_C6_LPC"))
# # Check clusters on  the first two PCs
# seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = 2)
# mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
# pca <- seur_obj[["pca"]]
# 
# # Get the total variance:
# total_variance <- sum(matrixStats::rowVars(mat))
# eigValues = (pca@stdev)^2  ## EigenValues
# varExplained = eigValues / total_variance
# 
# # Checking the first PC
# Stdev(object = seur_obj[["pca"]])[1]
# Idents(seur_obj) <- "ml_full_rej_collin"
# 
# pdf(paste0(resultsdir,'/PCplot_ML_omitted.pdf') ,width=6,height=5,paper='special')
# pc1_viz <-2
# pc2_viz <- 1
# y_label = "PC2"#paste0(paste0(paste0("PC",pc2_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
# x_label = "PC1"#paste0(paste0(paste0("PC",pc1_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
# PC_dimred <- DimPlot(seur_obj, reduction = "pca",label = F, dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
# print(PC_dimred)
# dev.off()
# 
# saveRDS(seur_obj,paste0(resultsdir,"/filtered_collins_ML_subset_26012022.rds"))
# 
# # Check Catala
# cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20220126/cornea_integrated_normalized_unannotated_rPCA_30dims_ML.rds")
# 
# cells.to.subset <- c('Ca_C10_WSE','Ca_C8_LSC','Ca_C6_LCE','Ca_C5_BCE','Ca_C7_BCE')
# seur_obj <- cornea.combined
# seur_obj <-subset(x=seur_obj, subset = (ml_full_catala == "Ca_C10_WSE"|ml_full_catala == "Ca_C8_LSC"|ml_full_catala == "Ca_C6_LCE"|ml_full_catala == "Ca_C5_BCE"|
#                                           ml_full_catala == "Ca_C7_BCE"))
# 
#  cluster_order <- c('Ca_C10_CLC','Ca_C8_LSC','Ca_C6_LCE','Ca_C5_BCE','Ca_C7_WSEL')
# # 
# # # Generate new cluster labeling:
#  current.cluster.ids <- c('Ca_C10_WSE','Ca_C8_LSC','Ca_C6_LCE','Ca_C5_BCE','Ca_C7_BCE')
# # 
#  new.cluster.ids <- c('Ca_C10_CLC','Ca_C8_LSC','Ca_C6_LCE','Ca_C5_BCE','Ca_C7_WSEL')
# # 
#  seur_obj$ml_full_catala <- plyr::mapvalues(x = seur_obj$ml_full_catala, from = current.cluster.ids, to = new.cluster.ids)
#  
#  cluster_order <- c('Ca_C10_CLC','Ca_C8_LSC','Ca_C6_LCE','Ca_C5_BCE','Ca_C7_WSEL','Unknown')
#  # 
#  # # Generate new cluster labeling:
#  current.cluster.ids <- c('Ca_C10_WSE','Ca_C8_LSC','Ca_C6_LCE','Ca_C5_BCE','Ca_C7_BCE','Unknown')
#  # 
#  new.cluster.ids <- c('Ca_C10_CLC','Ca_C8_LSC','Ca_C6_LCE','Ca_C5_BCE','Ca_C7_WSEL','Unknown')
#  # 
#  seur_obj$ml_full_rej_catala <- plyr::mapvalues(x = seur_obj$ml_full_rej_catala, from = current.cluster.ids, to = new.cluster.ids)
#  
# 
# # Check clusters on  the first two PCs
# seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = 2)
# mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
# pca <- seur_obj[["pca"]]
# 
# # Get the total variance:
# total_variance <- sum(matrixStats::rowVars(mat))
# eigValues = (pca@stdev)^2  ## EigenValues
# varExplained = eigValues / total_variance
# 
# # Checking the first PC
# Stdev(object = seur_obj[["pca"]])[1]
# 
# Idents(seur_obj) <- "ml_full_catala"
# pdf(paste0(resultsdir,'/PCplot_ML_non-omitted_ca.pdf') ,width=6,height=5,paper='special')
# pc1_viz <-2
# pc2_viz <- 1
# y_label = "PC2"#paste0(paste0(paste0("PC",pc2_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
# x_label = "PC1"#paste0(paste0(paste0("PC",pc1_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
# PC_dimred <- DimPlot(seur_obj, reduction = "pca",label = F, dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
# print(PC_dimred)
# Idents(seur_obj) <- "ml_full_rej_catala"
# PC_dimred <- DimPlot(seur_obj, reduction = "pca",label = F, dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
# print(PC_dimred)
# dev.off()
# 
# saveRDS(seur_obj,paste0(resultsdir,"/unfiltered_catalas_ML_subset_26012022.rds"))
# 
# # Checking the PCs after removing unknown
# seur_obj <-subset(x=seur_obj, subset = (ml_full_rej_catala == "Ca_C10_CLC"|ml_full_rej_catala == "Ca_C8_LSC"|ml_full_rej_catala == "Ca_C6_LCE"|ml_full_rej_catala == "Ca_C5_BCE"|
#                                           ml_full_rej_catala == "Ca_C7_WSEL"))
# 
# 
# # Check clusters on  the first two PCs
# seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = 2)
# mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
# pca <- seur_obj[["pca"]]
# 
# # Get the total variance:
# total_variance <- sum(matrixStats::rowVars(mat))
# eigValues = (pca@stdev)^2  ## EigenValues
# varExplained = eigValues / total_variance
# 
# # Checking the first PC
# Stdev(object = seur_obj[["pca"]])[1]
# Idents(seur_obj) <- "ml_full_rej_catala"
# 
# pdf(paste0(resultsdir,'/PCplot_ML_omitted_ca.pdf') ,width=6,height=5,paper='special')
# pc1_viz <-2
# pc2_viz <- 1
# y_label = "PC2"#paste0(paste0(paste0("PC",pc2_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
# x_label = "PC1"#paste0(paste0(paste0("PC",pc1_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
# PC_dimred <- DimPlot(seur_obj, reduction = "pca",label = F, dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
# print(PC_dimred)
# dev.off()
# 
# saveRDS(seur_obj,paste0(resultsdir,"/filtered_catalas_ML_subset_26012022.rds"))
# ###############################################################################
# # Continue with rpca log normalization 30 dims (looked best with integration)
# cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20220126/cornea_integrated_normalized_unannotated_rPCA_30dims_ML.rds")
# cornea.combined <-subset(x=cornea.combined, subset = (Condition == "Co"|Condition == "Ca"))
# 
# #cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20220114 cornea_integrated_sctransform_unannotated_rpca_30dims.rds")
# 
# DefaultAssay(cornea.combined) <- "RNA"
# #cornea.combined <- QCSeurObj(cornea.combined)
# #DefaultAssay(cornea.combined) <- "integrated"
# 
# # unbiasedly look at the PCs to selected the optimal one
# PCASeurObj(SeuratObject = cornea.combined,pcs = 30)
# 
# #cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20211223 cornea_integrated_sctransform_unannotated_rpca_15dims.rds")
# 
# #DefaultAssay(cornea.combined) <- "integrated"
# 
# # Run UMAP for two datasets with integration
# cornea.combined <- RunUMAP(cornea.combined, reduction = "pca", dims = 1:26)
# 
# # add the ml labels again
# vec1 <- c(reference$ml_labels_collin,query$ml_labels_collin,query2$ml_labels_collin)
# vec2 <- c(reference$ml_labels_rej_collin,query$ml_labels_rej_collin,query2$ml_labels_rej_collin)
# cornea.combined$ml_full_collin <- as.factor(vec1)
# cornea.combined$ml_full_rej_collin <- as.factor(vec2)
# 
# vec1 <- c(reference$ml_labels_catala,query$ml_labels_catala,query2$ml_labels_catala)
# vec2 <- c(reference$ml_labels_rej_catala,query$ml_labels_rej_catala,query2$ml_labels_rej_catala)
# cornea.combined$ml_full_catala <- as.factor(vec1)
# cornea.combined$ml_full_rej_catala <- as.factor(vec2)
# 
# # Plot the tree and the final clustering
# pdf(paste(resultsdir, "clustering_collins_catala_rpca_lognorm_26dims.pdf", sep = '/'), width = 6, height = 5)
# #p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
# p1 <- DimPlot(cornea.combined, reduction = 'umap', label=TRUE,label.size = 8)
# p2 <- DimPlot(cornea.combined, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p3 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
# p4 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
# p5 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'ml_labels_collin') + ggtitle("Louvain Clustering")
# p6 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'ml_labels_rej_collin') + ggtitle("Louvain Clustering")
# p7 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'ml_labels_catala') + ggtitle("Louvain Clustering")
# p8 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'ml_labels_rej_catala') + ggtitle("Louvain Clustering")
# #ape::plot.phylo(x = data.tree, direction = "downwards",)
# #
# p1
# p2
# p3
# p4
# p5
# p6
# p7
# p8
# #print(multiplot(p1, p2, p3, cols = 3))
# dev.off()
# 
# cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20220126/cornea_integrated_SCT_unannotated_rPCA_30dims_ML.rds")
# cornea.combined <-subset(x=cornea.combined, subset = (Condition == "Co"|Condition == "Ca"))
# 
# #cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20220114 cornea_integrated_sctransform_unannotated_rpca_30dims.rds")
# 
# DefaultAssay(cornea.combined) <- "RNA"
# #cornea.combined <- QCSeurObj(cornea.combined)
# #DefaultAssay(cornea.combined) <- "integrated"
# 
# # unbiasedly look at the PCs to selected the optimal one
# #PCASeurObj(SeuratObject = cornea.combined,pcs = 30)
# 
# #cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20211223 cornea_integrated_sctransform_unannotated_rpca_15dims.rds")
# 
# #DefaultAssay(cornea.combined) <- "integrated"
# 
# # Run UMAP for two datasets with integration
# cornea.combined <- RunUMAP(cornea.combined, reduction = "pca", dims = 1:26)
# 
# # add the ml labels again
# vec1 <- c(reference$ml_labels_collin,query$ml_labels_collin,query2$ml_labels_collin)
# vec2 <- c(reference$ml_labels_rej_collin,query$ml_labels_rej_collin,query2$ml_labels_rej_collin)
# cornea.combined$ml_full_collin <- as.factor(vec1)
# cornea.combined$ml_full_rej_collin <- as.factor(vec2)
# 
# vec1 <- c(reference$ml_labels_catala,query$ml_labels_catala,query2$ml_labels_catala)
# vec2 <- c(reference$ml_labels_rej_catala,query$ml_labels_rej_catala,query2$ml_labels_rej_catala)
# cornea.combined$ml_full_catala <- as.factor(vec1)
# cornea.combined$ml_full_rej_catala <- as.factor(vec2)
# 
# # Plot the tree and the final clustering
# pdf(paste(resultsdir, "clustering_collins_catala_rpca_sct_26dims.pdf", sep = '/'), width = 6, height = 5)
# #p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
# p1 <- DimPlot(cornea.combined, reduction = 'umap', label=TRUE,label.size = 8)
# p2 <- DimPlot(cornea.combined, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p3 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
# p4 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
# p5 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'ml_full_collin') + ggtitle("Louvain Clustering")
# p6 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'ml_full_rej_collin') + ggtitle("Louvain Clustering")
# p7 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'ml_full_catala') + ggtitle("Louvain Clustering")
# p8 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'ml_full_rej_catala') + ggtitle("Louvain Clustering")
# #ape::plot.phylo(x = data.tree, direction = "downwards",)
# #
# p1
# p2
# p3
# p4
# p5
# p6
# p7
# p8
# #print(multiplot(p1, p2, p3, cols = 3))
# dev.off()
# 
# cornea.combined2 <-subset(x=cornea.combined, subset = (ml_full_rej_collin != "Unknown"))
# droplevels(cornea.combined2$ml_full_rej_collin)
# 
# # Plot the tree and the final clustering
# pdf(paste(resultsdir, "clustering__filtered_collins_rpca_sct_26dims.pdf", sep = '/'), width = 6, height = 5)
# #p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
# p1 <- DimPlot(cornea.combined2, reduction = 'umap', label=TRUE,label.size = 8)
# p2 <- DimPlot(cornea.combined2, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p3 <- DimPlot(cornea.combined2, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
# p4 <- DimPlot(cornea.combined2, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
# p6 <- DimPlot(cornea.combined2, label=F,label.size = 5, group.by = 'ml_full_rej_collin') + ggtitle("Louvain Clustering")
# #ape::plot.phylo(x = data.tree, direction = "downwards",)
# #
# p1
# p2
# p3
# p4
# p6
# #print(multiplot(p1, p2, p3, cols = 3))
# dev.off()
# 
# cornea.combined2 <-subset(x=cornea.combined, subset = (ml_full_rej_catala != "Unknown"))
# droplevels(cornea.combined2$ml_full_rej_catala)
# 
# # Plot the tree and the final clustering
# pdf(paste(resultsdir, "clustering__filtered__catala_rpca_sct_26dims.pdf", sep = '/'), width = 6, height = 5)
# #p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
# p1 <- DimPlot(cornea.combined2, reduction = 'umap', label=TRUE,label.size = 8)
# p2 <- DimPlot(cornea.combined2, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p3 <- DimPlot(cornea.combined2, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
# p4 <- DimPlot(cornea.combined2, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
# p6 <- DimPlot(cornea.combined2, label=F,label.size = 5, group.by = 'ml_full_rej_catala') + ggtitle("Louvain Clustering")
# #ape::plot.phylo(x = data.tree, direction = "downwards",)
# #
# p1
# p2
# p3
# p4
# p6
# #print(multiplot(p1, p2, p3, cols = 3))
# dev.off()
# 
# cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20220126/cornea_integrated_SCT_unannotated_rPCA_30dims_ML.rds")
# 
# # add the ml labels again
# vec1 <- c(reference$ml_labels_collin,query$ml_labels_collin,query2$ml_labels_collin)
# vec2 <- c(reference$ml_labels_rej_collin,query$ml_labels_rej_collin,query2$ml_labels_rej_collin)
# cornea.combined$ml_full_collin <- as.factor(vec1)
# cornea.combined$ml_full_rej_collin <- as.factor(vec2)
# 
# vec1 <- c(reference$ml_labels_catala,query$ml_labels_catala,query2$ml_labels_catala)
# vec2 <- c(reference$ml_labels_rej_catala,query$ml_labels_rej_catala,query2$ml_labels_rej_catala)
# cornea.combined$ml_full_catala <- as.factor(vec1)
# cornea.combined$ml_full_rej_catala <- as.factor(vec2)
# 
# # subset for only corneal epithelium subpopulations
# cornea.combined.ca <-subset(x=cornea.combined, subset = (ml_full_catala == "Ca_C10_CLC"|ml_full_catala == "Ca_C8_LSC"|ml_full_catala == "Ca_C6_LCE"|ml_full_catala == "Ca_C5_BCE"|
#                                                         ml_full_catala == "Ca_C7_WSEL"))
# 
# DefaultAssay(cornea.combined.ca) <- "RNA"
# 
# # Run UMAP for two datasets with integration
# 
# droplevels(cornea.combined.ca$ml_full_catala)
# 
# 
# cornea.combined.ca <-subset(x=cornea.combined.ca, subset = (ml_full_rej_catala != "Ca_C1_ASK"|ml_full_rej_catala != "Unknown"))
# 
# droplevels(cornea.combined.ca$ml_full_rej_catala)
# 
# cornea.combined.ca<-FindVariableFeatures(cornea.combined.ca)
# cornea.combined.ca <- ScaleData(cornea.combined.ca)
# PCASeurObj(SeuratObject = cornea.combined.ca,pcs = 30)
# 
# # Generating UMAP based on optimal PCs
# cornea.combined.ca <- RunUMAP(cornea.combined.ca$ml_full_rej_catala, reduction = "pca", dims = 1:22)
# 
# # Plot the tree and the final clustering
# pdf(paste(resultsdir, "clustering_corn_filtered_catala_rpca_sct_22dims.pdf", sep = '/'), width = 6, height = 5)
# #p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
# p1 <- DimPlot(cornea.combined.ca, reduction = 'umap', label=TRUE,label.size = 8)
# p2 <- DimPlot(cornea.combined.ca, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p3 <- DimPlot(cornea.combined.ca, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
# p4 <- DimPlot(cornea.combined.ca, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
# p6 <- DimPlot(cornea.combined.ca, label=F,label.size = 5, group.by = 'ml_full_rej_catala') + ggtitle("Louvain Clustering")
# #ape::plot.phylo(x = data.tree, direction = "downwards",)
# #
# p1
# p2
# p3
# p4
# p6
# #print(multiplot(p1, p2, p3, cols = 3))
# dev.off()
# 
# # subset for only corneal epithelium subpopulations Collin
# cornea.combined.co <-subset(x=cornea.combined, subset = (ml_full_collin == "Co_C1_SCE"|ml_full_collin == "Co_C10_BCE"|ml_full_collin == "Co_C2_CjS"|ml_full_collin == "Co_C3_CjB"|
#                                                            ml_full_collin == "Co_C5_LNCP"|ml_full_collin == "Co_C6_LPC"))
# 
# DefaultAssay(cornea.combined.co) <- "RNA"
# 
# # Run UMAP for two datasets with integration
# 
# droplevels(cornea.combined.co$ml_full_collin)
# droplevels(cornea.combined.co$ml_full_rej_collin)
# 
# cornea.combined.co <-subset(x=cornea.combined.co, subset = (ml_full_rej_collin != "Unknown"))
# 
# #PCASeurObj(SeuratObject = cornea.combined.ca,pcs = 30)
# 
# #FindVariableFeatures(cornea.combined.ca)
# cornea.combined.co <- RunUMAP(cornea.combined.co, reduction = "pca", dims = 1:22)
# 
# # Plot the tree and the final clustering
# pdf(paste(resultsdir, "clustering_corn_filtered_collin_rpca_sct_26dims.pdf", sep = '/'), width = 6, height = 5)
# #p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
# p1 <- DimPlot(cornea.combined.co, reduction = 'umap', label=TRUE,label.size = 8)
# p2 <- DimPlot(cornea.combined.co, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering")
# p3 <- DimPlot(cornea.combined.co, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
# p4 <- DimPlot(cornea.combined.co, label=F,label.size = 5, group.by = 'Condition') + ggtitle("Louvain Clustering")
# p6 <- DimPlot(cornea.combined.co, label=F,label.size = 5, group.by = 'ml_full_rej_collin') + ggtitle("Louvain Clustering")
# #ape::plot.phylo(x = data.tree, direction = "downwards",)
# #
# p1
# p2
# p3
# p4
# p6
# #print(multiplot(p1, p2, p3, cols = 3))
# dev.off()
# 
# saveRDS(cornea.combined.co,paste0(resultsdir,"/filtered_collin_SCT_rPCA_22dims_ML_subset_28012022.rds"))
# 
# ###############################################################################
# # generating own clusters SCT rPCA cornea
# 
# cornea.combined <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20220126/cornea_integrated_SCT_unannotated_rPCA_30dims_ML.rds")
# 
# # Selecting corneal populations (based on collin)
# cornea.combined <-subset(x=cornea.combined, subset = (ml_full_collin == "Co_C1_SCE"|ml_full_collin == "Co_C10_BCE"|ml_full_collin == "Co_C2_CjS"|ml_full_collin == "Co_C3_CjB"|
#                                                            ml_full_collin == "Co_C5_LNCP"|ml_full_collin == "Co_C6_LPC"))
# 
# # Select 22 PCs for optimal cornea visualization
# ClusteringResSeurObjIntegrated(SeuratObject = cornea.combined,dimensions=22)
# 
# chosen_res <- 0.15
# 
# DefaultAssay(cornea.combined) <- "integrated"
# cornea.combined <- RunUMAP(cornea.combined, reduction = "pca", dims = 1:22)
# 
# seur_obj <- cornea.combined
# seur_obj <- FindNeighbors(seur_obj, dims = 1:22)
# res <- chosen_res  #cluster_variable_name <- paste0("RNA_snn_res.", res)
#   cluster_variable_name <- paste0("integrated_snn_res.", res)
#   seur_obj <- FindClusters(seur_obj, verbose = FALSE, resolution = res, graph.name = "integrated_snn")
#   seur_obj <- BuildClusterTree(seur_obj)
# 
# # Final chosen resolution
# pdf(paste(resultsdir,'clust_integrated_corn_015_22pcs.pdf',sep="/") ,width=6,height=5,paper='special')
# DimPlot(seur_obj, label=TRUE, label.size = 6, group.by = 'integrated_snn_res.0.15') + ggtitle("Louvain Clustering") + ggtitle(paste0("cluster resolution ", res))
# dev.off()
# 
# 
# #saveRDS(seur_obj,paste0(resultsdir,"/unfiltered_costum_SCT_rPCA_22dims_015res_subset_28012022.rds"))
# 
# ################################################################################
# cornea.combined <-readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20220128/unfiltered_costum_SCT_rPCA_22dims_015res_subset_28012022.rds")
# 
# # Continued with resolution: 
# current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16,17)
# 
# new.cluster.ids <- c('cell1','cell2','cell3','cell4','cell5','cell6','cell7','cell8','cell9','cell10','cell11','cell12','cell13','cell14','cell15','cell16','cell17','cell18')
# cornea.combined$costum_clustering <- plyr::mapvalues(x = cornea.combined$integrated_snn_res.0.15, from = current.cluster.ids, to = new.cluster.ids)
# cluster_order <- c('cell1','cell2','cell3','cell4','cell5','cell6','cell7','cell8','cell9','cell10','cell11','cell12','cell13','cell14','cell15','cell16','cell17','cell18') # change if nessesary
# 
# # Re-level object@ident
# cornea.combined$costum_clustering <- factor(x = cornea.combined$costum_clustering, levels = cluster_order)
# 
# # Add custom labels to the phylogenetic tree
# data.tree <- Tool(object = cornea.combined, slot = "BuildClusterTree")
# data.tree$tip.label <- levels(cornea.combined$costum_clustering)
# 
# # Generating a color matrix of all cell_types
# cells <- unique(levels(cornea.combined$costum_clustering))
# col_mat <- viridis(length(cells),alpha = 0.5,option = "H",begin = 0.1)
# dfcol <- as.data.frame(col_mat,cells)
# 
# pdf(paste(resultsdir, "20c.clustering_costum_tree.pdf", sep = '/'), width = 8, height = 8)
# #p1 <- DimPlot(cornea.combined, label=TRUE, group.by = 'cell_type')
# p1 <- DimPlot(cornea.combined, label=TRUE,label.size = 8)
# p2 <- DimPlot(cornea.combined, label=TRUE,label.size = 8, group.by = 'costum_clustering',cols = dfcol$col_mat) + ggtitle("Louvain Clustering")
# p3 <- DimPlot(cornea.combined, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
# ape::plot.phylo(x = data.tree, direction = "downwards",)
# #p1
# p2
# p3
# dev.off()

################################################################################
# Generating marker genes overview figures for all datasets and the meta_atlas

# Import the nessesary datasets
#seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/meta_all.rds")
seur_obj<- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20220128/filtered_collin_SCT_rPCA_22dims_ML_subset_28012022.rds")

#seur_obj <- cornea.combined

# See above
#reference <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated.rds")
#query <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_lapointe2021/20221201 laPointeRNAannotated.rds")
#query2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_Chen2021/20221222 ChenRNAannotated.rds")

marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/LSC_Marker_Genes.csv"
marker_genes_df <- read.table(marker_gene_file, header = TRUE, sep = ',',comment.char = "#", stringsAsFactors = F)

DefaultAssay(seur_obj) <- "RNA"
# Setting two cores as max on the system
SnowParam(workers = 2)

# # This vector will contain all marker genes across datsets
# fullmarkers2<-NULL
# 
# seur_obj$costum_clustering <- seur_obj$ml_full_rej_collin

# Make dotplots for all datasets with their cell populations and marker genes
for (paper in unique(marker_genes_df$abreviation)){
  # 
  # if (paper == "LaPointe2021"){
  #   seur_obj <- query
  # } else if (paper == "Lako2021"){
  #   seur_obj <- reference
  # } else if (paper == "Chen2021") {
  #   seur_obj <- query2}
  # 
  seur_obj@active.ident <- seur_obj$scVI_label
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes_rpca_int/'),paper)
  print(marker_dir)
  system(paste("mkdir -p ", marker_dir))
  sub_marker_df <- marker_genes_df[marker_genes_df$abreviation == paper,]
  
  fullmarkers <- NULL
  for (cell_type in unique(sub_marker_df$name_plot)){
    marker_genes <- sub_marker_df[sub_marker_df$name_plot == cell_type,]$marker_genes
    marker_genes <- unlist(strsplit(marker_genes, ','))
    fullmarkers <- c(fullmarkers,marker_genes)
    m <- length(fullmarkers)
    n <- length(unique(sub_marker_df$name_plot))
    if (m<30){
      w <- 10
    }else{
      w <- m/5
    }
  }
  fullmarkers <- unique(fullmarkers)
  fullmarkers2<-c(fullmarkers2,fullmarkers)
  pdf(paste0(paste(paste(marker_dir, paper, sep = '/'), sep = '_'),'full_markers.pdf'), width =w, height = 5)
  print(DotPlot(object = seur_obj, features = fullmarkers,cluster.idents = T,)+ RotatedAxis())
  dev.off()
}

# Getting all markers from all datasets in overview
#seur_obj <- cornea.combined
seur_obj<- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20220128/filtered_collin_SCT_rPCA_22dims_ML_subset_28012022.rds")
#seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/meta_all.rds")
DefaultAssay(seur_obj) <- "RNA"

seur_obj$costum_clustering <- seur_obj$ml_full_rej_collin
seur_obj@active.ident <- seur_obj$costum_clustering

# Marker genes for each population in the meta-atlas
for (paper in unique(marker_genes_df$abreviation)){
  
  seur_obj@active.ident <- seur_obj$costum_clustering
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes_split'))
  print(marker_dir)
  system(paste("mkdir -p ", marker_dir))
  sub_marker_df <- marker_genes_df[marker_genes_df$abreviation == paper,]
  
  fullmarkers <- NULL
  for (cell_type in unique(sub_marker_df$name_plot)){
    marker_genes <- sub_marker_df[sub_marker_df$name_plot == cell_type,]$marker_genes
    marker_genes <- unlist(strsplit(marker_genes, ','))
    #fullmarkers <- c(fullmarkers,marker_genes)
    m <- length(marker_genes)
    n <- length(unique(sub_marker_df$name_plot))
    if (m<30){
      w <- 10
    }else{
      w <- m/5
    }
    pdf(paste0(paste(paste(marker_dir,paper,cell_type,sep = '/')),".pdf"), width =w, height = 5)
    print(DotPlot(object = seur_obj, features = marker_genes,cluster.idents = T,)+ RotatedAxis())
    dev.off()
  }
}

fullmarkers2 <- unique(fullmarkers2)
m <- length(fullmarkers2)
pdf(paste0(paste(paste(resultsdir, sep = '/'), sep = '_'),'full_markers_all_datasets.pdf'), width =m/5, height = 5)
print(DotPlot(object = seur_obj, features = fullmarkers2,cluster.idents = T,)+ RotatedAxis())
dev.off()

fullmarkers3 <- c("PMEL","KRT4","ACKR1","MITF","KERA","UBE2C","CCL3","ACTA2","MMP3","COL4A3","ITGB4","TFAP2A","CPVL","GPHA2","TYRP1","S100A9","TP63")
fullmarkers4 <- c("CEBPD","CXCL14","TP63","CPVL","CASP14","GPHA2","BMI1","NOTCH1","BIRC5","KRT14","KRT15","KRT3","PAX6","KRT13","S100A8","S100A9","MUC4","KRT19")

pdf(paste0(paste(paste(resultsdir, sep = '/'), sep = '_'),'full_markers_all_datasets_sub.pdf'), width =7, height = 4)
print(DotPlot(object = seur_obj, features = fullmarkers3,cluster.idents = T,)+ RotatedAxis())
print(DotPlot(object = seur_obj, features = fullmarkers4,cluster.idents = T,)+ RotatedAxis())
dev.off()
#saveRDS(cornea.combined, file = paste(resultsdir,"meta_all_unannotated.rds"))

# Check clusters on  the first two PCs
seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = 2)
mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
pca <- seur_obj[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

# Checking the first PC
Stdev(object = seur_obj[["pca"]])[1]

Idents(seur_obj) <- "costum_clustering"
pdf(paste0(resultsdir,'/PCtest.pdf') ,width=12,height=6,paper='special')
pc1_viz <-2
pc2_viz <- 1
y_label = "PC2"#paste0(paste0(paste0("PC",pc2_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
x_label = "PC1"#paste0(paste0(paste0("PC",pc1_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
PC_dimred <- DimPlot(seur_obj, reduction = "pca",label = T, dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
print(PC_dimred)
dev.off()
  #heights=c(3,1))

seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/meta_subset_raw.rds")
#seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = 2)
mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
pca <- seur_obj[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

# Checking the first PC
Stdev(object = seur_obj[["pca"]])[1]

Idents(seur_obj) <- "costum_clustering"
pdf(paste0(resultsdir,'/PCtestsubset.pdf') ,width=12,height=12,paper='special')
pc1_viz <- 1
pc2_viz <- 2
y_label = paste0(paste0(paste0("PC",2), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
x_label = paste0(paste0(paste0("PC",1), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
PC_dimred <- DimPlot(seur_obj, reduction = "pca",label = T, dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
print(PC_dimred)
dev.off()
# 
# # Additional: interesting genes
# Markergenes <- c("IRF1", "IRF5", "FOXC1", "NR2F2", "PAX6", "ELF3", "TNF", "CXCL1", "CXCL2", "CXCL3", "PTGS2", "NFKBIA", "OTX1", "ELK3", "HES4","KRT14")
# 
# #remove integrated assays
# #cornea.combined@assays$SCT<- NULL
# #cornea.combined@assays$integrated<- NULL
# DefaultAssay(cornea.combined) <- "RNA"
# 
# pdf(paste(resultsdir,'interesting_rpca.pdf',sep= '/'), width = 20, height = 15)
# FeaturePlot(cornea.combined, features = Markergenes)
# dev.off()

Idents(cornea.combined) <- "costum_clustering"

#marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/LSC_Marker_Genes.csv" #lapointe
marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/LSC_Marker_Genes.csv"
marker_genes_df <- read.table(marker_gene_file, header = TRUE, sep = ',',comment.char = "#", stringsAsFactors = F)

# Setting two cores as max on the system
SnowParam(workers = 2)

for (paper in unique(marker_genes_df$abreviation)){
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes_rpca_int/'),paper)
  print(marker_dir)
  system(paste("mkdir -p ", marker_dir))
  sub_marker_df <- marker_genes_df[marker_genes_df$abreviation == paper,]

  for (cell_type in unique(sub_marker_df$name_plot)){
    #print(cell_type)
    marker_genes <- sub_marker_df[sub_marker_df$name_plot == cell_type,]$marker_genes
    marker_genes <- unlist(strsplit(marker_genes, ','))

    m <- length(unique(sub_marker_df$name_plot))
    n <- length(marker_genes)
    k <- 10 ## your LEN
    marker_gene_sets <- split(marker_genes, rep(1:ceiling(n/k), each=k)[1:n])

    pdf(paste0(paste(paste(marker_dir, cell_type, sep = '/'), sep = '_'),'stacked_markers.pdf'), width = 10, height = n*2)
    print(StackedVlnPlot(obj = cornea.combined, cols2 = dfcol$col_mat,features = unlist(marker_gene_sets)))
    dev.off()

    for (count_type in c('raw_counts', 'seurat_norm_counts','SAVER_imputed_counts')){        #loop  over the different type of normalized data
      #for (count_type in c('seurat_norm_counts')){        #loop  over the different type of normalized data

      print(count_type)

      pdf(paste0(paste(paste(marker_dir, cell_type, sep = '/'), count_type, sep = '_'),'.pdf'), width = 20, height = 15)

      for (genes in  marker_gene_sets){#loop over subsets of max 12 genes to vizualize
        print(genes)
        plot_list1 <- DimPlot(cornea.combined, label=TRUE, combine = F, pt.size = 0.1)
        #plot_list1 <- append(plot_list1, DimPlot(cornea.combined, label=TRUE, group.by = 'custom_clustering', combine = F, pt.size = 0.1))
        count_plots <- plot_list1
        count_plots_clusters <- plot_list1

        if (count_type == 'raw_counts'){
          count_plots <- append(count_plots, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA', slot = 'counts')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA', slot = 'counts')))}
        if(count_type == 'seurat_norm_counts'){
          count_plots <- append(count_plots, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA')))}

        if(count_type == 'SAVER_imputed_counts'){
          count_plots <- append(count_plots, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA_SAVER')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(cornea.combined, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA_SAVER')))}


        plot_list1 <- append(plot_list1,FeaturePlot(cornea.combined, features = genes, pt.size = 0.2, ncol = 4, combine=F))
        print(cowplot::plot_grid(plotlist = plot_list1))
        #print(FeaturePlot(cornea.combined, features = genes, pt.size = 0.2, ncol = 4))

        p1 <- cowplot::plot_grid(plotlist = count_plots)
        title <- ggdraw() + draw_label(paste0(count_type,"of marker genes per timepoint"), fontface = 'bold')
        print(cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1)))
        p1 <- cowplot::plot_grid(plotlist = count_plots_clusters)
        title <- ggdraw() + draw_label(paste0(count_type,"of marker genes per cluster"), fontface = 'bold')
        print(cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1)))
      }
      dev.off()}
  }
}

# Save the object
saveRDS(cornea.combined, file = paste(resultsdir,"cornea_integrated_rpca_22dims_026res.rds"))

seur_obj <- cornea.combined
#Lets find marker genes for each cluster
cluster.markers <- FindAllMarkers(seur_obj, only.pos = TRUE)

# making a nice heatmap of expression in each cluster not significant yet
n_genes <- 40
heatmap.markers <- cluster.markers %>% group_by(cluster) %>% top_n(n_genes, avg_log2FC)
DoHeatmap(seur_obj, features = heatmap.markers$gene) + NoLegend()

#filter on only sigi
cluster.markers$gene_name_shorter <- str_sub(cluster.markers$gene,end = -1)
cluster.markers_fc <- cluster.markers[cluster.markers$p_val_adj < 0.01,]

cluster.markers_fc <- cluster.markers_fc[(cluster.markers_fc$avg_log2FC > 0.58) ,]
# avg_logFC...
table(cluster.markers_fc$cluster) # number of significantly expressed genes with a high log2FC (>0.58)

# GO terms
plot_list <- list()
cluster_GOs <- c()
#expressed_genes <- row.names(counts(sce_qc)[rowSums(counts(sce_qc))>20,])
expressed_genes <- rownames(cluster.markers_fc)

unique(expressed_genes)

#saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated2.rds") # 10 PC
#saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated_4clus.rds")
#seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated2.rds") # 10 PC
#seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated_4clus.rds")

###########################################################################
# Appending the GO terms of the unannotated clusters
plot_list <- list()
cluster_GOs <- c()
expressed_genes <- rownames(cluster.markers_fc)

for (cluster in unique(cluster.markers_fc$cluster)){
  print(cluster)
  df_subset <- cluster.markers_fc[cluster.markers_fc$cluster == cluster,]
  cluster_genes <- unique(str_sub(df_subset$gene,end = -1))
  PC1_ego <- enrichGO(gene = cluster_genes,
                      #universe = expressed_genes,
                      OrgDb         = 'org.Hs.eg.db',
                      keyType       = "SYMBOL",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      #pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
  if(is.null(PC1_ego)){
    go_plot1 <- rectGrob(gp=gpar(col=NA))
  }else{
    go_plot1 <- barplot(PC1_ego, showCategory=10)
  }
  plot_list <- c(plot_list, list(go_plot1))
  cluster_GOs <- append(cluster_GOs, (paste('cluster of cell type ', cluster, sep = '')))
}
dev.off()

# plotting the GO-terms
grob2 <- ggarrange(plotlist= plot_list, nrow = 18, ncol = 1, labels = cluster_GOs) #, align = 'v')
pdf(paste(resultsdir,'/all_DEGS_heatmap.pdf',sep="") ,width=10,height=40,paper='special')
print(ggarrange(grob2, ncol =1 , nrow = 18, widths= c(6, 4)))
print(grob2)
dev.off()

# Heatmap cluster marking genes:
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

heatmap.markers <- cluster.markers_fc
cluster.markers_fc$log2FC <- cluster.markers_fc$avg_log2FC
top10 <- cluster.markers_fc %>%  group_by(cluster) %>%  top_n(n = 5)
top5_genes <- DoHeatmap(slot = "data",seur_obj, features = top10$gene) + NoLegend()
#top10 <- cluster.markers_fc %>% group_by(cluster) %>% top_n(n = 10, wt = 'avg_log2FC')
#marker_top5 <- DoHeatmap(slot= "data",seur_obj, features = top10$gene) + NoLegend()

n_genes <- 100
heatmap.markers <- cluster.markers_fc %>% group_by(cluster) #%>% top_n(n_genes, avg_logFC)

cell_type_plot <- DimPlot(seur_obj, label=TRUE,pt.size = 2) + ggtitle("Timepoint")
cluster_plot <- DimPlot(seur_obj, label=TRUE, group.by = 'costum_clustering', pt.size = 2) + ggtitle("Louvain Clustering")
grob_clustering = ggarrange(plotlist = list(cell_type_plot, cluster_plot) , ncol =1) #labels = cluster_GOs) #align = 'v')
nclust<- length(unique(as.numeric(cluster.markers_fc$cluster)))
RGB_colours_ggplot <- as.list(as.character(gg_color_hue(nclust)))
names <- as.numeric(unique(heatmap.markers$cluster))
col_fun = colorRamp2(names, RGB_colours_ggplot)#make sure the rows correspond to ggplot colour mapping

seur_obj2 <- seur_obj
seur_obj@assays$SCT <- NULL
seur_obj@assays$integrated <- NULL

mat <- as.matrix(seur_obj@assays$RNA[heatmap.markers$gene,]) # removed scale data
mat <- rbind(mat,seur_obj@meta.data$costum_clustering)
mat <- mat[,order(mat[nrow(mat),])]

column_ha <- HeatmapAnnotation(cluster = mat[nrow(mat),], col =list(cluster = col_fun))
breaks <- mat[nrow(mat),]
mat <- mat[-nrow(mat),]
row_ha = rowAnnotation(adj_p_vallue = heatmap.markers$p_val_adj)
clust_heatmap <- grid.grabExpr(draw(Heatmap(mat, column_split = breaks,row_split = as.numeric(heatmap.markers$cluster), cluster_columns = F, cluster_rows = F,show_row_names = F,show_column_names = F,row_names_gp = gpar(fontsize = 6), top_annotation = column_ha, left_annotation = row_ha, row_names_rot = -40)))

###########################################################################
# Plotting the GO terms of the un-annotated clusters and the heatmap
#plot_list <- list()
#cluster_GOs <- c()
#expressed_genes <- rownames(cluster.markers_fc)

pdf(paste0(resultsdir,'/9.GO_terms_clusters.pdf') ,width=12,height=120,paper='special')
#cowplot::plot_grid(plotlist =grob2, ncol = 1, nrow = 1)#, rel_heights = c(1, 4,2), rel_widths = c(3,2,2))
print(grob2)
dev.off()

pdf(paste0(resultsdir,'/9.complex_heatmap.pdf') ,width=40,height=24,paper='special')
cowplot::plot_grid(top5_genes,grob_clustering, ncol = 3, nrow = 1, rel_heights = c(1, 4,2), rel_widths = c(3,2,2))
dev.off()

# Write the markers as a table
write.table(cluster.markers_fc, file = paste0(resultsdir,'/9b.cluster_markers.csv'), sep = ',')

################################################################################
# Annotate the clusters
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20211224 cornea_integrated_rpca_22dims_026res.rds")

Markergenes <- c("CXCL3","RBP1", "CCN2", "MYH10")

#remove integrated assays
#cornea.combined@assays$SCT<- NULL
#cornea.combined@assays$integrated<- NULL
DefaultAssay(seur_obj) <- "RNA"

pdf(paste(resultsdir,'interesting_cell14.pdf',sep= '/'), width = 20, height = 15)
FeaturePlot(seur_obj, features = Markergenes)
dev.off()

cluster_order <- c('BCE','StC','CjC','WSCE','LSC','PLNSC','LNSC',
                   'nsKer','CSSC','nsKer2','BVCE','Mel','taLSC',
                   'CE','IC','ESD','CE')

# Generate new cluster labeling:
current.cluster.ids <- c('cell1','cell2','cell3','cell4','cell5','cell6','cell7',
                         'cell8','cell9','cell10','cell11','cell12','cell13',
                         'cell14','cell15','cell16','cell17')

new.cluster.ids <- c('BCE','StC','CjC','WSCE','LSC','PLNSC','LNSC',
                     'nsKer','CSSC','nsKer2','BVCE','Mel','taLSC',
                     'CE','IC','ESD','CE')

# Generating a color matrix of all cell_types
cells <- unique(levels(seur_obj$costum_clustering))
col_mat <- viridis(length(cells),alpha = 0.5,option = "H",begin = 0.1)
dfcol <- as.data.frame(col_mat,cells)

seur_obj$costum_clustering <- plyr::mapvalues(x = seur_obj$costum_clustering, from = current.cluster.ids, to = new.cluster.ids)

# Re-level object@ident
seur_obj$costum_clustering <- factor(x = seur_obj$costum_clustering, levels = cluster_order)
seur_obj$costum_clustering

# Add custom labels to the phylogenetic tree
data.tree <- Tool(object = seur_obj, slot = "BuildClusterTree")
data.tree$tip.label <- levels(seur_obj$costum_clustering)

# Plot the tree and the final clustering
pdf(paste(resultsdir, "13.clustering_costum_final.pdf", sep = '/'), width = 10, height = 8)
#p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
p1 <- DimPlot(seur_obj, label=TRUE,label.size = 8)
p2 <- DimPlot(seur_obj, label=TRUE,label.size = 8, group.by = 'costum_clustering',cols = dfcol$col_mat) + ggtitle("Louvain Clustering")
p3 <- DimPlot(seur_obj, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering")
ape::plot.phylo(x = data.tree, direction = "downwards",)
#
p1
p2
p3
#print(multiplot(p1, p2, p3, cols = 3))
dev.off()

# perform DEG and GO terms?
#....


####
# compare cell dist
cluster_order <- c('Co','Ca','Li')

# Generate new cluster labeling:
current.cluster.ids <- c('Lako','Lapointe','Chen')

new.cluster.ids <- c('Co','Ca','Li')

seur_obj$Condition <- plyr::mapvalues(x = seur_obj$Condition, from = current.cluster.ids, to = new.cluster.ids)

newdf <- NULL

for (xx in unique(seur_obj$costum_clustering)) {
  print(xx)
  for (j in unique(seur_obj$Condition)){
    val<-sum(seur_obj$costum_clustering==xx & seur_obj$Condition==j)
    col<- dfcol[rownames(dfcol)%in% xx,]
    df <- data.frame(values=val,Dataset=j,cell=xx,color=col)
    newdf <- rbind(newdf, df)
  }
}

pdf(paste(resultsdir,'Cell_dist_datasets.pdf',sep="/") ,width=10,height=8,paper='special')
# Stacked barchart showing the composition
print(ggplot(newdf, aes(fill=Dataset, y=values, x=cell)) +
        geom_bar(position="stack", stat="identity") + labs(x="Cell type",y="Number of cells"))
dev.off()

for (xx in unique(seur_obj$costum_clustering)) {
  print(xx)
  val<-sum(seur_obj$costum_clustering==xx)
  col<- dfcol[rownames(dfcol)%in% xx,]
  df <- data.frame(values=val,Dataset="meta",cell=xx,color=col)
  newdf <- rbind(newdf, df)
}


# Re-order that the colors make sense for the bar_plot
newdf <- newdf[order(newdf$cell),]
dfcol2 <- dfcol[order(rownames(dfcol)),]

pdf(paste(resultsdir,'Cell_dist_datasets_stacked.pdf',sep="/") ,width=5,height=5,paper='special')
print(ggplot(newdf, aes(fill=cell, y=values, x=Dataset)) +
  geom_bar(position="fill", stat="identity", width = 0.3) + labs(x="Dataset",y="Cell proportions")+scale_fill_manual(values = dfcol2))
dev.off()

# Generate the pseudobulk-table for correlation between datasets
sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

#pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc$sample <- seur_obj@meta.data$costum_clustering
sce_qc$sample

sce_qc$rep <- seur_obj@meta.data$orig.ident
sce_qc$rep

unique(seur_obj@meta.data$orig.ident)

sce_qc$reps <- "Collin"
sce_qc$reps[sce_qc$rep == "lapointeGSM5651509" | sce_qc$rep == "lapointeGSM5651511" | sce_qc$rep == "lapointeGSM5651513" |
              sce_qc$rep == "lapointeGSM56511515"| sce_qc$rep == "lapointeGSM5651117"| sce_qc$rep == "lapointeGSM5651119"|
              sce_qc$rep == "lapointeGSM5651510" | sce_qc$rep == "lapointeGSM5651512" | sce_qc$rep == "lapointeGSM5651514" |
              sce_qc$rep == "lapointeGSM56511516"| sce_qc$rep == "lapointeGSM5651118"|sce_qc$rep == "lapointeGSM5651120"] <- "Catala"
sce_qc$reps[sce_qc$rep == "ChenGSM4646295" | sce_qc$rep == "ChenGSM4646296"] <- "Li"

pseudobulk_df <- NULL

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))
pseudobulk_df[["Collin"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$reps == "Collin"]))
pseudobulk_df[["Catala"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$reps == "Catala"]))
pseudobulk_df[["Li"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$reps == "Li"]))

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

pdf(paste(resultsdir,'Correlation_datasets_raw_Collin_v_Catala.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$Collin, ylab = "Collin")
ggqqplot(pseudobulk_df$Catala, ylab = "Catala")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$Catala, pseudobulk_df$Collin,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=Catala, y=Collin)) +
  geom_point()+
  geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$Catala,na.rm=T)/10*1,y = max(pseudobulk_df$Collin,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()

pdf(paste(resultsdir,'Correlation_datasets_raw_Collin_v_Li.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$Collin, ylab = "Collin")
ggqqplot(pseudobulk_df$Li, ylab = "Li")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$Li, pseudobulk_df$Collin,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=Li,y=Collin)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$Li,na.rm=T)/10*1,y = max(pseudobulk_df$Collin,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()

pdf(paste(resultsdir,'Correlation_datasets_raw_Catala_v_Li.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$Catala, ylab = "Catala")
ggqqplot(pseudobulk_df$Li, ylab = "Li")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$Li, pseudobulk_df$Catala,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=Li,y=Catala)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$Li,na.rm=T)/10*1,y = max(pseudobulk_df$Catala,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()

# Save the annotated meta-atlas object with all cell populations
saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/meta_all.rds")
write.table(data.frame("ID"=rownames(dfcol),dfcol), file = paste0(resultsdir,'/dfcol_19_01_2022.tsv'), sep = '\t',quote = F, row.names = F)

################################################################################
# Subset for corneal epithelium and limbus only
cells.to.subset <- c('BCE','CjC','WSCE','LSC','PLNSC','LNSC','taLSC')
seur_obj2 <- seur_obj
seur_obj <-subset(x=seur_obj, subset = (costum_clustering == "BCE"|costum_clustering == "CjC"|costum_clustering == "WSCE"|costum_clustering == "LSC"|
                                          costum_clustering == "PLNSC"|costum_clustering == "LNSC"|costum_clustering == "taLSC"))

dfcol2 <- dfcol
dfcol <- dfcol[rownames(dfcol) %in% cells.to.subset,]
dfcol <- as.data.frame(dfcol,cells.to.subset)

####
# compare cell dist

newdf <- NULL

for (xx in unique(seur_obj$costum_clustering)) {
  print(xx)
  for (j in unique(seur_obj$Condition)){
    val<-sum(seur_obj$costum_clustering==xx & seur_obj$Condition==j)
    col<- dfcol[rownames(dfcol)%in% xx,]
    df <- data.frame(values=val,Dataset=j,cell=xx,color=col)
    newdf <- rbind(newdf, df)
  }
}

pdf(paste(resultsdir,'Cell_dist_datasets_subset.pdf',sep="/") ,width=10,height=8,paper='special')
# Stacked barchart showing the composition
print(ggplot(newdf, aes(fill=Dataset, y=values, x=cell)) +
        geom_bar(position="stack", stat="identity") + labs(x="Cell type",y="Number of cells"))
dev.off()

for (xx in unique(seur_obj$costum_clustering)) {
  print(xx)
  val<-sum(seur_obj$costum_clustering==xx)
  col<- dfcol[rownames(dfcol)%in% xx,]
  df <- data.frame(values=val,Dataset="meta",cell=xx,color=col)
  newdf <- rbind(newdf, df)
}


# Re-order that the colors make sense for the bar_plot
newdf <- newdf[order(newdf$cell),]
dfcol2 <- dfcol[order(rownames(dfcol)),]

pdf(paste(resultsdir,'Cell_dist_datasets_stacked_subset.pdf',sep="/") ,width=2*length(unique(newdf$Dataset)),height=5,paper='special')
print(ggplot(newdf, aes(fill=cell, y=values, x=Dataset)) +
        geom_bar(position="fill", stat="identity", width = 0.3) + labs(x="Dataset",y="Cell proportions")+scale_fill_manual(values = dfcol2))
dev.off()

# Generate the pseudobulk-table for correlation between datasets
sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

#pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc$sample <- seur_obj@meta.data$costum_clustering
sce_qc$sample

sce_qc$rep <- seur_obj@meta.data$orig.ident
sce_qc$rep

unique(seur_obj@meta.data$orig.ident)

sce_qc$reps <- "Collin"
sce_qc$reps[sce_qc$rep == "lapointeGSM5651509" | sce_qc$rep == "lapointeGSM5651511" | sce_qc$rep == "lapointeGSM5651513" |
              sce_qc$rep == "lapointeGSM56511515"| sce_qc$rep == "lapointeGSM5651117"| sce_qc$rep == "lapointeGSM5651119"|
              sce_qc$rep == "lapointeGSM5651510" | sce_qc$rep == "lapointeGSM5651512" | sce_qc$rep == "lapointeGSM5651514" |
              sce_qc$rep == "lapointeGSM56511516"| sce_qc$rep == "lapointeGSM5651118"|sce_qc$rep == "lapointeGSM5651120"] <- "Catala"
sce_qc$reps[sce_qc$rep == "ChenGSM4646295" | sce_qc$rep == "ChenGSM4646296"] <- "Li"

pseudobulk_df <- NULL

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))
pseudobulk_df[["Collin"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$reps == "Collin"]))
pseudobulk_df[["Catala"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$reps == "Catala"]))
pseudobulk_df[["Li"]] <- as.vector(rowSums(counts(sce_qc)[,sce_qc$reps == "Li"]))

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

pdf(paste(resultsdir,'Correlation_datasets_subset_Collin_v_Catala.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$Collin, ylab = "Collin")
ggqqplot(pseudobulk_df$Catala, ylab = "Catala")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$Catala, pseudobulk_df$Collin,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=Catala, y=Collin)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$Catala,na.rm=T)/10*1,y = max(pseudobulk_df$Collin,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()

pdf(paste(resultsdir,'Correlation_datasets_subset_Collin_v_Li.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$Collin, ylab = "Collin")
ggqqplot(pseudobulk_df$Li, ylab = "Li")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$Li, pseudobulk_df$Collin,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=Li,y=Collin)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$Li,na.rm=T)/10*1,y = max(pseudobulk_df$Collin,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()

pdf(paste(resultsdir,'Correlation_datasets_subset_Catala_v_Li.pdf',sep="/") ,width=5,height=5,paper='special')

# test for normality with qq-plots
ggqqplot(pseudobulk_df$Catala, ylab = "Catala")
ggqqplot(pseudobulk_df$Li, ylab = "Li")

# data not normally distributed thus spearman correlation
spearmancortest <- cor.test(pseudobulk_df$Li, pseudobulk_df$Catala,  method = "spearman", use = "complete.obs")
spearmancorval <- round(as.numeric(as.character(unname(spearmancortest['estimate'][[1]]))),2)

print(ggplot(pseudobulk_df, aes(x=Li,y=Catala)) +
        geom_point()+
        geom_smooth(method=lm))+
  scale_color_manual(aesthetics = "color")+
  geom_text(x= max(pseudobulk_df$Li,na.rm=T)/10*1,y = max(pseudobulk_df$Catala,na.rm = T)/10*9, label = paste0("Rho: ",spearmancorval), parse = TRUE,check_overlap = T)

dev.off()

# Save the annotated meta-atlas object with all cell populations

# Process the subset to raw
seur_obj@active.ident<-seur_obj$costum_clustering
seur_obj@assays$SCT <- NULL
seur_obj@assays$integrated <- NULL

saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/meta_subset_raw.rds")
write.table(data.frame("ID"=rownames(dfcol),dfcol), file = paste0(resultsdir,'/dfcol_subset_19_01_2022.tsv'), sep = '\t',quote = F, row.names = F)

varfeatures<-seur_obj2@assays$integrated@data@Dimnames[[1]]
write.table(varfeatures, file = paste0(resultsdir,'/varfeatures.txt'), sep = '\t',quote = F, row.names = F,col.names = F)


################################################################################
# volcanoplots
knots <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/pseudotime/20220104/earlyDERes_k15_30.tsv",sep='\t',row.names=1)
# Generate artificial replicates for each cell type based on extraction

lakovstcalc3 <- as.data.frame(lakovstcalc2@listData$log2FoldChange,row.names = rownames(lakovst2))
lakovstcalc3$padj <- lakovstcalc2@listData$padj

v2 <- c()
indx <- grepl('logFC', colnames(knots))

pdf(paste(resultsdir,'Volcanoplots_datasets.pdf',sep="/") ,width=6,height=6,paper='special')

# Perform mulitple testing correction upon the Pvalues from trade-seq
knots$padj <- p.adjust(knots$pvalue, method = "bonferroni", n = length(knots$pvalue))

for (i in colnames(knots)[indx]){
  print(i)
  lakovstcalc3 <- as.data.frame(knots[,colnames(knots)==i],row.names = rownames(knots))

  lakovstcalc3$padj <- knots$padj

  # remove true zero values
  lakovstcalc3 <- lakovstcalc3[lakovstcalc3$padj!=0,]
  colnames(lakovstcalc3) <- c("log2FoldChange","padj")
  print(head(lakovstcalc3))
  write.table(lakovstcalc3, file = paste0(resultsdir,'/', i, '_comp.tsv'),sep = "\t", quote = F)

  # printing the volcanoplots
  dataset1_vs_dataset2 <- lakovstcalc3

  dataset1_vs_dataset2$noms<-rownames(dataset1_vs_dataset2)
  # Create new categorical column ------------------------------------------------
  dataset1_vs_dataset2 <- dataset1_vs_dataset2 %>%
    mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                 log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                                 TRUE ~ "ns"))

  # Add colour, size and alpha (transparency) to volcano plot --------------------
  vircols<-viridis(n = 2)
  cols <- c("up" = vircols[2], "down" = "blue", "ns" = vircols[1])
  sizes <- c("up" = 2, "down" = 2, "ns" = 1)
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

  print(dataset1_vs_dataset2 %>%
          ggplot(aes(x = log2FoldChange,
                     y = -log10(padj),#-log10(padj),
                     fill = gene_type,
                     size = gene_type,
                     alpha = gene_type)) +
          geom_point(shape = 21, # Specify shape and colour as fixed local parameters
                     colour = "black") +
          geom_hline(yintercept = -log10(0.05),
                     linetype = "dashed") +
          geom_vline(xintercept = c(log2(0.5), log2(2)),
                     linetype = "dashed") +
          scale_fill_manual(values = cols) + # Modify point colour
          scale_size_manual(values = sizes) + # Modify point size
          scale_alpha_manual(values = alphas) + # Modify point transparency
          scale_x_continuous(breaks = c(seq(-10, 10, 2)),
                             limits = c(-10, 10)) +
          scale_y_continuous(limits = c(0, 10)) +
          geom_text_repel(
            data = subset(dataset1_vs_dataset2, padj < 0.05 & log2FoldChange>2.5),
            aes(label = noms),
            nudge_x = 8,
            segment.size  = 0.2,
            segment.color = "grey50",
            size = 3,
            box.padding = 0.5,force = 5,
            #box.padding = unit(0.35, "lines"),
            max.overlaps = 10)+#,
          geom_text_repel(
            data = subset(dataset1_vs_dataset2, padj < 0.05 & log2FoldChange< (-2.5)),
            aes(label = noms),
            nudge_x = -8,
            segment.size  = 0.2,
            segment.color = "grey50",
            size = 3,
            box.padding = 0.5,force = 5,
            #box.padding = unit(0.35, "lines"),
            max.overlaps = 10)+
            #point.padding = unit(0.3, "lines"))+
          labs(fill="Gene type color",size="Size",title=i))
}
dev.off()

###########################################################################
# Splitting the pseudobulk datasets unbiased into two replicates for Z-score calculation and over pseudotime
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/meta_all.rds")

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc$sample <- seur_obj@meta.data$costum_clustering

sce_qc$sample <- factor(sce_qc$sample)
sce_qc$sample

sce_qc$rep <- seur_obj@meta.data$orig.ident
sce_qc$rep

unique(seur_obj@meta.data$orig.ident)

sce_qc$reps <- "rep1"

# use rep2 for Catala and one subject from Li
sce_qc$reps[sce_qc$rep == "ChenGSM4646296"] <- "rep2"
sce_qc$reps[sce_qc$rep == "lapointeGSM5651509" | sce_qc$rep == "lapointeGSM5651511" | sce_qc$rep == "lapointeGSM5651513" | sce_qc$rep == "lapointeGSM56511515"| sce_qc$rep == "lapointeGSM5651117"| sce_qc$rep == "lapointeGSM5651119"] <- "rep2"
sce_qc$reps[sce_qc$rep == "lapointeGSM5651510" | sce_qc$rep == "lapointeGSM5651512" | sce_qc$rep == "lapointeGSM5651514" | sce_qc$rep == "lapointeGSM56511516"| sce_qc$rep == "lapointeGSM5651118"|sce_qc$rep == "lapointeGSM5651120"] <- "rep2"

# putting values of artificial rep2 and rep1 into the pseudobulkdataframecolumns
for (sample in unique(sce_qc$sample)){
  pseudobulk_df[ , paste0(sample,"_1")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep1"])
  pseudobulk_df[ , paste0(sample,"_2")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep2"])
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulk_reps_DE_datasets_markers_split2.tsv'), sep = '\t',quote = F, row.names = F)

# For the deseq2 matrix vs ESC
lakocountfile <- pseudobulk_df
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score_datasets"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

#collist<- c("LiCo","StCSC")

# for complex heatmap
#lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/pseudobulk.tsv', sep = '\t', header = TRUE, row.names = 1)
#lakovst <- as.matrix(lakocountfile,row.names="ID")

lakovst <- lakocountfile

# Generate coldata dataframe
coldata <- NULL

# conditions
j <- as.character(unname(unique(seur_obj$costum_clustering)))
b<-2 # Or some other number
j<-sapply(j, function (x) rep(x,b))
j<-as.vector(j)

# Change _1 & _2 to condition 1 and _3 & _4 to condition 2
c <- c("1","2")
d <- paste(j,c,sep="_")

# reps
cols <- colnames(pseudobulk_df)

# type
xx <- rep("paired-end",length(cols))

coldata <- data.frame(cols=cols,condition=d,condition2=j,type=xx)

#coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/col2.tsv', sep = '\t', header = TRUE, row.names = 1)
#coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210802/col.tsv', sep = '\t', header = TRUE, row.names = 1)

# for complex heatmap
#coldata <- coldata[1:22,]
rownames(coldata)<- coldata$cols
coldata <- coldata[,c("condition","condition2","type")]

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# checking the data
head(lakovst,2)

coldata

# setting the rowdata in coldata similar to the coldata in the count matrix
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(lakovst))

all(rownames(coldata) == colnames(lakovst))

lakovst <- lakovst[, rownames(coldata)]
all(rownames(coldata) == colnames(lakovst))

dds <- DESeqDataSetFromMatrix(countData = round(lakovst),
                              colData = coldata,
                              design = ~ condition)
dds
###################################################################
# for complex heatmap
#dds2 <- DESeq(dds)

# transpose and scale the matrix per row Z-score per gene per sample is scaled upon
vsd <- assay(vst(dds,blind = T))
Z <- t(scale(t(vsd)))
Z
Z_score <- as.data.frame(Z,row.names = rownames(lakovst))

# joining the columns on the means sequential
n <- 2
Z_joined <- t(rowMeans(t(Z_score), as.integer(gl(ncol(Z_score), n, ncol(Z_score)))) / n)

df %>% mutate(mean_all = rowMeans(.),
              mean_sel = rowMeans(select(., select_vars)))

#  generate list of factors
vec <- rownames(coldata)
x <- split(vec, ceiling(seq_along(vec)/2))

scoretable <- as.data.frame(do.call(cbind, lapply(x, function(i) rowMeans(Z_score[, i]))), row.names = rownames(lakovst))
colnames(scoretable) <- unique(coldata$condition2)
write.table(scoretable, file = paste0(resultsdir, '/Zscoretable_markers_split2.tsv'),sep = "\t", quote = F,row.names = T,col.names = T)
scoretable1 <- scoretable
################################################################################
# Subset over pseudotime
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/meta_subset.rds")

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc$sample <- seur_obj@meta.data$costum_clustering

sce_qc$sample <- factor(sce_qc$sample)
sce_qc$sample

sce_qc$rep <- seur_obj@meta.data$orig.ident
sce_qc$rep

unique(seur_obj@meta.data$orig.ident)

sce_qc$reps <- "rep1"

# use rep2 for lapointe and one subject from Li
sce_qc$reps[sce_qc$rep == "ChenGSM4646296"] <- "rep2"
sce_qc$reps[sce_qc$rep == "lapointeGSM5651509" | sce_qc$rep == "lapointeGSM5651511" | sce_qc$rep == "lapointeGSM5651513" | sce_qc$rep == "lapointeGSM56511515"| sce_qc$rep == "lapointeGSM5651117"| sce_qc$rep == "lapointeGSM5651119"] <- "rep2"
sce_qc$reps[sce_qc$rep == "lapointeGSM5651510" | sce_qc$rep == "lapointeGSM5651512" | sce_qc$rep == "lapointeGSM5651514" | sce_qc$rep == "lapointeGSM56511516"| sce_qc$rep == "lapointeGSM5651118"|sce_qc$rep == "lapointeGSM5651120"] <- "rep2"

# putting values of artificial rep2 and rep1 into the pseudobulkdataframecolumns
for (sample in unique(sce_qc$sample)){
  pseudobulk_df[ , paste0(sample,"_1")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep1"])
  pseudobulk_df[ , paste0(sample,"_2")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep2"])
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulk_reps_DE_datasets_pseudotime_split2.tsv'), sep = '\t',quote = F, row.names = F)

# For the deseq2 matrix vs ESC
lakocountfile <- pseudobulk_df
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score_datasets"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

#collist<- c("LiCo","StCSC")

# for complex heatmap
lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/pseudobulk.tsv', sep = '\t', header = TRUE, row.names = 1)
#lakovst <- as.matrix(lakocountfile,row.names="ID")

lakovst <- lakocountfile

# Generate coldata dataframe
coldata <- NULL

# conditions
j <- as.character(unname(unique(seur_obj$costum_clustering)))
b<-2 # Or some other number
j<-sapply(j, function (x) rep(x,b))
j<-as.vector(j)

# Change _1 & _2 to condition 1 and _3 & _4 to condition 2
c <- c("1","2")
d <- paste(j,c,sep="_")

# reps
cols <- colnames(pseudobulk_df)

# type
xx <- rep("paired-end",length(cols))

coldata <- data.frame(cols=cols,condition=d,condition2=j,type=xx)

#coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/col2.tsv', sep = '\t', header = TRUE, row.names = 1)
#coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210802/col.tsv', sep = '\t', header = TRUE, row.names = 1)

# for complex heatmap
#coldata <- coldata[1:22,]
rownames(coldata)<- coldata$cols
coldata <- coldata[,c("condition","condition2","type")]

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# checking the data
head(lakovst,2)

coldata

# setting the rowdata in coldata similar to the coldata in the count matrix
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(lakovst))

all(rownames(coldata) == colnames(lakovst))

lakovst <- lakovst[, rownames(coldata)]
all(rownames(coldata) == colnames(lakovst))

dds <- DESeqDataSetFromMatrix(countData = round(lakovst),
                              colData = coldata,
                              design = ~ condition)
dds
###################################################################
# for complex heatmap
#dds2 <- DESeq(dds)

# transpose and scale the matrix per row Z-score per gene per sample is scaled upon
vsd <- assay(vst(dds,blind = T))
Z <- t(scale(t(vsd)))
Z
Z_score <- as.data.frame(Z,row.names = rownames(lakovst))

# joining the columns on the means sequential
n <- 2
Z_joined <- t(rowMeans(t(Z_score), as.integer(gl(ncol(Z_score), n, ncol(Z_score)))) / n)

df %>% mutate(mean_all = rowMeans(.),
              mean_sel = rowMeans(select(., select_vars)))

#  generate list of factors
vec <- rownames(coldata)
x <- split(vec, ceiling(seq_along(vec)/2))

scoretable <- as.data.frame(do.call(cbind, lapply(x, function(i) rowMeans(Z_score[, i]))), row.names = rownames(lakovst))
colnames(scoretable) <- unique(coldata$condition2)
write.table(scoretable, file = paste0(resultsdir, '/Zscoretable_pseudotime_split2.tsv'),sep = "\t", quote = F,row.names = T,col.names = T)
scoretable2<-scoretable
################################################################################
# Plotting the Z-score heatmaps of the top 25 target genes for each
# cell type in the cornea
#lakorna <- scoretable1
lakorna <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score_datasets20220117/Zscoretable_markers_split2.tsv",sep="\t")

x <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20211224/9b.cluster_markers.csv", sep=",",row.names=1)
# Filter out genes with a true zero p-value and a p-adjusted value score below 0.05
# x <- x[x$p_val_adj!=0,]
# x <- x[x$p_val_adj<0.05,]


x$genes <- rownames(x)

z <- subset(lakorna, rownames(lakorna) %in% (x$genes))
matz <- as.matrix(z)

# Re-name unannotated cluster names to annotated cluster names
as.factor(x$cluster)

cluster_order <- c('BCE','StC','CjC','WSCE','LSC','PLNSC','LNSC',
                   'nsKer','CSSC','nsKer2','BVCE','Mel','taLSC',
                   'CE','IC','ESD','CE')

# Generate new cluster labeling:
current.cluster.ids <- c('cell1','cell2','cell3','cell4','cell5','cell6','cell7',
                         'cell8','cell9','cell10','cell11','cell12','cell13',
                         'cell14','cell15','cell16','cell17')

new.cluster.ids <- c('BCE','StC','CjC','WSCE','LSC','PLNSC','LNSC',
                     'nsKer','CSSC','nsKer2','BVCE','Mel','taLSC',
                     'CE','IC','ESD','CE')



x$cluster <- plyr::mapvalues(x = as.factor(x$cluster), from = current.cluster.ids, to = new.cluster.ids)
x$cluster <- as.factor(x$cluster)
vec2 <- NULL

for (i in colnames(z)) {
  matsel <- x[x$cluster==i,]
  vec <- rownames(matsel)[1:25]
  print(vec)
  vec2 <- c(vec2,vec)
}


vec2 <- unique(vec2)

z <- subset(lakorna, rownames(lakorna) %in% (vec2))

matz <- as.matrix(z)

# plot the heatmap with the top 25 marker genes Z-score for each cluster based
# on the log value

pdf(paste(resultsdir,'heatmap_Z_score_annotated_all_top25.pdf',sep="/") ,width=8,height=15,paper='special')
f1 = colorRamp2(c(-1, 0, 1), c("darkblue", "#EEEEEE", "darkred"), space = "RGB")

ht1 = Heatmap(matz, col = f1, cluster_columns = T,cluster_rows = T, name = "Z-score RNA count")
htlist1 <- ht1
print(draw(htlist1))
dev.off()

################################################################################
table(cluster.markers_fc$cluster) # number of significantly expressed genes with a high log2FC (>0.58)

# GO terms
plot_list <- list()
cluster_GOs <- c()
#expressed_genes <- row.names(counts(sce_qc)[rowSums(counts(sce_qc))>20,])
expressed_genes <- rownames(cluster.markers_fc)

unique(expressed_genes)

###########################################################################
# Appending the GO terms of the annotated clusters
cluster.markers_fc <- x
plot_list <- list()
cluster_GOs <- c()
expressed_genes <- rownames(cluster.markers_fc)

for (cluster in unique(cluster.markers_fc$cluster)){
  print(cluster)
  df_subset <- cluster.markers_fc[cluster.markers_fc$cluster == cluster,]
  cluster_genes <- unique(str_sub(df_subset$gene,end = -1))
  PC1_ego <- enrichGO(gene = cluster_genes,
                      #universe = expressed_genes,
                      OrgDb         = 'org.Hs.eg.db',
                      keyType       = "SYMBOL",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      #pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

  cc.simp <- simplify( PC1_ego, cutoff = 0.4, by = "p.adjust", select_fun = min )

  if(is.null(PC1_ego)){
    go_plot1 <- rectGrob(gp=gpar(col=NA))
  }else{
    go_plot1 <- barplot(PC1_ego, showCategory=10)
  }
  plot_list <- c(plot_list, list(go_plot1))
  cluster_GOs <- append(cluster_GOs, (paste('cluster of cell type ', cluster, sep = '')))

  # Sankey diagrams from genes
  df <- as.data.frame(row.names =cc.simp@result$Description,x=cc.simp@result$geneID)
  df$go <- rownames(df)

  colnames(df)<- c("genes","go")
  df$score <- cc.simp@result$p.adjust
  df2 <- df %>%
    mutate(genes = strsplit(as.character(genes), "/")) %>%
    unnest(genes)

  df2 <- df2[order(df2$score,decreasing = F),]

  df2.first <- df2[match(unique(df2$genes), df2$genes),]
  df2.first <- df2.first[df2.first$genes %in% cluster_genes[1:25],]

  pdf(paste(resultsdir,'/',cluster,'marker_GO_sankey.pdf',sep=""),width=10,height=length(df2.first$genes)*0.8,paper='special')
  print(ggplot(data = df2.first,
               aes(axis1 = genes, axis2 = go,
                   y = -log(score))) +
          scale_x_discrete(limits = c("genes", "go"))+#, expand = c(.2, .05)) +
          geom_alluvium() +
          geom_stratum() +
          geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
          #ggfittext::geom_fit_text(stat = "stratum", width = 1/4, min.size = 3) +
          theme_minimal())#+

  #ggtitle("passengers on the maiden voyage of the Titanic",
  #"stratified by demographics and survival"))

  #print(alluvial(df2.first[,1:2], freq = -log(df2.first$score)),right = T)#,col=extracted2$col,border = NA,blocks = TRUE)
  dev.off()
}

grob2 <- ggarrange(plotlist= plot_list, nrow = 5, ncol = 4, labels = cluster_GOs) #, align = 'v')
pdf(paste0(resultsdir,'/9.GO_terms_clusters.pdf') ,width=35,height=20,paper='special')
#cowplot::plot_grid(plotlist =grob2, ncol = 1, nrow = 1)#, rel_heights = c(1, 4,2), rel_widths = c(3,2,2))
print(grob2)
dev.off()

################################################################################
# Z-score markers for subsetted only
lakorna <- scoretable2

x <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_3_datasets/20211224/9b.cluster_markers.csv", sep=",",row.names=1)
# Filter out genes with a true zero p-value and a p-adjusted value score below 0.05
# x <- x[x$p_val_adj!=0,]
# x <- x[x$p_val_adj<0.05,]


x$genes <- rownames(x)

#z <- subset(lakorna, rownames(lakorna) %in% (x$genes))

#matz <- as.matrix(z)

# Re-name unannotated cluster names to annotated cluster names
as.factor(x$cluster)

cluster_order <- c('BCE','StC','CjC','WSCE','LSC','PLNSC','LNSC',
                   'nsKer','CSSC','nsKer2','BVCE','Mel','taLSC',
                   'CE','IC','ESD','CE')

# Generate new cluster labeling:
current.cluster.ids <- c('cell1','cell2','cell3','cell4','cell5','cell6','cell7',
                         'cell8','cell9','cell10','cell11','cell12','cell13',
                         'cell14','cell15','cell16','cell17')

new.cluster.ids <- c('BCE','StC','CjC','WSCE','LSC','PLNSC','LNSC',
                     'nsKer','CSSC','nsKer2','BVCE','Mel','taLSC',
                     'CE','IC','ESD','CE')



x$cluster <- plyr::mapvalues(x = as.factor(x$cluster), from = current.cluster.ids, to = new.cluster.ids)
x$cluster <- as.factor(x$cluster)
vec2 <- NULL

for (i in colnames(z)) {
  matsel <- x[x$cluster==i,]
  vec <- rownames(matsel)[1:25]
  print(vec)
  vec2 <- c(vec2,vec)
}
vec2 <- unique(vec2)

z <- subset(lakorna, rownames(lakorna) %in% (vec2))

matz <- as.matrix(z)

vec2 <- NULL

for (i in colnames(z)) {
  matsel <- x[x$cluster==i,]
  vec <- rownames(matsel)[1:25]
  print(vec)
  vec2 <- c(vec2,vec)
}
vec2 <- unique(vec2)

z <- subset(lakorna, rownames(lakorna) %in% (vec2))

matz <- as.matrix(z)

# plot the heatmap with the top 25 marker genes Z-score for each cluster based
# on the log value

pdf(paste(resultsdir,'heatmap_Z_score_annotated_subset_top25.pdf',sep="/") ,width=15,height=15,paper='special')
f1 = colorRamp2(c(-3, 0, 3), c("blue", "#EEEEEE", "red"), space = "RGB")

ht1 = Heatmap(matz, col = f1, cluster_columns = T,cluster_rows = T, name = "Z-score RNA count")
htlist1 <- ht1
print(draw(htlist1))
dev.off()

################################################################################
# Plotting the heatmaps over pseudotime trajectories
# subselect overlapping genes in the two datasets
# Select the top genes begin v end test across lineages from trade-seq

x <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/pseudotime/20220104/startres.tsv", sep="\t",row.names=1)
x$padj <- p.adjust(x$pvalue, method = "bonferroni", n = length(x$pvalue))

# Filter out genes with a true zero p-value and a p-adjusted value score below 0.05
x <- x[x$padj!=0,]
x <- x[x$padj<0.05,]

x$genes <- rownames(x)

z <- subset(lakorna, rownames(lakorna) %in% (x$genes))

list_of_cells1 = c("LSC","taLSC","BCE","WSCE")
list_of_cells3 = c("LSC","taLSC","CjC")
list_of_cells4 = c("LSC","taLSC","LNSC")
veclst <- list(list_of_cells1,list_of_cells3,list_of_cells4)#,

facs2 <- NULL

pdf(paste(resultsdir,'heatmaps_trajectories.pdf',sep="/") ,width=4,height=6,paper='special')

for (i in veclst){
  matz <- as.matrix(z)
  matz <- matz[,colnames(matz)%in%i]
  matz<-matz[,match(i, colnames(matz))]
  # Making the complex heatmap
  f1 = colorRamp2(c(-3, 0, 3), c("blue", "#EEEEEE", "red"), space = "RGB")

  ht1 = Heatmap(matz, col = f1, cluster_columns = F,cluster_rows = T, name = "Z-score RNA count",)
  htlist1 <- ht1
  print(draw(htlist1, column_title = "Pseudotime"))
}
dev.off()

################################################################################
# Generate the scANANSE barcode (bam/peak) files
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/meta_all.rds")

################################################################################
# write as function?
scANANSEdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/"

# Setting the scANANSE dir
path <- paste0(scANANSEdir,"RNA_peaks")
p
system(paste("mkdir -p ", path))

path <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/"
# Extracting metadata from sthe Seurat object
scRNA_metadata <- seur_obj@meta.data
# Generate barcode files for each cell population to be used in "bamslice" for each sample
scRNA_metadata$sample <- seur_obj$Condition
scRNA_metadata$sample

# Delete extra annotation to barcodes
vec <- rownames(scRNA_metadata)
logvec <- duplicated(sub(".*_", "", vec))

# Delete cells with a duplicate ID across datasets
scRNA_metadata <- scRNA_metadata[!logvec,]

vec <- rownames(scRNA_metadata)
vec <- sub(".*_", "", vec)
rownames(scRNA_metadata) <- vec

write.csv(scRNA_metadata, paste0(path, "/RNA_metadata.csv"), row.names=FALSE)

# Split the dataset only on condition and cell
for (z in unique(scRNA_metadata$sample)){
  for (i in unique(scRNA_metadata$costum_clustering)){
    df <- scRNA_metadata[i==scRNA_metadata$costum_clustering & scRNA_metadata$sample == z,]
    write.table(rownames(df), paste0(path, "/", i,"_", z,"_metadata_split.csv"),col.names = F, row.names=F, quote =F)
  }
}

# Generate the csv file with cells and barcodes
barcodes_csv<-c()
cells <- c()

for (z in unique(scRNA_metadata$sample)){
  for (i in unique(scRNA_metadata$costum_clustering)){
    cells <- c(cells,paste0(i,"_",z))
    barcodes_csv <- c(barcodes_csv,paste0(path, "/", i,"_",z,"_metadata_split.csv"))
  }
}
df <- as.data.frame(x=barcodes_csv,row.names = cells)
df$library_id <- cells

write.table(df, paste0(path, "/","RNA_barcodes_split.csv"),col.names = T, row.names=F, quote =F,sep = ",")

# Split the dataset only on cell meta-atlas
  for (i in unique(scRNA_metadata$costum_clustering)){
    df <- scRNA_metadata[i==scRNA_metadata$costum_clustering,]
    write.table(rownames(df), paste0(path, "/", i,"_metadata_sample.csv"),col.names = F, row.names=F, quote =F)
  }

# Generate the csv file with cells and barcodes
barcodes_csv<-c()
cells <- c()

for (i in unique(scRNA_metadata$costum_clustering)){
  cells <- c(cells,i)
  barcodes_csv <- c(barcodes_csv,paste0(path, "/", i,"_metadata_sample.csv"))
}

df <- as.data.frame(x=barcodes_csv,row.names = cells)
df$library_id <- cells

write.table(df, paste0(path, "/","RNA_barcodes.csv"),col.names = T, row.names=F, quote =F,sep = ",")

################################################################################
# Generating heatmaps from motif analysis output from SCEPIA and gimme motifs on
# scATAC-seq



##########
coldata2 <- NULL
dds2 <- NULL
#i<- "cell17"

# Calculation padj value for each cell with DEseq2 & plotting the volcanoplots with DEseq2
pdf(paste(resultsdir,'Volcanoplots_datasets.pdf',sep="/") ,width=8,height=8,paper='special')

for (i  in unique(coldata$condition2)){
  print(i)
  coldata2 <- coldata[coldata$condition2==i,]
  lakovst2 <- lakovst[, rownames(coldata2)]
  if (sum(colMeans(lakovst2[1:4])==0)==0){ #this makes cells where no values are found in one of each datasets to be skipped
    dds2 <- DESeqDataSetFromMatrix(countData = round(lakovst2),
                                  colData = coldata2,
                                  design = ~ condition)
    lakovstcalc2 <- DESeq(dds2)
    lakovstcalc2<-results(lakovstcalc2)
    lakovstcalc3 <- as.data.frame(lakovstcalc2@listData$log2FoldChange,row.names = rownames(lakovst2))
    lakovstcalc3$padj <- lakovstcalc2@listData$padj
    colnames(lakovstcalc3) <- c("log2FoldChange","padj")
    print(head(lakovstcalc3))
    write.table(lakovstcalc3, file = paste0(resultsdir,'/', i, '_comp.tsv'),sep = "\t", quote = F)

    # printing the volcanoplots
    dataset1_vs_dataset2 <- lakovstcalc3

    dataset1_vs_dataset2$noms<-rownames(dataset1_vs_dataset2)
    # Create new categorical column ------------------------------------------------
    dataset1_vs_dataset2 <- dataset1_vs_dataset2 %>%
      mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                   log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                                   TRUE ~ "ns"))

    # Add colour, size and alpha (transparency) to volcano plot --------------------
    vircols<-viridis(n = 2)
    cols <- c("up" = vircols[2], "down" = "blue", "ns" = vircols[1])
    sizes <- c("up" = 2, "down" = 2, "ns" = 1)
    alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

    print(dataset1_vs_dataset2 %>%
      ggplot(aes(x = log2FoldChange,
                 y = -log10(padj),
                 fill = gene_type,
                 size = gene_type,
                 alpha = gene_type)) +
      geom_point(shape = 21, # Specify shape and colour as fixed local parameters
                 colour = "black") +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed") +
      geom_vline(xintercept = c(log2(0.5), log2(2)),
                 linetype = "dashed") +
      scale_fill_manual(values = cols) + # Modify point colour
      scale_size_manual(values = sizes) + # Modify point size
      scale_alpha_manual(values = alphas) + # Modify point transparency
      scale_x_continuous(breaks = c(seq(-10, 10, 2)),
                         limits = c(-10, 10)) +
      scale_y_continuous(limits = c(0, 10)) +
      geom_text_repel(
      data = subset(dataset1_vs_dataset2, padj < 0.05 & abs(log2FoldChange)>5),
      aes(label = noms),
      size = 3,
      box.padding = unit(0.35, "lines"),max.overlaps = 15,
      point.padding = unit(0.3, "lines"))+
      labs(fill="Gene type color",size="Size",title=i))
  }
}
dev.off()


# develop of volcanoplot
#dataset1_vs_dataset2 <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score_datasets20211220/cell4_comp.tsv", sep = '\t', header = TRUE, row.names=1,comment.char = "")

#TEST #dataset1_vs_dataset2 <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710MECESCpseudobulkpadj.tsv', sep = '\t', header = TRUE, row.names=1,comment.char = "")

# Create new categorical column ------------------------------------------------
dataset1_vs_dataset2 <- dataset1_vs_dataset2 %>%
  mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                               log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                               TRUE ~ "ns"))

# Add colour, size and alpha (transparency) to volcano plot --------------------
vircols<-viridis(n = 2)
cols <- c("up" = vircols[2], "down" = "blue", "ns" = vircols[1])
sizes <- c("up" = 2, "down" = 2, "ns" = 1)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

top20_up <- dataset1_vs_dataset2[dataset1_vs_dataset2$gene_type=="up",]
top20_up<-top20_up[order(top20_up$log2FoldChange,decreasing = T),]
top20_up<-rownames(top20_up)[1:20]

dataset1_vs_dataset2$noms<-rownames(dataset1_vs_dataset2)
dataset1_vs_dataset2 %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = gene_type,
             size = gene_type,
             alpha = gene_type)) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),
                     limits = c(-10, 10)) +
  scale_y_continuous(limits = c(0, 10)) +
  geom_text_repel(
    data = subset(dataset1_vs_dataset2, padj < 0.05 & abs(log2FoldChange)>5),
    aes(label = noms),
    size = 3,
    box.padding = unit(0.35, "lines"),max.overlaps = 15,
    point.padding = unit(0.3, "lines"))+
  labs(fill="Gene type color",alpha="Transparency",size="Size")
