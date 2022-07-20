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
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera_0.7/"

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
# load in all datasets
path2 <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/collin_kera/"

list.files (path = path2, pattern = "SRR")
nam3 <- NULL
v2 <- c()
for (i in list.files (path = path2, pattern = "SRR")){
  print(i)
  v <- paste0(path2,i,"/outs/filtered_feature_bc_matrix/")
  nam <- paste("A", i, sep = "")
  nam <- Read10X(data.dir = v)
  nam2 <- paste("B", i, sep = "")
  nam2 <- CreateSeuratObject(counts = nam, project = paste0("kera",i), min.cells = 3, min.features = 200)
  assign(paste0(i), nam2)
  v2 <- c(v2,i)
}

Combined <- SRR12386367
remove(nam)
remove(nam2)
remove(SRR12386367)

# merge the first with the rest of the datasets
seur_obj_all <- merge(Combined, y = SRR12386368, add.cell.ids = v2, project = "kera2022")

# Remove all separate objects
rm(list = ls()[grep("SRR", ls())])

# Perform quality control
seur_obj <- QCSeurObj(SeuratObject = seur_obj_all,pct_MT = 30)
rm(seur_obj_all)

seur_obj <- PCASeurObj(SeuratObject = seur_obj)

seur_obj <- QCSeurObjClustering(SeuratObject = seur_obj,dimensions=14)

ClusteringResSeurObj(SeuratObject = seur_obj,dimensions=14)

chosen_res <- 0.1

# Final chosen resolution
pdf(paste(resultsdir,'8b.clusdendrogram_0.1_14dims.pdf',sep="/") ,width=8,height=8,paper='special')
seur_obj <- FindNeighbors(seur_obj, dims = 1:14)
res <- chosen_res
cluster_variable_name <- paste0("RNA_snn_res.", res)
seur_obj <- FindClusters(seur_obj, verbose = FALSE, resolution = res, graph.name = "RNA_snn")
seur_obj <- BuildClusterTree(seur_obj)
PlotClusterTree(seur_obj, direction = "downwards")
DimPlot(seur_obj, label=TRUE, label.size = 6, group.by = 'seurat_clusters') + ggtitle("Louvain Clustering") + ggtitle(paste0("cluster resolution ", res))
DimPlot(seur_obj, label=TRUE, label.size = 6, group.by = 'orig.ident') + ggtitle("Louvain Clustering") + ggtitle(paste0("cluster resolution ", res))
dev.off()

############################################################################
# Rename your cluster numbers to cellnumbers

# Name them according to the meta-atlas

# Export the raw data for SVM prediction done on 11-05-2022
#test_dat<-CreateSeuratObject(counts=seur_obj@assays$RNA@counts)
#saveRDS(test_dat,paste0(resultsdir,"/test_dat_kera.rds"))

#### Create 5- fold cross-validation data ####
#install.packages("rBayesianOptimization") # only run on first install
#library(rBayesianOptimization)
# 
#source("Cross_Validation.R")
#Cross_Validation("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/training_labels_meta.csv", OutputDir = paste0(resultsdir,"/meta_cv/"))

# Load in the data
### Load in the predicted data with labels back to seur_obj

# # predict cell types based on the machine learning algorithm
meta_labels <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/predictions_kera/SVM_Pred_Labels.csv")
seur_obj@meta.data$meta_labels <- meta_labels$X0

# adding the probability of the labels
meta_labels_prob <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/predictions_kera/SVMrej_Prob.csv")
seur_obj@meta.data$meta_labels_prob <- meta_labels_prob$X0

# Only the labels with a prob above 0.7 (quick glance -> mostly LESC (limbal epithelial stem cells))
meta_labels_0.7 <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/jupyter_notebooks/predictions_kera/SVMrej_Pred_Labels.csv")
seur_obj@meta.data$meta_labels_0.7 <- meta_labels_0.7$X0

# Plot the tree and the final clustering
colsordered<-c("#83FF52CC","#9AFE42CC" ,"#26EDA5CC","#3E3891CC","#455CD0CC","#F3C83ACC","#FABA39CC","#FEAA33CC","#FE992CCC","#36AAF9CC","#7A0403CC","#FABA39CC","#1AE4B6CC","#A41301CC","#900A01CC",NA)
colsordered2<-c("#9AFE42CC","#455CD0CC","#F3C83ACC","#FE992CCC","#36AAF9CC","#1AE4B6CC","#A41301CC",NA)

DimPlot(seur_obj, label=F,label.size = 5, group.by = 'meta_labels_0.7',cols = colsordered2,order = rev(c("LESC","CE","CSSC","CF","EC","IC","CDH19+","Unknown"))) 

pdf(paste(resultsdir, "13.clustering_costum_final.pdf", sep = '/'), width = 6, height = 5)
#p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
p1 <- DimPlot(seur_obj, label=F,label.size = 8, group.by = 'orig.ident')
p2 <- DimPlot(seur_obj, label=F,label.size = 8, group.by = 'meta_labels') + ggtitle("Louvain Clustering")
p3 <- DimPlot(seur_obj, label=F,label.size = 5, group.by = 'meta_labels_0.7',cols = colsordered2,order = rev(c("LESC","CE","CSSC","CF","EC","IC","CDH19+","Unknown"))) + ggtitle("Louvain Clustering")
#p4 <- DimPlot(seur_obj, label=F,label.size = 5, group.by = 'meta_labels_prob') + ggtitle("Louvain Clustering")+ NoLegend()
p4 <- FeaturePlot(seur_obj, features = 'meta_labels_prob', cols =c('white','dodgerblue3'))
#
p1
p2
p3
p4
#print(multiplot(p1, p2, p3, cols = 3))
dev.off()

seur_obj$Condition <- seur_obj@meta.data$orig.ident

seur_KCN <- seur_obj

newdf <- NULL

for (xx in unique(seur_obj@meta.data$meta_labels)) {
  print(xx)
  for (j in unique(seur_obj$Condition)){
    val<-sum(seur_obj@meta.data$meta_labels==xx & seur_obj$Condition==j)
    #col<- dfcol[rownames(dfcol)%in% xx,]
    df <- data.frame(values=val,Dataset=j,cell=xx)
    newdf <- rbind(newdf, df)
  }
}

pdf(paste(resultsdir,'Cell_dist_datasets_epi_full.pdf',sep="/") ,width=5,height=5,paper='special')
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

# stacked bar with transitions
library(ggalluvial)
#> Loading required package: ggplot2
library(tidyverse)

newdf <- NULL

for (xx in unique(seur_obj@meta.data$meta_labels_0.7)) {
  print(xx)
  for (j in unique(seur_obj$Condition)){
    val<-sum(seur_obj@meta.data$meta_labels_0.7==xx & seur_obj$Condition==j)
    perc <-val*100/sum(seur_obj$Condition==j)
    #col<- dfcol[rownames(dfcol)%in% xx,]
    df <- data.frame(values=val,Dataset=j,cell=xx,perc=perc)
    newdf <- rbind(newdf, df)
  }
}
newdf%>%
  ggplot(aes(y = perc, x = Dataset, fill = as.character(cell))) +
  geom_flow(aes(alluvium = cell), alpha= .5, color = "white",
            curve_type = "linear", 
            width = .5) +
  geom_col(width = .5, color = "white") +
  scale_fill_brewer(palette = "RdBu")+
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())



seur_obj2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")

newdf2 <- NULL

for (xx in unique(seur_obj2$scVI_label)) {
  print(xx)
  for (j in unique(seur_obj2$Condition)){
    val<-sum(seur_obj2$scVI_label==xx & seur_obj2$Condition==j)
    perc <-val*100/sum(seur_obj2$Condition==j)
    #col<- dfcol[rownames(dfcol)%in% xx,]
    df <- data.frame(values=val,Dataset=j,cell=xx,perc=perc)
    newdf2 <- rbind(newdf2, df)
  }
}
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 12))


newdf2%>%
  ggplot(aes(y = perc, x = Dataset, fill = as.character(cell))) +
  geom_flow(aes(alluvium = cell), alpha= .3, color = "white",
            curve_type = "linear", 
            width = .5) +
  geom_col(width = .5, color = "white") +
  scale_fill_manual(values = mycolors)+
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())

set.seed(0)
data_bar <- data.frame(
  stringsAsFactors = F,
  Sample = rep(c("A", "B"), each = 10),
  Percentage = runif(20),
  Taxon = rep(1:10, by = 2)
)

# full_df
newdf2 <- newdf2[newdf2$Dataset=="Co"|newdf2$Dataset=="Ca",]

full <- rbind(newdf, newdf2)

# Import the right color scheme
dfcol <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220510/unannotated_meta_colors.tsv")
# Add annotated colors
dfcol$V3 <- c("Cj","CE","LE","LESC","CSSC","CF","LSC","SK","LE","Cj","SK","TSK","CE","CF","Mel","EC","Ves","IC","CE",'CDH19+',"MF")
full$col <- dfcol$V2[match(full$cell,dfcol$V3)]

# WIP
vec <- c(6,3,2,8,0,9,1,12,18,4,7,10,11,5,13,15,16,14,17,20,19)
vec <- vec+1
ord <- unique(dfcol[vec,]$V3)
ord <- c(ord,"Unknown")

dfcol<-dfcol %>% arrange(factor(V3, levels = ord)) 
dfcol <- dfcol[!duplicated(dfcol$V3),]

dfcol2 <- rbind(dfcol, c(21,"grey",NA))

#####
full<-full  %>% arrange(factor(Dataset, levels = c("Ca","Co","keraSRR12386367","keraSRR12386368")))%>% arrange(factor(cell, levels = ord))
df$derma <- factor(df$derma, levels = df$derma)
rownames(dfcol) <- NULL

colsordered<-c("#83FF52CC","#9AFE42CC" ,"#26EDA5CC","#3E3891CC","#455CD0CC","#F3C83ACC","#FABA39CC","#FEAA33CC","#FE992CCC","#36AAF9CC","#7A0403CC","#FABA39CC","#1AE4B6CC","#A41301CC","#900A01CC",NA)

pdf(paste(resultsdir,'stacked_bar_KCN.pdf',sep="/") ,width=5,height=5,paper='special')
full%>%
  ggplot(aes(y = perc, x = Dataset, fill = factor(cell,levels = ord),alluvium = cell, 
             stratum = factor(cell,levels = ord))) +
  geom_flow(data = full,alpha= .3, color = "white", #, fill = factor(cell,levels = ord) aes(alluvium = cell), 
            curve_type = "linear", 
            width = .3) +
  geom_col(width = .3, color = "white") +
  scale_fill_manual(values = colsordered,name = "cells")+
  #scale_fill_manual(values = mycolors)+
  scale_y_continuous(NULL, expand = c(0,0)) +
  scale_x_discrete(labels=c("keraSRR12386367" = "KCN_1", "keraSRR12386368" = "KCN_2"
                             ))+
  cowplot::theme_minimal_hgrid() +
  labs(colour="cells")+
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
dev.off()

## first need to calcualte the actual percentage by group
data_bar %>%
  group_by(Sample) %>%
  mutate(perc = Percentage* 100/sum(Percentage)) %>%
  ggplot(aes(y = perc, x = Sample, fill = as.character(Taxon))) +
  geom_flow(aes(alluvium = Taxon), alpha= .5, color = "white",
            curve_type = "linear", 
            width = .5) +
  geom_col(width = .5, color = "white") +
  scale_fill_manual(values = dfcol$V2)+
  #scale_fill_brewer(palette = "RdBu")+
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())



# Pseudobulk table for kera comparisons
# Splitting the pseudobulk datasets unbiased into two replicates for Z-score calculation and over pseudotime
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")

seur_obj <-subset(x=seur_obj, subset = (batch == "Co"|batch == "Ca"))
seur_obj$batch<-droplevels(seur_obj$batch)
seur_obj$batch

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

# putting values of artificial rep2 and rep1 into the pseudobulkdataframecolumns
for (sample in unique(sce_qc$sample)){

  sce_qc$reps <- "rep1"
  
  # use rep2 for Catala
  sce_qc$reps[sce_qc$rep == "CatalaGSM5651509" | sce_qc$rep == "CatalaGSM5651511" | sce_qc$rep == "CatalaGSM5651513" | sce_qc$rep == "CatalaGSM56511515"| sce_qc$rep == "CatalaGSM5651117"| sce_qc$rep == "CatalaGSM5651119"] <- "rep2"
  sce_qc$reps[sce_qc$rep == "CatalaGSM5651510" | sce_qc$rep == "CatalaGSM5651512" | sce_qc$rep == "CatalaGSM5651514" | sce_qc$rep == "CatalaGSM56511516"| sce_qc$rep == "CatalaGSM5651118"|sce_qc$rep == "CatalaGSM5651120"] <- "rep2"
  
  pseudobulk_df[ , paste0(sample,"_1")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep1"])
  pseudobulk_df[ , paste0(sample,"_2")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep2"])
  # }
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulk_reps_DE_datasets_markers_split2_kera.tsv'), sep = '\t',quote = F, row.names = F)




# DEG expression in-between populations
pseudobulk_df<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera_0.7/20220519/pseudobulk_reps_DE_datasets_markers_split2_kera.tsv',header = T,row.names = 1)

# Selecting the pseudobulk_df and adding the kera conditions
pseudobulk_df <- pseudobulk_df[,]
pseudobulk_df_meta <- pseudobulk_df



sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc$sample <- seur_obj@meta.data$meta_labels_0.7

sce_qc$sample <- factor(sce_qc$sample)
sce_qc$sample

sce_qc$rep <- seur_obj@meta.data$orig.ident
sce_qc$rep

for (sample in unique(sce_qc$sample)){
  pseudobulk_df[ , paste0(sample)] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample])
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

# Select only the conditions with enough cells
vec <- c("LESC","CSSC","CE")
pseudobulk_df <- pseudobulk_df[,colnames(pseudobulk_df)%in%vec]

write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulkdf_full_kera.tsv'), sep = '\t',quote = F, row.names = F)

pseudobulk_df <- NULL
pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

# putting values of artificial rep2 and rep1 into the pseudobulkdataframecolumns
for (sample in unique(sce_qc$sample)){
  pseudobulk_df[ , paste0(sample,"_1")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$rep == "keraSRR12386367"])
  pseudobulk_df[ , paste0(sample,"_2")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$rep == "keraSRR12386368"])
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

# Select only the conditions with enough cells
vec <- c("LESC_1","CSSC_1","CE_1","LESC_2","CSSC_2","CE_2")
pseudobulk_df <- pseudobulk_df[,colnames(pseudobulk_df)%in%vec]
colnames(pseudobulk_df)

# Alter names for the kera condition
colnames(pseudobulk_df) <- c("Kera_CE_1","Kera_CE_2","Kera_LESC_1","Kera_LESC_2","Kera_CSSC_1","Kera_CSSC_2")

# Subselect the meta atlas pseudobulk for conditions in keratoconus populations
vec <- c("LESC_1","CSSC_1","CE_1","LESC_2","CSSC_2","CE_2")
pseudobulk_df_meta <- pseudobulk_df_meta[,colnames(pseudobulk_df_meta)%in%vec]

pseudobulk_df$Gene <- rownames(pseudobulk_df)
pseudobulk_df_meta$Gene <- rownames(pseudobulk_df_meta)

pseudo_joined <- pseudobulk_df %>% left_join(pseudobulk_df_meta, by="Gene")
rownames(pseudo_joined) <- pseudo_joined$Gene
pseudo_joined$Gene<- NULL
  
lakocountfile3 <- pseudo_joined
coldata <- NULL

# conditions
v2 <- colnames(lakocountfile3)
v2 <- v2 %>% str_replace("_[0-9]", "")
j <- unique(v2)
b<-2 # Or some other number
j<-sapply(j, function (x) rep(x,b))
j<-as.vector(j)

# Change _1 & _2 to condition 1 and _3 & _4 to condition 2
c <- c("1","2")
d <- paste(j,c,sep="_")

# reps
cols <- colnames(lakocountfile3)

# type
xx <- rep("paired-end",length(cols))

coldata <- data.frame(cols=cols,condition=j,condition2=d,type=xx)

rownames(coldata) <- coldata$cols
coldata <- coldata[,c("condition","type")]

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# checking the data
head(lakovst,2)

coldata

lakovst <- as.matrix(lakocountfile3)

# setting the rowdata in coldata similar to the coldata in the count matrix
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(lakovst))

all(rownames(coldata) == colnames(lakovst))

lakovst <- lakovst[, rownames(coldata)]
all(rownames(coldata) == colnames(lakovst))
lakovst <- na.omit(lakovst)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(lakovst),
                              colData = coldata,
                              design = ~ condition)
dds

##################################################################
dds2 <- DESeq(dds)
dds2

lakovstcalc <- results(dds2)
lakovstcalc3 <- as.data.frame(lakovstcalc@listData$log2FoldChange,row.names = rownames(lakovst))

for (i in unique(dds2$condition)){
  for (j in unique(dds2$condition)){
    if (j != i){
      print(j)
      lakovstcalc <- results(dds2,contrast = c("condition", i,j))
      lakovstcalc3 <- as.data.frame(lakovstcalc@listData$log2FoldChange,row.names = rownames(lakovst))
      lakovstcalc3$padj <- lakovstcalc@listData$padj
      colnames(lakovstcalc3) <- c("log2FoldChange","padj")
      write.table(data.frame("resid"=rownames(lakovstcalc3),lakovstcalc3), file = paste0(resultsdir,"/", i,"_",j,'_pseudobulkpadj.tsv'),sep = "\t", quote = F,row.names=F,col.names = T)
    }
  }
}

# transpose and scale the matrix per row Z-score per gene per sample is scaled upon
vsd <- assay(vst(dds,blind = T))
Z <- t(scale(t(vsd)))
Z
Z_score <- as.data.frame(Z,row.names = rownames(lakovst))

# joining the columns on the means sequential
vec <- rownames(coldata)
x <- split(vec, ceiling(seq_along(vec)/2))
scoretable <- as.data.frame(do.call(cbind, lapply(x, function(i) rowMeans(Z_score[, i]))), row.names = rownames(lakovst))
colnames(scoretable) <- unique(coldata$condition)
write.table(scoretable, file = paste0(resultsdir, '/Zscoretable_markers_split2.tsv'),sep = "\t", quote = F,row.names = T,col.names = T)
scoretable1 <- scoretable

# DEG in kera vs meta (reference)
genes_df_CSSC_kera <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera_0.7/20220516/Kera_CSSC_CSSC_pseudobulkpadj.tsv", header = TRUE, sep = '\t',comment.char = "#", stringsAsFactors = F)
#/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera/20220516/Kera_CSSC_CSSC_pseudobulkpadj.tsv
genes_df_CSSC_kera<-genes_df_CSSC_kera[!is.na(genes_df_CSSC_kera$padj),]
genes_df_CSSC_kera <- genes_df_CSSC_kera[genes_df_CSSC_kera$padj < 0.05,]

# DEG in kera vs meta (reference)
genes_df_CE_kera <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera_0.7/20220516/Kera_CE_CE_pseudobulkpadj.tsv", header = TRUE, sep = '\t',comment.char = "#", stringsAsFactors = F)
#/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera/20220516/Kera_CSSC_CSSC_pseudobulkpadj.tsv
genes_df_CE_kera<-genes_df_CE_kera[!is.na(genes_df_CE_kera$padj),]
genes_df_CE_kera <- genes_df_CE_kera[genes_df_CE_kera$padj < 0.05,]

# DEG in kera vs meta (reference)
genes_df_LESC_kera <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera_0.7/20220516/Kera_LESC_LESC_pseudobulkpadj.tsv", header = TRUE, sep = '\t',comment.char = "#", stringsAsFactors = F)
#/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera/20220516/Kera_CSSC_CSSC_pseudobulkpadj.tsv
genes_df_LESC_kera<-genes_df_LESC_kera[!is.na(genes_df_LESC_kera$padj),]
genes_df_LESC_kera <- genes_df_LESC_kera[genes_df_LESC_kera$padj < 0.05,]

# plot as heatmaps

rownames(genes_df_CSSC_kera) <- genes_df_CSSC_kera$resid
genes_df_CSSC_kera$resid <- NULL
mat_RNA <- as.matrix(genes_df_CSSC_kera)

pdf(paste(resultsdir,'/complexheatmap_all_DEG_kera_stromal.pdf',sep="/") ,width=13,height=3,paper='special')
ht1 = Heatmap(t(mat_RNA), cluster_columns = T,cluster_rows = T, name = "Relative expression CSSC")
print(ht1)
dev.off()

rownames(genes_df_LESC_kera) <- genes_df_LESC_kera$resid
genes_df_LESC_kera$resid <- NULL
mat_RNA2 <- as.matrix(genes_df_LESC_kera)

pdf(paste(resultsdir,'/complexheatmap_all_DEG_kera_LESC.pdf',sep="/") ,width=8.5,height=3,paper='special')
ht2 = Heatmap(t(mat_RNA2), cluster_columns = T,cluster_rows = T, name = "Relative expression LESC")
print(ht2)
dev.off()

mat_LSC <- mat_RNA

rownames(genes_df_CE_kera) <- genes_df_CE_kera$resid
genes_df_CE_kera$resid <- NULL
mat_RNA3 <- as.matrix(genes_df_CE_kera)

pdf(paste(resultsdir,'/complexheatmap_all_DEG_kera_CE.pdf',sep="/") ,width=7.5,height=3,paper='special')
ht3 = Heatmap(t(mat_RNA3), cluster_columns = T,cluster_rows = T, name = "Relative expression CE")
print(ht3)
dev.off()


pdf(paste(resultsdir,'/complexheatmap_all_DEG_kera_all.pdf',sep="/") ,width=27,height=3,paper='special')
htlist <- ht1 + ht2 + ht3
draw(htlist)
dev.off()
# Select union of genes from lists to show upon Z-score table
vec <- unique(c(genes_df_CSSC_kera$resid,genes_df_CE_kera$resid,genes_df_LESC_kera$resid))
scoretable2 <- scoretable1[rownames(scoretable1)%in%vec,]
scoretable2 <- scoretable2[,c("LESC","Kera_LESC","CE","Kera_CE","CSSC","Kera_CSSC")]

# Generate the heatmap
mat_RNA <- as.matrix(scoretable2)

f5 = colorRamp2(c(-2, 0, 2), c("blue", "#EEEEEE", "red"), space = "RGB")
ht1 = Heatmap(mat_RNA, col = f5, cluster_columns = F,cluster_rows = T, name = "Relative expression")

pdf(paste(resultsdir,'/complexheatmap_all_DEG_kera.pdf',sep="/") ,width=4,height=20,paper='special')
print(ht1)
dev.off()

# add identities with different lists
tfs <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera/tfs_ananse.txt", header = F, sep = '\t',comment.char = "#", stringsAsFactors = F)
tfs <- as.vector(tfs$V1)


mat_tf <- scoretable2[rownames(scoretable2)%in%tfs,]
ht1 = Heatmap(mat_tf, col = f5, cluster_columns = F,cluster_rows = T, name = "Relative expression")

pdf(paste(resultsdir,'/complexheatmap_tfs.pdf',sep="/") ,width=4,height=4,paper='special')
print(ht1)
dev.off()

#### analysis of epigenomics factors
epig_factors <- read.table(file = "genes.txt",sep = '\t',header = T)

epig_factors <- epig_factors %>% group_by(Function) %>% arrange(desc(Function),.by_group = T)

epig_factors2 <- epig_factors[epig_factors$Protein.complex!="-",]

mat_epig <- scoretable2[rownames(scoretable2)%in%epig_factors2$HGNC.approved.symbol,]

ht1 = Heatmap(mat_epig, col = f5, cluster_columns = F,cluster_rows = T, name = "Relative expression")

pdf(paste(resultsdir,'/complexheatmap_epig.pdf',sep="/") ,width=4,height=2.5,paper='special')
print(ht1)
dev.off()

# ligands and receptors human_lr_pair_celltalkDB.txt
lig_rec <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera/human_lr_pair_celltalkDB.txt", header = T, sep = '\t',comment.char = "#", stringsAsFactors = F)

# ligands
mat_lig <-scoretable2[rownames(scoretable2)%in%lig_rec$ligand_gene_symbol,]

ht1 = Heatmap(mat_lig, col = f5, cluster_columns = F,cluster_rows = T, name = "Relative expression")

pdf(paste(resultsdir,'/complexheatmap_lig.pdf',sep="/") ,width=4,height=3,paper='special')
print(ht1)
dev.off()
# receptors
mat_rec <- scoretable2[rownames(scoretable2)%in%lig_rec$receptor_gene_symbol,]

ht1 = Heatmap(mat_rec, col = f5, cluster_columns = F,cluster_rows = T, name = "Relative expression")

pdf(paste(resultsdir,'/complexheatmap_rec.pdf',sep="/") ,width=4,height=3.5,paper='special')
print(ht1)
dev.off()

# Load in all marker genes
# check DEG and marker genes CSSC
markers<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/cluster_markers_all.csv', sep = ',', header = TRUE, row.names = 1)
data_new2 <- markers %>%                                      # Top N highest values by group
  arrange(desc(avg_log2FC)) %>% 
  group_by(cluster) #%>%
  # slice(1:5)
data_new2    

data_new3 <- data_new2[data_new2$cluster=="CSSC",]

# Check marker genes unbiased
data_new3$gene[data_new3$gene %in% genes_df_CSSC_kera$resid]

data_new3 <- data_new2[data_new2$cluster=="CF",]

# Check marker genes unbiased
data_new3$gene[data_new3$gene %in% genes_df_CSSC_kera$resid]

# Volcanoplot
genes_df_CSSC_kera <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera_0.7/20220516/Kera_CSSC_CSSC_pseudobulkpadj.tsv", header = TRUE, sep = '\t',comment.char = "#", stringsAsFactors = F)
#/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera/20220516/Kera_CSSC_CSSC_pseudobulkpadj.tsv
genes_df_CSSC_kera<-genes_df_CSSC_kera[!is.na(genes_df_CSSC_kera$padj),]

lakovstcalc3<-genes_df_CSSC_kera
rownames(lakovstcalc3)<-lakovstcalc3$resid
lakovstcalc3$resid<-NULL
#colnames(lakovstcalc3) <- c("log2FoldChange","padj")
print(head(lakovstcalc3))
write.table(lakovstcalc3, file = paste0(resultsdir,'/', "CSSC", '_comp.tsv'),sep = "\t", quote = F)

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

pdf(paste(resultsdir,'Volcanoplots_datasets.pdf',sep="/") ,width=8,height=6,paper='special')

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
        scale_x_continuous(breaks = c(seq(-10, 3, 1)),
                           limits = c(-10, 3)) +
        scale_y_continuous(limits = c(0, 6)) +
        geom_text_repel(
          data = subset(dataset1_vs_dataset2, padj < 0.05 & log2FoldChange>2),
          aes(label = noms),
          nudge_x = 8,
          segment.size  = 0.2,
          segment.color = "grey50",
          size = 3,
          box.padding = 0.5,force = 5,
          #box.padding = unit(0.35, "lines"),
          max.overlaps = 10)+#,
        geom_text_repel(
          data = subset(dataset1_vs_dataset2, padj < 0.05 & log2FoldChange< (-2)),
          aes(label = noms),
          nudge_x = -4,
          segment.size  = 0.2,
          segment.color = "grey50",
          size = 3,
          box.padding = 0.5,force = 5,
          #box.padding = unit(0.35, "lines"),
          max.overlaps = 10)+
        #point.padding = unit(0.3, "lines"))+
        labs(fill="Gene type color",size="Size",title="test"))
dev.off()

markers<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/cluster_markers_all.csv', sep = ',', header = TRUE, row.names = 1)
markers <- markers[markers$cluster%in%c("CSSC","SK","TSK","CF"),]
markers <- unique(markers$gene)
stromal_markers <- markers

# printing the volcanoplots
dataset1_vs_dataset2 <- lakovstcalc3

dataset1_vs_dataset2$noms<-rownames(dataset1_vs_dataset2)
# Create new categorical column ------------------------------------------------
`%!in%` <- Negate(`%in%`)
dataset1_vs_dataset2$noms[which(dataset1_vs_dataset2$padj <= 0.05 &rownames(dataset1_vs_dataset2) %in% stromal_markers)]

dataset1_vs_dataset2 <- dataset1_vs_dataset2 %>%
  mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 &rownames(dataset1_vs_dataset2)%!in% stromal_markers~ "up",
                               log2FoldChange <= -1 & padj <= 0.05 &rownames(dataset1_vs_dataset2)%!in% stromal_markers~ "down",
                               padj >= 0.05&rownames(dataset1_vs_dataset2) %in% stromal_markers~ "stromal marker",
                               padj <= 0.05&rownames(dataset1_vs_dataset2) %in% stromal_markers~ "marker down",
                               padj <= 0.05 &rownames(dataset1_vs_dataset2) %in% stromal_markers~ "marker up",
                               TRUE ~ "ns"))

cols <- c("up" = vircols[2], "down" = "blue", "ns" = "grey",
          "stromal marker" = vircols[1],"marker up"= "orange","marker down"="darkblue")
sizes <- c("up" = 2, "down" = 2, "ns" = 1,
           "stromal marker" = 2,"marker up" = 2, "marker down" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.3,
            "stromal marker"=1,"marker up" = 1, "marker down" = 1)

pdf(paste(resultsdir,'Volcanoplots_datasets_CSSC_no_name&markers.pdf',sep="/") ,width=4.5,height=4,paper='special')

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
        scale_x_continuous(breaks = c(seq(-10, 3, 1)),
                           limits = c(-10, 3)) +
        scale_y_continuous(limits = c(0, 6)) +
        labs(fill="Gene type color",size="Size",title="CSSC")+theme_minimal())
dev.off()

# Load in the markers for the limbal-epithelium
markers<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/cluster_markers_all.csv', sep = ',', header = TRUE, row.names = 1)
markers <- markers[markers$cluster%in%c("LSC","LESC","LE","CE","Cj"),]
markers <- unique(markers$gene)
epi_markers <- markers

# Volcanoplot
genes_df_LESC_kera <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera_0.7/20220516/Kera_LESC_LESC_pseudobulkpadj.tsv", header = TRUE, sep = '\t',comment.char = "#", stringsAsFactors = F)
#/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera/20220516/Kera_CSSC_CSSC_pseudobulkpadj.tsv
genes_df_LESC_kera<-genes_df_LESC_kera[!is.na(genes_df_LESC_kera$padj),]

lakovstcalc3<-genes_df_LESC_kera
rownames(lakovstcalc3)<-lakovstcalc3$resid
lakovstcalc3$resid<-NULL

# printing the volcanoplots
dataset1_vs_dataset2 <- lakovstcalc3

dataset1_vs_dataset2$noms<-rownames(dataset1_vs_dataset2)

dataset1_vs_dataset2$noms[which(dataset1_vs_dataset2$padj <= 0.05 &rownames(dataset1_vs_dataset2) %in% epi_markers)]

# Create new categorical column ------------------------------------------------
dataset1_vs_dataset2 <- dataset1_vs_dataset2 %>%
  mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 &rownames(dataset1_vs_dataset2)%!in% epi_markers~ "up",
                               log2FoldChange <= -1 & padj <= 0.05 &rownames(dataset1_vs_dataset2)%!in% epi_markers~ "down",
                               padj >= 0.05&rownames(dataset1_vs_dataset2) %in% epi_markers~ "limbal-epithelial marker",
                               padj <= 0.05&rownames(dataset1_vs_dataset2) %in% epi_markers~ "marker down",
                               padj <= 0.05 &rownames(dataset1_vs_dataset2) %in% epi_markers~ "marker up",
                               TRUE ~ "ns"))

cols <- c("up" = vircols[2], "down" = "blue", "ns" = "grey",
          "limbal-epithelial marker" = vircols[1],"marker up"= "orange","marker down"="darkblue")
sizes <- c("up" = 2, "down" = 2, "ns" = 1,
           "limbal-epithelial marker" = 2,"marker up" = 2, "marker down" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.3,
            "limbal-epithelial marker"=1,"marker up" = 1, "marker down" = 1)

pdf(paste(resultsdir,'Volcanoplots_datasets_LESC_no_name&markers.pdf',sep="/") ,width=4.5,height=4,paper='special')

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
        scale_x_continuous(breaks = c(seq(-7, 6, 1)),
                           limits = c(-7, 6)) +
        scale_y_continuous(limits = c(0, 6)) +
        labs(fill="Gene type color",size="Size",title="LESC")+theme_minimal())
dev.off()

# Volcanoplot
genes_df_CE_kera <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera_0.7/20220516/Kera_CE_CE_pseudobulkpadj.tsv", header = TRUE, sep = '\t',comment.char = "#", stringsAsFactors = F)
#/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera/20220516/Kera_CSSC_CSSC_pseudobulkpadj.tsv
genes_df_CE_kera<-genes_df_CE_kera[!is.na(genes_df_CE_kera$padj),]

lakovstcalc3<-genes_df_CE_kera
rownames(lakovstcalc3)<-lakovstcalc3$resid
lakovstcalc3$resid<-NULL

# printing the volcanoplots
dataset1_vs_dataset2 <- lakovstcalc3

dataset1_vs_dataset2$noms<-rownames(dataset1_vs_dataset2)
# Create new categorical column ------------------------------------------------

dataset1_vs_dataset2$noms[which(dataset1_vs_dataset2$padj <= 0.05 &rownames(dataset1_vs_dataset2) %in% epi_markers)]

dataset1_vs_dataset2 <- dataset1_vs_dataset2 %>%
  mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 &rownames(dataset1_vs_dataset2)%!in% epi_markers~ "up",
                               log2FoldChange <= -1 & padj <= 0.05 &rownames(dataset1_vs_dataset2)%!in% epi_markers~ "down",
                               padj >= 0.05&rownames(dataset1_vs_dataset2) %in% epi_markers~ "limbal-epithelial marker",
                               padj <= 0.05&rownames(dataset1_vs_dataset2) %in% epi_markers~ "marker down",
                               padj <= 0.05 &rownames(dataset1_vs_dataset2) %in% epi_markers~ "marker up",
                               TRUE ~ "ns"))

cols <- c("up" = vircols[2], "down" = "blue", "ns" = "grey",
          "limbal-epithelial marker" = vircols[1],"marker up"= "orange","marker down"="darkblue")
sizes <- c("up" = 2, "down" = 2, "ns" = 1,
           "limbal-epithelial marker" = 2,"marker up" = 2, "marker down" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.3,
            "limbal-epithelial marker"=1,"marker up" = 1, "marker down" = 1)

pdf(paste(resultsdir,'Volcanoplots_datasets_CE_no_name&markers.pdf',sep="/") ,width=4.5,height=4,paper='special')

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
        scale_x_continuous(breaks = c(seq(-8, 5, 1)),
                           limits = c(-8, 5)) +
        scale_y_continuous(limits = c(0, 6)) +
        labs(fill="Gene type color",size="Size",title="CE")+theme_minimal())
dev.off()

#######################################################################
# Integration on scVI UMAP of corneal meta-atlas

# set condition2 for the names of the corneal meta-atlas and KCN populations
# Generate new cluster labeling for KCN datsets:
current.cluster.ids <- c('CSSC','CF','CDH19+','CE','EC','IC','LESC','Unknown')

new.cluster.ids <- c('CSSC_KCN','CF_KCN','CDH19+_KCN','CE_KCN','EC_KCN','IC_KCN','LESC_KCN','Unknown')

seur_obj$meta_labels_0.7 <- plyr::mapvalues(x = as.factor(seur_obj$meta_labels_0.7), from = current.cluster.ids, to = new.cluster.ids)
seur_obj$meta_labels_0.7 <- as.factor(seur_obj$meta_labels_0.7)
seur_obj$Condition2 <- seur_obj$meta_labels_0.7

seur_meta <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")
seur_meta$Condition2 <- seur_meta$scVI_label

# Integrate the meta-atlas with scVI labels with predicted labels SVM
fullcornea_kcn <- merge(seur_meta, y = seur_obj, add.cell.ids = c("meta-atlas", "KCN"), project = "fullcornea_KCN")
factor(fullcornea_kcn$Condition2)

# Setting the correct parameters for the anndata conversion
fullcornea_kcn$batch <-unname(fullcornea_kcn$Condition2)
fullcornea_kcn@meta.data$batch <- unname(fullcornea_kcn$Condition2)
# generate training data for the machine learning model for every group

# Save it as a separate object raw integration of the four datasets
saveRDS(fullcornea_kcn, file = paste0(resultsdir,"/fullcornea_kcn.rds"))

# Embed this object into scVI integration (see jupyter notebook)

#######################################################################
# Single cell stromal marker gene heatmap
markers<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/cluster_markers_all.csv', sep = ',', header = TRUE, row.names = 1)
markers <- markers[markers$cluster%in%c("CSSC","SK","TSK","CF"),]
markers <- unique(markers$gene)

seur_obj <- seur_KCN
seur_obj2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")

seur_obj <-subset(x=seur_obj, subset = (meta_labels_0.7 == "CSSC"|meta_labels_0.7 == "SK"|meta_labels_0.7 =="TSK"| meta_labels_0.7 =="CF"))

cluster_order <- c('CSSC_KCN','CF_KCN')

# Generate new cluster labeling:
current.cluster.ids <- c('CSSC','CF')

new.cluster.ids <- c('CSSC_KCN','CF_KCN')

seur_obj$meta_labels_0.7 <- plyr::mapvalues(x = as.factor(seur_obj$meta_labels_0.7), from = current.cluster.ids, to = new.cluster.ids)
seur_obj$meta_labels_0.7 <- as.factor(seur_obj$meta_labels_0.7)

DoHeatmap(seur_obj, features = markers,group.by = "meta_labels_0.7") + NoLegend()

seur_obj$comp <- seur_obj$meta_labels_0.7
# Load in the meta-atlas with the corneal populations
#seur_obj2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")
seur_obj2 <-subset(x=seur_obj2, subset = (Condition == "Co")) #|Condition == "Ca"
seur_obj2 <-subset(x=seur_obj2, subset = (scVI_label == "CSSC"|scVI_label == "SK"|scVI_label =="TSK"| scVI_label =="CF"))

seur_obj2 <-DietSeurat(seur_obj2, counts = TRUE, data = TRUE, scale.data = FALSE)
seur_obj2 <- ScaleData(seur_obj2,features = markers)

seur_obj2$scVI_label <- droplevels(seur_obj2$scVI_label)
DoHeatmap(seur_obj2, features = markers,group.by ="scVI_label") + NoLegend()
seur_obj2$comp <- seur_obj2$scVI_label

# Join the single-cells of KCN & corneal meta-atlas to get rescaled values
big <- merge(seur_obj, y = seur_obj2, add.cell.ids = c("KCN", "meta"), project = "comp")
factor(big$comp)
big <- ScaleData(big,features = markers)

big@active.ident <- factor(big@meta.data$comp)
big@active.ident <- factor(x = big@active.ident, levels = c("CSSC_KCN","CSSC","CF")) #,"SK","TSK"
big$comp <- big@active.ident

meta_dat <- big@meta.data

pdf(paste(resultsdir,'single_cell.pdf',sep="/") ,width=15,height=15,paper='special')
print(DoHeatmap(big, features = markers,group.by ="comp",label=FALSE) + NoLegend())
dev.off()

pdf(paste(resultsdir,'single_cell_long.pdf',sep="/") ,width=15,height=80,paper='special')
print(DoHeatmap(big, features = markers,group.by ="comp",label=FALSE) + NoLegend())
dev.off()

# Subselecting transcription factors
tfs <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/factors_stromal.txt", header = F, sep = '\t',comment.char = "#", stringsAsFactors = F)
tfs <- as.vector(tfs$V1)
markers <- tfs

big <- ScaleData(big,features = markers)

big@active.ident <- factor(big@meta.data$comp)
big@active.ident <- factor(x = big@active.ident, levels = c("CSSC_KCN","CSSC","CF")) #,"SK","TSK"
big$comp <- big@active.ident

pdf(paste(resultsdir,'single_cell_tfs.pdf',sep="/") ,width=15,height=15,paper='special')
print(DoHeatmap(big, features = markers,group.by ="comp",label=FALSE) + NoLegend())
dev.off()

pdf(paste(resultsdir,'single_cell_long_tfs.pdf',sep="/") ,width=15,height=120,paper='special')
print(DoHeatmap(big, features = markers,group.by ="comp",label=FALSE) + NoLegend())
dev.off()

#########################################################
# Single cell limbal-epithelial marker gene heatmap
markers<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/cluster_markers_all.csv', sep = ',', header = TRUE, row.names = 1)
markers <- markers[markers$cluster%in%c("LSC","LESC","LE","CE","Cj"),]
markers <- unique(markers$gene)

seur_obj <- seur_KCN
seur_obj2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")

seur_obj <-subset(x=seur_obj, subset = (meta_labels_0.7 == "LESC"|meta_labels_0.7 == "CE"))

cluster_order <- c('LESC_KCN','CE_KCN')

# Generate new cluster labeling:
current.cluster.ids <- c('LESC','CE')

new.cluster.ids <- c('LESC_KCN','CE_KCN')

seur_obj$meta_labels_0.7 <- plyr::mapvalues(x = as.factor(seur_obj$meta_labels_0.7), from = current.cluster.ids, to = new.cluster.ids)
seur_obj$meta_labels_0.7 <- as.factor(seur_obj$meta_labels_0.7)

DoHeatmap(seur_obj, features = markers,group.by = "meta_labels_0.7") + NoLegend()

seur_obj$comp <- seur_obj$meta_labels_0.7
# Load in the meta-atlas with the corneal populations
#seur_obj2 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")
seur_obj2 <-subset(x=seur_obj2, subset = (Condition == "Co")) #|Condition == "Ca"
seur_obj2 <-subset(x=seur_obj2, subset = (scVI_label == "LSC"|scVI_label == "LESC"|scVI_label =="LE"| scVI_label =="Cj"| scVI_label =="CE"))

seur_obj2 <-DietSeurat(seur_obj2, counts = TRUE, data = TRUE, scale.data = FALSE)
seur_obj2 <- ScaleData(seur_obj2,features = markers)

seur_obj2$scVI_label <- droplevels(seur_obj2$scVI_label)
DoHeatmap(seur_obj2, features = markers,group.by ="scVI_label") + NoLegend()
seur_obj2$comp <- seur_obj2$scVI_label

# Join the single-cells of KCN & corneal meta-atlas to get rescaled values
big <- merge(seur_obj, y = seur_obj2, add.cell.ids = c("KCN", "meta"), project = "comp")
factor(big$comp)
big <- ScaleData(big,features = markers)

big@active.ident <- factor(big@meta.data$comp)
big@active.ident <- factor(x = big@active.ident, levels = c("LESC_KCN","CE_KCN","LSC","LESC","LE","Cj","CE")) #,"SK","TSK"
big$comp <- big@active.ident

meta_dat <- big@meta.data

pdf(paste(resultsdir,'single_cell_epi.pdf',sep="/") ,width=25,height=15,paper='special')
print(DoHeatmap(big, features = markers,group.by ="comp",label=FALSE) + NoLegend())
dev.off()

pdf(paste(resultsdir,'single_cell_long_epi.pdf',sep="/") ,width=25,height=80,paper='special')
print(DoHeatmap(big, features = markers,group.by ="comp",label=FALSE) + NoLegend())
dev.off()

########################################################################################
# Differential gene expression of CSSC

# Pseudobulk table for kera comparisons
# Splitting the pseudobulk datasets unbiased into two replicates for Z-score calculation and over pseudotime
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220228/4datasets_annotated_joined.rds")

seur_obj <-subset(x=seur_obj, subset = (scVI_label == "CSSC"))

seur_obj <-subset(x=seur_obj, subset = (batch == "Co"))
seur_obj$batch<-droplevels(seur_obj$batch)
seur_obj$batch

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

write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulkdf_CSSC.tsv'), sep = '\t',quote = F, row.names = F)

pseudobulk_df <- NULL
pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

# putting values of artificial rep2 and rep1 into the pseudobulkdataframecolumns

sce_qc$reps <- NULL
unique(sce_qc$rep)
# use rep2 for CollinSRR12386368
sce_qc$reps[sce_qc$rep == "CollinSRR12386359"|sce_qc$rep == "CollinSRR12386361"|sce_qc$rep == "CollinSRR12386363"|sce_qc$rep == "CollinSRR12386368"] <- "rep1"
sce_qc$reps[sce_qc$rep == "CollinSRR12386360"|sce_qc$rep == "CollinSRR12386362"|sce_qc$rep == "CollinSRR12386367"] <- "rep2"
  
pseudobulk_df[ , paste0(sample,"_1")] <- rowSums(counts(sce_qc)[,sce_qc$reps == "rep1"])
pseudobulk_df[ , paste0(sample,"_2")] <- rowSums(counts(sce_qc)[, sce_qc$reps == "rep2"])

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL

write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulk_reps_DE_datasets_markers_split2_kera_collin.tsv'), sep = '\t',quote = F, row.names = F)

# DEG expression in-between populations
pseudobulk_df<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_collin_kera_0.7/20220623/pseudobulk_reps_DE_datasets_markers_split2_kera_collin.tsv',header = T,row.names = 1)

pseudobulk_df_meta2 <- pseudobulk_df

seur_obj <- seur_KCN

seur_obj <-subset(x=seur_obj, subset = (meta_labels_0.7 == "CSSC"))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

sce_qc$sample <- seur_obj$meta_labels_0.7

sce_qc$sample <- factor(sce_qc$sample)
sce_qc$sample

sce_qc$rep <- seur_obj@meta.data$orig.ident
sce_qc$rep

sce_qc$reps <- NULL
unique(sce_qc$rep)

pseudobulk_df <- NULL
pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

pseudobulk_df[ , paste0(sample,"_1")] <- rowSums(counts(sce_qc)[,sce_qc$rep == "keraSRR12386367"])
pseudobulk_df[ , paste0(sample,"_2")] <- rowSums(counts(sce_qc)[, sce_qc$rep == "keraSRR12386368"])

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL
colnames(pseudobulk_df) <- c("CSSC_KCN_1","CSSC_KCN_2")
pseudobulk_df$Gene<-rownames(pseudobulk_df)
pseudobulk_df_meta2$Gene<-rownames(pseudobulk_df_meta2)

#WIP
# Continue from here!
pseudo_joined<- NULL
pseudo_joined <- pseudobulk_df %>% left_join(pseudobulk_df_meta2, by="Gene")

rownames(pseudo_joined) <- pseudo_joined$Gene
pseudo_joined$Gene<- NULL

lakocountfile3 <- pseudo_joined[c("CSSC_KCN_1","CSSC_KCN_2","CSSC_1","CSSC_2")]
coldata <- NULL

# conditions
v2 <- colnames(lakocountfile3)
v2 <- v2 %>% str_replace("_[0-9]", "")
j <- unique(v2)
b<-2 # Or some other number
j<-sapply(j, function (x) rep(x,b))
j<-as.vector(j)

# Change _1 & _2 to condition 1 and _3 & _4 to condition 2
c <- c("1","2")
d <- paste(j,c,sep="_")

# reps
cols <- colnames(lakocountfile3)

# type
xx <- rep("paired-end",length(cols))

coldata <- data.frame(cols=cols,condition=j,condition2=d,type=xx)

rownames(coldata) <- coldata$cols
coldata <- coldata[,c("condition","type")]

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# checking the data
head(lakovst,2)

coldata

lakovst <- as.matrix(lakocountfile3)

# setting the rowdata in coldata similar to the coldata in the count matrix
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(lakovst))

all(rownames(coldata) == colnames(lakovst))

lakovst <- lakovst[, rownames(coldata)]
all(rownames(coldata) == colnames(lakovst))
lakovst <- na.omit(lakovst)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(lakovst),
                              colData = coldata,
                              design = ~ condition)
dds

##################################################################
dds2 <- DESeq(dds)
dds2

lakovstcalc <- results(dds2)
lakovstcalc3 <- as.data.frame(lakovstcalc@listData$log2FoldChange,row.names = rownames(lakovst))

for (i in unique(dds2$condition)){
  for (j in unique(dds2$condition)){
    if (j != i){
      print(j)
      lakovstcalc <- results(dds2,contrast = c("condition", i,j))
      lakovstcalc3 <- as.data.frame(lakovstcalc@listData$log2FoldChange,row.names = rownames(lakovst))
      lakovstcalc3$padj <- lakovstcalc@listData$padj
      colnames(lakovstcalc3) <- c("log2FoldChange","padj")
      write.table(data.frame("resid"=rownames(lakovstcalc3),lakovstcalc3), file = paste0(resultsdir,"/", i,"_",j,'_pseudobulkpadj.tsv'),sep = "\t", quote = F,row.names=F,col.names = T)
    }
  }
}

# transpose and scale the matrix per row Z-score per gene per sample is scaled upon
vsd <- assay(vst(dds,blind = T))
Z <- t(scale(t(vsd)))
Z
Z_score <- as.data.frame(Z,row.names = rownames(lakovst))

# joining the columns on the means sequential
vec <- rownames(coldata)
x <- split(vec, ceiling(seq_along(vec)/2))
scoretable <- as.data.frame(do.call(cbind, lapply(x, function(i) rowMeans(Z_score[, i]))), row.names = rownames(lakovst))
colnames(scoretable) <- unique(coldata$condition)

# DEG all general genes?
