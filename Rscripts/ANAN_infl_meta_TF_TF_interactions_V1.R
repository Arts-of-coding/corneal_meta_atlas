# Loading in important packages
library(chorddiag)
library(circlize)
library(viridis)
library(ggvenn)
#install.packages('venn')
library(venn)
library(dplyr)
library(tidyr)
library(tidyselect)
#devtools::install_github("hms-dbmi/UpSetR")
library(UpSetR)
library(ComplexHeatmap)
library(ggplot2)
library(reshape)
#install.packages("dtwclust")
library(dtwclust)
library(clustree)
#BiocManager::install("dorothea")
library(dorothea)

#####################################
# Storing results
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/ANANSE_TF_TF/"

# Setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

#First we load the necessary packages:
#  ```{r load, message=FALSE}
library(dorothea)
library(ggplot2)
library(dplyr)
#```

#Here is how to retrieve all regulons from human:
#  ```{r model}
net <- dorothea::dorothea_hs
head(net)
#```

n_genes <- net %>%
  group_by(tf) %>%
  summarize(n = n())
ggplot(data=n_genes, aes(x=n)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  xlab('Number of target genes') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")

n_edges <- net %>%
  group_by(confidence) %>%
  summarize(n = n())
ggplot(data=n_edges, aes(x=confidence, y=log10(n), color=confidence, fill=confidence)) +
  geom_bar(stat="identity") +
  theme(text = element_text(size=12)) +
  xlab('log10(Number of edges)') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")

prop <- net %>%
  group_by(tf, mor) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(mor == 1)
ggplot(data=prop, aes(x=freq)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  xlab('% of positive edges') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")

corneal_TFs <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/ANANSE_TF_TF/20220708/corneal_tfs.csv")
corneal_TFs <- corneal_TFs$x
stromal_TFs <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/ANANSE_TF_TF/20220708/stromal_tfs.csv")
stromal_TFs <- stromal_TFs$x

net_corneal <- net[net$tf %in% corneal_TFs,]
net_corneal <- net_corneal[net_corneal$target %in% corneal_TFs,]

write.table(net_corneal, file = paste0(resultsdir,'/net_epi_nofilt.csv'), sep = ',',row.names = F)

net_corneal <- net_corneal[net_corneal$confidence != "E" & net_corneal$confidence != "D",]

#Output it to a df for network analysis
write.table(net_corneal, file = paste0(resultsdir,'/net_epi.csv'), sep = ',',row.names = F)

net_stromal <- net[net$tf %in% stromal_TFs,]
net_stromal <- net_stromal[net_stromal$target %in% stromal_TFs,]


write.table(net_stromal, file = paste0(resultsdir,'/net_stromal_nofilt.csv'), sep = ',',row.names = F)

net_stromal <- net_stromal[net_stromal$confidence != "E" & net_stromal$confidence != "D",]

#Output it to a df for network analysis
write.table(net_stromal, file = paste0(resultsdir,'/net_stromal.csv'), sep = ',',row.names = F)

# single cell test
LSC_TFs <- fullheat$LSC$factor[fullheat$LSC$factor %in% corneal_TFs]
net_LSC <- net[net$tf %in% LSC_TFs,]
net_LSC <- net_LSC[net_LSC$target %in% LSC_TFs,]
net_LSC <- net_LSC[net_LSC$confidence != "E" & net_LSC$confidence != "D",]

#Output it to a df for network analysis
write.table(net_LSC, file = paste0(resultsdir,'/net_LSC.csv'), sep = ',',row.names = F)

# Load in networks with DEG ESC to stromal & ESC to epi 

















write.table(corneal_TFs, file = paste0(resultsdir,'/corneal_tfs.csv'), sep = ',',row.names = F)
write.table(stromal_TFs, file = paste0(resultsdir,'/stromal_tfs.csv'), sep = ',',row.names = F)




# Read in the files to analyse with chord
files <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE_31032022/ANANSE_files_ESCcomparison_0804_0.1_sel_250.csv",header = T, comment.char = '#') # meta_no_sel_union

# Read in transcription factor annotation files if they are available
expected <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/ANANSE/Expected_tfs.csv",header = T, comment.char = '#')
exgeneral <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/ANANSE/Expectedgeneral_tfs.csv",header = T, comment.char = '#')

# Generate the interaction dataframe
factorsheat1 <- NULL
fullheat <- list()
fullheat2 <- list()
for (i in unique(files$cell_id)){
  factorsheat1 <- read.table(files$network[files$cell_id==i],header = T)
  factorsheat2 <- factorsheat1[,c("factor","direct_targets","influence_score")]
  colnames(factorsheat2) <- c("factor","direct_targets",i)
  fullheat2[[i]] <- factorsheat2
  factorsheat2 <- factorsheat1[,c("factor","influence_score")]
  colnames(factorsheat2) <- c("factor",i)
  fullheat[[i]] <-factorsheat2
}

merged.data.frame = Reduce(function(...) merge(..., all=T), fullheat)
merged.data.frame[merged.data.frame == 0] <- NA

merged.data.frame2 <- merged.data.frame
merged.data.frame2[!is.na(merged.data.frame2)] <- 1
merged.data.frame2[is.na(merged.data.frame2)] <- 0
merged.data.frame2$factor <- merged.data.frame$factor
merged.data.frame2 <-merged.data.frame2[rowSums(merged.data.frame2[, -1])>0, ]

for(i in 2:ncol(merged.data.frame2)){ merged.data.frame2[ , i] <- as.integer(merged.data.frame2[ , i]) }
colnoms <- colnames(merged.data.frame2)[colnames(merged.data.frame2) != "factor"]

fill = viridis(length(files$cell_id),alpha = 0.8,option = "H")

pdf(paste(resultsdir,'/intersectionblob.pdf',sep="/") ,width=10,height=5,paper='special')
upset(merged.data.frame2,nintersects = NA,
      sets = colnoms, 
      order.by="freq", matrix.color="black", point.size=1,
      sets.bar.color=fill)
dev.off()

# add this!: set_order = c("a", "b", "c")

# Extracting unique factors corneal epithelium
lstnames <- c()
fullheatfacts= list()
for (i in unique(files$cell_id)[1:5]){
  lstnames <- c(lstnames,paste0(merged.data.frame2$i))
  factorsheat3 <- read.table(files$network[files$cell_id==i],header = T)
  factorsheat4 <- factorsheat3[,"factor"]
  fullheatfacts[[i]] <- factorsheat4
}

# Extracting the shared factors of epithelial
sharedallep <- Reduce(f = intersect, x = fullheatfacts)

# Extracting unique factors 
lstnames <- c()
fullheatfacts= list()
for (i in unique(files$cell_id)){
  lstnames <- c(lstnames,paste0(merged.data.frame2$i))
  factorsheat3 <- read.table(files$network[files$cell_id==i],header = T)
  factorsheat4 <- factorsheat3[,"factor"]
  fullheatfacts[[i]] <- factorsheat4
}

# Extracting the shared factors of all populations
sharedall <- Reduce(f = intersect, x = fullheatfacts)

###################
# Extracting unique factors for each cell population to list
uniq <- list()
for (v in unique(files$cell_id)){
  fullheatfacts2 <- fullheatfacts
  fullheatfacts2[[v]] <- NULL
  fac_interest <- fullheatfacts[[v]]
  `%notin%` <- Negate(`%in%`)
  uniq[[v]] <- fac_interest[fac_interest %notin% unlist(fullheatfacts2)]
}

###############################################
# Generating the complex heatmaps quantitation for all unique factors + shared factors and influence scores

## load in the z-score normalized dataset:
lakorna <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score_datasets20220225/Zscoretable_markers_split2.tsv", sep="\t",row.names=1)

lakorna <- na.omit(lakorna)

# for filtering
lakorna <- lakorna[, !names(lakorna) %in% c("IC", "Mel", "Ves")]

# Combine all lists into one giant dataframe
merged.data.frame = Reduce(function(...) merge(..., all=T), fullheat)

# heatdfshared <- subset(merged.data.frame, merged.data.frame$factor %in% sharedall)
# uniqtfs <- unlist(uniq)
# heatdfuniq <- subset(merged.data.frame, merged.data.frame$factor %in% uniqtfs)
# rownames(heatdfuniq) <- heatdfuniq$factor
# heatdfuniq$factor <- NULL
# 
# rownames(heatdfshared) <- heatdfshared$factor
# heatdfshared$factor <- NULL
# 
# # Say influence score is 0 for uniq factors that have NA
# heatdfuniq[is.na(heatdfuniq)] <- 0
# 
# matuniq <- as.matrix(heatdfuniq)
# matshared <- as.matrix(heatdfshared)
# 
# f1 = colorRamp2(c(0, 1), c("white","purple"), space = "RGB")
# f2 = colorRamp2(c(0, 1), c("white", "red"), space = "RGB")
# 
# ht1 = Heatmap(matuniq, col= f1,cluster_rows = T, cluster_columns = T, name = "influence score")
# ht2 = Heatmap(matshared, col= f1,cluster_rows = T, cluster_columns = T, name = "influence score")
# 
# dev.off()
# pdf(paste(resultsdir,'/complexheatmap.pdf',sep="/") ,width=15,height=8,paper='special')
# ht1
# ht2
# dev.off()

# Generating a dataframe from influence scores
heatdfall <- merged.data.frame
rownames(heatdfall) <- heatdfall$factor
heatdfall$factor <- NULL
heatdfall[is.na(heatdfall)] <- 0

# Selecting a cutoff for the influence values
maxvec <- unname(apply(heatdfall[,], 1, max))
heatdfall$maxscore <- maxvec
heatdfall= heatdfall[heatdfall$maxscore > 0.3,]
heatdfall$maxscore <- NULL

# Selecting tfs based on the cutoff in the Z-score table
z <- subset(lakorna, rownames(lakorna) %in% (merged.data.frame$factor))
z <- z[rownames(z) %in% rownames(heatdfall),]

# Generate the matrix 
matall <- as.matrix(heatdfall)
matz <- as.matrix(z)

ht4 = Heatmap(matall, col= viridis(100), cluster_rows = T, cluster_columns = F, name = "Influence score")

roword <- row_order(ht4)
matall <- matall[roword,]
colord <- column_order(ht4)
matall <- matall[, colord]
matz <- matz[rownames(matall),colnames(matall)]

ht4 = Heatmap(matall, col= viridis(100), cluster_rows = F, cluster_columns = F, name = "Influence score")

# 0.75 final score
####### factors that have a score >0.75 in at least one celpopulation corneal epithelium
# Selecting a cutoff for the influence values

heatdfall3 <- heatdfall
heatdfall3 <- heatdfall[,1:5]
maxvec <- unname(apply(heatdfall3[,], 1, max))
heatdfall3$maxscore <- maxvec
heatdfall3= heatdfall3[heatdfall3$maxscore > 0.75,]
heatdfall3$maxscore <- NULL

matall3 <- as.matrix(heatdfall3)
# if(length(vec5) != 0) {
#vec5 <- rownames(matz[rownames(matz) %in% vec5 ,])
#vec6 <- which(rownames(matz) %in% vec5, arr.ind = T)
# plotting the list

DE_epi_stromal <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/RNA_intra_cpm/20220421/epi_stromal_pseudobulkpadj.tsv", sep="\t",row.names=1)

pdf(paste(resultsdir,'/', paste0('score_sel_0.75_ep.pdf'),sep="/") ,width=4,height=20,paper='special')
ht30 = Heatmap(matall3, col= viridis(100), cluster_rows = T, cluster_columns = F, name = "Influence score")
htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('0.75'))

dev.off()

vircols <- viridis(100)

fvir = colorRamp2(c(0.7,0.8,0.9,1), c(vircols[1],vircols[33],vircols[66],vircols[100]), space = "RGB")
fvir2 = colorRamp2(c(0.5,0.625,0.75,0.875,1), c(vircols[1],vircols[25],vircols[50],vircols[75],vircols[100]), space = "RGB")


pdf(paste(resultsdir,'/', paste0('score_sel_0.75_epi_cols.pdf'),sep="/") ,width=4,height=18,paper='special')
ht30 = Heatmap(matall3, col= fvir, cluster_rows = T, cluster_columns = F, name = "Influence score")
htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('0.75'))

dev.off()


pdf(paste(resultsdir,'/', paste0('score_sel_0.75_epi_rot.pdf'),sep="/") ,width=20,height=3.5,paper='special')
matall4 <- t(matall3)

f4 = colorRamp2(c(-5, 0, 5), c("orange", "#EEEEEE", "darkolivegreen"), space = "RGB")

f5 = colorRamp2(c(1e-09, 0.05), c("red", "blue"), space = "RGB")

DE_epi_stromal <- t(DE_epi_stromal)

DE_epi_stromal <- DE_epi_stromal[,colnames(DE_epi_stromal)%in%colnames(matall4)]
DE_epi_stromal <- DE_epi_stromal[,colnames(matall4)]

DE_epi_stromal <- as.data.frame(DE_epi_stromal)

ha = HeatmapAnnotation(
  log2FC_epi_stromal = unlist(unname(DE_epi_stromal[1,])), 
  padj = unlist(unname(DE_epi_stromal[2,])),
  col = list(log2FC_epi_stromal = f4,
             padj = f5)
  ,
  annotation_name_side = "left"
)
ht30 = Heatmap(matall4, col= fvir, cluster_rows = F, cluster_columns = T, name = "Influence score",column_names_rot = 90, row_names_side = "left",top_annotation = ha)
htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('0.75'))

dev.off()

pdf(paste(resultsdir,'/', paste0('score_sel_0.75_epi_rot_0.5cut.pdf'),sep="/") ,width=19,height=3,paper='special')
matall4 <- t(matall3)
ht30 = Heatmap(matall4, col= fvir2, cluster_rows = F, cluster_columns = T, name = "Influence score",column_names_rot = 90, row_names_side = "left")
htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('0.75'))

dev.off()

DE_epi_stromal <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/RNA_intra_cpm/20220421/epi_stromal_pseudobulkpadj.tsv", sep="\t",row.names=1)

pdf(paste(resultsdir,'/', paste0('score_sel_0.75_epi_rot_0.5cut_DE.pdf'),sep="/") ,width=20,height=3.5,paper='special')
matall4 <- t(matall3)

f4 = colorRamp2(c(-5, 0, 5), c("orange", "#EEEEEE", "darkolivegreen"), space = "RGB")

f5 = colorRamp2(c(1e-09, 0.05), c("red", "blue"), space = "RGB")

DE_epi_stromal <- t(DE_epi_stromal)

DE_epi_stromal <- DE_epi_stromal[,colnames(DE_epi_stromal)%in%colnames(matall4)]
DE_epi_stromal <- DE_epi_stromal[,colnames(matall4)]

DE_epi_stromal <- as.data.frame(DE_epi_stromal)

ha = HeatmapAnnotation(
  log2FC_epi_stromal = unlist(unname(DE_epi_stromal[1,])), 
  padj = unlist(unname(DE_epi_stromal[2,])),
  col = list(log2FC_epi_stromal = f4,
             padj = f5)
  ,
  annotation_name_side = "left"
)
ht30 = Heatmap(matall4, col= fvir2, cluster_rows = F, cluster_columns = T, name = "Influence score",column_names_rot = 90, row_names_side = "left",top_annotation = ha)
htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('0.75'))

dev.off()

# Cluster again based on log2foldchange stromal
group = kmeans(t(matall4), centers = 10)$cluster

DE_epi_stromal <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/RNA_intra_cpm/20220421/stromal_epi_pseudobulkpadj.tsv", sep="\t",row.names=1)

f4 = colorRamp2(c(-5, 0, 5), c("purple", "#EEEEEE", "orange"), space = "RGB")

f5 = colorRamp2(c(1e-09, 0.05), c("red", "blue"), space = "RGB")


DE_epi_stromal <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/RNA_intra_cpm/20220421/epi_stromal_pseudobulkpadj.tsv", sep="\t",row.names=1)
DE_epi_stromal <- t(DE_epi_stromal)

DE_epi_stromal <- DE_epi_stromal[,colnames(DE_epi_stromal)%in%colnames(matall4)]
DE_epi_stromal <- DE_epi_stromal[,colnames(matall4)]

DE_epi_stromal <- as.data.frame(DE_epi_stromal)


ha = HeatmapAnnotation(
  log2FC_epi_stromal = unlist(unname(DE_epi_stromal[1,])), 
  padj = unlist(unname(DE_epi_stromal[2,])),
  col = list(log2FC_epi_stromal = f4,
             padj = f5)
  ,
  annotation_name_side = "left"
)

pdf(paste(resultsdir,'/', paste0('score_sel_0.75_epi_rot_0.5cut_DE_clusterwithin.pdf'),sep="/") ,width=17,height=3.5,paper='special')
Heatmap(matall4, col= fvir2, cluster_rows = F, cluster_columns = cluster_within_group(DE_epi_stromal, group),column_split = 10,
        name = "Influence score",column_names_rot = -90, row_names_side = "left", bottom_annotation =ha)
dev.off()

# add cell fates
ha_met <- rowAnnotation(foo = anno_block(gp = gpar(fill = "lightgrey"),
                                         labels = "Limbal-epithelial", labels_rot = 270,
                                         labels_gp = gpar(col = "black", fontsize = 12)))

pdf(paste(resultsdir,'/', paste0('score_sel_0.75_epi_rot_0.5cut_DE_clusterwithin.pdf'),sep="/") ,width=18.5,height=3.5,paper='special')
Heatmap(matall4, col= fvir2, cluster_rows = F, cluster_columns = cluster_within_group(DE_epi_stromal, group),column_split = 10,
        name = "Influence score",column_names_rot = -90, row_names_side = "left", bottom_annotation =ha,right_annotation = ha_met)
dev.off()

# Set TFs to subselect
corneal_TFs <-colnames(matall4)

### stromal >0.75
heatdfall3 <- heatdfall
heatdfall3 <- heatdfall[,6:7]
maxvec <- unname(apply(heatdfall3[,], 1, max))
heatdfall3$maxscore <- maxvec
heatdfall3= heatdfall3[heatdfall3$maxscore > 0.75,]
heatdfall3$maxscore <- NULL

matall3 <- as.matrix(heatdfall3)
# if(length(vec5) != 0) {
#vec5 <- rownames(matz[rownames(matz) %in% vec5 ,])
#vec6 <- which(rownames(matz) %in% vec5, arr.ind = T)
# plotting the list

pdf(paste(resultsdir,'/', paste0('score_sel_0.75_strom.pdf'),sep="/") ,width=3,height=20,paper='special')
ht30 = Heatmap(matall3, col= viridis(100), cluster_rows = T, cluster_columns = F, name = "Influence score")
htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('0.75'))

dev.off()

pdf(paste(resultsdir,'/', paste0('score_sel_0.75_strom_cols.pdf'),sep="/") ,width=3,height=20,paper='special')
ht30 = Heatmap(matall3, col= fvir, cluster_rows = T, cluster_columns = F, name = "Influence score")
htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('0.7'))

dev.off()


pdf(paste(resultsdir,'/', paste0('score_sel_0.75_str_rot.pdf'),sep="/") ,width=20,height=2.5,paper='special')
matall4 <- t(matall3)

DE_stromal_epi <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/RNA_intra_cpm/20220421/stromal_epi_pseudobulkpadj.tsv", sep="\t",row.names=1)

f4 = colorRamp2(c(-5, 0, 5), c("darkolivegreen", "#EEEEEE", "orange"), space = "RGB")

f5 = colorRamp2(c(1e-09, 0.05), c("red", "blue"), space = "RGB")

DE_stromal_epi <- t(DE_stromal_epi)

DE_stromal_epi <- DE_stromal_epi[,colnames(DE_stromal_epi)%in%colnames(matall4)]
DE_stromal_epi <- DE_stromal_epi[,colnames(matall4)]

DE_stromal_epi <- as.data.frame(DE_stromal_epi)

ha = HeatmapAnnotation(
  log2FC_stromal_epi = unlist(unname(DE_stromal_epi[1,])), 
  padj = unlist(unname(DE_stromal_epi[2,])),
  col = list(log2FC_stromal_epi = f4,
             padj = f5)
  ,
  annotation_name_side = "left"
)

ht30 = Heatmap(matall4, col= fvir, cluster_rows = F, cluster_columns = T, name = "Influence score",column_names_rot = 90, row_names_side = "left",top_annotation = ha)
htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('0.75'))

dev.off()


pdf(paste(resultsdir,'/', paste0('score_sel_0.75_str_rot_0.5cut.pdf'),sep="/") ,width=19,height=2,paper='special')
matall4 <- t(matall3)

DE_stromal_epi <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/RNA_intra_cpm/20220421/stromal_epi_pseudobulkpadj.tsv", sep="\t",row.names=1)

f4 = colorRamp2(c(-5, 0, 5), c("darkolivegreen", "#EEEEEE", "orange"), space = "RGB")

f5 = colorRamp2(c(1e-09, 0.05), c("red", "blue"), space = "RGB")

DE_stromal_epi <- t(DE_stromal_epi)

DE_stromal_epi <- DE_stromal_epi[,colnames(DE_stromal_epi)%in%colnames(matall4)]
DE_stromal_epi <- DE_stromal_epi[,colnames(matall4)]

DE_stromal_epi <- as.data.frame(DE_stromal_epi)

ha = HeatmapAnnotation(
  log2FC_stromal_epi = unlist(unname(DE_stromal_epi[1,])), 
  padj = unlist(unname(DE_stromal_epi[2,])),
  col = list(log2FC_stromal_epi = f4,
             padj = f5)
  ,
  annotation_name_side = "left"
)


ht30 = Heatmap(matall4, col= fvir2, cluster_rows = F, cluster_columns = T, name = "Influence score",column_names_rot = 90, row_names_side = "left")
htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('0.75'))

dev.off()

DE_stromal_epi <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/RNA_intra_cpm/20220421/stromal_epi_pseudobulkpadj.tsv", sep="\t",row.names=1)

f4 = colorRamp2(c(-5, 5), c("darkolivegreen", "#EEEEEE", "orange"), space = "RGB")

f5 = colorRamp2(c(1e-09, 0.05), c("red", "blue"), space = "RGB")

DE_stromal_epi <- t(DE_stromal_epi)

DE_stromal_epi <- DE_stromal_epi[,colnames(DE_stromal_epi)%in%colnames(matall4)]
DE_stromal_epi <- DE_stromal_epi[,colnames(matall4)]

DE_stromal_epi <- as.data.frame(DE_stromal_epi)

ha = HeatmapAnnotation(
  log2FC_stromal_epi = unlist(unname(DE_stromal_epi[1,])), 
  padj = unlist(unname(DE_stromal_epi[2,])),
  col = list(log2FC_stromal_epi = f4,
             padj = f5)
             ,
             annotation_name_side = "left"
)

pdf(paste(resultsdir,'/', paste0('score_sel_0.75_str_rot_0.5cut_DE.pdf'),sep="/") ,width=20,height=2.5,paper='special')
matall4 <- t(matall3)
ht30 = Heatmap(matall4, col= fvir2, cluster_rows = F, cluster_columns = T, name = "Influence score",column_names_rot = 90, row_names_side = "left", top_annotation = ha)
htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('0.70'))

dev.off()

# Cluster again based on log2foldchange stromal
group = kmeans(t(matall4), centers = 3)$cluster

DE_epi_stromal <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/RNA_intra_cpm/20220421/stromal_epi_pseudobulkpadj.tsv", sep="\t",row.names=1)

f4 = colorRamp2(c(-5, 0, 5), c("purple", "#EEEEEE", "orange"), space = "RGB")

f5 = colorRamp2(c(1e-09, 0.05), c("red", "blue"), space = "RGB")


DE_epi_stromal <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/RNA_intra_cpm/20220421/epi_stromal_pseudobulkpadj.tsv", sep="\t",row.names=1)
DE_epi_stromal <- t(DE_epi_stromal)

DE_epi_stromal <- DE_epi_stromal[,colnames(DE_epi_stromal)%in%colnames(matall4)]
DE_epi_stromal <- DE_epi_stromal[,colnames(matall4)]

DE_epi_stromal <- as.data.frame(DE_epi_stromal)


ha = HeatmapAnnotation(
  log2FC_epi_stromal = unlist(unname(DE_epi_stromal[1,])), 
  padj = unlist(unname(DE_epi_stromal[2,])),
  col = list(log2FC_epi_stromal = f4,
             padj = f5)
  ,
  annotation_name_side = "left"
)

pdf(paste(resultsdir,'/', paste0('score_sel_0.75_str_rot_0.5cut_DE_clusterwithin.pdf'),sep="/") ,width=16,height=2.5,paper='special')
Heatmap(matall4, col= fvir2, cluster_rows = F, cluster_columns = cluster_within_group(DE_epi_stromal, group),column_split = 4,
        name = "Influence score",column_names_rot = -90, row_names_side = "left", bottom_annotation =ha)
dev.off()

# add cell fates
ha_met <- rowAnnotation(foo = anno_block(gp = gpar(fill = "lightgrey"),
                                         labels = "Stromal", labels_rot = 270,
                                         labels_gp = gpar(col = "black", fontsize = 12)))


pdf(paste(resultsdir,'/', paste0('score_sel_0.75_str_rot_0.5cut_DE_clusterwithin.pdf'),sep="/") ,width=16,height=2.7,paper='special')
Heatmap(matall4, col= fvir2, cluster_rows = F, cluster_columns = cluster_within_group(DE_epi_stromal, group),column_split = 4,
        name = "Influence score",column_names_rot = -90, row_names_side = "left", bottom_annotation =ha,right_annotation = ha_met)
dev.off()

# pdf(paste(resultsdir,'/', paste0('score_sel_0.75_rot.pdf'),sep="/") ,width=15,height=5,paper='special')
# matall4 <- t(matall3)
# ht30 = Heatmap(matall4, col= viridis(100), cluster_rows = F, cluster_columns = T, name = "Influence score")
# htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
# draw(htlist2, column_title = paste0(i,'0.75'))
# 
# dev.off()

# generate intersection plot of factors with a score >0.70
heatdfall3 <- heatdfall

merged.data.frame = Reduce(function(...) merge(..., all=T), fullheat)
merged.data.frame[merged.data.frame == 0] <- NA

maxvec <- unname(apply(heatdfall3[,1:5], 1, max))
heatdfall3$maxscore <- maxvec
heatdfall3= heatdfall3[heatdfall3$maxscore > 0.75,1:5]

merged.data.frame <- merged.data.frame[merged.data.frame$factor%in%rownames(heatdfall3),]
merged.data.frame2 <- merged.data.frame
merged.data.frame2[!is.na(merged.data.frame2)] <- 1
merged.data.frame2[is.na(merged.data.frame2)] <- 0
merged.data.frame2$factor <- merged.data.frame$factor
merged.data.frame2 <-merged.data.frame2[rowSums(merged.data.frame2[, -1])>0, ]
merged.data.frame2 <-merged.data.frame2[,1:6]

for(i in 2:ncol(merged.data.frame2)){ merged.data.frame2[ , i] <- as.integer(merged.data.frame2[ , i]) }
colnoms <- colnames(merged.data.frame2)[colnames(merged.data.frame2) != "factor"]

fill = viridis(length(files$cell_id[1:5]),alpha = 0.8,option = "H")

pdf(paste(resultsdir,'/intersection_0.7_epi.pdf',sep="/") ,width=4,height=5,paper='special')
upset(merged.data.frame2,nintersects = NA,
      sets = rev(colnoms),keep.order=T,
      order.by="freq", matrix.color="black", point.size=1,
      sets.bar.color=fill)
dev.off()

heatdfall3 <- heatdfall

merged.data.frame = Reduce(function(...) merge(..., all=T), fullheat)
merged.data.frame[merged.data.frame == 0] <- NA

maxvec <- unname(apply(heatdfall3[,6:7], 1, max))
heatdfall3$maxscore <- maxvec
heatdfall3= heatdfall3[heatdfall3$maxscore > 0.75,]
merged.data.frame <- merged.data.frame[merged.data.frame$factor%in%rownames(heatdfall3),]
merged.data.frame2 <- merged.data.frame
merged.data.frame2[!is.na(merged.data.frame2)] <- 1
merged.data.frame2[is.na(merged.data.frame2)] <- 0
merged.data.frame2$factor <- merged.data.frame$factor
merged.data.frame2 <-merged.data.frame2[rowSums(merged.data.frame2[, -1])>0, ]
merged.data.frame2 <-merged.data.frame2[,c(1,6:7)]

for(i in 2:ncol(merged.data.frame2)){ merged.data.frame2[ , i] <- as.integer(merged.data.frame2[ , i]) }
colnoms <- colnames(merged.data.frame2)[colnames(merged.data.frame2) != "factor"]

fill = viridis(length(files$cell_id[6:7]),alpha = 0.8,option = "H")

pdf(paste(resultsdir,'/intersection_0.75_str.pdf',sep="/") ,width=2,height=5,paper='special')
upset(merged.data.frame2,nintersects = NA,
      sets = colnoms,keep.order=T,
      order.by="freq", matrix.color="black", point.size=1,
      sets.bar.color=fill)
dev.off()

# Set TFs to subselect
stromal_TFs <-colnames(matall4)

# TF-TF interaction networks
diffnetwork_full <- NULL
for (i in unique(files$cell_id)){
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  diffnetwork_full <- rbind(diffnetwork_full,diffnetwork)
}

names(diffnetwork_full) <- diffnetwork_full[1,]
diffnetwork_full <- diffnetwork_full[-1,]

diff_epi <- diffnetwork_full[diffnetwork_full$source%in%corneal_TFs & diffnetwork_full$target%in%corneal_TFs,]

epi_net <-read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE_22042022/outs_v1/epi/narrow/full_network_includeprom.txt")
names(epi_net) <- epi_net[1,]
epi_net <- epi_net[-1,]

epi_net2 <- epi_net[as.double(epi_net$activity)>0.8,]
head(as.double(epi_net$activity))
i <- "LSC"

library(igraph)
for (i in unique(files$cell_id)){
diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
factors <- read.table(files$network[files$cell_id==i],header = T)
factors <- factors[order(factors$influence_score,decreasing = T),]
# selecting the top 40 factors
#top <- factors$factor[1:25]
top <- corneal_TFs
top <- top[top %in% factors$factor[factors$influence_score>0.8]]
# fullheat2 gives the nodes stats
fullheat3 <- list()
fullheat4 <- list()

factorsheat3 <- read.table(files$diffnetwork[files$cell_id==i],header = T)
factorsheat4 <- factorsheat3[factorsheat3$source %in% top,]
factorsheat4 <- factorsheat4[factorsheat4$target %in% top,]

extracted2 <- factorsheat4[,c("source","target","weight")]
colnames(extracted2) <- c("from","to","weight")
edges_set <- extracted2[,c("from","to","weight")]

edges_set<-edges_set[!duplicated(edges_set), ]

# Count the number of degree for each node:

dirtar <- fullheat2[[i]]$direct_targets
noms <- fullheat2[[i]]$factor
degree <- setNames(as.numeric(dirtar), noms)

g = graph_from_data_frame(d=edges_set,directed=TRUE)#, vertices=nodes
deg <- degree(g, mode="all")


uni_all <- seq(min(degree), max(degree))

colors <- data.frame(color = heat.colors(length(uni_all), rev = T),
                      levels = uni_all)
degree2 <- degree[names(degree) %in% names(deg)]
degree2 <- degree2[names(deg)]

V(g)$degree <- degree2
# Use match to get the index of right color, matching on levels
V(g)$color <- colors$color[match(V(g)$degree, colors$levels)]

# use degree as labels
V(g)$label <- names(degree2)

#saveRDS(g,paste(resultsdir,'/',i, 'network.rds',sep="/"))

#g<-readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/ANANSE_influence_analysis_meta20220315/network_LSC.rds")
pdf(paste(resultsdir,"/",i,'full_network.pdf',sep=""),width=11,height=12,paper='special')
print(plot(g,vertex.size=log(deg+2)*3,layout=layout_with_graphopt,#,# layout = layout_with_graphopt, 
           edge.width=-1/log(edge.attributes(g)[["weight"]])*5,edge.arrow.size = 0.2))
k <- as_data_frame(g, what = "vertices")
k <- k[order(k$degree,decreasing = F),]
print(legend("bottomright",title = "Number of targets", legend=levels(as.factor(k$degree)), col =k$color , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, horiz = FALSE, inset = c(0.1, 0.1)))
dev.off()

#export network for nicer visualization in cytoscape
write_graph(graph = g,file = paste(resultsdir,"/",i,'full_network.graphml',sep=""),
            format = "graphml"
)

factorsheat4 <- factorsheat4[factorsheat4$weight>0.2,]
extracted2 <- factorsheat4[,c("source","target","weight")]
colnames(extracted2) <- c("from","to","weight")
edges_set <- extracted2[,c("from","to","weight")]

edges_set<-edges_set[!duplicated(edges_set), ]

# Count the number of degree for each node:

dirtar <- fullheat2[[i]]$direct_targets
noms <- fullheat2[[i]]$factor
degree <- setNames(as.numeric(dirtar), noms)

g = graph_from_data_frame(d=edges_set,directed=TRUE)#, vertices=nodes
deg <- degree(g, mode="all")


uni_all <- seq(min(degree), max(degree))

colors <- data.frame(color = heat.colors(length(uni_all), rev = T),
                     levels = uni_all)
degree2 <- degree[names(degree) %in% names(deg)]
degree2 <- degree2[names(deg)]

V(g)$degree <- degree2
# Use match to get the index of right color, matching on levels
V(g)$color <- colors$color[match(V(g)$degree, colors$levels)]

# use degree as labels
V(g)$label <- names(degree2)

#saveRDS(g,paste(resultsdir,'/',i, 'network.rds',sep="/"))

#g<-readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/ANANSE_influence_analysis_meta20220315/network_LSC.rds")
pdf(paste(resultsdir,"/",i,'full_network_0.2.pdf',sep=""),width=11,height=12,paper='special')
print(plot(g,vertex.size=log(deg+2)*3,layout=layout_with_graphopt,#,# layout = layout_with_graphopt, 
           edge.width=-1/log(edge.attributes(g)[["weight"]])*5,edge.arrow.size = 0.2))
k <- as_data_frame(g, what = "vertices")
k <- k[order(k$degree,decreasing = F),]
print(legend("bottomright",title = "Number of targets", legend=levels(as.factor(k$degree)), col =k$color , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, horiz = FALSE, inset = c(0.1, 0.1)))
dev.off()

#export network for nicer visualization in cytoscape
write_graph(graph = g,file = paste(resultsdir,"/",i,'full_network_0.2.graphml',sep=""),
            format = "graphml"
)
}

#################### scores of top 25 genes with the highest rowvar

# Generate the interaction dataframe
fullheat3 <- list()
fullheat4 <- list()
for (i in unique(files$cell_id)){
  factorsheat1 <- read.table(files$network[files$cell_id==i],header = T)
  factorsheat2 <- factorsheat1[,c("factor","factor_fc","influence_score")]
  colnames(factorsheat2) <- c("factor","factor_fc",i)
  fullheat3[[i]] <- factorsheat2
  factorsheat2 <- factorsheat1[,c("factor","factor_fc")]
  colnames(factorsheat2) <- c("factor",i)
  fullheat4[[i]] <-factorsheat2
}

merged.data.frame3 = Reduce(function(...) merge(..., all=T), fullheat4)

rownames(merged.data.frame3) <- merged.data.frame3$factor

merged.data.frame3[is.na(merged.data.frame3)] <- 0

merged.data.frame3$factor <- NULL
library(dplyr)
library(matrixStats)

var_fact <- merged.data.frame3 %>% replace(is.na(.), 0) %>% mutate(row_wise_var = rowVars(as.matrix(merged.data.frame3[,])))

var_fact$factor <- rownames(var_fact)
# Generating a dataframe from influence scores
# heatdfall <- merged.data.frame3
# rownames(heatdfall) <- heatdfall$factor
# heatdfall$factor <- NULL
# heatdfall[is.na(heatdfall)] <- 0

pdf(paste(resultsdir,'/', paste0('complexheatmap_25_var.pdf'),sep="/") ,width=5,height=12,paper='special')
top30 <- c()
df_merge <- NULL
for (i in unique(files$cell_id)){
  print(i)
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  factors <- read.table(files$network[files$cell_id==i],header = T)
  df_merge <- merge(var_fact,factors,by="factor")
  df_merge$varscaled <- df_merge$influence_score*abs(log10(df_merge$row_wise_var))
  df_merge <- factors[order(df_merge$varscaled,decreasing = T),]
  # selecting the top 30 factors
  df_merge <- df_merge[df_merge$influence_score>0.3,]
  top <- df_merge$factor[1:30]
  heatdfall2 <- heatdfall[rownames(heatdfall) %in% top,]
  matall2 <- as.matrix(heatdfall2)
  # if(length(vec5) != 0) {
  #vec5 <- rownames(matz[rownames(matz) %in% vec5 ,])
  #vec6 <- which(rownames(matz) %in% vec5, arr.ind = T)
  # plotting the list
  
  ht30 = Heatmap(matall2, col= viridis(100), cluster_rows = T, cluster_columns = F, name = "Influence score")
  htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
  draw(htlist2, column_title = paste0(i,' top 25var'))
  top30 <- c(top30,top)
}
dev.off()

# joined top 50 heatmap corneal epithelium
# Generate a heatmap of the top 20 found factors for each cell population
pdf(paste(resultsdir,'/', paste0('top25var_joined.pdf'),sep="/") ,width=5,height=16,paper='special')
top30 <- unique(top30)
ht30 = Heatmap(matall, col= viridis(100), cluster_rows = F, cluster_columns = F, name = "Influence score")
vec5 <- rownames(matall[rownames(matall) %in% top30 ,])
vec6 <- which(rownames(matall) %in% vec5, arr.ind = T)

htlist2 <- ht30 + rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('corneal top 25'))

matall20 <- matall
matall20 <- matall20[rownames(matall20) %in% top30,]

ht30 = Heatmap(matall20, col= viridis(100), cluster_rows = T, cluster_columns = F, name = "Influence score")
#vec5 <- rownames(matall[rownames(matall) %in% top30 ,])
#vec6 <- which(rownames(matall) %in% vec5, arr.ind = T)

htlist2 <- ht30 #+ rowAnnotation(link = anno_mark(at =  vec6,labels = vec5))#,labels_gp = gpar(col = fontcolors)
draw(htlist2, column_title = paste0('corneal 25var'))

dev.off()






###
# delete below?


# Order the rows and columns based on Z-score gene expression table
col_mat = viridis(3,alpha = 0.8)
f3 = colorRamp2(c(-3, 0, 3), c("blue", "#EEEEEE","red"), space = "RGB")
ht3 = Heatmap(matz, col= f3, cluster_rows = F, cluster_columns = F, name = "Z-score RNA count")

vec3 <- sharedall[sharedall %in% rownames(heatdfall)]
vec3 <- rownames(matz[rownames(matz) %in% vec3 ,])

vec4 <- unname(uniqtfs[uniqtfs %in% rownames(heatdfall)])
vec4 <- rownames(matz[rownames(matz) %in% vec4 ,])

vec1 <- which(rownames(matz) %in% vec3, arr.ind = T)
vec2 <- which(rownames(matz) %in% vec4, arr.ind = T)

# Setting the colors of the general transcription factors expected to be found
#neurexpected <- unlist(strsplit(exgeneral[exgeneral$cell_type == "neural",]$expected_tfs,","))
eyeexpected <- unlist(strsplit(exgeneral[exgeneral$cell_type == "eye",]$expected_tfs,","))
epiexpected <- unlist(strsplit(exgeneral[exgeneral$cell_type == "epidermal",]$expected_tfs,","))

if (!length(vec1) == 0) {
  # linking the celltype specific factors
  fontcolors <- rep('black', length(vec3))
  #row_idx <- which(vec3 %in% neurexpected)
  #fontcolors[row_idx] <- 'purple'
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% eyeexpected)
  fontcolors[row_idx] <- 'red'
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% epiexpected)
  fontcolors[row_idx] <- 'orange'
  
  # linking the expected factors as a blue color
  expected2 <- unlist(strsplit(expected[expected$cell_type == i,]$expected_tfs,","))
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% expected2)
  fontcolors[row_idx] <- 'blue'
  htlist1 <- ht4 + ht3 + rowAnnotation(link = anno_mark(at =  vec1,labels = vec3,labels_gp = gpar(col = fontcolors)))
}

# linking the expected factors as a blue color for the single cell populations
expected2 <- unlist(strsplit(expected$expected_tfs,","))

# linking the celltype specific factors
fontcolors <- rep('black', length(vec4))
#row_idx <- which(vec4 %in% neurexpected)
#fontcolors[row_idx] <- 'purple'

# linking the celltype specific factors
row_idx <- which(vec4 %in% eyeexpected)
fontcolors[row_idx] <- 'red'

# linking the celltype specific factors
row_idx <- which(vec4 %in% epiexpected)
fontcolors[row_idx] <- 'orange'

# linking the celltype specific factors
row_idx <- which(vec4 %in% expected2)
fontcolors[row_idx] <- 'blue'

htlist2 <- ht4 + ht3 + rowAnnotation(link = anno_mark(at =  vec2,labels = vec4,labels_gp = gpar(col = fontcolors)))

pdf(paste(resultsdir,'/complexheatmap_rna.pdf',sep="/") ,width=10,height=10,paper='special')
if (!length(vec1) == 0) {
  draw(htlist1, column_title = "shared_factors")
}
draw(htlist2, column_title = "unique_factors")
dev.off()

# Top 25 selection and intersection

for (i in unique(files$cell_id)){
  print(i)
  #diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  factors <- read.table(files$network[files$cell_id==i],header = T)
  factors <- factors[order(factors$influence_score,decreasing = T),]
  # selecting the top 25 factors
  top <- factors$factor[1:30]
  vec5 <- top
  vec5 <- vec5[vec5 %in% rownames(heatdfall)]
  print(vec5)
  vec25 <- c(vec5,vec5)
}

# Subselecting the matrix
matall25 <- matall[rownames(matall) %in% unique(vec25),]
ht7 = Heatmap(matall25, col= viridis(100), cluster_rows = T, cluster_columns = T, name = "Influence score")

# Order the rows and columns based on Z-score gene expression table
matz25 <- matz[rownames(matz) %in% unique(vec25),]

roword <- row_order(ht7)
matall25 <- matall25[roword,]
colord <- column_order(ht7)
matall25 <- matall25[, colord]
matz25 <- matz25[rownames(matall25),colnames(matall25)]

###
col_mat = viridis(3,alpha = 0.8)
f3 = colorRamp2(c(-3, 0, 3), c("blue", "#EEEEEE","red"), space = "RGB")
ht8 = Heatmap(matz25, col= f3, cluster_rows = F, cluster_columns = F, name = "Z-score RNA count")

vec6 <- vec25
vec5 <- which(vec6 %in% rownames(matz25) , arr.ind = T)

# linking the celltype specific factors
fontcolors <- rep('black', length(vec6))
#row_idx <- which(vec6 %in% neurexpected)
#fontcolors[row_idx] <- 'purple'

# linking the celltype specific factors
row_idx <- which(vec6 %in% eyeexpected)
fontcolors[row_idx] <- 'red'

# linking the celltype specific factors
row_idx <- which(vec6%in% epiexpected)
fontcolors[row_idx] <- 'orange'
htlist3 <- ht7 + ht8 + rowAnnotation(link = anno_mark(at =  vec5,labels = vec6,labels_gp = gpar(col = fontcolors)))

pdf(paste(resultsdir,'/complexheatmap_rna25.pdf',sep="/") ,width=10,height=5,paper='special')
draw(htlist3, column_title = "shared_factors_top25")
dev.off()

# Single the cell populations for a clear overview
pdf(paste(resultsdir,'/complexheatmap_rna_singled.pdf',sep="/") ,width=10,height=10,paper='special')
for (i in unique(files$cell_id)) {
  vec5 <- uniq[[i]]
  vec5 <- vec5[vec5 %in% rownames(heatdfall)]
  if(length(vec5) == 0) {
    next
  }
  vec5 <- rownames(matz[rownames(matz) %in% vec5 ,])
  vec6 <- which(rownames(matz) %in% vec5, arr.ind = T)
  
  # linking the expected factors as a blue color
  expected2 <- unlist(strsplit(expected[expected$cell_type == i,]$expected_tfs,","))
  
  # linking the celltype specific factors
  fontcolors <- rep('black', length(vec6))
  #row_idx <- which(vec6 %in% neurexpected)
  #fontcolors[row_idx] <- 'purple'
  
  # linking the celltype specific factors
  row_idx <- which(vec6 %in% eyeexpected)
  fontcolors[row_idx] <- 'red'
  
  # linking the celltype specific factors
  row_idx <- which(vec6 %in% epiexpected)
  fontcolors[row_idx] <- 'orange'
  
  # linking the celltype specific factors
  row_idx <- which(vec6 %in% expected2)
  fontcolors[row_idx] <- 'blue'
  
  htlist2 <- ht4 + ht3 + rowAnnotation(link = anno_mark(at =  vec6,labels = vec5,labels_gp = gpar(col = fontcolors)))
  draw(htlist2, column_title = i)
}
dev.off()





##############################
# Extracting unique factors shared across populations of interest
list_of_cells = c("LPCs","CSB","CjS")

list_of_vecs <- list()
fullheatfacts2 <- fullheatfacts
for (v in list_of_cells){
  list_of_vecs[[v]]<-fullheatfacts[[v]]
  fullheatfacts2[[v]] <- NULL
}

fac_interest <- Reduce(intersect, list_of_vecs)
`%notin%` <- Negate(`%in%`)

vecname <- paste("int_",paste(list_of_cells, collapse = ''),sep="")
assign(vecname,fac_interest[fac_interest %notin% unlist(fullheatfacts2)])
facs <- fac_interest[fac_interest %notin% unlist(fullheatfacts2)]

# Generate a heatmap of interesting factors
pdf(paste(resultsdir,'/', paste0(vecname,'complexheatmap_interesting.pdf'),sep="/") ,width=10,height=10,paper='special')
vec5 <- facs
vec5 <- vec5[vec5 %in% rownames(heatdfall)]
if(length(vec5) != 0) {
  vec5 <- rownames(matz[rownames(matz) %in% vec5 ,])
  vec6 <- which(rownames(matz) %in% vec5, arr.ind = T)
  
  # linking the expected factors as a blue color
  expected2 <- unlist(strsplit(expected$expected_tfs,","))
  
  # linking the celltype specific factors
  fontcolors <- rep('black', length(vec5))
  #row_idx <- which(vec5 %in% neurexpected)
  #fontcolors[row_idx] <- 'purple'
  
  # linking the celltype specific factors
  row_idx <- which(vec5 %in% eyeexpected)
  fontcolors[row_idx] <- 'red'
  
  # linking the celltype specific factors
  row_idx <- which(vec5 %in% epiexpected)
  fontcolors[row_idx] <- 'orange'
  
  # linking the celltype specific factors
  row_idx <- which(vec5 %in% expected2)
  fontcolors[row_idx] <- 'blue'
  htlist2 <- ht4 + ht3 + rowAnnotation(link = anno_mark(at =  vec6,labels = vec5,labels_gp = gpar(col = fontcolors)))
  draw(htlist2, column_title = vecname)
}
dev.off()

# Subselecting interesting cells
colsub <- c("LPCs","LNPCs","CjS","CSB","CB","CSSCs","StC")
list_of_cells1 = c("LPCs","LNPCs","CjS","CSB","CB")
list_of_cells4 = c("LPCs","CjS","CSB","CB")
list_of_cells5 = c("LNPCs","CjS","CSB","CB")
list_of_cells6 = c("LNPCs","CB")
list_of_cells3 = c("LPCs","CjS")
list_of_cells9 = c("LPCs","CB")
list_of_cells10 = c("LNPCs","CjS")
list_of_cells2 = c("CjS","CSB","CB")
list_of_cells7 = c("LNPCs","LPCs")
list_of_cells8 = c("LNPCs","LPCs","CB")
list_of_cells11 = c("LNPCs","LPCs","CjS")
list_of_cells12 = c("CSSCs","StC")
list_of_cells13 = c("LPCs","CSB","CjS")
list_of_cells14 = c("LNPCs","CB","CSB")

#list_of_cells1 = c("LiCo")
#list_of_cells2 = c("StCSC")
#list_of_cells3 = c("CSSCs","StC")

veclst <- list(list_of_cells1,list_of_cells2,list_of_cells3,list_of_cells4,list_of_cells5,list_of_cells6,
               list_of_cells7,list_of_cells8,list_of_cells9,list_of_cells10,list_of_cells11,list_of_cells12,list_of_cells13,list_of_cells14)#,

facs2 <- NULL

for (i in veclst){
  print(i)
  fullheatfacts2 <- fullheatfacts
  list_of_vecs <- list()
  for (v in i){
    print(v)
    list_of_vecs[[v]]<-fullheatfacts[[v]]
    fullheatfacts2[[v]] <- NULL
  }
  
  fac_interest <- Reduce(intersect, list_of_vecs)
  `%notin%` <- Negate(`%in%`)
  
  vecname <- paste("int_",paste(i, collapse = ''),sep="")
  assign(vecname,fac_interest[fac_interest %notin% unlist(fullheatfacts2)])
  facs <- fac_interest[fac_interest %notin% unlist(fullheatfacts2)]
  facs2 <- c(facs2,facs)
}

# Add all expected factors in the subselection
#expected2 <- unlist(strsplit(expected$expected_tfs,","))
#facs2 <- c(facs2,epiexpected,eyeexpected,neurexpected,expected2)

# Selecting a cutoff for the influence values
heatdfall2 <- heatdfall[rownames(heatdfall) %in% (facs2),]

# subselecting columns based on interesting cells
heatdfall2 <- heatdfall2[, colnames(heatdfall2) %in% colsub]


# Selecting tfs based on the cutoff in the Z-score table
#z <- subset(lakorna, rownames(lakorna) %in% (merged.data.frame$factor))
#zz <- z[rownames(z) %in% rownames(heatdfall2),]
zz <- z[rownames(z) %in% rownames(heatdfall2),]

# subselecting columns based on interesting cells
zz <- zz[, colnames(zz) %in% colsub]

# Generate the matrix 
#matall2 <- as.matrix(heatdfall2)
matall2 <- as.matrix(heatdfall2)
matz2 <- as.matrix(zz)

col.order <- c("LPCs","LNPCs","CB","CSB","CjS","CSSCs","StC")
matall2 <- matall2[,col.order]

ht5 = Heatmap(matall2, col= viridis(100), cluster_rows = T, cluster_columns = F, name = "Influence score")

# Order the rows and columns based on Z-score gene expression table
roword <- row_order(ht5)
matall2 <- matall2[roword,]
colord <- column_order(ht5)
matall2 <- matall2[, colord]
matz2 <- matz2[rownames(matall2),colnames(matall2)]

col_mat = viridis(3,alpha = 0.8)
f3 = colorRamp2(c(-3, 0, 3), c("blue", "#EEEEEE","red"), space = "RGB")
ht6 = Heatmap(matz2, col= f3, cluster_rows = F, cluster_columns = F, name = "Z-score RNA count")

vec3 <- rownames(matz2)
vec1 <- which(rownames(matz2) %in% vec3, arr.ind = T)

#pdf(paste(resultsdir,'/complexheatmap_rna_subselected.pdf',sep="/") ,width=10,height=13.5,paper='special')
pdf(paste(resultsdir,'/complexheatmap_rna_subselected3.pdf',sep="/") ,width=10,height=10,paper='special')
if (!length(vec1) == 0) {
  # linking the celltype specific factors
  fontcolors <- rep('black', length(vec3))
  #row_idx <- which(vec3 %in% neurexpected)
  #fontcolors[row_idx] <- 'purple'
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% eyeexpected)
  fontcolors[row_idx] <- 'red'
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% epiexpected)
  fontcolors[row_idx] <- 'orange'
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% expected2)
  fontcolors[row_idx] <- 'red'
  htlist1 <- ht5 + ht6 + rowAnnotation(link = anno_mark(at =  vec1,labels = vec3,labels_gp = gpar(col = fontcolors)))
}
draw(htlist1, column_title = "subselected")#, row_km = 2, cluster_rows = TRUE)
dev.off()

# Line plots of interesting factors together (WIP)
datalineplot <- NULL
# classifications of time series data

datalineplot2 <- NULL
datalineplot <-as.data.frame(as.table(matall2))
matall4 <- matall2
for (res in 4:20){
dtw_cluster = tsclust(heatdfall2, type="partitional",k=res,
                      distance="dtw_basic",centroid = "pam",seed=1234,trace=T,
                      args = tsclust_args(dist = list(window.size = 5)))


dfall3 <- NULL
matall3 <- as.data.frame(unname(dtw_cluster@cluster),rownames(matall2))

datalineplot[[paste0("dtw", res)]]<- rep(matall3[rownames(matall3)%in%datalineplot$Var1,], length(unique(datalineplot$Var2)))


}


pdf(paste(resultsdir,'/grouped_tfs_pseudotime_clustree.pdf',sep="/") ,width=10,height=20,paper='special')
datalineplot2 <- datalineplot[1:46,4:16]

print(clustree(datalineplot2,prefix = "dtw"))
print(clustree(
  datalineplot2,
  prefix = "dtw",
  #exprs = c("data", "counts", "scale.data"),
  assay = NULL,
  node_colour = "sc3_stability"
))
dev.off()

########### corresponding expression in scRNA-seq of the factors (violin plots)
# function for stacked plots was based on https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
library(Seurat)
library(patchwork)
library(ggplot2)

unique(datalineplot$Var1)
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.05, 0, -0.05, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.05, 0, -0.05, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

seur_obj <- readRDS('/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated.rds')
seur_obj2 <-subset(x=seur_obj, subset = (costum_clustering == "LPCs"|costum_clustering == "LNPCs"|costum_clustering == "CB"|costum_clustering == "CSB"|costum_clustering == "CjS"))
seur_obj2$costum_clustering
seur_obj2 <- droplevels(x = seur_obj2)
seur_obj2$costum_clustering <- droplevels(x = seur_obj2$costum_clustering)

features<- unique(datalineplot$Var1)

seur_obj2@meta.data$seurat_clusters <- seur_obj2$costum_clustering
seur_obj2@active.ident <- seur_obj2$costum_clustering

new_order <- c("LPCs","LNPCs","CB","CSB","CjS")
seur_obj2@active.ident <- factor(seur_obj2@active.ident, levels = new_order)

# Selected 13 clusters from clustree
pdf(paste(resultsdir,'/grouped_tfs_pseudotime.pdf',sep="/") ,width=6,height=8,paper='special')
for (i in 1:13){
datalineplot5 <- datalineplot[datalineplot$dtw13==i,]
print(ggplot(datalineplot5, aes(x=Var2, y=Freq, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point() +
  labs(color = "Transcription factors",y = "Influence score",x= "Ident"))
  print(StackedVlnPlot(obj = seur_obj2, features = unique(datalineplot5$Var1)))
}
dev.off()

#### interesting subselection according to clustering
int_nfkb <- c("NFKB1","RELA","ESRRA","HES4")

pdf(paste(resultsdir,'/grouped_tfs_NFKB.pdf',sep="/") ,width=6,height=8,paper='special')
datalineplot5 <- datalineplot[datalineplot$Var1 %in%int_nfkb,]
print(ggplot(datalineplot5, aes(x=Var2, y=Freq, group=Var1)) +
          geom_line(aes(color=Var1))+
          geom_point() +
          labs(color = "Transcription factors",y = "Influence score",x= "Ident"))
print(StackedVlnPlot(obj = seur_obj2, features = unique(datalineplot5$Var1)))
dev.off()


int_cornea <- c("ZBTB33","ELK4","PAX6","USF1","DPB","KLF5","FOSL2","EHF","ETV2","GRHL1","ELF3","GRHL2")

pdf(paste(resultsdir,'/grouped_tfs_cornea.pdf',sep="/") ,width=6,height=5,paper='special') 
datalineplot5 <- datalineplot[(datalineplot$Var1 %in%int_cornea)& (datalineplot$Var2 %notin% c("CSSCs","StC")),]
print(ggplot(datalineplot5, aes(x=Var2, y=Freq, group=Var1)) +
        geom_line(aes(color=Var1))+
        geom_point() +
        labs(color = "Transcription factors",y = "Influence score",x= "Ident"))
dev.off()

pdf(paste(resultsdir,'/grouped_tfs_cornea_expr.pdf',sep="/") ,width=6,height=12,paper='special') 
print(StackedVlnPlot(obj = seur_obj2, features = unique(datalineplot5$Var1)))
dev.off()

fontcolors <- rep('black', length(datalineplot$Var1))
#row_idx <- which(datalineplot$Var1 %in% neurexpected)
#fontcolors[row_idx] <- 'purple'

# linking the celltype specific factors
row_idx <- which(datalineplot$Var1 %in% eyeexpected)
fontcolors[row_idx] <- 'red'

# linking the celltype specific factors
row_idx <- which(datalineplot$Var1 %in% epiexpected)
fontcolors[row_idx] <- 'orange'

# linking the celltype specific factors
row_idx <- which(datalineplot$Var1 %in% expected2)
fontcolors[row_idx] <- 'blue'

datalineplot$type <- fontcolors

datalineplot1 <- datalineplot[datalineplot$type =="red",]
ggplot(datalineplot1, aes(x=Var2, y=Freq, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point() +
  labs(color = "Transcription factors")

datalineplot2 <- datalineplot[datalineplot$type =="blue",]
ggplot(datalineplot2, aes(x=Var2, y=Freq, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point() +
  labs(color = "Transcription factors")

datalineplot3 <- datalineplot[datalineplot$type =="black",]
ggplot(datalineplot3, aes(x=Var2, y=Freq, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point() +
  labs(color = "Transcription factors")

datalineplot4 <- datalineplot[datalineplot$type =="orange",]
ggplot(datalineplot4, aes(x=Var2, y=Freq, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point() +
  labs(color = "Transcription factors")




StackedVlnPlot(obj = seur_obj2, features = features)


#####################################
#rating chord diagrams of the top 10 factors
# note: perhaps based on the shared and unique factors?

########################## WIP
# subselecting factors based on the plot (by hand first)


circos.par(gap.degree=0.5)

myFun <- function(vector) {
  ind <- ave(rep(1, length(vector)), vector, FUN = length)
  print(ind)
  thresh <- max(ind)
  vector[ind > thresh - 1] ## added "+1" to match your terminology
}

top <- c("ZBTB33","ELK4","PAX6","FOXQ1","TFAP2C","OTX1","IRF6","E2F4","ELK1","TGIF1","RFX5","TP63","ASCL2")
length(top)
cells <- c("LPCs","LNPCs","CB","CjS","CSB")

# add a list in the code for the upset R plot
fullheat <- list()

pdf(paste(resultsdir,"/",cells,'chordsselected.pdf',sep=""),width=6,height=6,paper='special') 
for (i in cells){
  print(i)
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  #factors <- read.table(files$network[files$cell_id==i],header = T)
  #factors <- factors[order(factors$sumScaled,decreasing = T),]
  
  # selecting the top 20 factors
  #top <- factors$factor[1:20]
  
  colnames(diffnetwork) <- c("from","to","value")
  df <- data.frame(diffnetwork)
  
  # selecting
  df3 <- df[df$from %in% top,]
  df3 <- df3[order(df3$value,decreasing = T) & (df3$value >0.8),]
  vec <- df3$to
  vec2 <- unique(myFun(vec))
  df3 <- df3[(df3$to %in% vec2),]
  vec3 <- df3$from
  
  #ind <- ave(rep(1, length(vec)), vec, FUN = length)
  #thresh <- max(ind)
  
  vec4 <- unique(myFun(vec3))
  df3 <- df3[(df3$from %in% vec4),]
  #vec <- df3$to
  #vec2 <- unique(myFun(vec,length(top)))
  # for (z in unique(top)){
  #   print(z)
  #   df4 <- df3
  #   df5 <- df4[df4$from == z,]
  #   df6 <- df5[df5$to == head(df5$to,10),]
  #   df7 = NULL
    
    # Determine if one of the top 10 factors also influences the factor of interest
    # for (q in unique(head(df5$to,10))){
    #   print(q)
    #   df7 = rbind(df7, df4[df4$from == q & df4$to == z,])
    #   df7 = rbind(df7, df4[df4$from == q & df4$to %in% head(df5$to,10),])
    #   
    #   # If a factor is found influencing the factor of interest or another factor in the top 10 then add to original dataframe
    #   if (!is.null(df7)){
    #     df6 = rbind(df6, df7)#[df4$from == q & df4$to %in% head(df5$to,10) | df4$to == z,])
    #   }
    # }
    
    # making the dataframe consist of only unique rows
    df6 <- unique(df3)
    
    # Setting the right color scale
    u1 <- unique(unlist(df6[, c("from", "to")]))
    
    col_mat = viridis(length(u1),alpha = 0.8,option = "H")
    
    chordDiagram(df3, directional = 1, direction.type = c("diffHeight", "arrows"),
                 link.arr.type = "big.arrow", grid.col = col_mat)
    title(paste(i , "joint regulation factors"))
  #}
}

dev.off()

install.packages("alluvial")
library(alluvial)



df3$freq <- sum(df3$from)
alluvial(df3[,1:2], freq = df3$value,col=ifelse( df3$from == "PAX6" & df3$to == df3$to, "red", "gray") )#, freq=df3$from, border=NA)

## alluvial from top 10s
top <- c("ZBTB33","ELK4","PAX6","FOXQ1","TFAP2C","OTX1","IRF6","E2F4","ELK1","TGIF1","RFX5","TP63","ASCL2")
length(top)
cells <- c("LPCs","LNPCs","CB","CjS","CSB")

# add a list in the code for the upset R plot
fullheat <- list()

pdf(paste(resultsdir,"/",cells,'chordsselected.pdf',sep=""),width=6,height=6,paper='special') 
for (i in cells){
  print(i)
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  #factors <- read.table(files$network[files$cell_id==i],header = T)
  #factors <- factors[order(factors$sumScaled,decreasing = T),]
  
  # selecting the top 20 factors
  #top <- factors$factor[1:20]
  
  colnames(diffnetwork) <- c("from","to","value")
  df <- data.frame(diffnetwork)
  
  # selecting
  df3 <- df[df$from %in% top,]
  df3 <- df3[order(df3$value,decreasing = T) & (df3$value >0.8),]
  vec <- df3$to
  vec2 <- unique(myFun(vec))
  df3 <- df3[(df3$to %in% vec2),]
  vec3 <- df3$from
  
  #ind <- ave(rep(1, length(vec)), vec, FUN = length)
  #thresh <- max(ind)
  
  vec4 <- unique(myFun(vec3))
  df3 <- df3[(df3$from %in% vec4),]
  #vec <- df3$to
  #vec2 <- unique(myFun(vec,length(top)))
  # for (z in unique(top)){
  #   print(z)
  #   df4 <- df3
  #   df5 <- df4[df4$from == z,]
  #   df6 <- df5[df5$to == head(df5$to,10),]
  #   df7 = NULL
  
  # Determine if one of the top 10 factors also influences the factor of interest
  # for (q in unique(head(df5$to,10))){
  #   print(q)
  #   df7 = rbind(df7, df4[df4$from == q & df4$to == z,])
  #   df7 = rbind(df7, df4[df4$from == q & df4$to %in% head(df5$to,10),])
  #   
  #   # If a factor is found influencing the factor of interest or another factor in the top 10 then add to original dataframe
  #   if (!is.null(df7)){
  #     df6 = rbind(df6, df7)#[df4$from == q & df4$to %in% head(df5$to,10) | df4$to == z,])
  #   }
  # }
  
  # making the dataframe consist of only unique rows
  df6 <- unique(df3)
  
  # Setting the right color scale
  u1 <- unique(unlist(df6[, c("from", "to")]))
  
  col_mat = viridis(length(u1),alpha = 0.8,option = "H")
  
  chordDiagram(df3, directional = 1, direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow", grid.col = col_mat)
  title(paste(i , "joint regulation factors"))
  #}
}

chordDiagram(df3)

### UPSET intersect of factors in different cell types  WIP
factorsheat1 <- NULL
fullheat <- list()
for (i in unique(files$cell_id)){
  factorsheat1 <- read.table(files$network[files$cell_id==i],header = T)
  factorsheat2 <- factorsheat1[,c("factor","targetScaled")]
  colnames(factorsheat2) <- c("factor",i)
  fullheat[[i]] <- factorsheat2
}

merged.data.frame = Reduce(function(...) merge(..., all=T), fullheat)
merged.data.frame[merged.data.frame == 0] <- NA

merged.data.frame2 <- merged.data.frame
merged.data.frame2[!is.na(merged.data.frame2)] <- 1
merged.data.frame2[is.na(merged.data.frame2)] <- 0
merged.data.frame2$factor <- merged.data.frame$factor
merged.data.frame2 <-merged.data.frame2[rowSums(merged.data.frame2[, -1])>0, ]

for(i in 2:ncol(merged.data.frame2)){ merged.data.frame2[ , i] <- as.integer(merged.data.frame2[ , i]) }
colnoms <- colnames(merged.data.frame2)[colnames(merged.data.frame2) != "factor"]

fill = viridis(length(files$cell_id),alpha = 0.8,option = "H")

pdf(paste(resultsdir,'/intersectionblob.pdf',sep="/") ,width=10,height=5,paper='special')
upset(merged.data.frame2,nintersects = NA,
      sets = colnoms, 
      order.by="freq", matrix.color="black", point.size=1,
      sets.bar.color=fill)
dev.off()

########################################################3




df5 <- NULL
cells <- c("LPCs","LNPCs")#,"CB","CjS","CSB")
col_mat = viridis(length(cells),alpha = 0.8,option = "H")
dfcol <- as.data.frame(col_mat,cells)

for (i in cells){
  print(i)
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  factors <- read.table(files$network[files$cell_id==i],header = T)
  factors <- factors[order(factors$sumScaled,decreasing = T),]
  # selecting the top factors or iteresting factors
  #top <- factors$factor[1:20]
  top <- c("PAX6","TP63","TGIF1","SP3","ELK4","ZBTB33")
  colnames(diffnetwork) <- c("from","to","value")
  df <- data.frame(diffnetwork)
  df <- df[df$value>0.80,]
  df3 <- df[df$from %in% top,]
  df3 <- df3[order(df3$value,decreasing = T),]
  
  # Selecting the top 20 target genes
  top2 <- df3$to[1:20]
  df4 <- df3[df3$to %in% top2,]
  df4$cell <- i
  df4$col <- dfcol[rownames(dfcol) == i,]
  df5 <- rbind(df5, df4)

}

col_mat = viridis(length(cells),alpha = 0.8,option = "H")
ord <- list(NULL, with(df5, order(df5$cell,df5$to), NULL))

pdf(paste(resultsdir,"/",'top20top20sankeyblob.pdf',sep=""),width=8,height=10,paper='special') 
alluvial(df5[,1:2], freq = df5$value,col=df5$col,border = NA,blocks = TRUE,ordering=ord)#,col=ifelse( df3$from == "PAX6" & df3$to == df3$to, "red", "gray") )#, freq=df3$from, border=NA)
dev.off()




pdf(paste(resultsdir,"/",'chordssingleESC.pdf',sep=""),width=6,height=6,paper='special') 
for (i in unique(files$cell_id)){
  print(i)
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  factors <- read.table(files$network[files$cell_id==i],header = T)
  factors <- factors[order(factors$sumScaled,decreasing = T),]
  
  # selecting the top 20 factors
  top <- factors$factor[1:20]
  
  colnames(diffnetwork) <- c("from","to","value")
  df <- data.frame(diffnetwork)
  
  # selecting
  df3 <- df[df$from %in% top,]
  df3 <- df3[order(df3$value,decreasing = T),]

  for (z in unique(top)){
    print(z)
    df4 <- df3
    df5 <- df4[df4$from == z,]
    df6 <- df5[df5$to == head(df5$to,10),]
    df7 = NULL
    
    # Determine if one of the top 10 factors also influences the factor of interest
    for (q in unique(head(df5$to,10))){
      print(q)
      df7 = rbind(df7, df4[df4$from == q & df4$to == z,])
      df7 = rbind(df7, df4[df4$from == q & df4$to %in% head(df5$to,10),])
      
      # If a factor is found influencing the factor of interest or another factor in the top 10 then add to original dataframe
      if (!is.null(df7)){
        df6 = rbind(df6, df7)#[df4$from == q & df4$to %in% head(df5$to,10) | df4$to == z,])
      }
    }
    
    # making the dataframe consist of only unique rows
    df6 <- unique(df6)
    
    # Setting the right color scale
    u1 <- unique(unlist(df6[, c("from", "to")]))
    col_mat = viridis(length(u1),alpha = 0.8,option = "H")
  
    chordDiagram(df6, directional = 1, direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow", grid.col = col_mat)
    title(paste(i , "top 10 regulated genes by", z))
  }
}

dev.off()
