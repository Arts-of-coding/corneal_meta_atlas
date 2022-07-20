# Deseq2 significance calculation peaks
##based upon http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
library(dplyr)
library(stringr)
library(SingleCellExperiment)
library(Seurat)
# loading in the count matrix and coldata files
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/RNA_intra_cpm/"
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

pseudobulk_df<- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220309/pseudobulk_reps_DE_datasets_markers_split2.tsv',header = T,row.names = 1)


# For the deseq2 matrix vs ESC
#lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220309/pseudobulkdf_full.tsv', sep = '\t', header = TRUE, row.names = 1)

## setting up results directory
# If rowsums are smaller than one, this gene is not important for GRN analysis comparison with ESC
lakocountfile3 <- pseudobulk_df[rowSums(pseudobulk_df)>1,]

# for complex heatmap
#lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/pseudobulk.tsv', sep = '\t', header = TRUE, row.names = 1)
lakovst <- as.matrix(lakocountfile3)

# Generate coldata dataframe
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

# adjusting the conditions for stromal
# colnames(lakovst)<-gsub("CSSC", "stromal", colnames(lakovst))
# colnames(lakovst)<-gsub("CF", "stromal", colnames(lakovst))

coldata$condition<-gsub("CSSC", "stromal", coldata$condition)
coldata$condition<-gsub("CF", "stromal", coldata$condition)

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

#lakovstcalc <- results(dds2)
#lakovstcalc3 <- as.data.frame(lakovstcalc@listData$log2FoldChange,row.names = rownames(lakovst))

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

# Generate coldata dataframe
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

coldata$condition<-gsub("LSC", "epi", coldata$condition)
coldata$condition<-gsub("LESC", "epi", coldata$condition)
coldata$condition<-gsub("LE", "epi", coldata$condition)
coldata$condition<-gsub("Cj", "epi", coldata$condition)
coldata$condition<-gsub("CE", "epi", coldata$condition)

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

#lakovstcalc <- results(dds2)
#lakovstcalc3 <- as.data.frame(lakovstcalc@listData$log2FoldChange,row.names = rownames(lakovst))

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

# Now with stromal vs epi included
coldata$condition<-gsub("CSSC", "stromal", coldata$condition)
coldata$condition<-gsub("CF", "stromal", coldata$condition)

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

#lakovstcalc <- results(dds2)
#lakovstcalc3 <- as.data.frame(lakovstcalc@listData$log2FoldChange,row.names = rownames(lakovst))

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

# Generating TPM files
# For the deseq2 matrix vs ESC
lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220309/pseudobulkdf_full.tsv', sep = '\t', header = TRUE, row.names = 1)


# If rowsums are smaller than one, this gene is not important for GRN analysis comparison
lakocountfile2 <- lakocountfile[rowSums(lakocountfile)>1,]

lakocountfile3 <- lakocountfile2 %>% 
rowwise() %>% 
mutate(epi=sum(LESC, LSC, LE, Cj, CE, na.rm = F)) 

lakocountfile3 <- lakocountfile3 %>% 
rowwise() %>% 
mutate(stromal=sum(CSSC, CF, na.rm = F))


lakocountfile3 <- lakocountfile3[,c("epi","stromal")]
rownames(lakocountfile3) <- rownames(lakocountfile2)

write.table(data.frame("ID"=rownames(lakocountfile3),lakocountfile3),file = paste0(resultsdir,"/", 'intra_comp_pseudobulk.tsv'),sep = "\t", quote = F,row.names=F,col.names = T)
