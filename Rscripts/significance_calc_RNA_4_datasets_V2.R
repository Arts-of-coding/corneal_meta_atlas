# Deseq2 significance calculation peaks
##based upon http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
library(dplyr)
library(stringr)
# loading in the count matrix and coldata files

# For the deseq2 matrix vs ESC
lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220311_old/pseudobulk_reps_DE_datasets_markers_split2.tsv', sep = '\t', header = TRUE, row.names = 1)
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/scANANSE/RNA_DE/"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

#collist<- c("LiCo","StCSC")


# Add the ESC counts to the count file
ESCcountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/ESC_ANANSE/ESCcountssel.tsv', sep = '\t', header = TRUE, row.names = 1)
ESCcountfile$X <- rownames(ESCcountfile)
colnames(ESCcountfile)<- c("ESC_1","ESC_2","Gene.stable.ID")

# Convert gene stable ID to gene names
ESCnamefile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/ESC_ANANSE/mart_export.txt', sep = '\t', header = TRUE)

ESCjoined <- ESCcountfile %>% right_join(ESCnamefile, by="Gene.stable.ID")
ESCjoined <- ESCjoined %>% mutate_all(na_if,"")
ESCjoined <- ESCjoined[!is.na(ESCjoined$Gene.name),]
ESCjoined <- ESCjoined[!duplicated(ESCjoined$Gene.name),]
rownames(ESCjoined) <- ESCjoined$Gene.name

cornearna <- lakocountfile
lakocountfile$Gene.name<- rownames(lakocountfile)
zz <- merge(ESCjoined, lakocountfile, all = TRUE)
zz[is.na(zz)] <- 0

# Join the ocunts of ESC data
lakocountfile2 <- zz

lakocountfile2 <- lakocountfile2[!duplicated(lakocountfile2$Gene.name),]
lakocountfile2[is.na(lakocountfile2)] <- 0

rownames(lakocountfile2)<- lakocountfile2$Gene.name

lakocountfile2$Gene.stable.ID <- NULL
lakocountfile2$Gene.name<- NULL


# If rowsums are smaller than one, this gene is not important for GRN analysis comparison with ESC
lakocountfile3 <- lakocountfile2[rowSums(lakocountfile2)>1,]

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
