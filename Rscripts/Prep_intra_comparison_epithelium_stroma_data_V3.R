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




# Adjust condition if epithelial or stromal population and generate new dataframes from it
# stromal with joined number (due to CPM q-quantile normalization) 
# of reads from CSSC and CF (sub) and CSSC, CF, MF, SK, TSK (full)
lakocountfile3 <- lakocountfile






write.table(data.frame("ID"=rownames(lakocountfile3),lakocountfile3),file = paste0(resultsdir,"/", 'cornea_stromal_sub_pseudobulk.tsv'),sep = "\t", quote = F,row.names=F,col.names = T)

# corneal epithelial with joined number of reads from LSC, LESC, LE, Cj and CE
lakocountfile3 <- lakocountfile

write.table(data.frame("ID"=rownames(lakocountfile3),lakocountfile3),file = paste0(resultsdir,"/", 'cornea_epithelial_pseudobulk.tsv'),sep = "\t", quote = F,row.names=F,col.names = T)
