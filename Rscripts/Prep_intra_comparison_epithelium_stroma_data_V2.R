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

# For the deseq2 matrix vs ESC
lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq_integrated_cornea_4_datasets/20220309/pseudobulkdf_full.tsv', sep = '\t', header = TRUE, row.names = 1)

## setting up results directory




# Adjust condition if epithelial or stromal population and generate new dataframes from it
# stromal with joined number (due to CPM q-quantile normalization) 
# of reads from CSSC and CF (sub) and CSSC, CF, MF, SK, TSK (full)
lakocountfile3 <- lakocountfile

write.table(data.frame("ID"=rownames(lakocountfile3),lakocountfile3),file = paste0(resultsdir,"/", 'cornea_stromal_full_pseudobulk.tsv'),sep = "\t", quote = F,row.names=F,col.names = T)

write.table(data.frame("ID"=rownames(lakocountfile3),lakocountfile3),file = paste0(resultsdir,"/", 'cornea_stromal_sub_pseudobulk.tsv'),sep = "\t", quote = F,row.names=F,col.names = T)

# corneal epithelial with joined number of reads from LSC, LESC, LE, Cj and CE
lakocountfile3 <- lakocountfile

write.table(data.frame("ID"=rownames(lakocountfile3),lakocountfile3),file = paste0(resultsdir,"/", 'cornea_epithelial_pseudobulk.tsv'),sep = "\t", quote = F,row.names=F,col.names = T)
