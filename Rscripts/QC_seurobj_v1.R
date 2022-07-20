QCSeurObj <- function(SeuratObject,pct_MT=30) {
  ###########################################################################
  ## Filtering of the datasets
  seur_obj_all <- SeuratObject
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  seur_obj_all[["percent.mt"]] <- PercentageFeatureSet(seur_obj_all, pattern = "^MT-")
  
  # QC tresholds of counts, features and mt percentage
  seur_obj <- subset(seur_obj_all, subset = nCount_RNA > 2000 & nFeature_RNA > 1000 & percent.mt < pct_MT)
  
  ###########################################################################
  # Plotting percentage MTs before and after
  if (length(unique(unlist(seur_obj$orig.ident)))<6){
    width2 <- 6
  } else {
    width2 <- length(unique(unlist(seur_obj$orig.ident)))/2
  }
  
  pdf(paste(resultsdir,'3.QC_depth_ERCC_MT.pdf',sep="/") ,width=width2,height=6,paper='special')
  print(VlnPlot(object = seur_obj_all, features = ("nCount_RNA"), assay = 'RNA') + scale_y_continuous(limits = c(0,20000)) + ggtitle('total UMI counts all cells'))
  print(VlnPlot(object = seur_obj, features = c("nCount_RNA"), assay = 'RNA') + scale_y_continuous(limits = c(0,20000)) + ggtitle('total UMI counts filtered cells'))
  print(VlnPlot(object = seur_obj_all, features = c("nFeature_RNA"), assay = 'RNA') + scale_y_continuous(limits = c(0,5000)) + ggtitle('genes measured all cells'))
  print(VlnPlot(object = seur_obj, features = c("nFeature_RNA"), assay = 'RNA')  + scale_y_continuous(limits = c(0,5000)) + ggtitle('genes measured filtered cells'))
  print(VlnPlot(object = seur_obj_all, features = c("percent.mt"))+ scale_y_continuous(limits = c(0,40)) + ggtitle('% mitochondrial reads all cells'))
  print(VlnPlot(object = seur_obj, features = c("percent.mt"))+ scale_y_continuous(limits = c(0,40)) + ggtitle('% mitochondrial reads filtered cells'))
  print(FeatureScatter(seur_obj_all, feature1 = ("nCount_RNA"), feature2 = ("nFeature_RNA")) + ggtitle('% counts vs genes measured all cells'))
  print(FeatureScatter(seur_obj, feature1 = ("nCount_RNA"), feature2 = ("nFeature_RNA")) + ggtitle('% counts vs genes measured filtered cells'))
  dev.off()
  
  ###########################################################################
  # normalization of the data
  Idents(seur_obj) <- "cell_type"
  Idents(seur_obj_all) <- "cell_type"
  
  seur_obj <- NormalizeData(
    object = seur_obj, assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  
  ###########################################################################
  # finding variable features within the datasets and determine cell_cycle influence
  seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000)
  seur_obj <- ScaleData(seur_obj, features = rownames(seur_obj), assay = 'RNA')
  
  seur_obj <- CellCycleScoring(seur_obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
  seur_obj <- RunPCA(seur_obj, features = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))
  
  pdf(paste(resultsdir,'5.cell_cycle_markers.pdf',sep="/") ,width=12,height=6,paper='special')
  Idents(seur_obj) <- "Phase"
  print(PCAPlot(seur_obj))
  print(RidgePlot(seur_obj, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), assay= 'RNA', ncol = 2, slot = "data"))
  print(RidgePlot(seur_obj, features = c("S.Score", "G2M.Score")))
  dev.off()
  
  return(seur_obj)
}