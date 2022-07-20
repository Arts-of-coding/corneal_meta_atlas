QCSeurObjClustering <- function (SeuratObject,dimensions=30) {
  seur_obj <- SeuratObject
  seur_obj <- RunUMAP(seur_obj, dims = 1:dimensions)
  
  # plotting the quality control upon the UMAP to see if clustering is not driven by a specific factor
  pdf(paste(resultsdir,'7a.norm_umap.pdf',sep="/") ,width=8,height=8,paper='special')
  print(DimPlot(seur_obj, label=FALSE)+ ggtitle("original identity"))
  print(FeaturePlot(seur_obj, features = 'nCount_RNA', cols =c('white','purple')))
  print(FeaturePlot(seur_obj, features = 'nFeature_RNA', cols =c('white','purple')))
  print(FeaturePlot(seur_obj, features = 'percent.mt', cols =c('white','purple')))
  print(DimPlot(seur_obj, group.by = 'Phase'))
  print(FeaturePlot(seur_obj, features = 'S.Score', cols =c('white','dodgerblue3')))
  print(FeaturePlot(seur_obj, features = 'G2M.Score', cols =c('white','green4')))
  dev.off()
  return(seur_obj)
}