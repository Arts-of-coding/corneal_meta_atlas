PCASeurObj <- function (SeuratObject,
                        pcs=50) {
  seur_obj <- SeuratObject
  ###########################################################################
  # PCA plots for all genes
  seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = pcs)
  mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
  pca <- seur_obj[["pca"]]
  
  # Get the total variance:
  total_variance <- sum(matrixStats::rowVars(mat))
  eigValues = (pca@stdev)^2  ## EigenValues
  varExplained = eigValues / total_variance
  
  # Checking the first PC
  Stdev(object = seur_obj[["pca"]])[1]
  
  Idents(seur_obj) <- "cell_type"
  pdf(paste0(resultsdir,'/6.Principle_components.pdf') ,width=12,height=6,paper='special')
  print(ElbowPlot(seur_obj))
  for (pcs in c(1:15)){
    pc1_viz <- pcs*2-1
    pc2_viz <- pcs*2
    y_label = paste0(paste0(paste0("PC",pc2_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
    x_label = paste0(paste0(paste0("PC",pc1_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
    PC_dimred <- DimPlot(seur_obj, reduction = "pca", dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
    PC1_genes <- DimHeatmap(seur_obj, dims = c(pc1_viz), fast = FALSE)
    PC2_genes <- DimHeatmap(seur_obj, dims = c(pc2_viz), fast = FALSE)
    final_plot <- grid.arrange(as_grob(PC_dimred),as_grob(PC2_genes),as_grob(PC1_genes),
                               ncol=2,
                               as.table=TRUE)
    #heights=c(3,1))
    print(final_plot)}
  dev.off()
  return(seur_obj)
}