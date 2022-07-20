ClusteringResSeurObj<- function(SeuratObject,dimensions=24) {
  chosen_dims <- dimensions
  ###########################################################################
  # Perform clustering with resolution, to see which resolution fits the best
  pdf(paste(resultsdir, "8.clustering_resolution.pdf", sep = '/'), width = 15, height = 8)
  seur_obj <- SeuratObject
  seur_obj <- FindNeighbors(seur_obj, dims = 1:chosen_dims)
  for (i in 1:20){
    res <- i/20
    #cluster_variable_name <- paste0("RNA_snn_res.", res)
    cluster_variable_name <- paste0("RNA_snn_res.", res)
    seur_obj <- FindClusters(seur_obj, verbose = FALSE, resolution = res, graph.name = "RNA_snn")
    seur_obj <- BuildClusterTree(seur_obj)
    
    if (i != 1){
      
      #
      # seur_obj <- FindClusters(seur_obj, verbose = FALSE, resolution = res)
      # seur_obj <- BuildClusterTree(seur_obj)
      #p1 <- DimPlot(seur_obj,label=TRUE, group.by = 'seurat_clusters') + ggtitle("Louvain Clustering") + ggtitle(paste0("cluster resolution ", res))
      p2 <- ggplot(seur_obj@meta.data, aes(eval(parse(text= i))))+geom_bar(stat="count")
      print(clustree(
        seur_obj,
        prefix = "RNA_snn_res.",
        exprs = c("data", "counts", "scale.data"),
        assay = NULL,
        node_colour = "sc3_stability"
      ))
      print(p2)
    }
    #print(p1)
    
  }
  dev.off()
}