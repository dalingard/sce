reduceDimentions<-function(cells, hvgs, type="UMAP", colourby="cell_type"){
  if (colourby=="-"){
    colourby <- 'cell_type'
  }
  set.seed(100)
  cells <- runPCA(cells, exprs_values="batch_corrected", subset_row=hvgs, ncomponents=25)
  if (type=="UMAP"){
    cells <- runUMAP(cells, dimred = 'PCA', external_neighbors=TRUE)
    # Clustering
    g <- buildSNNGraph(cells, use.dimred = 'UMAP', k=28)
    colLabels(cells) <- factor(igraph::cluster_louvain(g)$membership)
    
    # Visualization.
    plotUMAP(cells, colour_by=colourby, by_exprs_values="batch_corrected")
  } else {
    # Clustering
    g <- buildSNNGraph(cells, use.dimred = 'PCA', k=28)
    colLabels(cells) <- factor(igraph::cluster_louvain(g)$membership)
    
    # Visualization.
    plotPCA(cells, colour_by=colourby, by_exprs_values="batch_corrected")
  }
}
