reduceDimentions<-function(cells, hvgs, geneOfInterest, type="UMAP", colourby="cell_type", norm_type="batch_corrected"){
  if (colourby=="Gene expression"){
    colourby <- geneOfInterest
  } else if (colourby=="Cell type" || colourby=="Cell-type based"){
    colourby <- "cell_type"
  } else if (colourby=="Batch") {
    colourby <- "batch"
  } else if (colourby=="Cluster-based"){
    colourby <- "label"
  }
  
  if (norm_type=="fpkm" || norm_type=="tpm"){
    counts <- assay(cells, norm_type)
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    assay(cells, "log2counts") <- log2(t(t(counts)/size.factors) + 1)
    norm_type <- "log2counts"
  }
  set.seed(100) 
  cells <- runPCA(cells, exprs_values="batch_corrected", subset_row=hvgs, ncomponents=25)
  if (type=="UMAP"){
    cells <- runUMAP(cells, dimred = 'PCA', external_neighbors=TRUE)
    # Clustering
    g <- buildSNNGraph(cells, use.dimred = 'UMAP', k=28)
    colLabels(cells) <- factor(igraph::cluster_louvain(g)$membership)
    
    # Visualization.
    plotUMAP(cells, colour_by=colourby, by_exprs_values=norm_type, point_alpha=1, point_size=3)
  } else {
    # Clustering
    g <- buildSNNGraph(cells, use.dimred = 'PCA', k=28)
    colLabels(cells) <- factor(igraph::cluster_louvain(g)$membership)
    
    # Visualization.
    plotPCA(cells, colour_by=colourby, by_exprs_values=norm_type, point_alpha=1, point_size=3)
  }
}
