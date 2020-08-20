reduceDimentions<-function(cells, geneOfInterest, type="UMAP", colourby="cell_type", norm_type="batch_corrected"){
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
    assay(cells, "log2counts") <- log2(assay(cells, norm_type) + 1)
    norm_type <- "log2counts"
  }
  
  # visualize
  if (type=="UMAP"){
    plotUMAP(cells, colour_by=colourby, by_exprs_values=norm_type, point_alpha=1, point_size=3)
  } else {
    plotPCA(cells, colour_by=colourby, by_exprs_values=norm_type, point_alpha=1, point_size=3)
  }
}

