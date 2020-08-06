# Expression pathways
paths<-function(cellType, pathway, sce, norm_type = "batch_corrected", avg="median")
{
  pathway <- substr(pathway, 4, 8)
  if (cellType=="Epiblast"){
    colour <- "#255836"
  } else if (cellType=="Trophectoderm"){
    colour <- "#1e4e78"
  } else if (cellType=="Primitive endoderm") {
    colour <- "#a3202b"
  } else if (cellType=="Morula"){
    colour <- "#108783"
  } else if (cellType=="t2iL+Go H9"){
    colour <- "#b8a904"
  } else if (cellType=="E8 H9"){
    colour <- "#b84304"
  }
  
  if (norm_type=="fpkm" || norm_type=="tpm"){
    counts <- assay(sce, norm_type)
    assay(sce, "log2counts") <- log2(assay(sce, norm_type) + 1)
    norm_type <- "log2counts"
  }
  
  cells <- assay(sce, norm_type)[,colData(sce)$cell_type==cellType]
  
  if (avg=="mean"){
    gene.data <- rowMeans(cells)
  } else {
    gene.data <- rowMedians(cells)
    names(gene.data) <- rownames(cells)
  }
  pv.out <- pathview(gene.data = gene.data, pathway.id = pathway, gene.idtype = "symbol", limit=ceiling(max(gene.data)), both.dirs = F, high = colour, mid = "white",
                     out.suffix = paste(avg,".",cellType, sep=""), node.sum = "max", kegg.dir = "./pathway_data")
  return(pv.out)
}
