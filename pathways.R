# Expression pathways
library(pathview)

paths<-function(cellType, pathway, sce, normalistion = "batch_corrected", avg="median")
{
  pathway <- substr(pathway, 4, 8)
  if (cellType=="Epiblast"){
    colour <- "#255836"
  } else if (cellType=="Trophectoderm"){
    colour <- "#1e4e78"
  } else {
    colour <- "#a3202b"
  }
  cells <- assay(sce, normalistion)[,!colData(sce)$cell_type!=cellType]
  if (avg=="mean"){
    gene.data <- rowMeans(cells)
  } else {
    gene.data <- rowMedians(cells)
    names(gene.data) <- rownames(cells)
  }
  pv.out <- pathview(gene.data = gene.data, pathway.id = pathway, gene.idtype = "symbol", limit=ceiling(max(gene.data)), both.dirs = F, high = colour, mid = "white",
                     out.suffix = paste(avg,".",cellType, sep=""), node.sum = "max")
  return(pv.out)
}

pathways<-list("04010","04012","04014","04015","04020","04022","04024","04064"
               ,"04066","04068","04071","04072","04150","04151","04152","04310"
               ,"04330","04340","04350","04370","04371","04390","04630","04668")
