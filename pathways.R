# Expression pathways
library(pathview)

paths<-function(cellType,colour,pathway,avg)
{
  cells <- assay(sce)[,!colData(sce)$cell_type!=cellType]
  if (avg=="mean"){
    gene.data <- rowMeans(cells)
  } else {
    gene.data <- rowMedians(cells)
    names(gene.data) <- rownames(cells)
  }
  pv.out <- pathview(gene.data = gene.data, pathway.id = pathway, gene.idtype = "symbol", limit=ceiling(max(gene.data)), both.dirs = F, high = colour, mid = "white",
                     out.suffix = paste(avg,".",cellType, sep=""), node.sum = "max")
  #pv.out$plot.data.gene
}

pathways<-list("04010","04012","04014","04015","04020","04022","04024","04064"
               ,"04066","04068","04071","04072","04150","04151","04152","04310"
               ,"04330","04340","04350","04370","04371","04390","04630","04668")

#single pathway
p <- "04010"
paths("Epiblast", "#255836", p, "mean")
paths("Trophectoderm", "#1e4e78", p, "mean")
paths("Primitive endoderm", "#a3202b", p, "mean")


#all pathways
for (p in pathways){
  paths("Epiblast", "#255836", p, "median")
  paths("Trophectoderm", "#1e4e78", p, "median")
  paths("Primitive endoderm", "#a3202b", p, "median")
}
