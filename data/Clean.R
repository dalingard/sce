library(scater)
library(scran)
library(batchelor)
library(patchwork)

load("data/sce_integrated_lateBlast.RData")

# Quality control Genes
# Filter genes without NA chr
sce <- sce[!is.na(rowData(sce)$chr),]

# Quality Control Cells
cellQC<-function(cells)
{
  # All genes in chromosome MT are considered mitochondrial
  is.mito <- rowData(cells)$chr == "MT"
  qcstats <- perCellQCMetrics(cells, subsets=list(Mito=is.mito))
  
  #filtering out high Mito%, Low no. of features and abnormal no. gene expressions
  filtered <- quickPerCellQC(qcstats, lib_size = "sum", n_features = "detected",
                             percent_subsets=c("subsets_Mito_percent"))
  
  
  # Diagnostic plots
  colData(cells) <- cbind(colData(cells), qcstats)
  cells$discard <- filtered$discard

  # Just to highlight how patchwork works but you can stick to gridExtra
  p_libsize <- plotColData(cells, x="batch", y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count")
  p_genes <- plotColData(cells, x="batch", y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features")
  p_mito <- plotColData(cells, x="batch", y="subsets_Mito_percent",
                        colour_by="discard") + ggtitle("Mito percent")
  
  cells <- cells[, !filtered$discard]
  return(list(flt_sce = cells, plots = (p_libsize + p_genes + p_mito)))
}

Blk_res <- cellQC(sce[, colData(sce)$batch == "Blk"])
Blk_res$plots
Blk <- Blk_res$flt_sce
Pet_res <- cellQC(sce[, colData(sce)$batch == "Pet"])
Pet_res$plots
Pet <- Pet_res$flt_sce
Yan_res <- cellQC(sce[, colData(sce)$batch == "Yan"])
Yan_res$plots
Yan <- Yan_res$flt_sce

# Get rid of no-show and lowly expressed genes
geneQC <- function(cells){
  cells <- cells[rowSums(counts(cells)) != 0, ]
  
  # Only keep genes with avg. expression across cells >= 1
  ave.counts <- rowMeans(counts(cells))
  keep <- ave.counts >= 1
  
  cells <- cells[keep, ]
  return(cells)
}
Blk <- geneQC(Blk)
Pet <- geneQC(Pet)
Yan <- geneQC(Yan)

# Subset all batches to the common “universe” of features
univ <- intersect(rownames(Blk), rownames(Pet))
univ <- intersect(univ, rownames(Yan))

Blk <- Blk[univ, ]
Pet <- Pet[univ, ]
Yan <- Yan[univ, ]

# Normalisation and batch correction
rescaled <- multiBatchNorm(Blk,Pet,Yan) #normalise

Blk <- rescaled[[1]]
Pet <- rescaled[[2]]
Yan <- rescaled[[3]]

# Identification of highly variable genes
Blk.dec <- modelGeneVar(Blk)
Pet.dec <- modelGeneVar(Pet)
Yan.dec <- modelGeneVar(Yan)
dec <- combineVar(Blk.dec, Pet.dec, Yan.dec)
hvgs <- dec$bio > 0

corrected <- rescaleBatches(Blk, Pet, Yan) #remove batch effect

#Update sce with normalised and batch-corrected assay
cells <- c(Blk$sample_accession,Pet$sample_accession,Yan$sample_accession)
sce <- sce[univ, cells]

assay(sce, "batch_corrected") <- assays(corrected)$corrected

# Dimensionality reduction
set.seed(9999)
sce <- runPCA(sce, ncomponents=25, 
              exprs_values="batch_corrected", 
              subset_row=hvgs)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)

# Clustering
g <- buildSNNGraph(sce, use.dimred = 'UMAP', k=28)
colLabels(sce) <- factor(igraph::cluster_louvain(g)$membership)


# Visualization.
plotUMAP(sce, colour_by="NANOG", by_exprs_values="batch_corrected")
plotUMAP(sce, colour_by="GATA3", by_exprs_values="batch_corrected")
plotUMAP(sce, colour_by="SOX17", by_exprs_values="batch_corrected")
plotUMAP(sce, colour_by="batch", by_exprs_values="batch_corrected")
plotUMAP(sce, colour_by="cell_type", by_exprs_values="batch_corrected")

# visualizing important markers for clusters
plotUMAP(sce, colour_by="label", by_exprs_values="batch_corrected") # show clusters

#find markers
markers <- multiMarkerStats(
  t=findMarkers(sce, direction="up"),
  wilcox=findMarkers(sce, test="wilcox", direction="up"),
  binom=findMarkers(sce, test="binom", direction="up")
)

interesting <- markers[[3]] 
interesting[1:10,1:9]

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


