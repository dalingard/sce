library(scater)
library(scran)
library(batchelor)

load("data/sce_integrated_lateBlast.RData")

# Quality control Genes
# Filter genes without NA chr
sce <- sce[!is.na(rowData(sce)$chr),]

# Filter out flat(ish) genes
row_sub = apply(assay(sce, "fpkm"), 1, function(row) all(row < 10))
sce <- sce[!row_sub[rowData(sce)$symbol],]

# Quality Control Cells
cellQC<-function(cells)
{
  is.mito <- grepl("^MT-", rownames(cells))
  qcstats <- perCellQCMetrics(cells, subsets=list(Mito=is.mito))
  #filtering out high Mito%, Low no. of features and abnormal no. gene expressions
  filtered <- quickPerCellQC(qcstats, percent_subsets=c("subsets_Mito_percent","sum","detected"))
  
  
  # # Diagnostic plots
  # colData(cells) <- cbind(colData(cells), qcstats)
  # cells$discard <- filtered$discard
  # 
  # gridExtra::grid.arrange(
  #   plotColData(cells, x="batch", y="sum", colour_by="discard") +
  #     scale_y_log10() + ggtitle("Total count"),
  #   plotColData(cells, x="batch", y="detected", colour_by="discard") +
  #     scale_y_log10() + ggtitle("Detected features"),
  #   plotColData(cells, x="batch", y="subsets_Mito_percent",
  #               colour_by="discard") + ggtitle("Mito percent"),
  #   ncol=1
  # )
  
  cells <- cells[, !filtered$discard]
  return(cells)
}

# Batch correction and Normalisation
Blk <- cellQC(sce[,!colData(sce)$batch!="Blk"])
Pet <- cellQC(sce[,!colData(sce)$batch!="Pet"])
Yan <- cellQC(sce[,!colData(sce)$batch!="Yan"])

rescaled <- multiBatchNorm(Blk,Pet,Yan) #normalise

Blk <- rescaled[[1]]
Pet <- rescaled[[2]]
Yan <- rescaled[[3]]

# Feature selection
Blk.dec <- modelGeneVar(Blk)
Pet.dec <- modelGeneVar(Pet)
Yan.dec <- modelGeneVar(Yan)
dec <- combineVar(Blk.dec, Pet.dec, Yan.dec)
hvgs <- dec$bio > 0

corrected <- rescaleBatches(Blk, Pet, Yan) #remove batch effect

#Update sce with normalised and batch-corrected assay
cells <- c(Blk$sample_accession,Pet$sample_accession,Yan$sample_accession)
sce <- sce[,cells]
assay(sce) <- NULL
assay(sce) <- NULL
assay(sce) <- NULL
assay(sce) <- NULL
assay(sce, "logcounts", withDimnames = FALSE) <- assay(corrected)

# Dimensionality reduction
set.seed(9999)
sce <- runPCA(sce, exprs_values="logcounts", subset_row=hvgs, ncomponents=3)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)

# Clustering
g <- buildSNNGraph(sce, use.dimred = 'UMAP', k=28)
colLabels(sce) <- factor(igraph::cluster_louvain(g)$membership)


# Visualization.
plotUMAP(sce, colour_by="NANOG")
plotUMAP(sce, colour_by="GATA3")
plotUMAP(sce, colour_by="SOX17")
plotUMAP(sce, colour_by="batch")
plotUMAP(sce, colour_by="cell_type")

# visualizing important markers for clusters
plotUMAP(sce, colour_by="label") # show clusters

#find markers
markers <- multiMarkerStats(
  t=findMarkers(sce, direction="up"),
  wilcox=findMarkers(sce, test="wilcox", direction="up"),
  binom=findMarkers(sce, test="binom", direction="up")
)

interesting <- markers[[3]] 
interesting[1:10,1:9]
plotUMAP(sce, colour_by="CNTFR")


# Expression pathways
library(pathview)

paths<-function(cellType,colour)
{
  cells <- assay(sce)[,!colData(sce)$cell_type!=cellType]
  gene.data <- rowMeans(cells)
  pv.out <- pathview(gene.data = gene.data, pathway.id = "04630", gene.idtype = "symbol", limit=ceiling(max(gene.data)), both.dirs = F, high = colour,
                     out.suffix = cellType, node.sum = "max")
  #pv.out$plot.data.gene
}

paths("Epiblast", "deepskyblue")
paths("Trophectoderm", "green3")
paths("Primitive endoderm", "orange")
