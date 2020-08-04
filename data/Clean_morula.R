library(scater)
library(scran)
library(batchelor)

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

# Get rid of no-show and lowly expressed genes
filterGenes<-function(cells)
{
  cells <- cells[rowSums(counts(cells)) != 0, ]
  # Only keep genes with avg. expression across cells > 0.04
  cells <- cells[rowMeans(assay(cells)) >= 1,]
  return(cells)
}

# START
load("data/sce_integrated_lateBlast_morula.RData")

# Quality control
# Filter genes without NA chr
sce <- sce[!is.na(rowData(sce)$chr),]
sce <- sce[,colData(sce)$cell_type == "Morula"]

Pet <- cellQC(sce[,colData(sce)$batch=="Pet"])
Yan <- cellQC(sce[,colData(sce)$batch=="Yan"])

hist(log10(assay(Blk)), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
abline(v=log10(1), col="blue", lwd=2, lty=2)

Pet <- filterGenes(Pet)
Yan <- filterGenes(Yan)

# Batch correction and Normalisation
univ <- intersect(rownames(Pet), rownames(Yan))

Pet <- Pet[univ, ]
Yan <- Yan[univ, ]

rescaled <- multiBatchNorm(Pet,Yan) #normalise

Pet <- rescaled[[1]]
Yan <- rescaled[[2]]

# Feature selection
Pet.dec <- modelGeneVar(Pet)
Yan.dec <- modelGeneVar(Yan)
dec <- combineVar(Pet.dec, Yan.dec)
hvgs <- dec$bio > 0 #highly variable genes

corrected <- rescaleBatches(Pet, Yan) #remove batch effect

#Update sce with normalised and batch-corrected assay
cells <- c(Pet$sample_accession,Yan$sample_accession)
sce <- sce[univ, cells]

assay(sce, "batch_corrected") <- assays(corrected)$corrected

# Dimensionality reduction
set.seed(100)
sce <- runPCA(sce, exprs_values="batch_corrected", subset_row=hvgs, ncomponents=50)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)

# Clustering
g <- buildSNNGraph(sce, use.dimred = 'UMAP', k=23)
colLabels(sce) <- factor(igraph::cluster_louvain(g)$membership)

# Visualization.
plotUMAP(sce, colour_by="NANOG", by_exprs_values="batch_corrected")
plotUMAP(sce, colour_by="GATA3", by_exprs_values="batch_corrected")
#plotUMAP(sce, colour_by="SOX17", by_exprs_values="batch_corrected")
plotUMAP(sce, colour_by="batch", by_exprs_values="batch_corrected")
plotUMAP(sce, colour_by="cell_type", by_exprs_values="batch_corrected")

# visualizing important markers for clusters
plotUMAP(sce, colour_by="label", by_exprs_values="batch_corrected") # show clusters

#find markers
markers <- multiMarkerStats(
  t=findMarkers(sce, direction="up", assay.type="batch_corrected"),
  wilcox=findMarkers(sce, test="wilcox", direction="up", assay.type="batch_corrected"),
  binom=findMarkers(sce, test="binom", direction="up", assay.type="batch_corrected")
)

interesting <- markers[[2]] 
interesting[1:10,1:9]
plotUMAP(sce, colour_by="HNRNPCL3", by_exprs_values="batch_corrected")

sce <- logNormCounts(sce)

morula_sce <- sce
morula_hvgs <- hvgs

save(morula_sce, morula_hvgs, file = "morula.RData")
