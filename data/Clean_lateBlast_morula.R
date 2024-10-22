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

Blk <- cellQC(sce[,colData(sce)$batch=="Blk"])
Pet <- cellQC(sce[,colData(sce)$batch=="Pet"])
Yan <- cellQC(sce[,colData(sce)$batch=="Yan"])

hist(log10(assay(Blk)), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
abline(v=log10(1), col="blue", lwd=2, lty=2)

Blk <- filterGenes(Blk)
Pet <- filterGenes(Pet)
Yan <- filterGenes(Yan)

# Batch correction and Normalisation
univ <- intersect(rownames(Blk), rownames(Pet))
univ <- intersect(univ, rownames(Yan))

Blk <- Blk[univ, ]
Pet <- Pet[univ, ]
Yan <- Yan[univ, ]

rescaled <- multiBatchNorm(Blk,Pet,Yan) #normalise

Blk <- rescaled[[1]]
Pet <- rescaled[[2]]
Yan <- rescaled[[3]]

# Feature selection
Blk.dec <- modelGeneVar(Blk)
Pet.dec <- modelGeneVar(Pet)
Yan.dec <- modelGeneVar(Yan)
dec <- combineVar(Blk.dec, Pet.dec, Yan.dec)
hvgs <- dec$bio > 0 #highly variable genes

corrected <- rescaleBatches(Blk, Pet, Yan) #remove batch effect

#Update sce with normalised and batch-corrected assay
cells <- c(Blk$sample_accession,Pet$sample_accession,Yan$sample_accession)
sce <- sce[univ, cells]

assay(sce, "batch_corrected") <- assays(corrected)$corrected

# Dimensionality reduction
set.seed(200)
sce <- runPCA(sce, exprs_values="batch_corrected", subset_row=hvgs, ncomponents=25)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)

# Clustering
g <- buildSNNGraph(sce, use.dimred = 'UMAP', k=23)
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
  t=findMarkers(sce, direction="up", assay.type="batch_corrected"),
  wilcox=findMarkers(sce, test="wilcox", direction="up", assay.type="batch_corrected"),
  binom=findMarkers(sce, test="binom", direction="up", assay.type="batch_corrected")
)

interesting <- markers[[4]] 
interesting[1:10,1:9]
plotUMAP(sce, colour_by="EXOSC3P1", by_exprs_values="batch_corrected")

sce <- logNormCounts(sce)

early_blastocyst_morula_sce <- sce
early_blastocyst_morula_hvgs <- hvgs

save(early_blastocyst_morula_sce, early_blastocyst_morula_hvgs, file = "early_blastocyst_morula.RData")
