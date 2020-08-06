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
  #cells <- cells[rowSums(counts(cells)) != 0, ]
  # Only keep genes with avg. expression across cells > 0.04
  #cells <- cells[rowMeans(assay(cells)) >= 1,]
  cells <- cells[rowSums(counts(cells) == 0) <= 0.98*length(colnames(counts(cells))), ]
  return(cells)
}

# START
load("data/sce_integrated_lateBlast_morula_hesc.RData")

# Quality control
# Filter genes without NA chr
sce <- sce[!is.na(rowData(sce)$chr),]
sce <- sce[,colData(sce)$cell_type != "Epiblast" & colData(sce)$cell_type != "Primitive endoderm" & colData(sce)$cell_type != "Trophectoderm"]

Pet <- cellQC(sce[,colData(sce)$batch=="Pet"])
Yan <- cellQC(sce[,colData(sce)$batch=="Yan"])
mes1 <- cellQC(sce[,colData(sce)$batch=="Mes_1"])
mes2 <- cellQC(sce[,colData(sce)$batch=="Mes_2"])
mes3 <- cellQC(sce[,colData(sce)$batch=="Mes_3"])

hist(log10(assay(mes3)), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
abline(v=log10(1), col="blue", lwd=2, lty=2)

Pet <- filterGenes(Pet)
Yan <- filterGenes(Yan)
mes1 <- filterGenes(mes1)
mes2 <- filterGenes(mes2)
mes3 <- filterGenes(mes3)

# Batch correction and Normalisation
univ <- intersect(rownames(Pet), rownames(Yan))
univ <- intersect(univ, rownames(mes1))
univ <- intersect(univ, rownames(mes2))
univ <- intersect(univ, rownames(mes3))

Pet <- Pet[univ, ]
Yan <- Yan[univ, ]
mes1 <- mes1[univ, ]
mes2 <- mes2[univ, ]
mes3 <- mes3[univ, ]

rescaled <- multiBatchNorm(Pet,Yan,mes1,mes2,mes3) #normalise

Pet <- rescaled[[1]]
Yan <- rescaled[[2]]
mes1 <- rescaled[[3]]
mes2 <- rescaled[[4]]
mes3 <- rescaled[[5]]

# Feature selection
Pet.dec <- modelGeneVar(Pet)
Yan.dec <- modelGeneVar(Yan)
mes1.dec <- modelGeneVar(mes1)
mes2.dec <- modelGeneVar(mes2)
mes3.dec <- modelGeneVar(mes3)
dec <- combineVar(Pet.dec, Yan.dec,mes1.dec,mes2.dec,mes3.dec)
hvgs <- dec$bio > 0 #highly variable genes

corrected <- rescaleBatches(Pet, Yan,mes1,mes2,mes3) #remove batch effect

#Update sce with normalised and batch-corrected assay
cells <- c(Pet$sample_accession,Yan$sample_accession,mes1$sample_accession,mes2$sample_accession,mes3$sample_accession)
sce <- sce[univ, cells]

assay(sce, "batch_corrected") <- assays(corrected)$corrected

# Dimensionality reduction
set.seed(100)
sce <- runPCA(sce, exprs_values="batch_corrected", subset_row=hvgs, ncomponents=14)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)

# Clustering
g <- buildSNNGraph(sce, use.dimred = 'UMAP', k=90)
colLabels(sce) <- factor(igraph::cluster_louvain(g)$membership)


# Visualization.
plotUMAP(sce, colour_by="NANOG", by_exprs_values="batch_corrected")
#plotUMAP(sce, colour_by="GATA3", by_exprs_values="batch_corrected")
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
plotUMAP(sce, colour_by="USP18", by_exprs_values="batch_corrected")

sce <- logNormCounts(sce)

morula_hesc_sce <- sce

saveRDS(morula_hesc_sce, file = "morula_hesc_sce.rds")
