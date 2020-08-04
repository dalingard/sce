generateBoxPlots <- function(toPlot, norm_type="batch_corrected", highlight=NULL){
  col_pal <- c("#4daf4a", "#e41a1c", "#377eb8", "#108783")
  if(!is.null(toPlot)){
    withProgress({
      p <- list()
      
      maxVal <- max(assay(toPlot, norm_type))
      minVal <- min(assay(toPlot, norm_type))
      
      if (norm_type=="fpkm" || norm_type =="tpm"){
        trans <- "log2"
        scale_y <- list(scale_y_continuous(
          breaks = scales::trans_breaks(trans, function(x) 2^x),
          labels = scales::trans_format(trans, scales::math_format(2^.x)),
          trans = trans
        ), annotation_logticks(base = 2, sides = "l"))
      } else {
        trans <- "identity"
        scale_y <- scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(floor(minVal), ceiling(maxVal)))
      }
      
      for (g in seq_along(rownames(toPlot))){
        tidyData <- as.data.frame(assay(toPlot, norm_type))[rownames(toPlot)[[g]],]
        tidyData <- cbind(t(tidyData), cell_type=colData(toPlot)$cell_type)
        colnames(tidyData) <- c("norm_type","cell_type")
        tidyData <- as.data.frame(tidyData)
        tidyData$cell_type <- factor(tidyData$cell_type)
        tidyData$norm_type <- as.numeric(tidyData$norm_type)
        
        if (!is.null(highlight) && rownames(toPlot)[[g]]==highlight){
          theme <- theme(legend.position = "none", panel.border = element_rect(colour = "red", fill=NA, size=2))
        } else {
          theme <- theme(legend.position = "none")
        }
        
        p[[g]] <- expr %>%
          ggplot(data=tidyData, mapping=aes(cell_type, norm_type, group=cell_type, fill=cell_type)) +
          geom_boxplot(outlier.colour = "red", outlier.shape = 8) +
          scale_fill_manual(breaks = c("Epiblast", "Primitive endoderm", "Trophectoderm", "Morula"), 
                            values=col_pal) + 
          labs(x = "", y = norm_type, title = rownames(toPlot)[[g]]) + theme_bw() +
          theme +
          scale_y +
          scale_x_discrete(labels=c("Epiblast" = "Epi", "Primitive endoderm" = "PE", "Trophectoderm" = "TE", "Morula" = "Morula"))
        
        
        incProgress()
      }
      return(plot_grid(plotlist = p, ncol = 2))
    }, message = "Making boxplots...")
  }else{
    return(NULL)
  }
}
