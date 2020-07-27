generateBoxPlots <- function(toPlot, norm="batch_corrected"){
  col_pal <- c("#4daf4a", "#e41a1c", "#377eb8")
  if(!is.null(toPlot)){
    withProgress({
      p <- list()
      for (g in seq_along(rownames(toPlot))){
        tidyData <- as.data.frame(assay(toPlot, norm))[rownames(toPlot)[[g]],]
        tidyData <- cbind(t(tidyData), cell_type=colData(toPlot)$cell_type)
        colnames(tidyData) <- c("norm","cell_type")
        tidyData <- as.data.frame(tidyData)
        tidyData$cell_type <- factor(tidyData$cell_type)
        tidyData$norm <- as.numeric(tidyData$norm)
        
        p[[g]] <- expr %>%
          ggplot(data=tidyData, mapping=aes(cell_type, norm, group=cell_type, fill=cell_type)) +
          geom_boxplot(outlier.colour = "red", outlier.shape = 8) +
          scale_fill_manual(values = col_pal) +
          labs(x = "", y = norm, title = rownames(toPlot)[[g]]) + theme_bw() +
          theme(legend.position = "none") +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
        
        incProgress()
      }
      return(plot_grid(plotlist = p, ncol = 2))
    }, message = "Making boxplots...")
  }else{
    return(NULL)
  }
}
