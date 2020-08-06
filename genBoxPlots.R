generateBoxPlots <- function(toPlot, norm_type="batch_corrected", highlight=NULL){
  col_pal <- c("#4daf4a", "#e41a1c", "#377eb8", "#108783", "#b8a904", "#b84304")
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
        scale_y <- scale_y_continuous(breaks = scales::pretty_breaks(n = 10), 
                                      trans = trans, 
                                      limits = c(minVal, maxVal))
      }
      
      for (g in seq_along(rownames(toPlot))){
        expr <- assay(toPlot, norm_type)[rownames(toPlot)[[g]],]
        tidyData <- tibble(cell_type = factor(toPlot$cell_type), 
                           norm_type = expr, gene = g)
        
        if(norm_type=="fpkm" || norm_type =="tpm"){
          tidyData$norm_type <- tidyData$norm_type + 1
        }
        
        if (!is.null(highlight) && rownames(toPlot)[[g]]==highlight){
          theme <- theme(legend.position = "none", 
                         panel.border = element_rect(colour = "red", 
                                                     fill=NA, 
                                                     size=2))
        } else {
          theme <- theme(legend.position = "none")
        }
        
        p[[g]] <- ggplot(data=tidyData, 
                         mapping=aes(cell_type, norm_type, 
                                     group=cell_type, fill=cell_type)) +
          geom_boxplot(outlier.colour = "red", outlier.shape = 8) +
          scale_fill_manual(breaks = c("Epiblast", "Primitive endoderm", "Trophectoderm", "Morula", "t2iL+Go H9", "E8 H9"), 
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
