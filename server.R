#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(cowplot)

pth_coords <- readRDS("www/data/pathway_gene_coords.rds")
expr <- readRDS("www/data/blastocyst_flt_rpkm.rds")
col_pal <- c("#4daf4a", "#e41a1c", "#377eb8")

# Determines genes inside pathway category based on click coordinates
get_genes <- function(pth, x, y){
  pth_flt <- pth_coords %>% 
    filter(pathway == pth) %>% 
    filter((x0 <= x & y0 <= y) & (x1 >= x & y1 >= y))
  
  if(nrow(pth_flt) > 0){
    return(pth_flt$genes)
  }else{
    return("")
  }
}

# Define server logic
shinyServer(function(input, output) {
  
  values <- reactiveValues(prows = 1)
   
  # Show pre-rendered pathway diagrams
  output$pthwy_img <- renderImage({
    if (is.null(input$spath))
      return(NULL)
    
    img_feats <- list(contentType = "image/png")
    img_feats$src <- paste0("www/figs/", input$spath, "_", input$ctype, ".png")
    img_feats$alt <- paste0(input$spath, " signalling")
    
    return(img_feats)
    
  }, deleteFile = FALSE)
  
  # Computes plot height based on number of rows
  plotHeight <- function(){
    return(values$prows * 200)
  }
  
  # Generate a boxplot based on double-clicked genes
  output$expr_boxplot <- renderPlot({
    if (is.null(input$image_click))
      return(NULL)
    genes <- get_genes(input$spath, input$image_click$x, input$image_click$y)
    if(genes != ""){
      withProgress({
        genes <- strsplit(genes, ",")[[1]]
        
        values$prows <- ceiling(length(genes)/4)
        
        # Generate one boxplot per gene
        p <- list()
        for(i in seq_along(genes)){
          p[[i]] <- expr %>%
            filter(Gene == genes[i]) %>% 
            ggplot(aes(cell_type, expression, fill = cell_type)) + 
            geom_boxplot(outlier.colour = "red", outlier.shape = 8) +
            scale_fill_manual(values = col_pal) +
            labs(x = "", y = "RPKM", title = genes[i]) + theme_bw() +
            theme(legend.position = "none")
          incProgress()
        }
        plot_grid(plotlist = p, ncol = 4)
      }, message = "Making boxplots...")
      
      
    }else{
      return(NULL)
    }
  })
  
  # Wrap plotOutput in renderUI
  output$ui_plot <- renderUI({
    plotOutput("expr_boxplot", height = plotHeight())
  })
  
})
