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
library(scater)
library(scran)

load("sce_batchCorrected.RData")
source('plotDimRed.R')
source('pathways.R')

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
shinyServer(function(input, output, session) {
  
  updateSelectInput(session, "goi", choices = rownames(sce))
  
  pathData <- reactiveValues(
    data = NULL
  )
  
  observe({
    if(is.null(input$pathway_click$x)) return(NULL)
    click <- c(input$pathway_click$x, input$pathway_click$y)
    print(click)
    print(pathData$data)
    #nearest_point <- which.min(apply(data[,1:2], 1, function(a) sum(((click-a)^2))))
    #id <- data$values[nearest_point]
  })
  
  output$pathway <- renderImage({
    
    width  <- (session$clientData$output_pathway_width*0.9)
    height <- session$clientData$output_pathway_height
    
    # For high-res displays, this will be greater than 1
    pixelratio <- session$clientData$pixelratio
    
    data <- paths(input$ctype, input$spath, sce, normalistion = input$pathwayNtype)
    pathData$data <- data
    #pathData$plot.data.gene
    
    filename <- normalizePath(file.path('./', paste(substr(input$spath,1,8), '.median.', input$ctype, '.png', sep='')))
    
    # Return a list containing the filename and alt text
    list(src = filename,
         width = width,
         alt = "Pathway")
    
  }, deleteFile = FALSE)
  
  output$dimred <- renderPlot({
    if (!(input$goi == '' && input$colourby == "Gene expression")){
      reduceDimentions(sce,hvgs,input$goi,input$redTech,input$colourby,input$ntype)
    }
  })
  
  observeEvent(input$colourby, {
    if (input$colourby == "Gene expression"){
      enable("goi")
      enable("ntype")
    } else {
      disable("goi")
      disable("ntype")
    }
  })
  
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
