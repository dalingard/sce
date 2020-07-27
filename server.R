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
library(png)

load("sce_batchCorrected.RData")
source('plotDimRed.R')
source('pathways.R')
source('genBoxPlots.R')

plotHeight <- function(numPlots){
  prows <- ceiling(numPlots/2)
  return(prows*250)
}

# Define server logic
shinyServer(function(input, output, session) {
  
  #provide options to selectize inputs
  updateSelectInput(session, "goi", choices = rownames(sce))
  updateSelectInput(session, "bp_goi", choices = rownames(sce))
  
  values <- reactiveValues(
    pathdata = NULL,
    normalisation = NULL,
    xsf = 1,
    ysf = 1,
    toPlot = NULL
  )
  
  #deactivate/activate inputs
  observeEvent(input$colourby, {
    if (input$colourby == "Gene expression"){
      enable("goi")
      enable("ntype")
    } else {
      disable("goi")
      disable("ntype")
    }
  })
  
  #on click of pathway get gene data to plot
  observeEvent(input$pathway_click,{
    values$xsf <- 1
    group <- filter(values$pathdata[[1]], input$pathway_click$x>=x*values$xsf-(0.5*width*values$xsf))
    group <- filter(group, input$pathway_click$x<=(x*values$xsf)+(0.5*width*values$xsf))
    group <- filter(group, input$pathway_click$y>=y*values$ysf-(0.5*height*values$ysf))
    group <- filter(group, input$pathway_click$y<=(y*values$ysf)+(0.5*height*values$ysf))
    entrezCells <- sce[!is.na(rowData(sce)$entrez), ]
    group <- unlist(strsplit(group$all.mapped, ","))
    group <- entrezCells[rowData(entrezCells)$entrez %in% group,]
    values$toPlot <- group
  })
  
  
  #plot pathway
  output$pathway <- renderImage({
    
    width  <- session$clientData$output_pathway_width*0.9
    height <- session$clientData$output_pathway_height
    pixelratio <- session$clientData$pixelratio
    
    data <- paths(input$ctype, input$spath, sce, normalistion = input$pathwayNtype)
    values$pathdata <- data
    
    filename <- normalizePath(file.path('./', paste(substr(input$spath,1,8), '.median.', input$ctype, '.png', sep='')))
    img <- readPNG(filename)
    size <- dim(img)
    values$xsf <- width/size[[1]]
    # Return a list containing the filename and alt text
    list(src = filename,
         #width = width,
         alt = "Pathway")
  }, deleteFile = FALSE)
  
  
  #plot PCA/UMAP
  output$dimred <- renderPlot({
    if (!(input$goi == '' && input$colourby == "Gene expression")){
      reduceDimentions(sce,hvgs,input$goi,input$redTech,input$colourby,input$ntype)
    }
  })
  
  #plot box plots on pathway tab
  output$expr_boxplot <- renderPlot({
    generateBoxPlots(values$toPlot, input$pathwayNtype)
  })
  
  #ui wrapper for boxplots on pathway tab
  output$ui_plot <- renderUI({
    numPlots <- length(rownames(values$toPlot))
    if (numPlots > 0){
      plotOutput("expr_boxplot", height=plotHeight(numPlots))
    }
  })
  
  #plot box plots on boxplots tab
  output$g_boxplot <- renderPlot({
    generateBoxPlots(values$toPlot, input$bp_ntype)
  })
  
  #ui wrapper for boxplots on boxplots tab
  output$bp_ui_plot <- renderUI({
    numPlots <- length(input$bp_goi)
    if (numPlots > 0){
      values$toPlot <- sce[rowData(sce)$symbol %in% input$bp_goi,]
      plotOutput("g_boxplot", height=plotHeight(numPlots))
    }
  })
  
})
