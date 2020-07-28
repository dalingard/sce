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
source('genBoxPlots.R')

plotHeight <- function(numPlots){
  prows <- ceiling(numPlots/2)
  return(prows*250)
}

getGenes <- function(input, values){
  values$toHighlight <- NULL
  group <- filter(values$pathdata[[1]], input$pathway_click$coords_img$x>=x-(0.5*width))
  group <- filter(group, input$pathway_click$coords_img$x<=(x)+(0.5*width))
  group <- filter(group, input$pathway_click$coords_img$y>=y-(0.5*height))
  group <- filter(group, input$pathway_click$coords_img$y<=(y)+(0.5*height))
  values$selectedGroup <- group
  return(values)
}

getCellsToPlot <- function(input, values){
  entrezCells <- sce[!is.na(rowData(sce)$entrez), ]
  group <- unlist(strsplit(values$selectedGroup$all.mapped, ","))
  group <- entrezCells[rowData(entrezCells)$entrez %in% group,]
  values$toPlot <- group
  return(values)
}
  
getGeneToHighlight <- function(input, values){
  if (length(rownames(values$toPlot))>1){
    cells <- assay(values$toPlot, input$pathwayNtype)[,!colData(values$toPlot)$cell_type!=input$ctype]
    gene.data <- rowMedians(cells)
    names(gene.data) <- rownames(cells)
    gene.data <- sort(gene.data, TRUE)
    values$toHighlight <- names(gene.data)[[1]]
    return(values)
  }
}

# Define server logic
shinyServer(function(input, output, session) {
  
  #provide options to selectize inputs
  updateSelectInput(session, "goi", choices = rownames(sce))
  updateSelectInput(session, "bp_goi", choices = rownames(sce))
  
  values <- reactiveValues(
    pathdata = NULL,
    normalisation = NULL,
    selectedGroup = NULL,
    toPlot = NULL,
    genesToPlot = NULL,
    toHighlight = NULL
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
  observeEvent(input$pathway_click, {
    values <- getGenes(input, values)
  })
  
  #plot pathway
  output$pathway <- renderImage({
    width  <- session$clientData$output_pathway_width*0.9
    height <- session$clientData$output_pathway_height
    pixelratio <- session$clientData$pixelratio
    
    data <- paths(input$ctype, input$spath, sce, normalistion = input$pathwayNtype)
    values$pathdata <- data
    
    filename <- normalizePath(file.path('./', paste(substr(input$spath,1,8), '.median.', input$ctype, '.png', sep='')))
    list(src = filename,
         width = width,
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
    generateBoxPlots(values$toPlot, input$pathwayNtype, values$toHighlight)
  })
  
  #ui wrapper for boxplots on pathway tab
  output$ui_plot <- renderUI({
    if (!is.null(values$selectedGroup)){
      numPlots <- length(unlist(strsplit(values$selectedGroup$all.mapped,",")))
      if (numPlots > 0){
        values <- getCellsToPlot(input, values)
        values <- getGeneToHighlight(input, values)
        plotOutput("expr_boxplot", height=plotHeight(numPlots))
      }
    }
  })
  
  #plot box plots on boxplots tab
  output$g_boxplot <- renderPlot({
    generateBoxPlots(values$genesToPlot, input$bp_ntype)
  })
  
  #ui wrapper for boxplots on boxplots tab
  output$bp_ui_plot <- renderUI({
    numPlots <- length(input$bp_goi)
    if (numPlots > 0){
      values$genesToPlot <- sce[rowData(sce)$symbol %in% input$bp_goi,]
      plotOutput("g_boxplot", height=plotHeight(numPlots))
    }
  })
  
})
