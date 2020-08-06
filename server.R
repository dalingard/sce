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
library(pathview)

source('plotDimRed.R')
source('pathways.R')
source('genBoxPlots.R')

#box plot height
plotHeight <- function(numPlots){
  prows <- ceiling(numPlots/2)
  return(prows*250)
}

#get clicked group
getGenes <- function(input, values){
  values$toHighlight <- NULL
  group <- filter(values$pathdata[[1]], input$pathway_click$coords_img$x>=x-(0.5*width)
                  & input$pathway_click$coords_img$x<=(x)+(0.5*width)
                  & input$pathway_click$coords_img$y>=y-(0.5*height)
                  & input$pathway_click$coords_img$y<=(y)+(0.5*height))
  values$selectedGroup <- group
  return(values)
}

#get genes from the clicked group
getGenesToPlot <- function(all.mapped, entrezCells){
  group <- unlist(strsplit(all.mapped, ","))
  group <- entrezCells[rowData(entrezCells)$entrez %in% group,]
  return(group)
}

#gets the gene with the highest median expression in a particular cell type group
getGeneToHighlight <- function(input, toPlot){
  if (length(rownames(toPlot))>1){
    cells <- assay(toPlot, convertNames(input$pathwayNtype))[,!colData(toPlot)$cell_type!=input$ctype]
    gene.data <- rowMedians(cells)
    names(gene.data) <- rownames(cells)
    gene.data <- sort(gene.data, TRUE)
    toHighlight <- names(gene.data)[[1]]
    return(toHighlight)
  }
}

#convert norm type labels to those used in sce
convertNames <- function(norm_type){
  if (norm_type == "Batch Corrected"){
    norm_type <- "batch_corrected"
  } else if (norm_type == "FPKM" | norm_type == "TPM"){
    norm_type <- tolower(norm_type)
  } else {
    norm_type <- "logcounts"
  }
  return(norm_type)
}

# Define server logic
shinyServer(function(input, output, session) {
  
  values <- reactiveValues(
    dataset = NULL,
    currently_loaded = c("None"),
    pathdata = NULL,
    normalisation = NULL,
    selectedGroup = NULL,
    toPlot = NULL,
    genesToPlot = NULL,
    toHighlight = NULL,
    markers = NULL,
    gs_toPlot = NULL
  )
  
  hideTab("navigation", "Pathways")
  hideTab("navigation", "Dimensionality Reduction")
  hideTab("navigation", "Box Plots")
  hideTab("navigation", "Gene Signatures")
  hideTab("navigation", "Differential Gene Expression")
  
  observeEvent(input$load_data, {
    withProgress({
    if (!is.null(input$dataset)){
      if ("Early Blastocyst" %in% input$dataset & "Morula" %in% input$dataset & "hESCs" %in% input$dataset){
        values$dataset <- readRDS(file = "sce/early_blastocyst_morula_hesc_sce.rds")
        values$currently_loaded <- c("Early Blastocyst", "Morula", "hESCs")
      } else if ("Early Blastocyst" %in% input$dataset & "Morula" %in% input$dataset){
        values$dataset <- readRDS(file = "sce/early_blastocyst_morula_sce.rds")
        values$currently_loaded <- c("Early Blastocyst", "Morula")
      } else if ("Early Blastocyst" %in% input$dataset & "hESCs" %in% input$dataset){
        values$dataset <- readRDS(file = "sce/early_blastocyst_hesc_sce.rds")
        values$currently_loaded <- c("Early Blastocyst", "hESCs")
      } else if ("Morula" %in% input$dataset & "hESCs" %in% input$dataset){
        values$dataset <- readRDS(file = "sce/morula_hesc_sce.rds")
        values$currently_loaded <- c("Morula", "hESCs")
      } else if ("Early Blastocyst" %in% input$dataset){
        values$dataset <- readRDS(file = "sce/early_blastocyst_sce.rds")
        values$currently_loaded <- c("Early Blastocyst")
      } else if ("Morula" %in% input$dataset){
        values$dataset <- readRDS(file = "sce/morula_sce.rds")
        values$currently_loaded <- c("Morula")
      } else if ("hESCs" %in% input$dataset){
        values$dataset <- readRDS(file = "sce/hesc_sce.rds")
        values$currently_loaded <- c("hESCs")
      }
      updateSelectizeInput(session, "ctype", choices = unique(values$dataset$cell_type), server = TRUE)
      incProgress()
      updateSelectizeInput(session, "goi", choices = rownames(values$dataset), server = TRUE)
      incProgress()
      updateSelectizeInput(session, "bp_goi", choices = rownames(values$dataset), server = TRUE)
      incProgress()
      updateSelectizeInput(session, "gs_custom_signature", choices = rownames(values$dataset), server = TRUE)
      incProgress()
      showTab("navigation", "Pathways")
      showTab("navigation", "Dimensionality Reduction")
      showTab("navigation", "Box Plots")
      showTab("navigation", "Gene Signatures")
      showTab("navigation", "Differential Gene Expression")
    }
    }, message = "Loading Data")
  })
  
  output$currently_loaded <- renderText({
    loaded <- paste(values$currently_loaded, collapse=', ' )
    paste("Currently Loaded: ",loaded)
  })
  
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
  
  #On tab change
  # observeEvent(input$navigation, {
  # 
  # })
  
  #on a change to either of the "compare with" inputs on diff gene expr page, fetch both sets of markers
  observeEvent({
    input$compare_group1
    input$compare_group2
  }, {
    if (input$compare_group1 != input$compare_group2){
      if (input$compare_by=="Cluster-based"){
        cells <- values$dataset[,(colData(values$dataset)$label == input$compare_group1 | colData(values$dataset)$label == input$compare_group2)]
        g <- factor(colData(cells)$label, levels=unique(colData(cells)$label))
      } else {
        cells <- values$dataset[,(colData(values$dataset)$cell_type == input$compare_group1 | colData(values$dataset)$cell_type == input$compare_group2)]
        g<- as.factor(colData(cells)$cell_type)
      }
      
      # For multiple combined tests
      # markers <- multiMarkerStats(
      #   t=findMarkers(cells, groups=g, direction="up", assay.type="batch_corrected"),
      #   wilcox=findMarkers(cells, groups=g, test="wilcox", direction="up", assay.type="batch_corrected"),
      #   binom=findMarkers(cells, groups=g, test="binom", direction="up", assay.type="batch_corrected")
      # )
      
      markers <- findMarkers(cells, groups = g, pval.type="any") 
      
      markers[[1]] <- cbind(names = rownames(markers[[1]]), markers[[1]])
      markers[[2]] <- cbind(names = rownames(markers[[2]]), markers[[2]])
      values$markers <- markers
    } else {
      values$markers <- NULL
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
    
    data <- paths(input$ctype, input$spath, values$dataset, norm_type = convertNames(input$pathwayNtype))
    values$pathdata <- data
    
    filename <- normalizePath(file.path('./', paste(substr(input$spath,1,8), '.median.', input$ctype, '.png', sep='')))
    list(src = filename,
         width = width,
         alt = "Pathway Unavailable")
  }, deleteFile = FALSE)
  
  
  #plot PCA/UMAP
  output$dimred <- renderPlot({
    if (!(input$goi == '' && input$colourby == "Gene expression")){
      reduceDimentions(values$dataset,input$goi,input$redTech,input$colourby,convertNames(input$ntype))
    }
  })
  
  #plot box plots on pathway tab
  output$expr_boxplot <- renderPlot({
    generateBoxPlots(values$toPlot, convertNames(input$pathwayNtype), values$toHighlight)
  })
  
  #ui wrapper for boxplots on pathway tab
  output$ui_plot <- renderUI({
    if (!is.null(values$selectedGroup)){
      numPlots <- length(unlist(strsplit(values$selectedGroup$all.mapped,",")))
      if (numPlots > 0){
        entrezCells <- values$dataset[!is.na(rowData(values$dataset)$entrez), ]
        values$toPlot <- getGenesToPlot(values$selectedGroup$all.mapped, entrezCells)
        values$toHighlight <- getGeneToHighlight(input, values$toPlot)
        plotOutput("expr_boxplot", height=plotHeight(numPlots))
      }
    }
  })
  
  #plot box plots on boxplots tab
  output$g_boxplot <- renderPlot({
    generateBoxPlots(values$genesToPlot, convertNames(input$bp_ntype))
  })
  
  #ui wrapper for boxplots on boxplots tab
  output$bp_ui_plot <- renderUI({
    numPlots <- length(input$bp_goi)
    if (numPlots > 0){
      values$genesToPlot <- values$dataset[rowData(values$dataset)$symbol %in% input$bp_goi,]
      plotOutput("g_boxplot", height=plotHeight(numPlots))
    }
  })
  
  #diff gene expr umap plot
  output$dge_Plot <- renderPlot({
    if (input$compare_by=="Cluster-based"){
      cs <- colData(values$dataset)$label
    } else {
      cs <- colData(values$dataset)$cell_type
    }
    cs <- sort(cs)
    updateSelectInput(session, "compare_group1", choices = cs)
    updateSelectInput(session, "compare_group2", choices = cs, selected = cs[[2]])
    reduceDimentions(values$dataset,"","UMAP",input$compare_by)
  })
  
  #diff gene expr gene table title
  output$comparison_table_title <- renderText({ 
    if (!is.null(values$markers)){
      "Potential Markers"
    }
  })
  
  #diff gene expr gene table 
  output$comparison_table <- renderDataTable({
    if (!is.null(values$markers)){
      if (names(values$markers)[[1]] == input$compare_group1){
        values$markers[[1]]
      } else {
        values$markers[[2]]
      }
    }
  }, options=list(pageLength=10))
  
  #on clicking visualize get all genes in pathway or custom signature
  observeEvent(input$gs_vis, {
    withProgress({
    genes <- vector()
    if (!is.null(input$gs_custom_signature)){
      genes <- input$gs_custom_signature
    } else {
      cells <- assay(values$dataset, convertNames(input$gs_norm))
      if (input$gs_based_on=="Mean"){
        gene.data <- rowMeans(cells)
      } else {
        gene.data <- rowMedians(cells)
        names(gene.data) <- rownames(cells)
      }
      pathway <- substr(input$gs_predefined_signatures, 4, 8)
      pv.out <- pathview(gene.data = gene.data, pathway.id = pathway, gene.idtype = "symbol", limit=ceiling(max(gene.data)), both.dirs = F, high = "purple", mid = "white",
                         out.suffix = paste(input$gs_based_on,".","all", sep=""), node.sum = "max", kegg.dir = "./pathway_data")
      
      entrezCells <- values$dataset[!is.na(rowData(values$dataset)$entrez), ]
      for (g in 1:length(pv.out$plot.data.gene$all.mapped)){
        group <- pv.out$plot.data.gene[g,]
        newGenes <- getGenesToPlot(group$all.mapped, entrezCells)
        genes <- union(genes,rownames(newGenes))
        incProgress()
      }
      
    }
    
    cells <- values$dataset[rowData(values$dataset)$symbol %in% genes,]
    cells <- assay(cells, convertNames(input$gs_norm))
    
    if (input$gs_based_on=="Mean"){
      average.gene.expression <- colMeans(cells)
    } else {
      average.gene.expression <- colMedians(cells)
      names(average.gene.expression) <- colnames(cells)
    }
    cells <- values$dataset
    colData(cells) <- cbind(colData(cells), average.gene.expression)
    values$gs_toPlot <- cells
    }, message='Generating Plot')
  })
  
  #plot PCA/UMAP on gene signature page
  output$gs_plot <- renderPlot({
    if (!is.null(values$gs_toPlot)){
      reduceDimentions(values$gs_toPlot,"",input$gs_dimred,"average.gene.expression",convertNames(input$gs_norm))
    }
  })
})
