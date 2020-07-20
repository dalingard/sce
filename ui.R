#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

pth <- c("hsa04010.MAPK", "hsa04012.ErbB", "hsa04014.Ras", "hsa04015.Rap1",
         "hsa04020.Calcium", "hsa04022.cGMP-PKG", "hsa04024.cAMP", "hsa04064.NFKB", 
         "hsa04066.HIF-1", "hsa04068.FoxO", "hsa04071.Sphingolipid", 
         "hsa04072.Phospholipase D", "hsa04150.mTOR", "hsa04151.PI3K-AKT", 
         "hsa04152.AMPK", "hsa04310.WNT", "hsa04330.Notch", "hsa04340.Hedgehog", 
         "hsa04350.TGF-beta", "hsa04370.VEGF", "hsa04371.Apelin", 
         "hsa04390.Hippo", "hsa04630.JAK-STAT", "hsa04668.TNF")
cell <- c("Epi", "PE", "TE")

# Define UI for application that draws a histogram
shinyUI(fluidPage(theme = "styles.css", title = "Embryo signalling",
  
  # Application title
  titlePanel("Signalling in the human early embryo"),
  
  # Application details
  fluidRow(column(12, 
                  p(span("This Shiny App colours genes on KEGG pathway "),
                    span("diagrams according to their expression in the "),
                    span("different cell types of the human blastocyst. "),
                    span("Signal transduction pathways can be chosen from "),
                    span("the dropdown menu on the left and blastocyst cell "),
                    span("types from the one on the right. Epi corresponds "),
                    span("to the epiblast, PE to the primitive endoderm "),
                    span("and TE to the trophectoderm. Gene expression is "),
                    HTML(paste0(span("shown as a colour range in log"), 
                                span("2", style = "vertical-align: sub;font-size: smaller;"), 
                                span("(RPKM +1) units "))),
                    span("and corresponds to single-cell RNA-seq data from "),
                    HTML(paste0(a(href = "https://www.nature.com/articles/nsmb.2660", "Yan et al."), " and ", 
                                a(href = "https://dev.biologists.org/content/142/18/3151", "Blakeley et al. "))),
                    span("As a result, rectangles on the KEGG diagram are "),
                    span("coloured with the median expression of the gene "),
                    span("across single cells of the same type. When rectangles "),
                    span("represent more than one gene, the maximum median is "),
                    span("represented. Clicking on a rectangle generates "),
                    span("boxplots with the expression distribution of the gene "),
                    span("or genes in the three blastocyst cell types.")
                    ), 
                  br())),
  
  fluidRow(
    column(12, wellPanel(
      fluidRow(
        column(6, selectInput('spath', 'Signalling pathway:', pth)),
        column(6, selectInput('ctype', 'Cell type:', cell))
      )
    ))
  ),
  fluidRow(
    column(12, div(id = "pthwy", imageOutput("pthwy_img", click = clickOpts("image_click", clip = FALSE))))
  ),
  fluidRow(
    column(12, br(), br())
  ),
  fluidRow(
    column(12, div(id = "bplot", uiOutput("ui_plot")))
  )
))
