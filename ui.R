#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)

pth <- c("hsa04010.MAPK", "hsa04012.ErbB", "hsa04014.Ras", "hsa04015.Rap1",
         "hsa04020.Calcium", "hsa04022.cGMP-PKG", "hsa04024.cAMP", "hsa04064.NFKB", 
         "hsa04066.HIF-1", "hsa04068.FoxO", "hsa04071.Sphingolipid", 
         "hsa04072.Phospholipase D", "hsa04150.mTOR", "hsa04151.PI3K-AKT", 
         "hsa04152.AMPK", "hsa04310.WNT", "hsa04330.Notch", "hsa04340.Hedgehog", 
         "hsa04350.TGF-beta", "hsa04370.VEGF", "hsa04371.Apelin", 
         "hsa04390.Hippo", "hsa04630.JAK-STAT", "hsa04668.TNF")

cell <- c("Epiblast", "Primitive endoderm", "Trophectoderm")
#genes <- rownames(sce)

# Define UI for application that draws a histogram
shinyUI(fluidPage(theme = "styles.css", title = "Embryo signalling", useShinyjs(),
                  
  navbarPage(title = "Signalling in the human early embryo",
    
    tabPanel("Pathways", 
             # Application details
             fluidRow(column(12, 
                             p(span("Description")
                             ), 
                             br())),
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 selectInput('spath', 'Signalling pathway:', pth),
                 selectInput('ctype', 'Cell type:', cell),
                 selectInput('pathwayNtype', 'Normalization:', c('batch_corrected','fpkm','tpm','logcounts'))
               ),
               mainPanel = mainPanel(plotOutput("pathway", click = "pathway_click")),
               position = "right"
             )
             ),
    tabPanel("Dimentionality Reduction", # Application details
             fluidRow(p(span("Description")), 
                             br()),
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 selectInput('redTech', 'Dimensionality reduction method:', c('UMAP', 'PCA')),
                 selectInput('colourby', 'Colour by:', c("Cell type", "Batch", "Gene expression")),
                 disabled(selectInput('ntype', 'Normalization:', c('batch_corrected','fpkm','tpm','logcounts'))),
                 disabled(
                 selectInput('goi', 'Gene of interest:',"")
                 )
               ),
               mainPanel = mainPanel(
                 div(id = "dimredPlot", plotOutput(outputId = "dimred", height = "800"))
               ),
               position = "right"
             )
    ) 
  )
))
