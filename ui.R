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

pth <-
  c(
    "hsa04010.MAPK",
    "hsa04012.ErbB",
    "hsa04014.Ras",
    "hsa04015.Rap1",
    "hsa04020.Calcium",
    "hsa04022.cGMP-PKG",
    "hsa04024.cAMP",
    "hsa04064.NFKB",
    "hsa04066.HIF-1",
    "hsa04068.FoxO",
    "hsa04071.Sphingolipid",
    "hsa04072.Phospholipase D",
    "hsa04150.mTOR",
    "hsa04151.PI3K-AKT",
    "hsa04152.AMPK",
    "hsa04310.WNT",
    "hsa04330.Notch",
    "hsa04340.Hedgehog",
    "hsa04350.TGF-beta",
    "hsa04370.VEGF",
    "hsa04371.Apelin",
    "hsa04390.Hippo",
    "hsa04630.JAK-STAT",
    "hsa04668.TNF"
  )

cell <- c("Epiblast", "Primitive endoderm", "Trophectoderm")

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  theme = "styles.css",
  title = "Embryo signalling",
  useShinyjs(),
  
  navbarPage(
    title = "Signalling in the human early embryo",
    
    tabPanel(
      "Pathways",
      # Application details
      fluidRow(column(12,
                      p(
                        span("Description")
                      ),
                      br())),
      sidebarLayout(
        sidebarPanel = sidebarPanel(
          selectInput('spath', 'Signalling pathway:', pth),
          selectInput('ctype', 'Cell type:', cell),
          selectInput(
            'pathwayNtype',
            'Normalization:',
            c('batch_corrected', 'fpkm', 'tpm', 'logcounts')
          ),
          uiOutput("ui_plot")
        ),
        mainPanel = mainPanel(imageOutput(
          "pathway",
          width = "100%",
          click = clickOpts("pathway_click", FALSE)
        )),
        position = "right"
      ),
      fluidRow()
    ),
    tabPanel(
      "Dimentionality Reduction",
      # Application details
      fluidRow(p(span("Description")),
               br()),
      sidebarLayout(
        sidebarPanel = sidebarPanel(
          selectInput('redTech', 'Dimensionality reduction method:', c('UMAP', 'PCA')),
          selectInput(
            'colourby',
            'Colour by:',
            c("Cell type", "Batch", "Gene expression")
          ),
          disabled(selectInput(
            'ntype',
            'Normalization:',
            c('batch_corrected', 'fpkm', 'tpm', 'logcounts')
          )),
          disabled(selectInput('goi', 'Gene of interest:', ""))
        ),
        mainPanel = mainPanel(div(
          id = "dimredPlot", plotOutput(outputId = "dimred", height = "500")
        )),
        position = "right"
      )
    ),
    tabPanel(
      "Box Plots",
      # Application details
      fluidRow(column(12,
                      p(
                        span("Description")
                      ),
                      br())),
      sidebarLayout(
        sidebarPanel = sidebarPanel(
          selectInput(
            'bp_ntype',
            'Normalization:',
            c('batch_corrected', 'fpkm', 'tpm', 'logcounts')
          ),
          selectizeInput(
            'bp_goi',
            'Genes of interest:',
            "",
            multiple = TRUE,
            options = list(maxItems = 8)
          )
        ),
        mainPanel = mainPanel(uiOutput("bp_ui_plot")),
        position = "right"
      )
    ),
    tabPanel("Gene Signatures",
             fluidRow(p(span("Description")),
                      br()),
             fluidRow(
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   selectInput('gs_dimred', 'Dimensionality reduction method:', c('UMAP', 'PCA')),
                   selectInput('gs_norm','Normalization:',c('batch_corrected', 'fpkm', 'tpm', 'logcounts')),
                   selectInput('gs_based_on','Based On:',c('Mean', 'Median')),
                   selectInput('gs_predefined_signatures','Predefined Signatures:',pth),
                   selectizeInput('gs_custom_signature','Custom Signature:',"",multiple = TRUE),
                   actionButton("gs_vis", "Visualize")
                 ),
                 mainPanel = mainPanel(plotOutput(outputId = "gs_plot", height = "500")),
                 position = "right"
               )
             )
    ),
    tabPanel(
      "Differential Gene Expression",
      fluidRow(p(span("Description")),
               br()),
      fluidRow(
        sidebarLayout(
          sidebarPanel = sidebarPanel(
            radioButtons(
              "compare_by",
              "Compare by:",
              choices = c("Cluster-based", "Cell-type based")
            ),
            selectInput('compare_group1', 'Compare:', ""),
            selectInput('compare_group2', 'With:', "")
          ),
          mainPanel = mainPanel(plotOutput(outputId = "dge_Plot", height = "500")),
          position = "right"
        )
      ),
      hr(),
      fluidRow(h1(textOutput(
        "comparison_table_title"
      ))),
      fluidRow(dataTableOutput("comparison_table")),
      hr(),
      fluidRow(h1(textOutput(
        "comparison_table_title2"
      ))),
      fluidRow(dataTableOutput("comparison_table2"))
    )
    
  )
))
