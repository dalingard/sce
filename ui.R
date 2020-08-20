library(shiny)
library(shinyjs)

# to include new pathways just add to pth in the format: hsa[5 digit id].[name]
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

# available normalisation methods to be displayed
normTypes <- c('Batch Corrected', 'FPKM', 'TPM', 'Log Norm Counts')

# available datasets to choose from and combine
dataSets <- c("Morula", "Blastocyst", "hESCs")

# User Interface
shinyUI(fluidPage(
  theme = "styles.css",
  title = "Embryo signalling",
  useShinyjs(),
  # Tabs
  navbarPage(
    id = "navigation",
    title = "Signalling in the Human Embryo",
    # Home Page
    tabPanel(
      "Home",
      fluidRow(p(
        span("This app is designed to aid in the analysis of single-cell RNA-seq data collected from the following studies; "),
        HTML(paste0(a(href = "https://www.nature.com/articles/nsmb.2660", "Yan et al."), ", ", 
                    a(href = "https://dev.biologists.org/content/142/18/3151", "Blakeley et al. "), ", ",
                    a(href = "", "Pet et al."), " and ",
                    a(href = "https://www.sciencedirect.com/science/article/pii/S2211124718320746", "Messmer et al.")
        )),
        span("Select which datasets you would like to inspect and click load data to begin. "),
        span("When more than one set is selected they will be integrated together."), align = "center"),
               br()),
      fluidRow(
        checkboxGroupInput(
          "dataset",
          "Select dataset to analyse",
          choices = dataSets,
          selected = NULL,
          inline = TRUE
        ),
        actionButton("load_data", "Load Data"),
        align = "center"
      ),
      br(),
      fluidRow(textOutput("currently_loaded"), align = "center")
    ),
    # Pathways Page
    tabPanel(
      "Pathways",
      fluidRow(column(12,
                      p(
                        span("Genes on KEGG pathway diagrams are coloured according to their expression in the different cell types. 
                             As a result, rectangles on the KEGG diagram are coloured with the median expression of the gene across single cells of the same type. 
                             When rectangles represent more than one gene, the maximum median is represented. 
                             Clicking on a rectangle generates boxplots with the expression distribution of the gene or genes in the cell types present in the loaded data.
                             The box plots representing the gene being used to colour the group is highlighted with a red border.")
                      ),
                      br())),
      sidebarLayout(
        sidebarPanel = sidebarPanel(
          selectizeInput('spath', 'Signalling pathway:', pth),
          selectizeInput('ctype', 'Cell type:', NULL),
          selectizeInput('pathwayNtype',
                      'Normalization:',
                      normTypes),
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
    # Dimensionality Reduction Page
    tabPanel(
      "Dimensionality Reduction",
      fluidRow(p(span("This page allows you to plot the single-cell data using either PCA or UMAP. Each point on the plot represents a cell in the data.
                      Choose to colour by gene expression and select a gene to colour each cell by how much it expressed that particular gene.
                      You can select genes by clicking them in the dropdown or by typing in the next gene name and selecting it.")),
               br()),
      sidebarLayout(
        sidebarPanel = sidebarPanel(
          selectizeInput('redTech', 'Dimensionality reduction method:', c('UMAP', 'PCA')),
          selectizeInput('colourby', 'Colour by:', c("Cell type", "Batch", "Gene expression")),
          disabled(selectizeInput('ntype', 'Normalization:', normTypes)),
          disabled(selectizeInput('goi', 'Gene of interest:', NULL))
        ),
        mainPanel = mainPanel(div(
          id = "dimredPlot", plotOutput(outputId = "dimred", height = "500")
        )),
        position = "right"
      )
    ),
    # Box plot page
    tabPanel(
      "Box Plots",
      fluidRow(column(12,
                      p(
                        span("Select up to 8 genes to view their expression distribution in the cell types present in the loaded data.
                             You can select genes by clicking them in the dropdown or by typing in the next gene name and selecting it.")
                      ),
                      br())),
      sidebarLayout(
        sidebarPanel = sidebarPanel(
          selectizeInput('bp_ntype',
                      'Normalization:',
                      normTypes),
          selectizeInput(
            'bp_goi',
            'Genes of interest:',
            NULL,
            multiple = TRUE,
            options = list(maxItems = 8)
          )
        ),
        mainPanel = mainPanel(uiOutput("bp_ui_plot")),
        position = "right"
      )
    ),
    # Gene Signatures Page
    tabPanel("Gene Signatures",
             fluidRow(p(
               span("This page produces a UMAP or PCA plot that colours each cell by their average expression of genes from a particular signalling pathway.
                    To use the pathway genes please ensure the custom signature box is empty. You can choose the average used using the 'Based On' dropdown. 
                    You can instead provide a list of genes in the custom signature box to colour the cells by their average expression of those genes.
                    The list must be comma seperated but is case insensitive. Genes not present in the data will be ignored. 
                    Click visualize to see the plot.
                    ")
             ),
             br()),
             fluidRow(
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   selectizeInput(
                     'gs_dimred',
                     'Dimensionality reduction method:',
                     c('UMAP', 'PCA')
                   ),
                   selectizeInput('gs_norm', 'Normalization:', normTypes),
                   selectizeInput('gs_based_on', 'Based On:', c('Mean', 'Median')),
                   selectizeInput('gs_predefined_signatures', 'Predefined Signatures:', pth),
                   textInput('gs_custom_signature', 'Custom Signature:', value = "", placeholder = "GENE1, GENE2..."),
                   actionButton("gs_vis", "Visualize")
                 ),
                 mainPanel = mainPanel(plotOutput(outputId = "gs_plot", height = "500")),
                 position = "right"
               )
             )),
    # Differential Gene Expression page
    tabPanel(
      "Differential Gene Expression",
      fluidRow(p(span("This page will produce a table of candidate marker genes that may differentiate the clusters or cell types being compared. 
                      When 'Compare' is different to 'With' the table will be generated.")),
               br()),
      fluidRow(
        sidebarLayout(
          sidebarPanel = sidebarPanel(
            radioButtons(
              "compare_by",
              "Compare by:",
              choices = c("Cluster-based", "Cell-type based")
            ),
            selectizeInput('compare_group1', 'Compare:', NULL),
            selectizeInput('compare_group2', 'With:', NULL)
          ),
          mainPanel = mainPanel(plotOutput(outputId = "dge_Plot", height = "500")),
          position = "right"
        )
      ),
      hr(),
      fluidRow(h1(textOutput("comparison_table_title"))),
      fluidRow(dataTableOutput("comparison_table"))
    )
    
  )
))
