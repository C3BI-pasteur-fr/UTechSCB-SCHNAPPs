suppressMessages(library(magrittr))

# list of menu Items
menuList <- list(
  shinydashboard::menuItem("Data Exploration",
           # id="dataExplorationID",
    tabName = "expore", icon = icon("wpexplorer"), startExpanded = FALSE,
    shinydashboard::menuSubItem("Expression", tabName = "DE_expression"),
    shinydashboard::menuSubItem("Panel plot", tabName = "DE_panelPlot")
    # shinydashboard::menuSubItem("Sorted plot", tabName = "DE_sortedPl")
  )
)

source(paste0(packagePath,  "/modulesUI.R"), local = TRUE)
# list of tab Items
tabList <- list(
  expressionTab = shinydashboard::tabItem(
    "DE_expression",
    shiny::fluidRow(div(
      htmltools::p(strong("\tInformation:")),
      htmltools::tags$ul(
        tags$li(
          strong("Clustering"),
          ":Clustering was performed with t-SNE followed by identification using DBSCAN"
        ),
        tags$li(
          strong("Cluster 0"),
          ":Cells that cannot be assigned to any cluster"
        ),
        tags$li(
          strong("3D Plot"),
          ":Enter gene name to visualize expression in a single cell"
        ),
        tags$li(
          strong("2D Plot"),
          ":Pick a cluster, highlight cells of interest to download gene expression matrix"
        )
      )
    )),
    br(),
    br(),
    fluidRow(
      column(
        6,
        fluidRow(
          column(
            4,
            textInput("DE_gene_id", "Enter gene", value = defaultValueSingleGene)
          )
        ),
        jqui_resizable(plotly::plotlyOutput("DE_tsne_plt"))
      ), column(
        6,
        clusterUI("DE_expclusters")
      )
    ),
    br(),
    fluidRow(column(
      12,
      jqui_resizable( plotOutput("DE_gene_vio_plot") %>% shinycssloaders::withSpinner())
    ))
  ),

  DE_panelPlotTab = shinydashboard::tabItem(
    "DE_panelPlot",
    tags$ul(
      tags$li(
        strong("Panel plot"),
        ":Select a cluster. Enter",
        strong("ONE"),
        "or",
        strong("MULTIPLE"),
        "gene ids to visualize expression in all clusters"
      ),
      tags$li("If the x-axis is a categorical value and the y-axis is UMI.counts the y-axis related to the count for that gene. Otherwise, all genes are used")
    ),
    fluidRow(
      column(
        2,
        uiOutput("DE_clusterSelectionPanelPlot")
      ),
      column(
        2,
        selectInput(
          "DE_dim_x",
          label = "X",
          choices = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne1"
        )
      ),
      column(
        2,
        selectInput(
          "DE_dim_y",
          label = "Y",
          choices = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne2"
        )
      ),
      column(
        2,

        textInput("DE_panelplotids", "Comma seperated gene names", value = defaultValueMultiGenes)
      )
    ),
    fluidRow(column(
      12,
      jqui_resizable(plotOutput("DE_panelPlot") )
    ))
  ),
  # DE_sortedPlot
  # shinydashboard::tabItem(
  #   "DE_sortedPl",
  #   tags$ul(
  #     tags$li(
  #       strong("Sorted plots")
  #     ),
  #     fluidRow(
  #       column(4,
  #              offset = 1,
  #              textInput("DE_geneIds", "comma separated list of genes for UmiCountPerGenes", value = "")),
  #       column(
  #       4,offset = 1,
  #       selectInput(
  #         "DE_dimension_y",
  #         label = "Y",
  #         choices = c("tsne1", "tsne2", "before.filter"),
  #         selected = "before.filter"
  #       )
  #     )
  #     ),fluidRow(
  #       column(
  #         10,
  #         offset = 1,
  #         jqui_resizable(plotly::plotlyOutput("DE_sortedPlot"))
  #         # imageOutput("DE_sortedPlot")
  #       )
  #     ),
  #     fluidRow(
  #       verbatimTextOutput("DE_SelectionText")
  #     )
  #   )
  # ),
  DE_scaterQCTab = shinydashboard::tabItem(
    "DE_scaterQC",
    tags$ul(
      tags$li(
        strong("Scater QC plots")
      ),
      checkboxInput("runScater", "run scater plot", FALSE),
      fluidRow(
        column(
          10,
          offset = 1,
          imageOutput("DE_scaterQC") %>% shinycssloaders::withSpinner() # PNG output with temp file
        )
      )
    )
  )
)
