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

source(paste0(packagePath, "/modulesUI.R"), local = TRUE)
# list of tab Items
tabList <- list(
  # DE_expression ----
  expressionTab = shinydashboard::tabItem(
    "DE_expression",
    shinydashboard::box(
      title = "Expression overview", solidHeader = TRUE, width = 12, status = "primary",
      footer = div(
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
      ),
      fluidRow(div(tags$h3("Expression based on subset of cells"), align = "center")),
      tags$p("Similar to Co-expression - selection, but with a focus on subsets of cells and genes."),
      tags$p("Limit the cells to visualize by cluster or any other factorial."),
      tags$p("Visualize specific genes in 3D and over the current clustering (not changeable)."),
      tags$p("Projections will be not be effected by this sub selection."),
      br(),
      fluidRow(
        column( offset = 3,
                width = 12, 
                cellSelectionUI("DE_Exp_dataInput")
        )),
      br(),
      fluidRow(
        column( offset = 3,
                width = 12, 
                textInput("DE_gene_id", "Enter gene(s) of interest", value = defaultValue("DE_gene_id", defaultValueSingleGene))
        )),
      
      br(),
      fluidRow(
        column(
          width = 12,
          clusterUI("DE_expclusters")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 3,
          selectInput("DE_expclusters_x",
                      label = "X",
                      choices = c(defaultValue("DE_expclusters_x", "tsne1"), "tsne2", "tsne3"),
                      selected = defaultValue("DE_expclusters_x", "tsne1")
          )),
          column(
            width = 3,
            selectInput("DE_expclusters_y",
                        label = "Y",
                        choices = c(defaultValue("DE_expclusters_y", "tsne1"), "tsne2", "tsne3"),
                        selected = defaultValue("DE_expclusters_y", "tsne2")
            )),
            column(
              width = 3,
              selectInput("DE_expclusters_z",
                          label = "Z",
                          choices = c(defaultValue("DE_expclusters_z", "tsne1"), "tsne2", "tsne3"),
                          selected = defaultValue("DE_expclusters_z", "tsne3")
              ))
        # ,
        # column(
        #   width = 3,
        #   selectInput("DE_expclusters_col",
        #               label = "color",
        #               choices = defaultValue("DE_expclusters_col", "sampleName"),
        #               selected = defaultValue("DE_expclusters_col", "sampleName")
        #   ))

        )
      ),
      fluidRow(
        column(
          width = 12,
          jqui_resizable(plotly::plotlyOutput("DE_tsne_plt"))
        )
      ),
      br(),
      fluidRow(
        column(
          width = 12,
          jqui_resizable(plotOutput("DE_gene_vio_plot") %>% withSpinner())
        )
      
    )
  ),
  # DE_panelPlot ----
  DE_panelPlotTab = shinydashboard::tabItem(
    "DE_panelPlot",
    shinydashboard::box(
      title = "Panel Plot", solidHeader = TRUE, width = 12, status = "primary",
      fluidRow(
        column(
          width = 12,
          
          footer = tags$ul(
            tags$li(
              strong("Panel plot"),
              ":Select a cluster. Enter",
              strong("ONE"),
              "or",
              strong("MULTIPLE"),
              "gene ids to visualize expression in all clusters"
            ),
            tags$li("If the x-axis is a categorical value and the y-axis is UMI.counts the y-axis related to the count for that gene. Otherwise, all genes are used. Only in this case the check box 'same scale' is used.")
          )
        )
      ),
      fluidRow(
        column(
          width = 12, offset = 1,
          actionButton("updatePanelPlot", "apply changes", width = "80%")
        )
      ),
      br(),
      cellSelectionUI("DE_PanelPlotCellSelection"),
      
      # fluidRow(
      #   column(width = 3,
      #          # uiOutput("DE_clusterSelectionPanelPlot")
      #          selectInput(inputId = "DE_clusterPP", label = "Clusters/Factor to use",
      #                      choices = c("dbCluster", "sampleNames"),
      #                      selected = "dbCluster")
      #   ),
      #   column(width = 3,
      #          selectInput(inputId = "DE_PPGrp", label = "Values to use",
      #                      choices = c("1","2"), selected = "1", multiple = TRUE)
      #   )),
      shinydashboard::box(width = 6, 
          fluidRow(
            column(
              width = 6,
              selectInput("DE_dim_x",
                          label = "X",
                          choices = c(defaultValue("DE_dim_x", "tsne1"), "tsne2", "tsne3"),
                          selected = defaultValue("DE_dim_x", "tsne1")
              )
            ),
            column(
              width = 6,
              selectInput("DE_dim_y",
                          label = "Y",
                          choices = c("tsne1", defaultValue("DE_dim_y", "tsne2"), "tsne3"),
                          selected = defaultValue("DE_dim_y", "tsne2")
              )
            )
          ),
          fluidRow(
            column(
              width = 4,
              checkboxInput("DE_panelplotSameScale", "same scale", value = defaultValue("DE_panelplotSameScale", TRUE))
            ),
            column(
              width = 8,
              selectInput("DE_nCol",
                          label = "number of columns for plot",
                          choices = c(1:10),
                          selected = defaultValue("DE_nCol", 4)
              )
            ), align = 'right'
          )
      ),
      fluidRow(
        column(
          width = 12,
          textInput("DE_panelplotids", "Comma seperated gene names", value = defaultValue("DE_panelplotids", defaultValueMultiGenes))
        ),
      ),
      fluidRow(column(
        12,
        jqui_resizable(plotOutput("DE_panelPlot"))
      )),
      br(),
      actionButton("save2HistPanel", "save to history")
    )
  ),
  # DE_scaterQC ----
  DE_scaterQCTab = shinydashboard::tabItem(
    "DE_scaterQC",
    shinydashboard::box(
      title = "Quality control plot from the scater package", solidHeader = TRUE, 
      width = 12, status = "primary", height = "1627px",
      fluidRow(
        column(12, offset = 1,
               actionButton("runScater", "apply changes", width = "80%")
        )
      ),
      br(),
      fluidRow(
        column(
          12,
          offset = 0, 
          div(style = "height:672px;",
              imageOutput("DE_scaterQC") 
          ) # PNG output with temp file
        )
      ),
      br(),
      actionButton("save2HistScater", "save to history")
    )
    
    
  )
)
