suppressMessages(library(magrittr))

source(paste0(packagePath, "/modulesUI.R"), local = TRUE)

# list of menu Items
menuList <- list(
  shinydashboard::menuItem("Co-expression",
    icon = icon("dashboard"),
    # id="coexpressionID",
    tabName = "coexpression", startExpanded = FALSE,
    shinydashboard::menuSubItem("All clusters", tabName = "coexpressionAll"),
    shinydashboard::menuSubItem("Selected", tabName = "coexpressionSelected"),
    shinydashboard::menuSubItem("Violin plot", tabName = "CoExpressionViolin"),
    shinydashboard::menuSubItem("SOM cluster", tabName = "SOMcluster")
  )
)



# list of tab Items
tabList <- list(
  # coexpressionAll ----
  coexpressionAllTab = shinydashboard::tabItem(
    "coexpressionAll",
    box(
      title = "Heat map of all cells", solidHeader = TRUE, width = 12, status = "primary",
      fluidRow(
        column(
          width = 12,
          textInput("coE_heatmap_geneids", "Comma seperated gene names", value = defaultValueMultiGenes)
        )
      ),
      fluidRow(
        column(
          width = 12,
          pHeatMapUI("coExpHeatmapModule")
        )
      )
    )
  ),

  # coexpression selected -----
  coexpressionSelectedTab = shinydashboard::tabItem(
    "coexpressionSelected",
    box(
      title = "Interrogate cells", solidHeader = TRUE, width = 12, status = "primary",
      footer = tags$p("2D plot for selecting cells"),
      fluidRow(
        column(
          width = 12,
          clusterUI("coE_selected")
        )
      ),
      box(
        title = "Detailed information on genes", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = FALSE, collapsed = TRUE,
        footer = "",
        fluidRow(column(
          width = 12,
          textInput("coE_heatmapselected_geneids", "Comma seperated gene names", value = defaultValueMultiGenes)
        )),
        br(),
        # fluidRow(checkboxInput(inputId = "coE_heatmapSelectedModuleShow", label = "calc heatmap", value = FALSE)),
        box(
          title = "Heatmap of selected cells and genes", solidHeader = TRUE, width = 12, status = "primary",
          collapsible = FALSE, collapsed = TRUE,
          fluidRow(
            column(
              width = 12,
              pHeatMapUI("coE_heatmapSelectedModule")
            )
          )
        ),
        br(),
        box(
          title = "Table with coefficient of variance", solidHeader = TRUE, width = 12, status = "primary",
          collapsible = FALSE, collapsed = TRUE,
          footer = div(
            tags$ul(
              tags$li(
                "for each gene we divice the standard diviation by the mean (coefficient of variance, cv)\n"
              )
            )
          ),
          # fluidRow(checkboxInput(inputId = "coEtgMinExprShow", label = "calc min expressing genes", value = FALSE)),
          fluidRow(
            column(
              width = 3,
              numericInput("coEtgMinExpr", "min UMI count per gene:",
                1,
                min = 0, max = 100000
              )
            ), column(
              width = 3,
              numericInput("coEtgPerc", "min percentage of cells expressing a genes:",
                60,
                min = 1, max = 100
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              tableSelectionUi("coE_topExpGenes"),
            )
          )
        ),
        box(
          title = "Table with correlation coefficients", solidHeader = TRUE, width = 12, status = "primary",
          collapsible = FALSE, collapsed = TRUE,
          footer = div(
            tags$ul(
              tags$li(
                "using Hmisc we calculate the correlation coefficient and p-values for the selected genes and cells\n"
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              # checkboxInput(inputId = "coE_topCCGenesShow", label = "calc correlations", value = FALSE)),
              tableSelectionUi("coE_topCCGenes")
            )
          )
        )
      )
    )
  ),
  # CoExpressionViolin ----
  expressionTab = shinydashboard::tabItem(
    "CoExpressionViolin",
    box(
      title = "Violin plot", solidHeader = TRUE, width = 12, status = "primary",
      footer = "for each cell we count how many of the genes specified have an expression larger or equal than the minimum exprssion.\nThese counts are then divided up for any variable that can be used as a factor (has less than 20 levels).",

      fluidRow(
        column(
          width = 12,
          checkboxInput("coE_showPermutations", "show Permutations", FALSE)
        )
      ),
      fluidRow(
        column(
          width = 4,
          textInput("coE_geneGrpVioIds", "Comma seperated gene names", value = defaultValueMultiGenes)
        ),
        column(
          width = 4,
          selectInput(
            "coE_dimension_xVioiGrp",
            label = "X",
            choices = c("dbCluster", "sampleName", "tsne3"),
            selected = "dbCluster"
          )
        ),
        column(
          width = 4,
          numericInput("coEminExpr", "min expression of genes:",
            1,
            min = 1, max = 100000
          )
        )
      ),
      br(),
      fluidRow(column(width = 12,
                      # jqui_resizable(plotly::plotlyOutput("coE_geneGrp_vio_plot") )
                      jqui_resizable(plotOutput("coE_geneGrp_vio_plot") )
      )),
      br(),
      actionButton("save2HistVio", "save to history")
    )
  ),
  # SOMcluster ----
  tabList = shinydashboard::tabItem(
    "SOMcluster",
    box(
      title = "Self organizing map (SOM)", solidHeader = TRUE, width = 12, status = "primary",
      footer = "Here, we calculate a SOM on all genes using the information from all cells. Then we ask, which other genes are in the same cluster as the gene of intereset.",
      fluidRow(
        column(
          width = 12, offset = 1,
          actionButton("updateSOMParameters", "apply changes", width = "80%")
        )
      ),
      fluidRow(
        cellSelectionUI("coE_SOM_dataInput"),
        box(
          fluidRow(
            column(width = 3,
                   numericInput("coE_dimSOM", "number of nodes per dimension",
                                20,
                                min = 2, max = 100
                   )
            ), 
            column(width = 3,
                   textInput("coE_geneSOM", "Gene of interest", value = defaultValueSingleGene)
            )
          )
        ),
        
        # column(width = 3,
        #        selectInput(inputId = "coE_clusterSOM", label = "Clusters/Factor to use", 
        #                    choices = c("dbCluster", "sampleName"),
        #                    selected = "dbCluster")
        # ),
        # column(width = 3,
        #        selectInput(inputId = "coE_clusterValSOM", label = "Values to use",
        #                    choices = c("1","2"), selected = "1", multiple = TRUE)
        # )
      ),
      br(),
      fluidRow(
        column(
          width = 12,
          pHeatMapUI("coE_heatmapSOM")
        )
      ),
      br(),
      fluidRow(column(
        width = 12,
        verbatimTextOutput("coE_somGenes")
      ))
    )
  )
)
