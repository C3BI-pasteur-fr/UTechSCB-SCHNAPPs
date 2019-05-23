library(magrittr)

source(paste0(packagePath,  "/modulesUI.R"))

# list of menu Items
menuList <- list(
  shinydashboard::menuItem("Co-expression", icon = icon("dashboard"),
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
  coexpressionAllTab = shinydashboard::tabItem(
    "coexpressionAll",
    fluidRow(
      column(
        7,
        offset = 1,

        textInput("coE_heatmap_geneids", "Comma seperated gene names", value = defaultValueMultiGenes)
      )
    ),
    fluidRow(column(
      10,
      offset = 1,
      pHeatMapUI("coExpHeatmapModule") %>% shinycssloaders::withSpinner()
    ))
  ),

  coexpressionSelectedTab = shinydashboard::tabItem(
    "coexpressionSelected",
    tags$p(
        strong("2D plot for selecting cells"),
        "some explaination.... (Vale?)"
    ),
    fluidRow(
      # column(
      #   2,
      #   textInput("coE_gene_id_sch", "Enter gene", value = defaultValueSingleGene)
      # ),
      column(
        12,
        clusterUI("coE_selected")
      )
    ),
    fluidRow(column(
      6,
      offset = 1,

      textInput("coE_heatmapselected_geneids", "Comma seperated gene names", value = defaultValueMultiGenes)
    )),
    # fluidRow(column(
    #   12,
    #   offset = 1,
    #   uiOutput("coE_heatmapNull")
    # )),
    fluidRow(column(
      12,
      offset = 0,
      pHeatMapUI("coE_heatmapSelectedModule")
    )),
    # br(),
    fluidRow(
      column(
        3,
        numericInput("coEtgMinExpr", "min UMI count per gene:",
          1,
          min = 0, max = 100000
        )
      ), column(
        3,
        numericInput("coEtgPerc", "min percentage of cells expressing a genes:",
          60,
          min = 1, max = 100
        )
      )
    ),
    tableSelectionUi("coE_topExpGenes") %>% shinycssloaders::withSpinner()
  ),
  expressionTab = shinydashboard::tabItem(
    "CoExpressionViolin",
    fluidRow(div(
      p(strong("\tInformation:")),
      tags$ul(
        tags$li(
          strong("Violin plot"),
          "for each cell we count how many of the genes specified have an expression larger or equal than the minimum exprssion.\nThese counts are then divided up for any variable that can be used as a factor (has less than 20 levels)."
        )
      )
    )),
    br(),
    br(),
    shinyBS::tipify(
      checkboxInput("coE_showPermutations", "show Permutations", FALSE),
      "check this if you are working on the cell/gene selection to avoid certain calculations"
    ),
    fluidRow(
      column(
        3,
       textInput("coE_geneGrpVioIds", "Comma seperated gene names", value = defaultValueMultiGenes)
      ),
      column(
        3,
        selectInput(
          "coE_dimension_xVioiGrp",
          label = "X",
          choices = c("dbCluster", "sampleName", "tsne3"),
          selected = "dbCluster"
        )
      ),
      column(
        3,
        numericInput("coEminExpr", "min expression of genes:",
          1,
          min = 1, max = 100000
        )
      )
    ),
    br(),
    fluidRow(column(
      12,
      # jqui_resizable(plotly::plotlyOutput("coE_geneGrp_vio_plot") )
      jqui_resizable(plotOutput("coE_geneGrp_vio_plot") )
    ))
  ),
  tabList = shinydashboard::tabItem(
    "SOMcluster",
    fluidRow(div(p(strong("\t Self organizing map and heatmap of cluster with gene of interest")))),
    br(),
    fluidRow(column(
      3,
      numericInput("coE_dimSOM", "number of nodes per dimension",
        20,
        min = 2, max = 100
      ),
      textInput("coE_geneSOM", "Gene of interest", value = defaultValueSingleGene)
    )),
    br(),
    fluidRow(column(
      10,
      offset = 1,
      pHeatMapUI("coE_heatmapSOM") %>% shinycssloaders::withSpinner()
    )),
    br(),
    fluidRow(column(
      10,
      offset = 1,
      verbatimTextOutput("coE_somGenes")
    ))
  )
)
