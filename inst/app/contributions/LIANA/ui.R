menuList <- list(
  shinydashboard::menuItem("Liana",
                           icon = icon("dashboard"),
                           tabName = "Liana", startExpanded = FALSE,
                           shinydashboard::menuSubItem("Liana basic", tabName = "LianaBasic")
  )
)

# SOMcluster ----
tabList = list(
  shinydashboard::tabItem(
    "LianaBasic",
    box(
      title = "Ligand - receptor analysis using Liana", solidHeader = TRUE, width = 12, status = "primary",
      footer = "LIANA enables the use of any combination of ligand-receptor methods and resources, and their consensus.",
      fluidRow(
        column(
          width = 12, offset = 1,
          actionButton("updateLianaParameters", "apply changes", width = "80%")
        )
      ),
      br(),
      fluidRow(
        cellSelectionUI("Liana_dataInput"),
        box(
          fluidRow(
            column(width = 3, 
                   sc_selectizeInput(inputId = "Liana_idents_col", label = "the cell labels to be used.", 
                                     choices = c("dbCluster"), selected = defaultValue("Liana_idents_col", "dbCluster"))),
            column(width = 3, 
                   sc_selectizeInput(inputId = "Liana_resource", label = "resource(s) to be used by the methods", 
                                     choices = c("OmniPath"), selected = defaultValue("Liana_resource", "OmniPath"))),
            column(width = 3, 
                   sc_selectizeInput(inputId = "Liana_method", label = "method(s) to be run via liana",
                                     multiple = T,
                                     choices = c("natmi", "connectome", "logfc", "sca", "cellphonedb"), selected = defaultValue("Liana_method", "natmi"))),
            column(width = 3,
                   sc_numericInput("Liana_min_cells", "minimum cell per cell identity to be considered for analysis",
                                   defaultValue("Liana_min_cells", 5),
                                   min = 2, max = 10000
                   )
            )
            # parallelize
            # whether to parallelize cellphonedb-like
            # 
            # workers
            # number of workers to be called
          )
        ),
        
      ),
      br(),
      fluidRow(
        column(
          width = 10,
          tableSelectionUi("Liana_raw_TableMod")
        )
      ),
      br(),
      fluidRow(column(
        12,
        jqui_resizable(plotOutput("Liana_dotPlot")) # %>% withSpinner()
      )),
      br(),
      actionButton("save2Hist_Liana_dotPlot", "save to history"),
      br(),
      fluidRow(column(
        12,
        jqui_resizable(plotOutput("Liana_Heatmap")) # %>% withSpinner()
      )),
      br(),
      actionButton("save2Hist_Liana_Heatmap", "save to history"),
      br()
    )
    
  )
)