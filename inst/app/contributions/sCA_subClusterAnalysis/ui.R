suppressMessages(library(magrittr))
source(paste0(packagePath,  "/modulesUI.R"), local = TRUE)

menuList <- list(
  shinydashboard::menuItem("Subcluster analysis", 
                           icon = icon("bar-chart-o"),
                           # id = 'subclusterID',
                           tabName = "subcluster", startExpanded = FALSE,
                           shinydashboard::menuSubItem("DGE analysis", tabName = "sCA_dge")
  )
)
# sCA_dge ----
tabList <- list(
  dgeTab = shinydashboard::tabItem(
    "sCA_dge",
    box(
      title = "Differential expression", solidHeader = TRUE, width = 12, status = 'primary', 
      footer = {
        tags$ul(
          tags$li(
            paste(
              ":Select a group of cells in plot1 and a different group of cells in plot2",
              "for identifying differential features between these subclusters"
            )
          ),
          tags$li(
            strong("colors:"),
            paste("colored by cluster identity")
          ),
          tags$li(
            strong("selection hint:"),
            paste("you can also slect by groups you have defined in other plots.")
          ),
          tags$li(
            strong("selection hint:"),
            paste('also check out "Gene.count" to verify that number genes per cell.')
          )
        )}
      ,
      
      fluidRow(
        column(width = 4,
          uiOutput("sCA_dgeClustersSelection")
        ),
        column(width = 4,
          selectInput(
            "sCA_subscluster_x1",
            label = "X",
            choice = c("tsne1", "tsne2", "tsne3"),
            selected = "tsne1"
          )
        ),
        column(width = 4,
          selectInput(
            "sCA_subscluster_y1",
            label = "Y",
            choice = c("tsne1", "tsne2", "tsne3"),
            selected = "tsne2"
          )
        )
      ),
      br(),
      fluidRow(
        column(width = 6,
          plotOutput("sCA_dge_plot1", brush = brushOpts(
            id = "db1"
          )) %>% withSpinner()
        ),
        column(width = 6,
          plotOutput("sCA_dge_plot2", brush = brushOpts(
            id = "db2"
          )) %>% withSpinner()
        )
      )
    ),
    # shinydashboard::tabItem(
    #   "diffExpMethod",
    #   list(
    #     tags$h3("Method to use for differential gene expression analysis"),
    box(
      title = "DGE method", solidHeader = TRUE, width = 12, status = 'primary', 
      collapsible = TRUE, collapsed = FALSE,
      fluidRow(
        column(width = 10, offset = 1,
               radioButtons(
                 inputId = "sCA_dgeRadioButton",
                 label = "Method to use",
                 choices = "dgeChoices",
                 selected = "DE_logNormalization",
                 width = "100%"
               )
        )
      )
    ),
    
    box(
      title = "Volcano plot", solidHeader = TRUE, width = 12, status = 'primary', 
      collapsible = FALSE, collapsed = FALSE,
      fluidRow(
        column(width = 12,
               verbatimTextOutput("sCA_volc_selected")
        )
      ),
      fluidRow(
        column(width = 12,
               jqui_resizable(plotly::plotlyOutput("sCA_volcanoPlot"))
        )
      )
    ),
    box(
      title = "Differentially Expressed Genes", solidHeader = TRUE, width = 12, status = 'primary', 
      collapsible = FALSE, collapsed = TRUE,
      fluidRow(
        column(12,
               tableSelectionUi("sCA_dgeTable")
        )
      )
    )
  )
)
