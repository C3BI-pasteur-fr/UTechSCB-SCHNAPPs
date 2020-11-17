suppressMessages(library(magrittr))
source(paste0(packagePath, "/modulesUI.R"), local = TRUE)

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
    tabBox(
      title = "Differential expression", width = 12, 
      tabPanel(title = "Select data", 
               width = 12, id = "modProj",
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
               cellSelectionUI("sCA_dataInput"),
               box(width = 6,
                   title = "Dimensions for plot",
                   fluidRow(
                     column(width = 6,
                            selectInput(
                              "sCA_subscluster_x1",
                              label = "X",
                              choices = c(defaultValue("sCA_subscluster_x1", "tsne1"), "tsne2", "tsne3"),
                              selected = defaultValue("sCA_subscluster_x1", "tsne1")
                            )
                     ),
                     column(width = 6,
                            selectInput(
                              "sCA_subscluster_y1",
                              label = "Y",
                              choices = c("tsne1", defaultValue("sCA_subscluster_y1", "tsne2"), "tsne3"),
                              selected = defaultValue("sCA_subscluster_y1", "tsne2")
                            )
                     )
                   )
               ),
               fluidRow(
                 column(
                   width = 10, offset = 1,
                   selectizeInput(
                     inputId = "sCA_dgeRadioButton",
                     label = "Method to use",
                     choices = defaultValue("sCA_dgeRadioButton", "DE_logNormalization"),
                     selected = defaultValue("sCA_dgeRadioButton", "DE_logNormalization"),
                     width = "100%"
                   )
                 )
               ),
               fluidRow(
                 column(width = 12, offset = 1,
                        actionButton("updateDGEParameters", "apply changes", width = "80%")
                 )
               ),
               
               
               br(),
               box(width = 12,
                   fluidRow(
                     column(width = 6,
                            # plotly::plotlyOutput("sCA_dge_plot1")
                            plotOutput("sCA_dge_plot1", brush = brushOpts(
                              id = "db1"
                            )) %>% withSpinner()
                     )
                     ,
                     column(width = 6,
                            # plotly::plotlyOutput("sCA_dge_plot2")
                            plotOutput("sCA_dge_plot2", brush = brushOpts(
                              id = "db2"
                            )) %>% withSpinner()
                     )
                   )
               ),
      ),
      # ,
      # box(
      #   title = "DGE method", solidHeader = TRUE, width = 12, status = "primary",
      #   collapsible = TRUE, collapsed = FALSE,
      #   
      # )
      
      br(),
      # box(width = 12,
      #     ),
      # br(),
      
      if ("manhattanly" %in% rownames(installed.packages()))
      tabPanel(
        title = "Volcano plot",  width = 12, 
        collapsible = FALSE, collapsed = FALSE,
        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = "sCA_volc_effectLimit", label = "x-axis threshold", value = defaultValue("sCA_volc_effectLimit", 1), min = 0.0, max = 10000,
              step = 0.1
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = "sCA_volc_pval", label = "y-axis threshold", value = defaultValue("sCA_volc_pval", 5), min = 0.0, max = 110000,
              step = 0.1
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            verbatimTextOutput("sCA_volc_selected")
          )
        ),
        fluidRow(
          column(
            width = 12,
            if ("manhattanly" %in% rownames(installed.packages()))
            jqui_resizable(plotly::plotlyOutput("sCA_volcanoPlot"))
          )
        ),
        br(),
        actionButton("save2HistVolc", "save to history")
      ),
      tabPanel(
        title = "Differentially Expressed Genes",  width = 12, 
        collapsible = FALSE, collapsed = TRUE,
        fluidRow(
          column(
            12,
            tableSelectionUi("sCA_dgeTable")
          )
        )
      )
    )
  )
)
