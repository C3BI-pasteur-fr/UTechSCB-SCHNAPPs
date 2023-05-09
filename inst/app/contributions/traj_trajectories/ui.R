require(shiny)
require(shinydashboard)
require(shinyBS)
require(shinycssloaders)
require(plotly)
# menu list
# defines the main entry

if (!is.null(.schnappsEnv$enableTrajectories)) {
  
  
  if (.schnappsEnv$enableTrajectories ){
    menuList <- list(
      menuItem("Trajectories",
               # id="trajectoryID",
               tabName = "TrajectoryList", startExpanded = FALSE,
               if("SCORPIUS" %in% installed.packages())
                 menuSubItem("Scorpius", tabName = "scorpiusTab")
               else {
                 cat(file = stderr(), "\n\nPlease install SCORPIUS:
                 devtools::install_github(\"rcannood/SCORPIUS\")\n\n")
               },
               if("ElPiGraph.R" %in% installed.packages()) 
                 menuSubItem("ELPIGraph", tabName = "elpiGraphTab")
               else {
                 cat(file = stderr(), "\n\nPlease install ElpiGraph:
                 devtools::install_github(\"Albluca/distutils\")  
                 devtools::install_github(\"Albluca/ElPiGraph.R\")\n\n")
               },
               if("Tempora" %in% installed.packages()) 
                 menuSubItem("Tempora", tabName = "temporaTab")
               else {
                 cat(file = stderr(), "\n\nPlease install ElpiGraph:
                 devtools::install_github(\"BaderLab/Tempora\") \n\n")
               }
      )
    )
  }
  
  
  # tabs with the actual content
  tabList <- list(
    tabItem(
      "scorpiusTab",
      shinydashboardPlus::box(
        title = "Scorpius trajectory inference", solidHeader = TRUE, width = 12, status = "primary",
        
        fluidRow(
          column( offset = 3,
                  width = 12, 
                  cellSelectionUI("Scorpius_dataInput")
          )),
        br(),
        fluidRow(
          column(
            width = 4,
            sc_selectInput("dimScorpiusX",
                           label = "Component 1",
                           choices = defaultValue("dimScorpiusX", "notyet"),
                           selected = defaultValue("dimScorpiusX", "notyet")
            ),
            sc_numericInput("scorpRepeat",
                            label = "number of permutations for random forrest",
                            min = 1, max = 100, step = 1,
                            value = defaultValue("scorpRepeat", 3)
            )
          ),
          column(
            width = 4,
            sc_selectInput("dimScorpiusY",
                           label = "Component 2",
                           choices = defaultValue("dimScorpiusY", "notyet"),
                           selected = defaultValue("dimScorpiusY", "notyet")
            ),
            sc_numericInput("scorpMaxGenes",
                            label = "max number of Genes",
                            min = 200, max = 20000, step = 10,
                            value = defaultValue("scorpMaxGenes", 500)
            )
          ),
          column(
            width = 4,
            sc_selectInput("dimScorpiusCol",
                           label = "Color by",
                           choices = defaultValue("dimScorpiusCol", "notyet"),
                           selected = defaultValue("dimScorpiusCol", "notyet")
            ),
            fileInput("trajInputFile",
                      "Choose .csv file with trajectory informaiton",
                      accept = c(
                        ".csv",
                        "text/comma-separated-values",
                        "text/tab-separated-values",
                        "text/plain",
                        ".csv",
                        ".tsv"
                      )
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            plotOutput("scropius_trajectory_plot", height = "672px") %>% jqui_resizable()
          )
        ),
        br(),
        fluidRow(
          column(
            width = 12, offset = 0,
            actionButton("updatetScorpiusParameters", "apply changes", width = "80%")
          )
        ),
        br(),
        # tags$h3("Heatmap "),
        fluidRow(
          column(
            width = 12,
            pHeatMapUI("scorpiusHeatmapPlotModule") # %>% withSpinner()
            # imageOutput('scorpiusHeatmapPlotModule', height = '672px')
          )
        ),
        fluidRow(
          column(
            width = 10,
            tableSelectionUi("scorpiusTableMod")
          )
        )
      )
    ),
    tabItem(
      "elpiGraphTab",
      shinydashboardPlus::box(
        title = "ElpiGraph trajectory inference", solidHeader = TRUE, width = 12, status = "primary",
        # tags$h3("trajectory by ElpiGraph"),
        # fluidRow(column(12,
        #                 offset = 1,
        #                 checkboxInput("elpiCalc", "calculate", FALSE)
        # )),
        fluidRow(
          column( offset = 3,
                  width = 12, 
                  cellSelectionUI("Elpi_dataInput")
          )),
        br(),
        fluidRow(
          column(
            3,
            sc_selectInput(
              "dimElpiX",
              label = "Component 1",
              choices = defaultValue("dimElpiX", "notyet"),
              selected = defaultValue("dimElpiX", "notyet")
            )
          ),
          column(
            3,
            sc_selectInput(
              "dimElpiY",
              label = "Component 2",
              choices = defaultValue("dimElpiY", "notyet"),
              selected = defaultValue("dimElpiY", "notyet")
            )
          ),
          column(
            3,
            sc_selectInput(
              "dimElpiCol",
              label = "Color by",
              choices = defaultValue("dimElpiCol", "notyet"),
              selected = defaultValue("dimElpiCol", "notyet")
            )
          ),
          column(
            3,
            sc_numericInput(
              inputId = "elpiSeed",
              label = "Seed",
              value = defaultValue("elpiSeed", 9),
              min = 1, max = 1000, step = 1
            )
          )
        ),
        fluidRow(
          column(
            3,
            sc_selectInput(
              "dimElpi",
              label = "Dimensions to use",
              choices = c("elpiPCA", "components"),
              selected = defaultValue("dimElpi", "components")
            )
          ),
          column(
            3,
            sc_selectInput(
              "ElpiMethod",
              label = "Method to use",
              choices = c(
                "computeElasticPrincipalCurve",
                "computeElasticPrincipalTree",
                "computeElasticPrincipalCircle"
              ),
              selected = defaultValue("ElpiMethod", "computeElasticPrincipalTree")
            )
          ),
          column(
            2,
            sc_numericInput(
              inputId = "elpiNumNodes",
              label = "Number of nodes",
              value = defaultValue("elpiNumNodes", 20),
              min = 10, max = 100, step = 1
            )
          ),
          column(
            2,
            sc_numericInput(
              inputId = "elpinReps",
              label = "Number of repeats",
              value = defaultValue("elpinReps", 1),
              min = 1, max = 50, step = 1
            )
          ),
          column(
            2,
            sc_numericInput(
              inputId = "elpiProbPoint",
              label = "probability of inclusing of a single point for each computation",
              value = defaultValue("elpiProbPoint", 0.6),
              min = 0.1, max = 1, step = 0.1
            )
          )
        ),
        fluidRow(column(
          12,
          plotOutput("elpi_plot", height = "672px") %>% jqui_resizable() # 
        )),
        br(),
        fluidRow(
          column(
            width = 12, offset = 0,
            actionButton("elpiCalc", "apply changes", width = "80%")
          )
        ),
        br(),
        fluidRow(
          column(
            2,
            selectInput(
              inputId = "elpiStartNode",
              label = "start node of trajectory analysis",
              choices = defaultValue("elpiStartNode", 1),
              selected = defaultValue("elpiStartNode", 1)
            )
          ),
          column(
            2,
            selectInput(
              inputId = "elpiEndNode",
              label = "end node of trajectory analysis",
              choices = defaultValue("elpiEndNode", 2),
              selected = defaultValue("elpiEndNode", 2)
            )
          ),
          column(
            2,
            sc_numericInput(
              inputId = "elpi_num_permutations",
              label = "elpi_num_permutations",
              value = defaultValue("elpi_num_permutations", 3)
            )
          ),
          column(
            2,
            sc_numericInput(
              inputId = "elpi_ntree",
              label = "elpi_ntree",
              value = defaultValue("elpi_ntree", 10000)
            )
          ),
          column(
            2,
            sc_numericInput(
              inputId = "elpi_ntree_perm",
              label = "elpi_ntree_perm",
              value = defaultValue("elpi_ntree_perm", 1000)
            )
          ),
          column(
            2,
            sc_numericInput(
              inputId = "elpi_nGenes",
              label = "number of output genes",
              value = defaultValue("elpi_nGenes", 50)
            )
          )
        ),
        shinydashboardPlus::box(
          which = "plot", width = 12,
          fluidRow(column(
            12,
            # plotOutput("elpi_heatmap", height = "672px")
            pHeatMapUI("elpiHeatmapPlotModule")
          )),
          br(),
          tags$h3("table"),
          fluidRow(column(
            12,
            offset = 0,
            tableSelectionUi("elpiTableMod")
          )),
          br(),
          # tags$h3("Heatmap "),
          fluidRow(column(
            12,
            plotOutput("elpi_histo", height = "672px") # %>% withSpinner()
          ))
        )
      )
    ),
    tabItem("temporaTab",
            tabBox(title = "Tempora", width = 12, id = "temporaTab",
                   tabPanel(
                     title = "Parameters for Tempora trajectory inference", solidHeader = TRUE, width = 12, 
                     value = "temporaParameters",
                     # The id lets us use input$tabset1 on the server to find the current tab
                     id = "tabsetTempora",
                     fluidRow(column(
                       width = 12, offset = 0,
                       br("uses transformed data"))),
                     br(),
                     fluidRow(
                       column( offset = 3,
                               width = 12, 
                               cellSelectionUI("Tempora_dataInput")
                       )),
                     br(),
                     fluidRow(
                       column(
                         width = 12, offset = 0,
                         actionButton("updatetTemporaParameters", "apply changes", width = "80%")
                       )
                     ),
                     fluidRow(
                       column(4,
                              offset = 0,
                              sc_selectInput("temporaCluster", "cluster points to be used",
                                             choices = c("dbCluster"),
                                             selected = defaultValue("temporaCluster", "dbCluster")
                              )
                       ),
                       column(4,
                              offset = 0,
                              sc_selectInput("temporaFactor", "time variable",
                                             choices = c("sampleNames", "dbCluster"),
                                             selected = defaultValue("temporaFactor", "sampleNames")
                              )),
                       column(4,
                              offset = 0,
                              sc_selectInput("temporaLevels", "Ordered time points",
                                             choices = c("AVC","MVE16.5","MV1","MV2"),
                                             multiple = TRUE,
                                             selected = defaultValue("temporaLevels", c("AVC","MVE16.5","MV1","MV2"))
                              )
                              
                              
                       )
                     ),
                     fluidRow(
                       column(
                         6,
                         offset = 0,
                         fileInput(
                           "temporaGMTFile",
                           "GMT file to use",
                           accept = c(".gmt"),
                           placeholder = "no file selected",
                           multiple = FALSE,
                         ) %>% setId(id="temporaGMTFile")
                       ),
                       column(
                         3,
                         sc_numericInput(
                           inputId = "temporaMinSz",
                           label = "min size of gene sets",
                           value = defaultValue("temporaMinSz", 5)
                         )
                       ),
                       column(
                         3,
                         sc_numericInput(
                           inputId = "temporaMaxSz",
                           label = "max size of gene sets",
                           value = defaultValue("temporaMaxSz", 200)
                         )
                       )
                     ),
                     fluidRow(column(
                       4,
                       sc_numericInput(
                         inputId = "temporaNPCs",
                         label = "N PCs to use",
                         value = defaultValue("temporaNPCs", 12)
                       )
                     ),
                     column(
                       4,
                       sc_numericInput(
                         inputId = "temporaDiff_thresh",
                         label = "Percent of permissible difference",
                         value = defaultValue("temporaDiff_thresh", 0.01)
                       )
                     ),
                     column(
                       4,
                       sc_numericInput(
                         inputId = "temporaPval_thresh",
                         label = "max p-value",
                         value = defaultValue("temporaPval_thresh", 0.99)
                       )
                     )
                     ),
                     fluidRow(column(
                       12,
                       plotOutput("tempora_screeplot", height = "672px") %>% jqui_resizable()
                     )),
                     br(),
                     actionButton("save2Hist_tempora_screeplot", "save to history"),
                     br(),
                     fluidRow(column(
                       12,
                       plotOutput("tempora_plot", height = "672px") %>% jqui_resizable() 
                     )),
                     br(),
                     actionButton("save2Hist_tempora_plot", "save to history"),
                     br(),
                     checkbsTT("temporaCluster"),
                     checkbsTT("temporaFactor"),
                     checkbsTT("temporaLevels"), 
                     checkbsTT("temporaGMTFile"),
                     checkbsTT("temporaMinSz"),
                     checkbsTT("temporaMaxSz"),
                     checkbsTT("temporaNPCs"),
                     checkbsTT("temporaDiff_thresh"),
                     checkbsTT("temporaPval_thresh"),
                     # n_pcs = 12
                     # difference_threshold = 0.01
                     # pval_threshold = 0.5
                     
                   ),
                   tabPanel(
                     title = "p-Values of GO terms", solidHeader = TRUE, width = 12, 
                     value = "temporaTable",
                     id = "tabsetTemporaGOTable",
                     fluidRow(column(
                       12,
                       offset = 0,
                       tableSelectionUi("temporaGOTableMod")
                     )),
                     br(),
                     fluidRow(column(
                       12,
                       plotOutput("temporaSelectedGOs", height = "672px") %>% jqui_resizable() 
                     )),
                     br(),
                     fluidRow(column(12,
                                     verbatimTextOutput("coE_temporaPWgenes"))),
                     actionButton("save2Hist_temporaSelectedGOs", "save to history")
                   ),
                   tabPanel(
                     title = "Trajectory inquiry", solidHeader = TRUE, width = 12, 
                     value = "temporaScorpius",
                     id = "tabsetTemporaScorpius",
                     fluidRow(column(
                       width = 4,
                       sc_selectInput("dimTemporaX",
                                      label = "Component 1",
                                      choices = defaultValue("dimTemporaX", "notyet"),
                                      selected = defaultValue("dimTemporaX", "notyet")
                       )),
                       column(
                         width = 4,
                         sc_selectInput("dimTemporaY",
                                        label = "Component 2",
                                        choices = defaultValue("dimTemporaY", "notyet"),
                                        selected = defaultValue("dimTemporaY", "notyet")
                         ))
                     ),
                     br(),
                     fluidRow(column(
                       12,
                       offset = 0,
                       plotly::plotlyOutput("tempora2dPlot", height = "672px") %>% jqui_resizable()
                     )),
                     br(),
                     actionButton("save2Hist_tempora2dPlot", "save to history"),
                     ######
                     br()
                   )
                   
                   
            )
    )
  )
}
# declare heavy calculations
# myHeavyCalculations <- NULL
