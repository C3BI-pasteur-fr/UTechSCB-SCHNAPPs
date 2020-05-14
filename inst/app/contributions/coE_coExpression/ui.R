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
                           shinydashboard::menuSubItem("alluvialTab", tabName = "alluvialTab")
                           # shinydashboard::menuSubItem("SOM cluster", tabName = "SOMcluster")
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
          textInput("coE_heatmap_geneids", "Comma seperated gene names", value = defaultValue("coE_heatmap_geneids", defaultValueMultiGenes))
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
    tabBox(title = "Interogate Cells", width = 12, id = "interogCells",
           tabPanel(
             title = "2D Plot", width = 12, value = "selected2DPlot",
             footer = tags$p("2D plot for selecting cells"),
             fluidRow(
               column(
                 width = 12,
                 clusterUI("coE_selected")
               )
             )
           ),
           tabPanel(
             title = "Detailed information on genes", 
             fluidRow(width = 12,
                      column( width = 12,
                              fluidRow(column(
                                width = 12,
                                textInput("coE_heatmapselected_geneids", "Comma seperated gene names", value = defaultValue("coE_heatmapselected_geneids", defaultValueMultiGenes))
                              )),
                              br(),
                              # fluidRow(checkboxInput(inputId = "coE_heatmapSelectedModuleShow", label = "calc heatmap", value = FALSE)),
                              tabBox(title = "cell group info", width = 12, id = "cellGroupInfo",
                                     tabPanel(
                                       title = "Heatmap of selected cells and genes", width = 12, 
                                       fluidRow(
                                         column(
                                           width = 12, offset = 1,
                                           actionButton("updateHeatMapSelectedParameters", "apply changes", width = "80%")
                                         )
                                       ),
                                       br(),
                                       fluidRow(
                                         column(
                                           width = 12,
                                           pHeatMapUI("coE_heatmapSelectedModule")
                                         )
                                       )
                                     ),
                                     
                                     tabPanel(
                                       title = "Table with coefficient of variance", width = 12, 
                                       footer = div(
                                         tags$ul(
                                           tags$li(
                                             "for each gene we divice the standard diviation by the mean (coefficient of variance, cv)\n"
                                           )
                                         )
                                       ),
                                       fluidRow(
                                         column(
                                           width = 12, offset = 1,
                                           actionButton("updateMinExprSelectedParameters", "apply changes", width = "80%")
                                         )
                                       ),
                                       br(),
                                       fluidRow(
                                         column(
                                           width = 3,
                                           numericInput("coEtgMinExpr", "min UMI count per gene:",
                                                        defaultValue("coEtgMinExpr", 1),
                                                        min = 0, max = 100000
                                           )
                                         ), column(
                                           width = 3,
                                           numericInput("coEtgPerc", "min percentage of cells expressing a genes:",
                                                        defaultValue("coEtgPerc", 60),
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
                                     tabPanel(
                                       title = "Table with correlation coefficients", width = 12,
                                       footer = div(
                                         tags$ul(
                                           tags$li(
                                             "using Hmisc we calculate the correlation coefficient and p-values for the selected genes and cells\n"
                                           )
                                         )
                                       ),
                                       fluidRow(
                                         column(
                                           width = 12, offset = 1,
                                           actionButton("updatetopCCGenesSelectedParameters", "apply changes", width = "80%")
                                         )
                                       ),
                                       br(),
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
             )
           )
    )
  ),
  # CoExpressionViolin ----
  expressionTab = shinydashboard::tabItem(
    "CoExpressionViolin",
    tabBox(title = "Violin plots", width = 12, id = "violinPlots",
           tabPanel(
             title = "permutated", value = "permViol",
             footer = "for each cell we count how many of the genes specified have an expression larger or equal than the minimum exprssion.\nThese counts are then divided up for any variable that can be used as a factor (has less than 20 levels).",
             
             fluidRow(
               column(
                 width = 12,
                 checkboxInput("coE_showPermutations", "show Permutations", value = defaultValue("coE_showPermutations",FALSE))
               )
             ),
             fluidRow(
               column(
                 width = 4,
                 textInput("coE_geneGrpVioIds", "Comma seperated gene names", value = defaultValue("coE_geneGrpVioIds", defaultValueMultiGenes))
               ),
               column(
                 width = 4,
                 selectInput(
                   "coE_dimension_xVioiGrp",
                   label = "X",
                   choices = c(defaultValue("coE_dimension_xVioiGrp", "dbCluster"), "sampleName", "tsne3"),
                   selected = defaultValue("coE_dimension_xVioiGrp", "dbCluster")
                 )
               ),
               column(
                 width = 4,
                 numericInput("coEminExpr", "min expression of genes:",
                              value = defaultValue("coEminExpr", 1),
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
           ),
           tabPanel(
             title = "grouped", value = "grpViol",
             fluidRow(
               column(
                 width = 4,
                 textInput("coE_geneGrpVioIds2", "Comma seperated gene names", value = defaultValue("coE_geneGrpVioIds2", defaultValueMultiGenes))
               ),
               column(
                 width = 4,
                 selectizeInput(
                   "coE_dimension_xVioiGrp2",
                   label = "X",
                   choices = c(defaultValue("coE_dimension_xVioiGrp2", "dbCluster"), "sampleName", "tsne3"),
                   selected = defaultValue("coE_dimension_xVioiGrp2", "dbCluster"), options = list(maxItems = 2)
                 )
               ),
               column(
                 width = 4,
                 numericInput("coEminExpr2", "min expression of genes:",
                              defaultValue("coEminExpr2", 1),
                              min = 1, max = 100000
                 )
               ),
             ),
             br(),
             fluidRow(column(width = 12,
                             jqui_resizable(plotly::plotlyOutput("coE_geneGrp_vio_plot2") )
                             # jqui_resizable(plotOutput("coE_geneGrp_vio_plot") )
             )),
             br(),
             actionButton("save2HistVio2", "save to history")
           )
           
    )
  ),
  tabItem("alluvialTab",
          box(
            title = "alluvial plot", solidHeader = TRUE, width = 12, status = 'primary', 
            fluidRow(
              column(width = 6, 
                     selectInput("alluiv1", "1st axsis", choices = defaultValue("alluiv1", "notyet"), selected = defaultValue("alluiv1", "notyet"))),
              column(width = 6, 
                     selectInput("alluiv2", "2nd axsis", choices = defaultValue("alluiv2", "notyet"), selected = defaultValue("alluiv2", "notyet")))
            ),
            fluidRow(
              column(width = 12, 
                     plotOutput("alluvial_plot") # %>% withSpinner()
              )
            )
          )
          
  )
)
