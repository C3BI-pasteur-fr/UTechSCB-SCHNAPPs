# inst/app/contributions/coE_coExpression/ui.R
# 
suppressMessages(library(magrittr))

source(paste0(packagePath, "/modulesUI.R"), local = TRUE)

# list of menu Items
menuList <- list(
  shinydashboard::menuItem("Co-expression",
                           icon = icon("tachometer-alt"),
                           expandedName = "coexpressionID",
                           tabName = "coexpression", startExpanded = FALSE,
                           shinydashboard::menuSubItem("All clusters", tabName = "coexpressionAll"),
                           shinydashboard::menuSubItem("Gene sets", tabName = "geneSets"),
                           shinydashboard::menuSubItem("Selected", tabName = "coexpressionSelected"),
                           shinydashboard::menuSubItem("Violin plot", tabName = "CoExpressionViolin"),
                           if(("ggalluvial" %in% rownames(installed.packages()))){
                             shinydashboard::menuSubItem("alluvialTab", tabName = "alluvialTab")
                           }else {
                             cat(file = stderr(), "Please install ggalluvial: install.packages('ggalluvial')")
                           }
                           # shinydashboard::menuSubItem("SOM cluster", tabName = "SOMcluster")
  )
)



# list of tab Items
tabList <- list(
  # coexpressionAll ----
  coexpressionAllTab = shinydashboard::tabItem(
    "coexpressionAll",
    tabBox(title = "Heatmap of all cells", width = 12, id = "heatMapTabBox",
      # title = "Heatmap of all cells", solidHeader = TRUE, width = 12, status = "primary",
      tabPanel(
        title = "Heat map",
        value = "heatmapPanel",
        fluidRow(
        column(
          width = 12,
          sc_textInput("coE_heatmap_geneids", "Comma separated gene names", 
                       width = '100%',
                       value = defaultValue("coE_heatmap_geneids", defaultValueMultiGenes)),
          fluidRow(column(
            width = 6,
            sc_selectInput("coE_subSampleFactor", "factor to subsample", selected = defaultValue("coE_subSampleFactor", "dbCluster"),
                           choices = c("dbCluster", "sample"))),
            column(
              width = 6,
              sc_numericInput("coE_nSubsample", "max number of cells per group", value = defaultValue("coE_nSubsample", 1000), min=10)
            )
          ),
        )
      ),
      fluidRow(
        column(
          width = 12,
          pHeatMapUI("coExpHeatmapModule")
        )
      )
      ),
      tabPanel(
        title = "Scran find markers", value = "scranFindMarkersPanel",
        fluidRow(
          column(
            width = 3,
            sc_numericInput("coE_nFindMarker", "number of markerGenes to plot (if gene names are empty)", 10, min=2),
            verbatimTextOutput("coE_objSize") 
          ),
          column(
            width = 3,
            sc_numericInput("coE_lfc", "min log fold change", 2, min=0.1),
            sc_numericInput("coE_nCPU", "number of CPUs to use", 2, min=1, max=(parallel::detectCores()-1) )
          ),
          column(
            width = 3,
            sc_selectInput("coE_direction", "direction", "up", choices = c("up", "down", "any")),
            sc_selectInput(inputId = "coE_scranFactor", label = "select factor for findmarkers", 
                           selected = defaultValue("coE_scranFactor", "dbCluster"), 
                           choices = c("please run normalization","dbCluster"), multiple = F)
          ),
          column(width = 3,
                 actionButton(inputId = "scranFindMarkerApply", label = "apply changes", width = "80%"))
          ),
        fluidRow(
          column(width = 12,
                 verbatimTextOutput("scranFindMarkersSelected") %>% jqui_resizable()
          )
        ),
        fluidRow(
          column(width = 3,
                 sc_selectInput(inputId = "coE_scranFindMarkerCluster", label = "select cluster", 
                                 selected = defaultValue("coE_scranFindMarkerCluster", "1"), 
                                 choices = c("please run scran.findclusters","1"), multiple = F)
                 )
        ),
        fluidRow(
          column(width = 12,
                 tableSelectionUi("coE_scranFindMarkerTable"))
        )
      )
    )
  ),
  # geneSetsTab ----
  geneSetsTab = shinydashboard::tabItem(
    "geneSets",
    tabBox(title = "Gene set plots", width = 12, id = "geneSetPlotsTabBox",
           tabPanel(
             title = "Dot plot",
             value = "dotPlotGeneSet",
             shinydashboardPlus::box(
               title = "Dot plot", solidHeader = TRUE, width = 12, status = "primary",
               fluidRow(
                 column(
                   width = 4,
                   sc_selectizeInput("coE_dotPlot_geneSets", label = "Gene sets for y axis", 
                                  choices = c("please load GMT file"),
                                  selected = .schnappsEnv$coE_dotPlot_geneSets,
                                  multiple = T),
                   sc_numericInput("coE_dotPlot_col.min", "Minimum scaled average expression threshold",
                                   defaultValue("coE_dotPlot_col.min", -2.5),
                                   min = -10, max = 10, step = 0.1
                   ),
                   sc_numericInput("coE_dotPlot_col.max", "Maximum scaled average expression threshold",
                                   defaultValue("coE_dotPlot_col.max", 2.5),
                                   min = -10, max = 10, step = 0.1
                   )
                 ),column(
                   width = 4,
                   sc_selectInput(
                     "coE_dimension_ydotPlotClusters",
                     label = "Y",
                     choices = c(defaultValue("coE_dimension_ydotPlotClusters", "dbCluster"), "sampleNames", "tsne3"),
                     selected = defaultValue("coE_dimension_ydotPlotClusters", "dbCluster")
                   ),
                   sc_numericInput("coE_dotPlot_dot.min", "The fraction of cells at which to draw the smallest dot",
                                   defaultValue("coE_dotPlot_dot.min", 0),
                                   min = 0, max = 1, step = 0.1
                   ),
                   sc_numericInput("coE_dotPlot_dot.scale", "Scale the size of the points",
                                   defaultValue("coE_dotPlot_dot.scale", 6),
                                   min = 1, max = 20, step = 0.1
                   )
                 ),column(
                   width = 4,
                   sc_selectInput(
                     "coE_dotPlot_col",
                     label = "col",
                     choices = rownames(x = brewer.pal.info),
                     selected = "RdBu"
                   ),
                   sc_selectInput(
                     "coE_dotPlot_scale.by",
                     label = "Scale the size of the points",
                     choices = c("size", "radius"),
                     selected = "radius"
                   )
                   # not needed? :
                   #' @param scale.by Scale the size of the points by 'size' or by 'radius'
                   #' @param scale.min Set lower limit for scaling, use NA for default
                   #' @param scale.max Set upper limit for scaling, use NA for default
                   
                 )
               ),
               
               fluidRow(
                 column(
                   width = 12,
                   # jqui_resizable(plotly::plotlyOutput("coE_dotPlot_GeneSets"))
                   plotlyOutput("coE_dotPlot_GeneSets")%>% jqui_resizable()
                 )
               ),
               br(),
               actionButton("save2histDotPlot", "save to history")
             )
           ),
           tabPanel(
             title = "Dot plot with Module Score",
             value = "dotPlotMscore",
             shinydashboardPlus::box(
               title = "Dot plot with Module Score", solidHeader = TRUE, width = 12, status = "primary",
               fluidRow(
                 column(
                   width = 4,
                   sc_selectizeInput("coE_dotPlotModuleScore_geneSets", label = "Gene sets for y axis", 
                                  choices = c("please load GMT file"),
                                  selected = .schnappsEnv$coE_dotPlotModuleScore_geneSets,
                                  multiple = T),
                   sc_numericInput("coE_dotPlotModuleScore_col.min", "Minimum scaled average expression threshold",
                                   defaultValue("coE_dotPlotModuleScore_col.min", -2.5),
                                   min = -10, max = 10, step = 0.1
                   ),
                   sc_numericInput("coE_dotPlotModuleScore_col.max", "Maximum scaled average expression threshold",
                                   defaultValue("coE_dotPlotModuleScore_col.max", 2.5),
                                   min = -10, max = 10, step = 0.1
                   )
                 ),column(
                   width = 4,
                   sc_selectInput(
                     "coE_dimension_ydotPlotModuleScoreClusters",
                     label = "Y",
                     choices = c(defaultValue("coE_dimension_ydotPlotModuleScoreClusters", "dbCluster"), "sampleNames", "tsne3"),
                     selected = defaultValue("coE_dimension_ydotPlotModuleScoreClusters", "dbCluster")
                   ),
                   sc_numericInput("coE_dotPlotModuleScore_dot.min", "The fraction of cells at which to draw the smallest dot",
                                   defaultValue("coE_dotPlotModuleScore_dot.min", 0),
                                   min = 0, max = 1, step = 0.1
                   ),
                   sc_numericInput("coE_dotPlotModuleScore_dot.scale", "Scale the size of the points",
                                   defaultValue("coE_dotPlotModuleScore_dot.scale", 6),
                                   min = 1, max = 20, step = 0.1
                   )
                 ),column(
                   width = 4,
                   sc_selectInput(
                     "coE_dotPlotModuleScore_col",
                     label = "col",
                     choices = rownames(x = brewer.pal.info),
                     selected = "RdBu"
                   ),
                   sc_selectInput(
                     "coE_dotPlotModuleScore_scale.by",
                     label = "Scale the size of the points",
                     choices = c("size", "radius"),
                     selected = "radius"
                   )
                   # not needed? :
                   #' @param scale.by Scale the size of the points by 'size' or by 'radius'
                   #' @param scale.min Set lower limit for scaling, use NA for default
                   #' @param scale.max Set upper limit for scaling, use NA for default
                   
                 )
               ),
               
               fluidRow(
                 column(
                   width = 12,
                   # jqui_resizable(plotly::plotlyOutput("coE_dotPlot_GeneSets"))
                   plotlyOutput("coE_dotPlot_GeneSetsModuleScore") %>% jqui_resizable()
                 )
               ),
               br(),
               actionButton("save2histDotPlotModuleScore", "save to history")
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
                                sc_textInput("coE_heatmapselected_geneids", "Comma separated gene names", 
                                             value = defaultValue("coE_heatmapselected_geneids", defaultValueMultiGenes))
                              )),
                              br(),
                              # fluidRow(sc_checkboxInput(inputId = "coE_heatmapSelectedModuleShow", label = "calc heatmap", value = FALSE)),
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
                                           sc_numericInput("coEtgMinExpr", "min UMI count per gene:",
                                                           defaultValue("coEtgMinExpr", 1),
                                                           min = 0, max = 100000
                                           )
                                         ), column(
                                           width = 3,
                                           sc_numericInput("coEtgPerc", "min percentage of cells expressing a genes:",
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
                                           # sc_checkboxInput(inputId = "coE_topCCGenesShow", label = "calc correlations", value = FALSE)),
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
             title = "combinations", value = "permViol",
             footer = "for each cell we count how many of the genes specified have an expression larger or equal than the minimum exprssion.\nThese counts are then divided up for any variable that can be used as a factor (has less than 20 levels).",
             fluidRow(
               column(
                 width = 12,
                 div(
                   p(
                     "For each cell we count how many of the genes specified have an expression within the specified range.\nThese counts are then divided up for any variable that can be used as a factor (has less than 20 levels)."
                   ),p(
                     "show expression values allows showing a regular violin plot."
                   ),
                   align = "left"
                 )
               )
             ),
             fluidRow(
               column(
                 width = 4,
                 sc_checkboxInput("coE_showPermutations", "show combinations", value = defaultValue("coE_showPermutations",FALSE))
               ),
               column(
                 width = 4,
                 sc_checkboxInput("coE_showExpression", "show expression values rather than count cells", value = defaultValue("coE_showExpression",FALSE))
               ),
               column(
                 width = 4,
                 sc_selectInput(
                   "coE_scale",
                   label = "Scale",
                   choices = c("area", "count", "width"),
                   selected = defaultValue("coE_scale", "count")
                 )
               )
             ),
             fluidRow(
               column(
                 width = 4,
                 sc_textInput("coE_geneGrpVioIds", "Comma separated gene names", value = defaultValue("coE_geneGrpVioIds", defaultValueMultiGenes))
               ),
               column(
                 width = 4,
                 sc_selectInput(
                   "coE_dimension_xVioiGrp",
                   label = "X",
                   choices = c(defaultValue("coE_dimension_xVioiGrp", "dbCluster"), "sampleNames", "tsne3"),
                   selected = defaultValue("coE_dimension_xVioiGrp", "dbCluster")
                 )
               ),
               column(
                 width = 4,
                 sliderInput(
                   "coEminMaxExpr",
                   label = "min/max value for heatmap",
                   min = -10000000,
                   max = 100000000,
                   value = defaultValue("coEminMaxExprValue", c(-1,1))
                 )
                 # sc_numericInput("coEminExpr", "min expression of genes:",
                 #              value = defaultValue("coEminExpr", 1),
                 #              min = -1000, max = 100000
                 # )
               )
             ),
             br(),
             fluidRow(column(width = 12,
                             # jqui_resizable(plotly::plotlyOutput("coE_geneGrp_vio_plot") )
                             plotOutput("coE_geneGrp_vio_plot")  %>% jqui_resizable( )
             )),
             br(),
             actionButton("save2HistVio", "save to history")
           ),
           tabPanel(
             title = "grouped", value = "grpViol",
             fluidRow(
               column(
                 width = 4,
                 sc_textInput("coE_geneGrpVioIds2", "Comma separated gene names", value = defaultValue("coE_geneGrpVioIds2", defaultValueMultiGenes))
               ),
               column(
                 width = 4,
                 sc_selectizeInput(
                   "coE_dimension_xVioiGrp2",
                   label = "X",
                   multiple = TRUE,
                   choices = c(defaultValue("coE_dimension_xVioiGrp2", "dbCluster"), "sampleNames", "tsne3"),
                   selected = defaultValue("coE_dimension_xVioiGrp2", "dbCluster"), 
                   options = list(maxItems = 2)
                 )
               ),
               column(
                 width = 4,
                 sliderInput(
                   "coEminMaxExpr2",
                   label = "min/max value for heatmap",
                   min = -10000000,
                   max = 100000000,
                   value = defaultValue("coEminMaxExpr2", c(-1,1))
                 )# sc_numericInput("coEminExpr2", "min expression of genes:",
                 #              defaultValue("coEminExpr2", 1),
                 #              min = -1000, max = 100000
                 # )
               ),
             ),
             br(),
             fluidRow(column(width = 12,
                             plotly::plotlyOutput("coE_geneGrp_vio_plot2")  %>% jqui_resizable()
                             # jqui_resizable(plotOutput("coE_geneGrp_vio_plot") )
             )),
             br(),
             actionButton("save2HistVio2", "save to history")
           )
           
    )
  ),
  tabItem("alluvialTab",
          shinydashboardPlus::box(
            title = "alluvial plot", solidHeader = TRUE, width = 12, status = 'primary', 
            fluidRow(
              column(width = 6, 
                     sc_selectInput("alluiv1", "1st axis", choices = defaultValue("alluiv1", "notyet"), selected = defaultValue("alluiv1", "notyet"))),
              column(width = 6, 
                     sc_selectInput("alluiv2", "2nd axis", choices = defaultValue("alluiv2", "notyet"), selected = defaultValue("alluiv2", "notyet")))
            ),
            fluidRow(
              column(width = 12, 
                     plotOutput("alluvial_plot")  #%>% jqui_resizable()
              )
            ),
            br(),
            actionButton("save2HistAlluvial", "save to history")
            
          )
          
  )
)
