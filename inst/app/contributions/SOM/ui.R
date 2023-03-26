# menuList <- list(
#   shinydashboard::menuItem("SOM",
#                            icon = icon("dashboard"),
#                            # id="coexpressionID",
#                            tabName = "som", startExpanded = FALSE,
#                            shinydashboard::menuSubItem("SOM cluster", tabName = "SOMcluster")
#   )
# )

# SOM /metaCell clustering ----
tabList = list(
  shinydashboard::tabItem(
    "metaCell",
    tabBox(title = "MetaCell", width = 12, id = "metaCellTabBox",
           tabPanel(
             title = "setup", width = 12, value = "mcSetup",
             box(
               title = "Meta cell", solidHeader = TRUE, width = 12, status = "primary",
               footer = "metaCell [REF] adapted to SCHNAPPs",
               fluidRow(
                 column(
                   width = 5, offset = 1,
                   actionButton("updateMetaCellParameters", "apply changes", width = "80%")
                 ),
                 column(
                   width = 5, offset = 1,
                   sc_checkboxInput("mC_recalc", "recalculate all",
                                    defaultValue("mC_recalc", TRUE))
                 )
               )
             ),
             br(),
             if(is.null(.schnappsEnv$historyPath)){
               span(("Please use the historyPath argument when calling schnapps"), style="color:red")
             }else {
               list(
                 span("cleaning of data is done with the tools supplied in SCHNAPPs"),
                 br(),
                 span("We a new gene set with all genes for which the scaled variance is T_vm and higher. 
          Then restrict this gene set to genes with at least T_tot UMIs across the entire dataset, 
          and also requires selected genes to have at least three cells for more than T_top3 UMIs were recorded."),
                 br(),
                 fluidRow(
                   cellSelectionUI("mC_metaCell_dataInput")
                 ),
                 br(),
                 fluidRow(
                   column(width = 3,
                          sc_numericInput("mC_bin_for_cutoff", "bin size for histogram used to auto find min_umis_cutoff",
                                          defaultValue("mC_bin_for_cutoff", 5),
                                          min = 2, max = 10000
                          )),
                   column(width = 3,
                          sc_numericInput("mC_T_vm", "T_vm: the threshold on normalized var/mean [0.2]",
                                          defaultValue("mC_T_vm", 0.2),
                                          min = 0, max = 1)
                   ),
                   column(width = 3,
                          sc_numericInput("mC_T_tot", "T_tot: total down sampled coverage threhsold",
                                          defaultValue("mC_T_tot", 100),
                                          min = 0, max = 1000000)
                   ),
                   column(width = 3,
                          sc_numericInput("mC_T_top3", "T_top3: threshold value for the third highest umi count for the gene (only genes with top3>T_top3 are used)",
                                          defaultValue("mC_T_top3", 2),
                                          min = 0, max = 1000000)
                   )
                 ),
                 fluidRow(
                   column(width = 3,
                          sc_numericInput("mC_K", "K: the guideline Knn parameter. The balanced will be constructed aiming at K edges per cell",
                                          defaultValue("mC_K", 100),
                                          min = 10, max = 1000000)
                   ),
                   column(width = 3,
                          sc_checkboxInput("mC_dsamp", "downsample?",
                                           defaultValue("mC_dsamp", TRUE))),
                   column(width = 3,
                          sc_numericInput("mC_min_mc_size", "target minimum metacell size. This is only an approximation and smaller MC may be returned by the algorithm",
                                          defaultValue("mC_min_mc_size", 20),
                                          min = 5, max = 1000000)),
                   column(width = 3,
                          sc_numericInput("mC_resamp_n", "number of resampling iterations",
                                          defaultValue("mC_resamp_n", 500),
                                          min = 5, max = 1000000)
                   )
                 ),
                 fluidRow(
                   column(width = 3,
                          sc_numericInput("mc_K2", "K2 determines the number of neighbors we wish to minimally associate with each cell.",
                                          defaultValue("mc_K2", 10),
                                          min = 2, max = 10000)),
                   column(width = 3,
                          sc_numericInput("mc_min_mc_size2","minimum mc size for graph cov",
                                          defaultValue("mc_min_mc_size2", 30),
                                          min = 2, max = 10000)),
                   column(width = 3,
                          sc_numericInput("mc_alpha","alpha: the threshold for filtering edges by their coclust weight is alpha * (Kth highest coclust on either node1 or node2)",
                                          defaultValue("mc_alpha", 2),
                                          min = 0, max = 10000)),
                   column(width = 3,
                          sc_numericInput("mc_T_lfc","T_lfc: log fold change for heatmap plot.",
                                          defaultValue("mc_T_lfc", 2),
                                          min = 0, max = 100))
                   
                 ),
                 fluidRow(
                   column(width = 3,
                          sc_numericInput("mc_T_gap", "T_gap: the minimal branch length for defining a supper metacell structure",
                                          defaultValue("mc_T_gap", 0.04),
                                          min = 0, max = 100))
                 )
               )
               
             }
           ),
           tabPanel(
             title = "QC", width = 12, value = "mc_QCtab",
             fluidRow(
               column(width = 12,
                      plotOutput("mc_metaDataDb_2dproj.2d_graph_proj") %>% jqui_resizable()
               )
             )
           ),
           tabPanel(
             title = "outliers", width = 12, value = "mc_outliertab",
             fluidRow(
               column(
                 width = 5, offset = 1,
                 actionButton("updateMetaCellOutliers", "update figures", width = "80%")
               )),
             br(),
             fluidRow(
               column(width = 12,
                      uiOutput("mc_outliers") %>% jqui_resizable()
               )
             )
           )
    )
  ),
  
  shinydashboard::tabItem(
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
      br(),
      fluidRow(
        cellSelectionUI("coE_SOM_dataInput"),
        box(
          fluidRow(
            column(width = 3,
                   sc_numericInput("coE_dimSOM", "number of nodes per dimension",
                                   defaultValue("coE_dimSOM", 20),
                                   min = 2, max = 100
                   )
            ), 
            column(width = 3,
                   sc_textInput("coE_geneSOM", "Gene of interest", value = defaultValue("coE_geneSOM", "notyet"))
            ),
            column(width = 3,
                   sc_selectInput("coE_distSOM",
                                  label = "Distance",
                                  choices = c("raw", "Spearman", "standardized"),
                                  selected = defaultValue("coE_distSOM","raw")
                   )
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
      fluidRow(column(
        width = 12,
        pHeatMapUI("coE_heatmapSOM")
      )
      ),
      br(),
      fluidRow(column(
        width = 12,
        verbatimTextOutput("coE_somGenes")
      )),
      br(),
      fluidRow(column(
        width = 12,
        plotOutput("coE_SOMcodebook")
      )),
      br(),
      fluidRow(column(
        width = 12,
        plotOutput("coE_SOMcomponents")
      )),
      br(),
      fluidRow(column(
        width = 12,
        plotOutput("coE_SOMuMat")
      )),
      br(),
      fluidRow(column(
        width = 6,
        sc_numericInput("coE_dimSOMX", "row",
                        defaultValue("coE_dimSOMX",1),
                        min = 1, max = 100
        )
      ),
      column(
        width = 6,
        sc_numericInput("coE_dimSOMY", "column",
                        defaultValue("coE_dimSOMY",1),
                        min = 1, max = 100
        )
      ),
      br(),
      fluidRow(column(
        width = 12,
        verbatimTextOutput("coE_somInfo")
      )),
      br(),
      fluidRow(column(
        width = 12,
        verbatimTextOutput("coE_somInfoSymbol")
      ))
      )
    )
  )
)