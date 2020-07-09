suppressMessages(require(magrittr))
source(paste0(packagePath, "/modulesUI.R"), local = TRUE)
menuList <- list(
  shinydashboard::menuItem("General QC",
                           # id="generalQCID",
                           tabName = "generalQC", icon = icon("thumbs-up"), startExpanded = FALSE,
                           shinydashboard::menuSubItem("UMI histogram", tabName = "gQC_umiHist"),
                           shinydashboard::menuSubItem("Sample histogram", tabName = "gQC_sampleHist"),
                           shinydashboard::menuSubItem("PC variance", tabName = "gQC_variancePC"),
                           shinydashboard::menuSubItem("Scater QC", tabName = "DE_scaterQC")
                           # ,
                           # shinydashboard::menuSubItem("TSNE plot", tabName = "gQC_tsnePlot"),
                           # shinydashboard::menuSubItem("Umap", tabName = "gQC_umapPlot")
  )
)



# modTab ----
modTab <- shinydashboard::tabItem(
  tabName = "modifyProj",
  fluidRow(div(h3("work with projections"), align = "center")),
  br(),
  tabBox(title = "modify projections", width = 12, id = "modProj",
         tabPanel(
           title = "Rename projections", solidHeader = TRUE, width = 12, value = "renameProj",
           fluidRow(
             column(
               width = 6,
               selectInput("oldPrj", "projections to copy + rename", choices = defaultValue("oldPrj", "notyet"), selected = defaultValue("oldPrj", "notyet"))
             ),
             column(
               width = 6,
               fluidRow(
                 column(
                   width = 8,
                   textInput("newPrj", "new name of Projection", value = defaultValue("newPrj", ""))
                 ),
                 column(
                   width = 4,
                   actionButton("updatePrjsButton", "rename")
                 )
               ),
               fluidRow(
                 column(
                   width = 8,
                   selectInput("delPrj", "projections to delete", choices = defaultValue("delPrj", "notyet"), selected = defaultValue("delPrj", "notyet"))
                 ),
                 column(
                   width = 4,
                   actionButton("delPrjsButton", "delete")
                 ),
                 tags$style(type = "text/css", "#updatePrjsButton { width:100%; margin-top: 25px;}"),
                 tags$style(type = "text/css", "#delPrjsButton { width:100%; margin-top: 25px;}")
               )
             )
           ),
           checkbsTT(item = "oldPrj"),
           checkbsTT(item = "newPrj"),
           checkbsTT(item = "updatePrjsButton"),
           checkbsTT(item = "delPrj"),
           checkbsTT(item = "delPrjsButton")
         ),
         tabPanel(
           title = "combine projections", solidHeader = TRUE, width = 12, value = "gQC_combProj",
           tags$p("Two factors can be combined by pasting the values per row and re-leveling"),
           br(),
           fluidRow(
             column(
               width = 6,
               selectInput("gQC_combPrj1", "1 st Projection", choices = c("notyet"), selected = "notyet")
             ),
             column(
               width = 6,
               
               selectInput("gQC_combPrj2", "2nd Projections", choices = c("notyet"), selected = "notyet")
             )),
           fluidRow(
             column(width = 6, 
                    textInput("gQC_newCombPrj", "name of new Projection", value = "")),
             column(
               width = 6,
               actionButton("gQC_updateCombPrjsButton", "apply")
             )
           ),
           fluidRow(column(
             width = 12,
             tableSelectionUi("gQC_projCombTableMod")
           )
           )
         ),
         tabPanel(
           title = "rename levels", width = 12, value = "gQC_renameLev",
           tags$p("rename the levels of a factor"),
           br(),
           
           fluidRow(
             column(
               width = 6,
               selectInput("gQC_rnProj", "Projection to modify", choices = c("notyet"), selected = "notyet")
             ),
             column(width = 6, 
                     textInput("gQC_newRnPrj", "name of new Projection", value = ""))
           ),
           
           fluidRow(
             column(
               width = 12,
               tags$p(tags$b("original values")),
               textOutput("gQC_orgLevels"),
               br()
             )
           ),
           fluidRow(
             column(
               width = 12,
               textAreaInput("gQC_renameLev", "new levels")
             )
           ),
           fluidRow(
             column(
               width = 6,
               actionButton("gQC_renameLevButton", "apply")
             )
           )
         ),
         tabPanel(
           title = "WIND", width = 12, value = "gQC_wind",
           tags$p("from WIND package compute hierachy of projection"),
           br(),
           fluidRow(
             column(
               width = 6,
               selectInput("gQC_windProj", "Projection to modify", choices = c("notyet"), selected = defaultValue("gQC_windProj", "notyet"))
             )
           ),
           fluidRow(
             column(
               width = 12,
               plotOutput("gQC_windHC")
             )
           )
         )
  )
)

tabList <- list(
  shinydashboard::tabItem(
    "gQC_umiHist",
    tags$h3("Histogram of UMI counts"),
    numericInput("gQC_binSize", "number of bins", 20, min=10),
    fluidRow(column(
      10,
      offset = 0,
      jqui_resizable(plotly::plotlyOutput("gQC_plotUmiHist"))
    )),
    br(),
    actionButton("save2Histumi", "save to history")
    
  ),
  
  shinydashboard::tabItem(
    "gQC_sampleHist",
    tags$h3("Histogram of cells per sample"),
    fluidRow(column(
      10,
      offset = 0,
      jqui_resizable(plotly::plotlyOutput("gQC_plotSampleHist"))
    )),
    br(),
    actionButton("save2HistSample", "save to history")
  ),
  
  shinydashboard::tabItem(
    "gQC_variancePC",
    tags$h3("Variance of PCs"),
    fluidRow(column(
      10,
      offset = 0,
      plotOutput("gQC_variancePCA") %>% withSpinner()
    )),
    br(),
    actionButton("save2Histvar", "save to history")
  ),
  
  tsnePlotTab = shinydashboard::tabItem(
    tabName = "gQC_tsnePlot",
    shinyjs::useShinyjs(),
    fluidRow(div(h3("TSNE Plot"), align = "center")),
    br(),
    box(
      title = "tSNE  parameters", solidHeader = TRUE, width = 12, status = "primary",
      fluidRow(
        column(
          width = 6,
          numericInput("gQC_tsneDim", "Tsne dimensions", 3, min = 3, max = 3)
        ),
        column(
          width = 6,
          numericInput("gQC_tsnePerplexity", "Perplexity", defaultValue("gQC_tsnePerplexity", 30), min = 1, max = 100)
        )
      ),
      box(
        title = "tSNE additional parameters", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        column(
          width = 6,
          numericInput("gQC_tsneTheta", "Theta", defaultValue("gQC_tsneTheta", 0.5), min = 0.0, max = 1, step = 0.1)
        ),
        column(
          width = 6,
          numericInput("gQC_tsneSeed", "Seed", defaultValue("gQC_tsneSeed", 1), min = 1, max = 10000)
        )
      ),
      fluidRow(
        column(
          width = 12, offset = 1,
          # uiOutput("updatetsneParametersButton")
          actionButton("updatetsneParameters", "apply changes", width = "80%")
          # ,
          #              style = "color: #fff; background-color: #A00272; border-color: #2e6da4")
        )
      )
    ),
    box(
      title = "3D plot", solidHeader = TRUE, width = 12, status = "primary",
      collapsible = TRUE, collapsed = FALSE,
      fluidRow(
        column(
          width = 3,
          selectInput("gQC_dim3D_x",
                      label = "X",
                      choices = c(defaultValue("gQC_dim3D_x", "tsne1"), "tsne2", "tsne3"),
                      selected = defaultValue("gQC_dim3D_x", "tsne1")
          )
        ),
        column(
          width = 3,
          selectInput("gQC_dim3D_y",
                      label = "Y",
                      choices = c("tsne1", defaultValue("gQC_dim3D_y", "tsne2"), "tsne3"),
                      selected = defaultValue("gQC_dim3D_y", "tsne2")
          )
        ),
        column(
          width = 3,
          selectInput("gQC_dim3D_z",
                      label = "Z",
                      choices = c("tsne1", "tsne2", defaultValue("gQC_dim3D_z", "tsne3")),
                      selected = defaultValue("gQC_dim3D_z", "tsne3")
          )
        ),
        column(
          width = 3,
          selectInput("gQC_col3D",
                      label = "colored by",
                      choices = defaultValue("gQC_col3D", "sampleNames"),
                      selected = defaultValue("gQC_col3D", "sampleNames")
          )
        )
      ),
      fluidRow(column(
        width = 12,
        jqui_resizable(plotly::plotlyOutput("gQC_tsne_main"))
      ))
    ),
    br(),
    box(
      title = "Table with all projections", solidHeader = TRUE, width = 12, status = "primary",
      collapsible = FALSE, collapsed = FALSE,
      fluidRow(column(
        width = 12,
        tableSelectionUi("gQC_projectionTableMod")
      ))
    )
  ),
  umapTab <- shinydashboard::tabItem(
    tabName = "gQC_umapPlot",
    box(
      title = "UMAP parameters", solidHeader = TRUE, width = 12, status = "primary",
      fluidRow(
        column(
          width = 12, offset = 1,
          shinyWidgets::actionBttn("activateUMAP", "apply changes", style = "bordered")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 3,
            selectInput("gQC_um_n_neighbors",
                      label = "N Neighbors",
                      choices = c(2:100), selected = defaultValue("gQC_um_n_neighbors", "15")
          )
        ),
        column(
          width = 3,
          selectInput("gQC_um_n_components",
                      label = "N components",
                      choices = c(2:20), selected = defaultValue("gQC_um_n_components", "2")
          )
        ),
        column(
          width = 3,
          selectInput("gQC_um_spread",
                      label = "spread",
                      choices = c(1:10), selected = defaultValue("gQC_um_spread", "1")
          )
        ),
        column(
          width = 3,
          selectInput("gQC_um_local_connectivity",
                      label = "local connectivity",
                      choices = 1:20, selected = defaultValue("gQC_um_local_connectivity", "1")
          )
        )
      ),
      box(
        title = "Addition UMAP options", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        fluidRow(
          column(
            width = 3,
            selectInput("gQC_um_randSeed",
                        label = "random seed",
                        choices = c(1:100), selected = defaultValue("gQC_um_randSeed", "1")
            ),
            selectInput("gQC_um_init",
                        label = "init",
                        choices = c("spectral", "random"), selected = defaultValue("gQC_um_init", "spectral")
            )
          ),
          column(
            width = 3,
            selectInput(
              "gQC_um_negative_sample_rate",
              label = "negative sample rate",
              choices = c(1:50), selected = defaultValue("gQC_um_negative_sample_rate", "5")
            ),
            selectInput("gQC_um_min_dist",
                        label = "min dist",
                        choices = seq(0.05, 0.5, 0.01), selected = defaultValue("gQC_um_min_dist", "0.01")
            )
          ),
          column(
            width = 3,
            selectInput("gQC_um_metric",
                        label = "metric",
                        choices = c("euclidean", "manhattan", "cosine", "hamming"),
                        selected = defaultValue("gQC_um_metric", "euclidean")
            ),
            selectInput("gQC_um_set_op_mix_ratio",
                        label = "set op mix ratio",
                        choices = seq(0, 1, 0.1), selected = defaultValue("gQC_um_set_op_mix_ratio", "1")
            )
          ),
          column(
            width = 3,
            selectInput("gQC_um_n_epochs",
                        label = "epochs",
                        choices = c(1:1000), selected = defaultValue("gQC_um_n_epochs", "200")
            ),
            selectInput(
              "gQC_um_bandwidth",
              label = "bandwidth",
              choices = c(1:20), selected = defaultValue("gQC_um_bandwidth", "1")
            )
          ),
        )
      ), # additional options box
    ),
    box(
      width = 12,
      fluidRow(column(
        width = 12,
        clusterUI("gQC_umap_main")
      ))
    )
  ),
  modTab
  
)
