# gQC_generalQC/ui.R

suppressMessages(require(magrittr))
source(paste0(packagePath, "/modulesUI.R"), local = TRUE)
menuList <- list(
  shinydashboard::menuItem("General QC",
                           # id="generalQCID",
                           tabName = "generalQC", icon = icon("thumbs-up"), startExpanded = FALSE,
                           shinydashboard::menuSubItem("UMI histogram", tabName = "gQC_umiHist"),
                           shinydashboard::menuSubItem("Sample histogram", tabName = "gQC_sampleHist"),
                           shinydashboard::menuSubItem("PC variance", tabName = "gQC_variancePC"),
                           
                           if("DoubletFinder" %in% installed.packages()){
                             shinydashboard::menuSubItem("Doublet Finder", tabName = "gQC_doubletFinder")},
                           shinydashboard::menuSubItem("Scater QC", tabName = "DE_scaterQC")
                           # shinydashboard::menuSubItem("TSNE plot", tabName = "gQC_tsnePlot"),
                           # shinydashboard::menuSubItem("Umap", tabName = "gQC_umapPlot")
  )
)

# gQC_doubletFinder ----

if("DoubletFinder" %in% installed.packages()){
  
  DoubletFinderTab <- shinydashboard::tabItem(
    tabName = "gQC_doubletFinder",
    fluidRow(div(h3("DoubletFinder"), align = "center")),
    br(),
    tabBox(title = "", width = 12, id = "doubletFinderBox",
           tabPanel(
             title = "Doublet Finder", solidHeader = TRUE, width = 12, value = "doubletFinderTab",
             id = "doublet.Finder.Tab",
             fluidRow(
               column(
                 width = 6,
                 sc_numericInput("GS_DF_dims", "PC imensions to use", min = 0, max = 200, step = 1, value = defaultValue("GS_DF_dims", 20))
               ),
               column(
                 width = 6,
                 sc_numericInput("GS_DF_nRecover", "number of cells recovered by CellRanger (before QC)", min = 100, max = 20000000, step = 1000, value = defaultValue("GS_DF_nRecover", 10000))
               )
             ),
             fluidRow(
               column(
                 width = 6,
                 sc_numericInput("GS_DF_pk", "neighborhood size to estimate likelihood of being doublet \n
                             (expressed as fraction of merged real-artifical data)", min = 0, max = Inf, step = 1, value = defaultValue("GS_DF_pk", 20))
               ),
               column(
                 width = 6,
                 sc_numericInput("GS_DF_pN", "number of artificial doublets generated \n
                             (expressed as fraction of merged real-artifical data)", min = 0, max = 1, step = 0.01, value = defaultValue("GS_DF_pN", 0.10))
               )
             ),
             br(),
             fluidRow(
               column(
                 width = 12, offset = 1,
                 actionButton("GS_DF_button", "apply changes", width = "80%")
               )
             ),
             fluidRow(
               column(
                 width = 12,
                 clusterUI("GS_DF_plot")
               )
             )
           )
    )
  )
}
# geneSetModTab -----

geneSetModTab <- shinydashboard::tabItem(
  tabName = "modifyGeneSets",
  fluidRow(div(h3("work with gene sets"), align = "center")),
  br(),
  tabBox(title = "modify gene sets", width = 12, id = "modGS",
         tabPanel(
           title = "Rename gene set", solidHeader = TRUE, width = 12, value = "renameGS",
           id = "rename.GS.Tab",
           fluidRow(
             column(
               width = 6,
               sc_selectInput("oldGS", "gene sets to copy + rename", choices = defaultValue("oldGS", "notyet"), selected = defaultValue("oldGS", "notyet")),
               verbatimTextOutput("gQC_renameGenes") %>% jqui_resizable()
             ),
             column(
               width = 6,
               fluidRow(
                 column(
                   width = 8,
                   sc_textInput("newGS", "new name of gene set", value = defaultValue("newGS", ""))
                 ),
                 column(
                   width = 4,
                   actionButton("updateGSButton", "rename")
                 )
               ),
               fluidRow(
                 column(
                   width = 8,
                   sc_selectInput("delGS", "gene set to delete", choices = defaultValue("delGS", "notyet"), selected = defaultValue("delGS", "notyet"))
                 ),
                 column(
                   width = 4,
                   actionButton("delGSButton", "delete")
                 ),
                 tags$style(type = "text/css", "#updateGSButton { width:100%; margin-top: 25px;}"),
                 tags$style(type = "text/css", "#delGSButton { width:100%; margin-top: 25px;}")
               )
             )
           )
         ),
         tabPanel(
           title = "edit gene set", solidHeader = TRUE, width = 12, value = "gQC_editGS",
           id = "gQC_editGS.Tab",
           tags$p("Select a gene set and then edit it"),
           fluidRow(
             column(
               width = 12,
               sc_selectizeInput(inputId="gQC_geneSetModifyInput",
                                 label = "Gene set to modify",
                                 multiple = FALSE,
                                 choices = c(defaultValue("gQC_geneSetModifyInput", "dummy")),
                                 selected = defaultValue("gQC_geneSetModifyInput", "dummy"), 
                                 width = "100%",
                                 options = list(maxItems = 1))
               ,
               verbatimTextOutput("gQC_geneSetModifyInputGL") %>% jqui_resizable(),
               sc_textInput(inputId= "gQC_geneSetModifynName", label="name of new gene set", value="newGeneSet"),
               sc_textInput(inputId= "gQC_geneSetModifynDesc", label="description of new gene set", value="newGeneSet"),
               sc_textInput(inputId= "gQC_geneSetModifyGenes", label="comma separated list of genes", value="") %>% jqui_resizable() ,
               actionButton("geneSetModifyButton", "rename")
             ))
         ),
         tabPanel(
           title = "search gene sets", solidHeader = TRUE, width = 12, value = "gQC_search",
           id = "gQC_search.Tab",
           fluidRow(
             column(
               width = 6,
               sc_textInput(inputId= "gQC_genesets_search", label="comma separated list of genes", value="") %>% jqui_resizable(),
               
               #,
               # actionButton("button_GSEA", "run GSEA")
             ),
             column(width = 6,
                    sc_numericInput(inputId = "gQC_genesets_maxPrint", label="max items to print",
                                    value=100, min=10, max=1000, step=10))),
           fluidRow(
             column(
               width=12,
               verbatimTextOutput("gQC_geneSetsearchOutput") %>% jqui_resizable()
             )
           )
         )
  )
)


# modTab ----
modTab <- 
  shinydashboard::tabItem(
    "modifyProj",
    fluidRow(div(h3("work with projections"), align = "center")),
    br(),
    tabBox(title = "modify projections", width = 12, id = "modProj",
           tabPanel(
             title = "Rename projections", solidHeader = TRUE, width = 12, value = "renameProj",
             id = "rename.Proj.Tab",
             fluidRow(
               column(
                 width = 6,
                 sc_selectInput("oldPrj", "projections to copy + rename", choices = defaultValue("oldPrj", "notyet"), selected = defaultValue("oldPrj", "notyet"))
               ),
               column(
                 width = 6,
                 fluidRow(
                   column(
                     width = 8,
                     sc_textInput("newPrj", "new name of Projection", value = defaultValue("newPrj", ""))
                   ),
                   column(
                     width = 4,
                     actionButton("updatePrjsButton", "rename")
                   )
                 ),
                 fluidRow(
                   column(
                     width = 8,
                     sc_selectInput("delPrj", "projections to delete", choices = defaultValue("delPrj", "notyet"), selected = defaultValue("delPrj", "notyet"))
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
             id = "combine.Proj.Tab",
             tags$p("Two factors can be combined by pasting the values per row and re-leveling"),
             br(),
             fluidRow(
               column(
                 width = 6,
                 sc_selectInput("gQC_combPrj1", "1 st Projection", choices = c("notyet"), selected = "notyet")
               ),
               column(
                 width = 6,
                 
                 sc_selectInput("gQC_combPrj2", "2nd Projections", choices = c("notyet"), selected = "notyet")
               )),
             fluidRow(
               column(width = 6,
                      sc_textInput("gQC_newCombPrj", "name of new Projection", value = "")),
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
             id = "rename.Levels.Tab",
             tags$p("rename the levels of a factor"),
             br(),
             
             fluidRow(
               column(
                 width = 6,
                 sc_selectInput("gQC_rnProj", "Projection to modify", choices = c("notyet"), selected = "notyet")
               ),
               column(width = 6,
                      sc_textInput("gQC_newRnPrj", "name of new Projection", value = ""))
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
             br(),
             fluidRow(
               column(
                 width = 6,
                 actionButton("gQC_renameLevButton", "apply")
               )
             )
           ),
           
           
           tabPanel(
             title = "rearrange levels", width = 12, value = "gQC_rearrangeLev",
             id = "rearrange.Levels.Tab",
             tags$p("rearrange the levels of a factor"),
             br(),
             fluidRow(
               column(
                 width = 6,
                 sc_selectInput("gQC_raProj", "Projection to modify", choices = c("notyet"), selected = defaultValue("gQC_raProj", "notyet"))
               ),
               column(width = 6,
                      sc_textInput("gQC_newRaPrj", "name of new Projection", value = ""))
             ),
             fluidRow(
               column(width = 12,
                      shinyjqui::orderInput('gQC_newRaLev', 'Rearrange levels', items = NULL, placeholder = 'choose projection...'),
               )),
             br(),
             fluidRow(
               column(
                 width = 6,
                 actionButton("gQC_rearrangeLevButton", "apply")
               )
             )
           ),
           
           
           tabPanel(
             title = "WIND", width = 12, value = "gQC_wind",
             id = "wind.Tab",
             tags$p("from WIND package compute hierachy of projection (https://github.com/haowulab/Wind)"),
             tags$p("to compute weighted normalized mutual information (wNMI) and weighted Rand index (wRI) to evaluate the clustering results by comparing a clustering output with a reference which has a hierarchical structure."),
             br(),
             fluidRow(
               column(
                 width = 6,
                 sc_selectInput("gQC_windProj", "Projection to use", choices = c("notyet"), selected = defaultValue("gQC_windProj", "notyet"))
               )
             ),
             fluidRow(
               column(
                 width = 12,
                 plotOutput("gQC_windHC") %>% jqui_resizable()
               )
             )
           )
    )
  )

tabList <- list(
  shinydashboard::tabItem(
    "gQC_umiHist",
    tags$h3("Histogram of UMI counts"),
    sc_numericInput("gQC_binSize", "number of bins", defaultValue("gQC_binSize",200), min=10),
    fluidRow(column(
      10,
      offset = 0,
      plotly::plotlyOutput("gQC_plotUmiHist") %>% jqui_resizable()
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
      plotly::plotlyOutput("gQC_plotSampleHist") %>% jqui_resizable()
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
      plotOutput("gQC_variancePCA") %>% jqui_resizable()
    )),
    br(),
    actionButton("save2Histvar", "save to history")
  ),
  
  shinydashboard::tabItem(
    tabName = "gQC_tsnePlot",
    # id = "gQC_tsneTab",
    # shinyjs::useShinyjs(),
    fluidRow(div(h3("tSNE Plot"), align = "center")),
    br(),
    shinydashboardPlus::box(
      title = "tSNE  parameters", solidHeader = TRUE, width = 12, status = "primary",
      fluidRow(
        column(
          width = 6,
          sc_numericInput("gQC_tsneDim", "Tsne dimensions", 3, min = 3, max = 3)
        ),
        column(
          width = 6,
          sc_numericInput("gQC_tsnePerplexity", "Perplexity", defaultValue("gQC_tsnePerplexity", 10), min = 1, max = 100)
        )
      ),
      shinydashboardPlus::box(
        title = "tSNE additional parameters", solidHeader = TRUE, width = 12, status = "primary",
        id = "tSNEParametersToggle",
        collapsible = TRUE, collapsed = TRUE,
        column(
          width = 6,
          sc_numericInput("gQC_tsneTheta", "Theta", defaultValue("gQC_tsneTheta", 0.5), min = 0.0, max = 1, step = 0.1)
        ),
        column(
          width = 6,
          sc_numericInput("gQC_tsneSeed", "Seed", defaultValue("gQC_tsneSeed", 1), min = 1, max = 10000)
        )
      ),
      checkbsTT(item = "gQC_tsneTheta"),
      checkbsTT(item = "gQC_tsnePerplexity"),
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
    shinydashboardPlus::box(
      title = "3D plot", solidHeader = TRUE, width = 12, status = "primary",
      collapsible = TRUE, collapsed = FALSE, id = "tsne3D",
      fluidRow(
        column(
          width = 3,
          sc_selectInput("gQC_dim3D_x",
                         label = "X",
                         choices = c(defaultValue("gQC_dim3D_x", "tsne1"), "tsne2", "tsne3"),
                         selected = defaultValue("gQC_dim3D_x", "tsne1")
          )
        ),
        column(
          width = 3,
          sc_selectInput("gQC_dim3D_y",
                         label = "Y",
                         choices = c("tsne1", defaultValue("gQC_dim3D_y", "tsne2"), "tsne3"),
                         selected = defaultValue("gQC_dim3D_y", "tsne2")
          )
        ),
        column(
          width = 3,
          sc_selectInput("gQC_dim3D_z",
                         label = "Z",
                         choices = c("tsne1", "tsne2", defaultValue("gQC_dim3D_z", "tsne3")),
                         selected = defaultValue("gQC_dim3D_z", "tsne3")
          )
        ),
        column(
          width = 3,
          sc_selectInput("gQC_col3D",
                         label = "color",
                         choices = defaultValue("gQC_col3D", "sampleNames"),
                         selected = defaultValue("gQC_col3D", "sampleNames")
          )
        )
      ),
      fluidRow(column(
        width = 12,
        plotly::plotlyOutput("gQC_tsne_main") %>% jqui_resizable()
      ))
    ),
    checkbsTT("gQC_tsnePerplexity"),
    checkbsTT("normalizationRadioButtonValue"),
    checkbsTT("updateNormalization"),
    
    br(),
    shinydashboard::box(
      title = "Table with all projections", solidHeader = TRUE, width = 12, status = "primary",
      collapsible = FALSE, collapsed = FALSE,
      fluidRow(column(
        width = 12,
        tableSelectionUi("gQC_projectionTableMod")
      ))
    )
  ),
  shinydashboard::tabItem(
    tabName = "gQC_umapPlot",
    
    shinydashboardPlus::box(
      title = "UMAP parameters", solidHeader = TRUE, width = 12, status = "primary",
      # id = "gQC_umapTab",
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
          sc_selectInput("gQC_um_n_neighbors",
                         label = "N Neighbors",
                         choices = c(2:100), selected = defaultValue("gQC_um_n_neighbors", "10")
          )
        ),
        column(
          width = 3,
          sc_selectInput("gQC_um_n_components",
                         label = "N components",
                         choices = c(2:20), selected = defaultValue("gQC_um_n_components", "3")
          )
        ),
        column(
          width = 3,
          sc_selectInput("gQC_um_spread",
                         label = "spread",
                         choices = c(1:10), selected = defaultValue("gQC_um_spread", "10")
          )
        ),
        column(
          width = 3,
          sc_selectInput("gQC_um_local_connectivity",
                         label = "local connectivity",
                         choices = 1:20, selected = defaultValue("gQC_um_local_connectivity", "5")
          )
        )
      ),
      shinydashboardPlus::box(
        title = "Addition UMAP options", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        id="addUMAPoptions",
        fluidRow(
          column(
            width = 3,
            sc_selectInput("gQC_um_randSeed",
                           label = "random seed",
                           choices = c(1:100), selected = defaultValue("gQC_um_randSeed", "1")
            ),
            sc_selectInput("gQC_um_init",
                           label = "init",
                           choices = c("spectral", "random"), selected = defaultValue("gQC_um_init", "spectral")
            )
          ),
          column(
            width = 3,
            sc_selectInput(
              "gQC_um_negative_sample_rate",
              label = "negative sample rate",
              choices = c(1:50), selected = defaultValue("gQC_um_negative_sample_rate", "5")
            ),
            sc_selectInput("gQC_um_min_dist",
                           label = "min dist",
                           choices = seq(0.05, 0.5, 0.01), selected = defaultValue("gQC_um_min_dist", "0.01")
            )
          ),
          column(
            width = 3,
            sc_selectInput("gQC_um_metric",
                           label = "metric",
                           choices = c("euclidean", "manhattan", "cosine", "hamming"),
                           selected = defaultValue("gQC_um_metric", "euclidean")
            ),
            sc_selectInput("gQC_um_set_op_mix_ratio",
                           label = "set op mix ratio",
                           choices = seq(0, 1, 0.1), selected = defaultValue("gQC_um_set_op_mix_ratio", "1")
            )
          ),
          column(
            width = 3,
            sc_numericInput(inputId = "gQC_um_n_epochs", label = "epochs:", 
                            value = as.numeric(defaultValue("gQC_um_n_epochs", 200)), min = 1, max = 1000),
            
            # selectizeInput("gQC_um_n_epochs",
            #             label = "epochs",
            #             choices = c(1:1000), selected = defaultValue("gQC_um_n_epochs", "200")
            # ),
            sc_selectInput(
              "gQC_um_bandwidth",
              label = "bandwidth",
              choices = c(1:20), selected = defaultValue("gQC_um_bandwidth", "1")
            ),
            checkbsTT("gQC_um_n_neighbors"),
            checkbsTT("gQC_um_spread"),
            checkbsTT("gQC_um_local_connectivity"),
            checkbsTT("gQC_um_init"),
            checkbsTT("gQC_um_negative_sample_rate"),
            checkbsTT("gQC_um_min_dist"),
            checkbsTT("gQC_um_metric"),
            checkbsTT("gQC_um_set_op_mix_ratio"),
            checkbsTT("gQC_um_bandwidth")
          )
        )
      ) # additional options box
    ),
    shinydashboardPlus::box(
      width = 12,
      fluidRow(column(
        width = 12,
        clusterUI("gQC_umap_main")
      ))
    )
  ),
  modTab,
  geneSetModTab,
  
  if("DoubletFinder" %in% installed.packages()){
    DoubletFinderTab
    }
  
)
