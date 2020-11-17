suppressMessages(require(shiny))
# require(shinyMCE)

source(paste0(packagePath, "/modulesUI.R"), local = TRUE)
# this is where the general tabs are defined:

# localContributionDir <- get(".SCHNAPPs_locContributionDir", envir = .schnappsEnv)
# defaultValueSingleGene <- get(".SCHNAPPs_defaultValueSingleGene", envir = .schnappsEnv)
# defaultValueMultiGenes <- get(".SCHNAPPs_defaultValueMultiGenes", envir = .schnappsEnv)
# defaultValueRegExGene <- get(".SCHNAPPs_defaultValueRegExGene", envir = .schnappsEnv)
# DEBUG <- get(".SCHNAPPs_DEBUG", envir = .schnappsEnv)
# DEBUGSAVE <- get(".SCHNAPPs_DEBUGSAVE", envir = .schnappsEnv)

source(paste0(packagePath, "/toolTips.R"), local = TRUE)


# inputTab ----
inputTab <- function() {
  shinydashboard::tabItem(
    tabName = "input",
    box(
      width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
      fluidRow(
        div(h3("SCHNAPPs Input"), align = "center")
      ),
      br(),
      fluidRow(div(
        h4(
          "Single Cell sHiNy APP(s)"
        ),
        align = "center"
      )),
      box(
        width = 12, solidHeader = FALSE, collapsible = FALSE, collapsed = FALSE,
        fluidRow(
          column(
            width = 12,
            div(
              p(
                "Shiny app for the exploration and analysis of single cell RNAseq data as it comes from 10X or MARSseq technologies or other. It is currently being developed based on user requests of the Cytometry and Biomarkers UTechS at the Institut Pasteur, Paris. The goal is to enable the users of our platform to explore their data, select cells they would like to work with and then perform the final analysis together with the bioinformatics support at Pasteur. We hope you might find it helpful as well."
              ),
              align = "left"
            )
          )
        )
      )
    ),
    br(),
    boxPlus(
      title = "Input files",
      # helpID = ,
      dropdown_icon = NULL,
      closable = FALSE,
      enable_dropdown = T,
      dropdown_menu = actionButton(inputId = "inputHelp", label = "", icon = icon("fas fa-question")),
      status = "primary", solidHeader = TRUE, width = 12,
      footer = div(
        "Choose one or more .RData/.Rds file with singleCellExperiment object OR one .txt/.csv file with count data to upload",
        br(),
        "Multiple files can be selected when using RData files with SingleCellExpression objects.",
        br(),
        "RData files are R data files generated using base::save()."
      ),
      
      fluidRow(
        column(
          6,
          offset = 3,
          fileInput(
            "file1",
            "Count data upload",
            accept = c(
              ".Rds", ".RData", ".Rdata", ".txt", ".csv"
            ),
            placeholder = "no file selected",
            multiple = TRUE,
          ) %>% setId(id="fileInput"), checkbsTT("fileInput")
        ),
        boxPlus(
          title = "Additional annotations", status = "primary", solidHeader = TRUE, width = 12,
          # helpID = "inputHelpAdd",
          closable = FALSE,
          dropdown_icon = NULL,
          enable_dropdown = T,
          dropdown_menu = actionButton(inputId = "inputHelpAdd", label = "", icon = icon("fas fa-question")),
          footer = "(Not required): Choose .CSV file with annotation to upload. This can be projections or gene based features.",
          collapsible = TRUE, collapsed = TRUE,
          fluidRow(
            column(
              6,
              offset = 3,
              fileInput(
                "annoFile",
                "Annotations to add",
                accept = c(
                  ".txt", ".csv", ".mtx"
                ),
                placeholder = "no file selected",
                multiple = TRUE
              ) %>% setId(id="fileInputAnnotation"), checkbsTT("annoFile")
            )
          )
          
        )
      )
    ),
    
    br(),
    box(
      title = "input options", solidHeader = TRUE, width = 12, status = "primary",
      fluidRow(
        column(
          6,
          
          checkboxInput("sampleInput", label = "sub sample", value = defaultValue("sampleInput", TRUE)),
          
          numericInput("subsampleNum",
                       label = "max number of cells",
                       min = 500, max = 10000, step = 100, value = defaultValue("subsampleNum", 1000)
          )
        ),
        
        checkbsTT("sampleInput"),
        column(
          6,
          
          radioButtons("whichscLog",
                       label = "Compute normalizations?",
                       choices = c(
                         "disable log" = "disablescEx_log",
                         "use scEx from loaded data" = "useLog",
                         "calculate normalization here" = "calcLog"
                       ),
                       selected = defaultValue("whichscLog", "disablescEx_log")
          )
          # checkboxInput("disablescEx_log", label = "disable Normalization", value = TRUE)
        ),
        checkbsTT("disablescEx_log")
      )
    ),
    
    br(),
    box(
      title = "before-filter counts", width = 12, solidHeader = TRUE, status = "primary",
      footer = "This regular expression will be used before filtering out genes. It is meant to keep track of genes that were removed from gene filtering. This will generate a projection called 'before.filter'.",
      fluidRow(column(
        6,
        
        textInput("beforeFilterRegEx", "regular expression to count genes/cell", value = "^MT-")
      )), checkbsTT("beforeFilterRegEx")
    )
  )
}

# geneSelectionTab ----
geneSelectionTab <- function() {
  shinydashboard::tabItem(
    tabName = "geneSelection",
    fluidRow(div(h3("Gene selection"), align = "center")),
    fluidRow(
      column(
        width = 12, offset = 1,
        actionButton("updateGeneSelectionParameters", "apply changes", width = "80%")
      )
    ), checkbsTT("minGenesGS"),
    checkbsTT("updateGeneSelectionParameters"),
    br(),
    box(
      title = "Gene selection parameters", solidHeader = TRUE, width = 12, status = "primary",
      fluidRow(
        column(
          width = 4,
          textInput("selectIds", "regular expression for selection of genes to be removed", value = defaultValue("selectIds", "^MT-|^RP|^MRP"))
        ), checkbsTT("selectIds"),
        # column(
        #   width = 4,
        #   h4("GeneList Selection"),
        #   shinyTree::shinyTree("geneListSelection", checkbox = TRUE)
        # ), checkbsTT("geneListSelection"),
        column(
          width = 4,
          h4("Min expression over all cells"),
          numericInput("minGenesGS", "Min # of UMIs over all cells", defaultValue("minGenesGS", 2), min = 1, max = 1000000)
        )
      ), checkbsTT("minGenesGS"),
      fluidRow(
        column(
          width = 8, align = "center", offset = 2,
          textInput("genesKeep", "genes to keep", value = defaultValue("genesKeep", ""))
        )
      )
    ), checkbsTT("genesKeep"),
    
    br(),
    fluidRow(
      tabBox(
        title = "Gene selection tables", width = 12, id = "geneselectiontb",
        tabPanel("Genes kept",
                 height = "250px", width = 12, value = "Genes kept",
                 # box(
                 #   title = "Genes kept, with mean Expression, and number of cells expressing min 1", solidHeader = TRUE, width = 12, status = "primary",
                 #   collapsible = FALSE, collapsed = TRUE,
                 fluidRow(
                   column(
                     12,
                     tableSelectionUi("gsSelectedGenesMod")
                   )
                   # ), checkbsTT("gsSelectedGenesMod")
                 )
        ),
        tabPanel("genes removed",
                 height = "250px", value = "genes removed",
                 # box(
                 #   title = "Genes removed, with mean Expression, and number of cells expressing min 1", solidHeader = TRUE, width = 12, status = "primary",
                 #   collapsible = FALSE, collapsed = TRUE,
                 fluidRow(
                   column(
                     12,
                     tableSelectionUi("gsRMGenesMod")
                   ), checkbsTT("gsRMGenesMod")
                 )
                 # )
        )
      )
    )
  )
}

# cellSelectionTab ----
cellSelectionTab <- function() {
  shinydashboard::tabItem(
    tabName = "cellSelection",
    fluidRow(div(h3("Cell selection"), align = "center")),
    br(),
    fluidRow(
      column(
        width = 12, offset = 1,
        actionButton("updateCellSelectionParameters", "apply changes", width = "80%")
      )
    ),
    checkbsTT("updateCellSelectionParameters"),
    br(),
    box(
      title = "Cell selection parameters", solidHeader = TRUE, width = 12, status = "primary",
      fluidRow(
        column(
          width = 6,
          textInput("minExpGenes", "List of genes with minimal expression", value = defaultValue("minExpGenes", defaultValueRegExGene)), # tool tip: '^CD7$|^KIT$
          textInput("minNonExpGenes", "List of genes that should not be expressed", value = defaultValue("minNonExpGenes",""))
        ),
        column(
          width = 6,
          numericInput("minGenes", "Min # of UMIs", defaultValue("minGenes", 2), min = 2, max = 1000000),
          numericInput("maxGenes", "Max # of UMIs", defaultValue("maxGenes", 1000000), min = 10, max = 1000000)
        )
      )
    ),
    box(
      title = "addition parameters", solidHeader = TRUE, width = 12, status = "primary",
      collapsible = TRUE, collapsed = TRUE,
      fluidRow(
        column(
          width = 12,
          textInput("cellPatternRM", "cells to be filtered out by pattern", value = defaultValue("cellPatternRM", "")),
          textInput("cellKeep", "cells to keep", value = defaultValue("cellKeep", "")),
          textInput("cellKeepOnly", "cells to keep (remove others)", value = defaultValue("cellKeepOnly", "")),
          textInput("cellsFiltersOut", "Cells to be removed", width = "100%", value = defaultValue("cellsFiltersOut", "")),
          textInput("cellSelectionComment", "Comment for selection of cells", value = defaultValue("cellSelectionComment", ""))
        )
      )
    ),
    checkbsTT("minExpGenes"),
    checkbsTT("minNonExpGenes"),
    checkbsTT("minGenes"),
    checkbsTT("maxGenes"),
    checkbsTT("cellSelectionComment"),
    checkbsTT("cellPatternRM"),
    checkbsTT("cellKeep"),
    checkbsTT("cellKeepOnly"),
    checkbsTT("cellsFiltersOut"),
    br(),
    box(
      title = "cell table", solidHeader = TRUE, width = 12, status = "primary",
      collapsible = FALSE, collapsed = TRUE,
      tableSelectionUi("cellSelectionMod")
    ),
    checkbsTT("cellSelectionMod"), br()
  )
}

# parse all parameters.R files under contributions to include in application
# allTabs holds all tabs regardsless of their location in the GUI
getparameterContributions <- function() {
  parameterContributions <- list()
  # localContributionDir <- .SCHNAPPs_locContributionDir
  parFiles <- dir(path = c(paste0(packagePath, "/contributions"), localContributionDir), pattern = "parameters.R", full.names = TRUE, recursive = TRUE)
  for (fp in parFiles) {
    myPparameters <- list()
    source(fp, local = TRUE)
    if (length(myPparameters) > 0) {
      for (li in myPparameters) {
        if (length(li) > 0) {
          if (DEBUG) cat(file = stderr(), paste(li$children[[1]], "\n"))
          parameterContributions[[length(parameterContributions) + 1]] <- li
        }
      }
    }
  }
  return(parameterContributions)
}

# submenu items for the paramters main tab

parameterItems <- function() {
  list(
    shinydashboard::menuSubItem("Normalization", tabName = "normalizations"),
    getparameterContributions(),
    shinydashboard::menuSubItem("General Parameters", tabName = "generalParameters"),
    shinydashboard::menuSubItem("TSNE plot", tabName = "gQC_tsnePlot"),
    shinydashboard::menuSubItem("Umap", tabName = "gQC_umapPlot"),
    shinydashboard::menuSubItem("Projections", tabName = "modifyProj")
  )
}

# generalParametersTab ----
generalParametersTab <- function() {
  shinydashboard::tabItem(
    "generalParameters",
    fluidRow(div(h2("General parameters"), align = "center")),
    br(),
    fluidRow(
      tabBox(title = "PCA", width = 12, id = "modPCA",
             tabPanel(
               title = "Parameters for PCA", solidHeader = TRUE, width = 12, value = "PCAparameters",
               # The id lets us use input$tabset1 on the server to find the current tab
               id = "tabsetPCA",
               fluidRow(
                 column(4,
                        offset = 0,
                        numericInput("pcaRank", "Number of components", defaultValue("pcaRank", 50), min = 2),
                        checkboxInput("pcaCenter", "center data", TRUE)
                 ),
                 column(4,
                        offset = 0,
                        numericInput("pcaN", "Number of variable genes to be used", defaultValue("pcaN", 500), min = 50),
                        checkboxInput("pcaScale", "scale data", defaultValue("pcaScale", TRUE))
                 ),
                 column(4,
                        offset = 0,
                        selectInput("hvgSelection","How to select highly variable genes.", 
                                    choices = c("getTopHVGs","vst", "mvp", "disp"),
                                    selected = defaultValue("hvgSelection", "vst")),
                        checkboxInput("useSeuratPCA", "use Seurat::RunPCA", defaultValue("useSeuratPCA", TRUE))
                 )
               ),
               checkbsTT(item = "pcaRank"),
               checkbsTT("pcaN"),
               checkbsTT("pcaCenter"),
               checkbsTT("pcaScale"),
               fluidRow(
                 column(12,
                        offset = 0,
                        textInput("genes4PCA", "Genes to be used for PCA", width = "100%", value = defaultValue("genes4PCA",""))
                 ), checkbsTT("genes4PCA")
               ),
               fluidRow(
                 column(12,
                        offset = 0,
                        textInput("genesRMPCA", "Genes NOT to be used for PCA", width = "100%", value = defaultValue("genesRMPCA", ""))
                 ), checkbsTT("genesRMPCA")
               ),
               fluidRow(
                 column(
                   width = 12, offset = 1,
                   actionButton("updatePCAParameters", "apply changes", width = "80%")
                 )
               )
             ),
             checkbsTT(item = "tabsetPCA"),
             
             tabPanel(
               title = "DimPlot for PCA", solidHeader = TRUE, width = 12, value = "dimPlotPCA",
               fluidRow(
                 column(12,
                        offset = 1,
                        actionButton("updateDimPlot", "generate plot",
                                     width = "80%",
                                     style = "color: #fff; background-color: #A00272; border-color: #2e6da4"
                        )
                 ),
               ),
               fluidRow(
                 column(
                   width = 12,
                   jqui_resizable(plotOutput("dimPlotPCA", height = "1400px"))
                 )
               ),
               checkbsTT(item = "dimPlotPCA"),
             ),
             tabPanel(
               title = "Loadings", solidHeader = TRUE, width = 12, value = "loadingsPCA",
               fluidRow(
                 column(12,
                        tableSelectionUi("PCAloadingsMod")
                 )
               )
             )
      ),
    ),
    fluidRow(
      tabBox(
        title = "Parameters for clustering", width = 12,
        id = "tabsetCluster",
        selected = defaultValue("tabsetCluster", "seurat_Clustering"),
        tabPanel("Seurat clustering",
                 width = 12,
                 value = "seurat_Clustering",
                 fluidRow(
                   column(
                     width = 6,
                     numericInput("seurClustDims", "Dimensions of PCA to use", min = 5, value = defaultValue("seurClustDims", 15), width = "100%"),
                     numericInput("seurClustk.param", "K for k-nearest neighbor algorithm", min = 20, value = defaultValue("seurClustk.param",15), width = "100%")
                   ),
                   column(
                     width = 6,
                     numericInput("seurClustresolution", "Value of the resolution parameter (below 1 -> smaller communities)", value = defaultValue("seurClustresolution", 0.5), min = 0.1, width = "100%"),
                     # TOD implement more options
                   )
                 ),
        ), # seurat clustering
        tabPanel("Quickcluster",
                 value = "scran_Cluster", # name of reactive to be used
                 width = 12,
                 fluidRow(
                   column(
                     width = 6,
                     selectInput("clusterSource", "use raw counts or normalized data?", choices = c("counts", "logcounts"), selected = defaultValue("clusterSource", "logcounts"), width = "100%"),
                     selectInput("clusterMethod", "clustering method to use", choices = c("hclust", "igraph"), selected = defaultValue("clusterMethod", "igraph"), width = "100%")
                   ),
                   column(
                     width = 6,
                     numericInput("minClusterSize", "minimum size of each cluster.", defaultValue("minClusterSize", 2), min = 2, width = "100%"),
                     selectInput("useRanks", "use ranks?\n", choices = c("TRUE", "FALSE"), selected = defaultValue("useRanks", "TRUE"), width = "100%")
                   )
                 ),
                 checkbsTT(item = "clusterSource"),
                 checkbsTT(item = "minClusterSize"),
                 checkbsTT(item = "clusterMethod"),
                 checkbsTT(item = "useRanks"),
                 fluidRow(
                   column(
                     12,
                     textInput("geneSelectionClustering", "Genes to be used for clustering", width = "100%", value = defaultValue("geneSelectionClustering", ""))
                   )
                 ),
                 checkbsTT(item = "geneSelectionClustering")
                 
        ), # quickclustering tab Panel
        tabPanel("SNNGraph",
                 value = "snnGraph",
                 width = 12,
                 fluidRow(
                   # column(
                   #   width = 4,
                   #   selectInput("snnClusterSource", "use raw counts or normalized data?", choices = c("counts", "logcounts"), selected = defaultValue("snnClusterSource", "logcounts"), width = "100%"),
                   # ),
                   column(
                     width = 4,
                     selectInput("snnType", "type to use", 
                                 selected = defaultValue("snnType", "rank"),
                                 choices = c("rank", "number", "jaccard"))
                   )
                 )),
        tabPanel("SMILR",
                 value = "simlrFunc",
                 width = 12,
                 fluidRow(
                   # column(
                   #   width = 4,
                   #   selectInput("snnClusterSource", "use raw counts or normalized data?", choices = c("counts", "logcounts"), selected = defaultValue("snnClusterSource", "logcounts"), width = "100%"),
                   # ),
                   column(
                     width = 4,
                     numericInput("smilr_nClust", "number of clusters (0 = estimate)", 
                                  value = defaultValue("smilr_nClust", 10),
                                  min = 0, max = 1000)
                   ),
                   column(
                     width = 4,
                     numericInput("smilr_maxClust", "max number of clusters when estimating)", 
                                  value = defaultValue("smilr_maxClust", 20),
                                  min = 2, max = 1000)
                   )
                 )),
        fluidRow(
          column(12, offset = 0, textOutput("Nclusters"))
        ),
        # checkbsTT(item = ""),
        fluidRow(
          column(
            width = 12, offset = 1,
            actionButton("updateClusteringParameters", "apply changes", width = "80%")
          )
        )
      )
    ),
    # fluidRow(div(h3("Parameters for clustering"), align = "left")),
    
    br(),
    # box(
    #   title = "Comments", solidHeader = TRUE, width = 12, status = "primary",
    #   collapsible = TRUE, collapsed = TRUE,
    #   if ("shinyMCE" %in% rownames(installed.packages())) {
    #     shinyMCE::tinyMCE(
    #       "descriptionOfWork",
    #       "Please describe your work. This will be included in the report."
    #     )
    #   } else {
    #     textInput("descriptionOfWork", "Please describe your work. This will be included in the report.")
    #   }
    # ),
    # checkbsTT(item = "descriptionOfWork"),
    # br(),
    box(
      title = "Colors", solidHeader = TRUE, width = 12, status = "primary",
      collapsible = TRUE, collapsed = TRUE,
      fluidRow(column(
        width = 12, offset = 1,
        actionButton("updateColors", "apply changes", width = "80%")
      )),
      br(),
      fluidRow(
        column(
          width = 6,
          uiOutput("sampleColorSelection")
        ),
        column(
          width = 6,
          uiOutput("clusterColorSelection")
        )
      )
    ),
    checkbsTT(item = "updateColors"),
    checkbsTT(item = "sampleColorSelection"),
    checkbsTT(item = "clusterColorSelection")
  )
}

# # renameTab
# renameTab <- function() {
#   shinydashboard::tabItem(
#     tabName = "renameProj",
#     fluidRow(div(h3("rename projections"), align = "center")),
#     br(),
#     box(
#       title = "Rename projections", solidHeader = TRUE, width = 12, status = "primary",
#       fluidRow(
#         column(
#           width = 6,
#           selectInput("oldPrj", "projections to copy + rename", choices = c("notyet"), selected = "notyet")
#         ),
#         column(
#           width = 6,
#           fluidRow(
#             column(
#               width = 8,
#               textInput("newPrj", "new name of Projection", value = "")
#             ),
#             column(
#               width = 4,
#               actionButton("updatePrjsButton", "rename")
#             )
#           ),
#           fluidRow(
#             column(
#               width = 8,
#               selectInput("delPrj", "projections to delete", choices = c("notyet"), selected = "notyet")
#             ),
#             column(
#               width = 4,
#               actionButton("delPrjsButton", "delete")
#             ),
#             tags$style(type = "text/css", "#updatePrjsButton { width:100%; margin-top: 25px;}"),
#             tags$style(type = "text/css", "#delPrjsButton { width:100%; margin-top: 25px;}")
#           )
#         )
#       ),
#       checkbsTT(item = "oldPrj"),
#       checkbsTT(item = "newPrj"),
#       checkbsTT(item = "updatePrjsButton"),
#       checkbsTT(item = "delPrj"),
#       checkbsTT(item = "delPrjsButton")
#     ),
# 
# 
#     checkbsTT(item = "updateColors"),
#     checkbsTT(item = "sampleColorSelection"),
#     checkbsTT(item = "clusterColorSelection")
#   )
# }

# # link to the content of the
# parametersTab  = tabItem(tabName = "normalizations",
#                             fluidRow(div(h3('Cell selection'), align = 'center')),
#                             br()
# )
if (DEBUG) {
  cat(file = stderr(), paste("end: tabs.R\n"))
}
