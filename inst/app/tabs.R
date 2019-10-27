suppressMessages(require(shiny))
# require(shinyMCE)

source(paste0(packagePath,  "/modulesUI.R"), local = TRUE)
# this is where the general tabs are defined:

# localContributionDir <- get(".SCHNAPPs_locContributionDir", envir = .schnappsEnv)
# defaultValueSingleGene <- get(".SCHNAPPs_defaultValueSingleGene", envir = .schnappsEnv)
# defaultValueMultiGenes <- get(".SCHNAPPs_defaultValueMultiGenes", envir = .schnappsEnv)
# defaultValueRegExGene <- get(".SCHNAPPs_defaultValueRegExGene", envir = .schnappsEnv)
# DEBUG <- get(".SCHNAPPs_DEBUG", envir = .schnappsEnv)
# DEBUGSAVE <- get(".SCHNAPPs_DEBUGSAVE", envir = .schnappsEnv)

source(paste0(packagePath,  "/toolTips.R"), local = TRUE)



# inputTab ----
inputTab <- shinydashboard::tabItem(
  tabName = "input",
  fluidRow(div(h3("SCHNAPPs Input"), align = "center")),
  br(),
  fluidRow(div(
    h5(
      "This app is designed for exploratory data analysis of processed RNA-Seq data of single cell experiments.<br/>
      Multiple files can be selected when using RData files with SingleCellExpression objects.<br/>
      RData files are R data files generated using base::save()."
    ),
    align = "center"
  )),
  fluidRow(column(
    8,
    offset = 1,
    fileInput(
      "annoFile",
      "(Not required): Choose .CSV file with annotation to upload",
      accept = c(
        ".txt",".csv", ".mtx"
      ),
      multiple = TRUE,
      width = '50%'
    ), checkbsTT("annoFile")
  )),
  br(),
  fluidRow(column(
    8,
    offset = 1,
    fileInput(
      "file1",
      "Choose one or more .RData/.Rds file with singleCellExperiment object OR one .txt/.csv file with count data to upload",
      accept = c(
        ".Rds",".RData", ".Rdata", ".txt", ".csv"
      ),
      multiple = TRUE,
      width = '50%'
    ), checkbsTT("file1")
  )),
  br(),
  fluidRow(column(
    4,
    offset = 1,
    checkboxInput("sampleInput", label = "sub sample", value = TRUE)),
    column(4,numericInput("subsampleNum", label = "max number of cells", 
                          min = 500, max = 10000, step = 100, value = 1000)
    )), checkbsTT("sampleInput"),
  fluidRow(column(
    4,
    offset = 1,
    radioButtons("whichscLog", label = "Compute normalizations?",
                choices = c("disable log" = "disablescEx_log" , 
                            "use scEx from loaded data" = "useLog" ,
                            "calculate normalization here" = "calcLog"),
                selected = "disablescEx_log")
    # checkboxInput("disablescEx_log", label = "disable Normalization", value = TRUE)
  )),checkbsTT("disablescEx_log"),
  
  br(),
  fluidRow(column(8,offset = 1,
                  textInput("beforeFilterRegEx", "regular expression to count genes/cell", value = "^MT-")
  )),checkbsTT("beforeFilterRegEx"),
  fluidRow(column(8, offset = 1, 
                  tags$div(
                    
                    tags$p("This regular expression will be used before filtering out genes.
               It is meant to keep track of genes that were removed from gene filtering. This will generate a projection
                           called 'before.filter'.")
                  )
  ))
)

# geneSelectionTab ----
geneSelectionTab <- shinydashboard::tabItem(
  tabName = "geneSelection",
  fluidRow(div(h3("Gene selection"), align = "center")),
  br(),
  fluidRow(div(
    h4(
      "Here we filter out genes"
    ),
    align = "center"
  )),
  fluidRow(
    column(3, offset = 1, 
           actionButton("updateGeneSelectionParameters", "apply changes"))),
  checkbsTT("updateGeneSelectionParameters"),
  fluidRow(
    column(3,
           offset = 1,
           textInput("selectIds", "regular expression for selection of genes to be removed", value = "^MT-|^RP|^MRP")
    ),checkbsTT("selectIds"),
    column(
      5,
      h4("GeneList Selection"),
      shinyTree::shinyTree("geneListSelection", checkbox = TRUE)
    ),checkbsTT("geneListSelection"),
    column(
      2,
      h4("Min expression over all cells"),
      numericInput("minGenesGS", "Min # of UMIs over all cells", 2, min = 2, max = 1000000)
    )
  ),checkbsTT("minGenesGS"),
  fluidRow(
    column(6,
           offset = 1,
           textInput("genesKeep", "genes to keep")
    )
  ),checkbsTT("genesKeep"),
  br(),
  fluidRow(
    h3("Genes kept, with mean Expression, and number of cells expressing min 1", align = "center"),
    br(),
    h4("Selected genes"),
    column(12,
           offset = 0,
           tableSelectionUi("gsSelectedGenesMod")
           # textOutput("gsSelectedGenes", inline = FALSE)
    )
  ),checkbsTT("gsSelectedGenesMod"),
  br(),
  fluidRow(
    h3("Genes removed, with mean Expression, and number of cells expressing min 1", align = "center"),
    h4("Selected genes"),
    br(),
    textOutput("gsrmGenes", inline = FALSE)
  ), checkbsTT("gsrmGenes"), br(),
  fluidRow(
    column(11,
           offset = 1,
           tableSelectionUi("gsRMGenesMod")
    ), checkbsTT("gsRMGenesMod")
  )
)


# cellSelectionTab ----
cellSelectionTab <- shinydashboard::tabItem(
  tabName = "cellSelection",
  fluidRow(div(h3("Cell selection"), align = "center")),
  br(),
  fluidRow(div(
    h4(
      "Here we filter out cells"
    ),
    align = "center"
  )),
  fluidRow(
    column(3, offset = 1, 
           actionButton("updateCellSelectionParameters", "apply changes"))),
  checkbsTT("updateCellSelectionParameters"),
  fluidRow(
    column(6,
           offset = 1,
           textInput("minExpGenes", "List of genes with minimal expression", value = defaultValueRegExGene) # tool tip: '^CD7$|^KIT$
    )
  ),checkbsTT("minExpGenes"),
  fluidRow(
    column(6,
           offset = 1,
           textInput("minNonExpGenes", "List of genes that should not be expressed", value = "")
    )
  ),checkbsTT("minNonExpGenes"),
  fluidRow(
    column(5,
           offset = 1,
           numericInput("minGenes", "Min # of UMIs", 2, min = 2, max = 1000000)
    ),
    column(
      5,
      numericInput("maxGenes", "Max # of UMIs", 1000000, min = 10, max = 1000000)
    )
  ), br(),checkbsTT("minGenes"),checkbsTT("maxGenes"),
  fluidRow(
    column(11,
           offset = 1,
           textInput("cellSelectionComment", "Comment for selection of cells")
    )
  ),checkbsTT("cellSelectionComment"),
  fluidRow(
    column(5,
           offset = 1,
           textInput("cellPatternRM", "cells to be filtered out by pattern"),
    ),checkbsTT("cellPatternRM"),
    column(5,
           offset = 0,
           textInput("cellKeep", "cells to keep")
    )
  ), br(),
  checkbsTT("cellKeep"),
  fluidRow(
    column(10,
           offset = 1,
           textInput("cellKeepOnly", "cells to keep (remove others)")
    )
  ),checkbsTT("cellKeepOnly"),
  fluidRow(
    column(10,
           offset = 1,
           textInput("cellsFiltersOut", "Cells to be removed", width = "100%")
    )
  ), br(),checkbsTT("cellsFiltersOut"),
  fluidRow(column(
    11,
    offset = 1,
    tableSelectionUi("cellSelectionMod")
  )), checkbsTT("cellSelectionMod"), br()
)


# parse all parameters.R files under contributions to include in application
# allTabs holds all tabs regardsless of their location in the GUI
parameterContributions <- list()
# localContributionDir <- .SCHNAPPs_locContributionDir
parFiles <- dir(path = c(paste0(packagePath,  "/contributions"), localContributionDir), pattern = "parameters.R", full.names = TRUE, recursive = TRUE)
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


# submenu items for the paramters main tab
parameterItems <- list(
  shinydashboard::menuSubItem("Normalization", tabName = "normalizations"),
  parameterContributions,
  shinydashboard::menuSubItem("General Parameters", tabName = "generalParameters"),
  shinydashboard::menuSubItem("TSNE plot", tabName = "gQC_tsnePlot"),
  shinydashboard::menuSubItem("Umap", tabName = "gQC_umapPlot")
)


# generalParametersTab ----
generalParametersTab <- shinydashboard::tabItem(
  "generalParameters",
  fluidRow(div(h2("General parameters"), align = "center")),
  br(),
  fluidRow(
    box(
      title = "Parameters for PCA", width = 12,
      # The id lets us use input$tabset1 on the server to find the current tab
      id = "tabsetPCA", 
      fluidRow(
        column(6, offset = 0,
               numericInput("pcaRank", "Number of components", 50, min = 2)),
        column(6, offset = 0,
               numericInput("pcaN", "Number of variable genes to be used", 500, min = 50)),
      ),
      checkbsTT(item="pcaRank"),
      checkbsTT("pcaN"),
      fluidRow(
        column(6, offset = 0,
               checkboxInput("pcaCenter", "center data", TRUE)
        ),
        column(6, offset = 0,
               checkboxInput("pcaScale", "scale data", TRUE)
        ),
        checkbsTT("pcaCenter"),
        checkbsTT("pcaScale")
      ),
      fluidRow(
        column(12,offset = 0,
               textInput("genes4PCA", "Genes to be used for PCA", width = "100%")
        ), checkbsTT("genes4PCA")
      )
      
    ),
    checkbsTT(item="tabsetPCA"),
  ),
  fluidRow(
    tabBox( title = "Parameters for clustering", width = 12,
            id = "tabsetCluster",
            tabPanel("Scran clustering", width = 12,
                     fluidRow(
                       column(6, offset = 0,
                              selectInput("clusterSource", "use raw counts or normalized data?", choices = c("PCA", "counts", "logcounts"), selected = "PCA", width = "100%")),
                       column(6, 
                              numericInput("minClusterSize", "minimum size of each cluster.", 2, min = 2, width = "100%"))
                     ),
                     checkbsTT(item="clusterSource"),
                     checkbsTT(item="minClusterSize"),
                     fluidRow(
                       column(6, 
                              selectInput("clusterMethod", "clustering method to use", choices = c("hclust", "igraph"), selected = "igraph", width = "100%")),
                       column(6, 
                              selectInput("useRanks", "use ranks?\n", choices = c("TRUE", "FALSE"), selected = "TRUE", width = "100%"))
                       
                     ),
                     checkbsTT(item=""),
                     checkbsTT(item=""),
                     
                     fluidRow(
                       column(12, offset = 0,
                              textInput("geneSelectionClustering", "Genes to be used for clustering", width = "100%")
                       )
                     ),
                     checkbsTT(item="geneSelectionClustering"),
                     fluidRow(
                       column(12, offset = 0, textOutput("Nclusters"))
                     ),
                     checkbsTT(item="")
            )
    )),
  # fluidRow(div(h3("Parameters for clustering"), align = "left")),
  br(),
  br(),
  fluidRow(div(h3("Comments"), align = "left")),
  fluidRow(
    if ("shinyMCE" %in% rownames(installed.packages()))
      shinyMCE::tinyMCE(
        "descriptionOfWork",
        "Please describe your work. This will be included in the report."
      ) else 
        textInput("descriptionOfWork", "Please describe your work. This will be included in the report.")
    
  ),
  checkbsTT(item="descriptionOfWork"),
  br(),
  fluidRow(div(h3("Colors"), align = "left")),
  fluidRow(
    actionButton("updateColors", "Update colours", icon = icon("update"))
  ),
  checkbsTT(item="updateColors"),
  fluidRow(column(4,offset = 1,
                  uiOutput('sampleColorSelection')
  ),
  checkbsTT(item="sampleColorSelection"),
  column(4,offset = 1,
         uiOutput('clusterColorSelection')
  ),
  checkbsTT(item="clusterColorSelection")
  )
  # ,
  # fluidRow(
  #   column(11,offset = 1,
  #          textOutput("descriptOfWorkOutput", inline = TRUE))
  # )
)

# renameTab ----
renameTab <- shinydashboard::tabItem(
  tabName = "renameProj",
  fluidRow(div(h3("rename projections"), align = "center")),
  br(),
  fluidRow(
    column(4, offset = 1,
           selectInput("oldPrj", "projections to copy + rename", choices = c("notyet"), selected = "notyet")),
    column(
      4,
      offset = 1,
      textInput("newPrj", "new name of Projection", value = "")
    ),
    column(2, offset = 0, 
           actionButton("updatePrjsButton", "rename")),
    tags$style(type='text/css', "#updatePrjsButton { width:100%; margin-top: 25px;}")
  ),
  checkbsTT(item="oldPrj"),
  checkbsTT(item="newPrj"),
  checkbsTT(item="updatePrjsButton"),
  fluidRow(
    column(4, offset = 1,
           selectInput("delPrj", "projections to delete", choices = c("notyet"), selected = "notyet")),
    column(2, offset = 0, 
           actionButton("delPrjsButton", "delete")),
    tags$style(type='text/css', "#delPrjsButton { width:100%; margin-top: 25px;}")
  ),
  checkbsTT(item="delPrj"),
  checkbsTT(item="delPrjsButton")
)


# # link to the content of the
# parametersTab  = tabItem(tabName = "normalizations",
#                             fluidRow(div(h3('Cell selection'), align = 'center')),
#                             br()
# )
if (DEBUG) {
  cat(file = stderr(), paste("end: tabs.R\n"))
}
