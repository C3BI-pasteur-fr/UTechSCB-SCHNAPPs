# # LIBRARIES -----------------------------------------------------------------
require(shiny)
require(reactlog)
require(shinyTree)
require(tibble)
require(plotly)
require(shinythemes)
require(ggplot2)
require(DT)
require(pheatmap)
require(threejs)
require(sm)
require(RColorBrewer)
require(mclust)
require(reshape2)
require(ggplot2)
require(knitr)
require(kableExtra)
require(shinyWidgets)
require(scater)
require(shinyMCE)
require(kohonen)
require(Rsomoclu)
require(gtools)
require(SingleCellExperiment)
require(Matrix)
require(colourpicker)
require(shinytest)
require(scran)
require(callr)
require(debugme)
require(BiocSingular)
require(shinyBS)

if (exists("devscShinyApp")) {
  if (devscShinyApp) {
    packagePath <- "inst/app"
  }
} else {
  packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE)
  packagePath <- paste0(packagePath, "/app/")
}
if (exists(".SCHNAPPs_locContributionDir")) {
  localContributionDir <- .SCHNAPPs_locContributionDir
} else {
  localContributionDir <- NULL
}
if (exists(".SCHNAPPs_defaultValueSingleGene")) {
  defaultValueSingleGene <- .SCHNAPPs_defaultValueSingleGene
} else {
  defaultValueSingleGene <- "CD52"
}
if (exists(".SCHNAPPs_defaultValueMultiGenes")) {
  defaultValueMultiGenes <- .SCHNAPPs_defaultValueMultiGenes
} else {
  defaultValueMultiGenes <- "CD52, S100A4, S100A9, S100A8"
}
if (exists(".SCHNAPPs_defaultValueRegExGene")) {
  defaultValueRegExGene <- .SCHNAPPs_defaultValueRegExGene
} else {
  defaultValueRegExGene <- "^MT-|^RP|^MRP"
}
if (exists(".SCHNAPPs_DEBUG")) {
  DEBUG <- .SCHNAPPs_DEBUG
} else {
  DEBUG <- FALSE
}
if (exists(".SCHNAPPs_DEBUGSAVE")) {
  DEBUGSAVE <- .SCHNAPPs_DEBUGSAVE
} else {
  DEBUGSAVE <- FALSE
}

# list available colors for samples and clusters, other colors are defined independantly.
if (!exists("allowedColors")) {
  allowedColors <- unique(c(
    "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e", "#c51b7d", "#e9a3c9",
    "#fde0ef", "#e6f5d0", "#a1d76a", "#4d9221", "#762a83", "#af8dc3", "#e7d4e8", "#d9f0d3",
    "#7fbf7b", "#1b7837", "#b35806", "#f1a340", "#fee0b6", "#d8daeb", "#998ec3", "#542788",
    "#b2182b", "#ef8a62", "#fddbc7", "#d1e5f0", "#67a9cf", "#2166ac", "#b2182b", "#ef8a62",
    "#fddbc7", "#e0e0e0", "#999999", "#4d4d4d"
  ))
}
Sys.setenv(DEBUGME = ".")
base::source(paste0(packagePath, "/serverFunctions.R"))


# enableBookmarking(store = "server")

# global variable with directory where to store files to be included in reports
reportTempDir <<- base::tempdir()

scShinyServer <- shinyServer(function(input, output, session) {
  session$onSessionEnded(stopApp)
  # TODO needs to be an option
  seed <- 2
  localContributionDir <- .SCHNAPPs_locContributionDir
  # cat(file = stderr(), paste("bernd2", str(environment())), "\n")
  # cat(file = stderr(), paste("bernd2", str(parent.env(environment()))), "\n")
  # cat(file = stderr(), paste("bernd2", str(parent.env(parent.env(environment())))), "\n")
  # cat(file = stderr(), paste("bernd2", localContributionDir, "\n"))
  base::set.seed(seed)
  # check that directory is availabl, otherwise create it
  if (DEBUG) {
    if (!dir.exists("~/SCHNAPPsDebug")) {
      base::dir.create("~/SCHNAPPsDebug")
    }
    # TODO ??? clean directory??
  }

  # files to be included in report
  # developers can add in outputs.R a variable called "myZippedReportFiles"
  zippedReportFiles <- c("Readme.txt", "report.html", "sessionData.RData", 
                         "normalizedCounts.csv", "variables.used.txt")

  base::options(shiny.maxRequestSize = 2000 * 1024^2)

  # TODO check if file exists
  # TODO have this as an option to load other files
  if (file.exists(paste0(packagePath, "/geneLists.RData"))) {
    base::load(file = paste0(packagePath, "/geneLists.RData"))
  } else {
    if (!exists("geneLists")) {
      geneLists <- list(emtpy = list())
    }
  }

  if (DEBUG) base::cat(file = stderr(), "ShinyServer running\n")
  # base calculations that are quite expensive to calculate
  # display name, reactive name to be executed
  heavyCalculations <- list(
    c("pca", "pca"),
    c("scran_Cluster", "scran_Cluster"),
    c("projections", "projections")
  )

  # base projections
  # display name, reactive to calculate projections
  projectionFunctions <<- list(
    c("sampleNames", "sample"),
    c("Gene count", "geneCount"),
    c("UMI count", "umiCount"),
    c("before filter", "beforeFilterPrj")
  )

  # differential expression functions
  # used in subcluster analysis
  diffExpFunctions <<- list()

  # load global reactives, modules, etc ----
  base::source(paste0(packagePath, "/reactives.R"), local = TRUE)
  base::source(paste0(packagePath, "/outputs.R"), local = TRUE)
  base::source(paste0(packagePath, "/modulesUI.R"), local = TRUE)
  base::source(paste0(packagePath, "/moduleServer.R"), local = TRUE)

  # bookmarking ----
  # couldn't get bookmarking to work, esp. with the input file
  # setBookmarkExclude(c("bookmark1"))
  # observeEvent(input$bookmark1, {
  #   if (DEBUG) cat(file = stderr(), paste("bookmarking: \n"))
  #   if (DEBUG) cat(file = stderr(), paste(names(input), collapse = "\n"))
  #
  #   session$doBookmark()
  #   if (DEBUG) cat(file = stderr(), paste("bookmarking: DONE\n"))
  # })
  # Need to exclude the buttons from themselves being bookmarked

  # load contribution reactives ----
  # parse all reactives.R files under contributions to include in application
  uiFiles <- base::dir(
    path = c(paste0(packagePath, "/contributions"), localContributionDir), pattern = "reactives.R",
    full.names = TRUE, recursive = TRUE
  )
  for (fp in uiFiles) {
    if (DEBUG) base::cat(file = stderr(), paste("loading: ", fp, "\n"))
    myHeavyCalculations <- NULL
    myProjections <- NULL
    myDiffExpFunctions <- NULL
    base::source(fp, local = TRUE)

    heavyCalculations <- append2list(myHeavyCalculations, heavyCalculations)
    projectionFunctions <- append2list(myProjections, projectionFunctions)
    diffExpFunctions <- append2list(myDiffExpFunctions, diffExpFunctions)
  }

  # update diffExpression radiobutton
  dgeChoices <- c()
  if (length(diffExpFunctions) > 0) {
    for (li in 1:length(diffExpFunctions)) {
      liVal <- diffExpFunctions[[li]]
      if (length(liVal) == 2) {
        dgeChoices <- c(dgeChoices, liVal[1])
      } else {
        # shouldn't happen
        error("number of values for normalization function is not 2\n")
      }
    }
  }
  updateRadioButtons(
    session = session, inputId = "sCA_dgeRadioButton",
    choices = dgeChoices
  )
  # make variable global
  diffExpFunctions <<- diffExpFunctions

  # load contribution outputs ----
  # parse all outputs.R files under contributions to include in application
  uiFiles <- base::dir(
    path = c(paste0(packagePath, "/contributions"), localContributionDir), pattern = "outputs.R",
    full.names = TRUE, recursive = TRUE
  )
  for (fp in uiFiles) {
    if (DEBUG) cat(file = stderr(), paste("loading: ", fp, "\n"))
    myHeavyCalculations <- NULL
    myProjections <- NULL
    myZippedReportFiles <- c()
    base::source(fp, local = TRUE)
    heavyCalculations <- append2list(myHeavyCalculations, heavyCalculations)
    projectionFunctions <<- append2list(myProjections, projectionFunctions)
    zippedReportFiles <- c(zippedReportFiles, myZippedReportFiles)
  }
}) # END SERVER

# shiny::showReactLog()

# enableBookmarking(store = "server")
