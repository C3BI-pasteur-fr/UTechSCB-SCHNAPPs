# # LIBRARIES -----------------------------------------------------------------
suppressMessages(require(shiny))
suppressMessages(require(shinyBS))
suppressMessages(require(reactlog))
suppressMessages(require(shinyTree))
suppressMessages(require(tibble))
suppressMessages(require(plotly))
suppressMessages(require(shinythemes))
suppressMessages(require(ggplot2))
suppressMessages(require(DT))
suppressMessages(require(pheatmap))
suppressMessages(require(threejs))
# suppressMessages(require(sm))
suppressMessages(require(RColorBrewer))
suppressMessages(require(mclust))
suppressMessages(require(reshape2))
suppressMessages(require(ggplot2))
suppressMessages(require(knitr))
suppressMessages(require(kableExtra))
suppressMessages(require(shinyWidgets))
suppressMessages(require(scater))
suppressMessages(require(kohonen))
suppressMessages(require(Rsomoclu))
suppressMessages(require(gtools))
suppressMessages(require(SingleCellExperiment))
suppressMessages(require(Matrix))
suppressMessages(require(colourpicker))
suppressMessages(require(shinytest))
suppressMessages(require(scran))
suppressMessages(require(callr))
suppressMessages(require(debugme))
suppressMessages(require(BiocSingular))

if (!exists(".schnappsEnv")) {
  .schnappsEnv <- new.env(parent = emptyenv())
}

if (exists("devscShinyApp")) {
  if (devscShinyApp) {
    packagePath <- "inst/app"
  }
} else {
  packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/")
}
if (exists(".SCHNAPPs_locContributionDir", envir = .schnappsEnv)) {
  localContributionDir <- get(".SCHNAPPs_locContributionDir", envir = .schnappsEnv)
} else {
  localContributionDir <- NULL
}
if (exists(".SCHNAPPs_defaultValueSingleGene", envir = .schnappsEnv)) {
  defaultValueSingleGene <- get(".SCHNAPPs_defaultValueSingleGene", envir = .schnappsEnv)
} else {
  defaultValueSingleGene <- "CD52"
}
if (exists(".SCHNAPPs_defaultValueMultiGenes", envir = .schnappsEnv)) {
  defaultValueMultiGenes <- get(".SCHNAPPs_defaultValueMultiGenes", envir = .schnappsEnv)
} else {
  defaultValueMultiGenes <- "CD52, S100A4, S100A9, S100A8"
}
if (exists(".SCHNAPPs_defaultValueRegExGene", envir = .schnappsEnv)) {
  defaultValueRegExGene <- get(".SCHNAPPs_defaultValueRegExGene", envir = .schnappsEnv)
} else {
  defaultValueRegExGene <- "^MT-|^RP|^MRP"
}
if (exists(".SCHNAPPs_DEBUG", envir = .schnappsEnv)) {
  DEBUG <- get(".SCHNAPPs_DEBUG", envir = .schnappsEnv)
} else {
  DEBUG <- FALSE
}
if (exists(".SCHNAPPs_DEBUGSAVE", envir = .schnappsEnv)) {
  DEBUGSAVE <- get(".SCHNAPPs_DEBUGSAVE", envir = .schnappsEnv)
} else {
  DEBUGSAVE <- FALSE
}

# list available colors for samples and clusters, other colors are defined independantly.
# Sys.setenv(DEBUGME = ".")
base::source(paste0(packagePath, "/defaultValues.R"), local = TRUE)
if (!exists("allowedColors")) {
  allowedColors <- unique(c(
    "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e", "#c51b7d", "#e9a3c9",
    "#fde0ef", "#e6f5d0", "#a1d76a", "#4d9221", "#762a83", "#af8dc3", "#e7d4e8", "#d9f0d3",
    "#7fbf7b", "#1b7837", "#b35806", "#f1a340", "#fee0b6", "#d8daeb", "#998ec3", "#542788",
    "#b2182b", "#ef8a62", "#fddbc7", "#d1e5f0", "#67a9cf", "#2166ac", "#b2182b", "#ef8a62",
    "#fddbc7", "#e0e0e0", "#999999", "#4d4d4d"
  ))
}
# Sys.setenv(DEBUGME = ".")
base::source(paste0(packagePath, "/serverFunctions.R"), local = TRUE)


# enableBookmarking(store = "server")

# global variable with directory where to store files to be included in reports
.schnappsEnv$reportTempDir <- base::tempdir()

scShinyServer <- shinyServer(function(input, output, session) {
  session$onSessionEnded(stopApp)
  # TODO needs to be an option
  seed <- 2
  # localContributionDir <- .SCHNAPPs_locContributionDir
  base::set.seed(seed)
  # check that directory is availabl, otherwise create it
  if (DEBUG) {
    if (!dir.exists("~/SCHNAPPsDebug")) {
      base::dir.create("~/SCHNAPPsDebug")
    }
    # TODO ??? clean directory??
  }

  if (exists("devscShinyApp")) {
    if (devscShinyApp) {
      packagePath <- "inst/app"
    }
  } else {
    packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/")
  }
  # files to be included in report
  # developers can add in outputs.R a variable called "myZippedReportFiles"
  zippedReportFiles <- c(
    "Readme.txt", "report.html", "sessionData.RData",
    "normalizedCounts.csv", "variables.used.txt"
  )

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
  projectionFunctions <- list(
    # c("sampleNames", "sample"),
    c("Gene count", "geneCount"),
    c("UMI count", "umiCount"),
    c("before filter", "beforeFilterPrj")
  )
  .schnappsEnv$projectionFunctions <- projectionFunctions

  # differential expression functions
  # used in subcluster analysis
  .schnappsEnv$diffExpFunctions <- list()
  diffExpFunctions <- list()

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
        stop("number of values for normalization function is not 2\n")
      }
    }
  }
  updateRadioButtons(
    session = session, inputId = "sCA_dgeRadioButton",
    choices = dgeChoices
  )
  .schnappsEnv$diffExpFunctions <- diffExpFunctions

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
    projectionFunctions <- append2list(myProjections, projectionFunctions)
    zippedReportFiles <- c(zippedReportFiles, myZippedReportFiles)
  }
  .schnappsEnv$projectionFunctions <- projectionFunctions
}) # END SERVER

# shiny::showReactLog()

# enableBookmarking(store = "server")
