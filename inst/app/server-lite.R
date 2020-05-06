
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
  
  ### history setup
  if (exists("historyPath", envir = .schnappsEnv)) {
    if (!is.null(x = .schnappsEnv$historyPath)) {
      .schnappsEnv$historyPath = paste0(.schnappsEnv$historyPath, "/hist_",format(Sys.time(), "%Y-%b-%d.%H.%M"))
      if (!dir.exists(.schnappsEnv$historyPath)){
        dir.create(.schnappsEnv$historyPath, recursive = T)
      }  
      if (!exists("historyFile", envir = .schnappsEnv)) {
        .schnappsEnv$historyFile = paste0("history.",format(Sys.time(), "%Y-%b-%d.%H.%M"),".Rmd")
      }
      if (is.null(.schnappsEnv$historyFile)) {
        .schnappsEnv$historyFile = "history2.Rmd"
      }
      .schnappsEnv$historyFile <- paste0(.schnappsEnv$historyPath,"/", basename(.schnappsEnv$historyFile))
      line=paste0("---\ntitle: \"history\"\noutput: html_document\n---\n\n```{r setup, include=FALSE}\nknitr::opts_chunk$set(echo = TRUE)
      suppressMessages(require(shiny))
      suppressMessages(require(shinyTree))
      suppressMessages(require(tibble))
      suppressMessages(require(plotly))
      suppressMessages(require(shinythemes))
      suppressMessages(require(ggplot2))
      suppressMessages(require(DT))
      suppressMessages(require(pheatmap))
      suppressMessages(require(threejs))
      suppressMessages(require(RColorBrewer))
      suppressMessages(require(mclust))
      suppressMessages(require(reshape2))
      suppressMessages(require(ggplot2))
      suppressMessages(require(knitr))
      suppressMessages(require(shinyWidgets))
      suppressMessages(require(scater))
      suppressMessages(require(kohonen))
      suppressMessages(require(Rsomoclu))
      suppressMessages(require(SingleCellExperiment))
      suppressMessages(require(Matrix))
      suppressMessages(require(colourpicker))
      # suppressMessages(require(shinytest))
      suppressMessages(require(scran))
      suppressMessages(require(BiocSingular))

      if (\"debugme\" %in% rownames(installed.packages())) {
        suppressMessages(require(debugme))
      }
      if (\"gtools\" %in% rownames(installed.packages())) {
        suppressMessages(require(gtools))
      }
      if (\"kableExtra\" %in% rownames(installed.packages())) {
        suppressMessages(require(kableExtra))
      }
      if (\"reactlog\" %in% rownames(installed.packages())) {
        suppressMessages(require(reactlog))
      }
      suppressMessages(require(Seurat))
      \n```\n" )
      write(line,file=.schnappsEnv$historyFile,append=FALSE)
      
    } else {
      rm("historyPath", envir = .schnappsEnv)
    }
  }
  
  
  # TODO check if file exists
  # TODO as parameter to load user specified information
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
  # TODO do we still need this?
  # heavyCalculations <- list(
  #   c("pca", "pca"),
  #   c("scran_Cluster", "scran_Cluster"),
  #   c("projections", "projections")
  # )
  
  # base projections
  # display name, reactive to calculate projections
  projectionFunctions <- list(
    # c("sampleNames", "sample"),
    # c("Gene count", "geneCount"),
    # c("UMI count", "umiCount"),
    # c("before filter", "beforeFilterPrj")
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
    # myHeavyCalculations <- NULL
    myProjections <- NULL
    myDiffExpFunctions <- NULL
    base::source(fp, local = TRUE)
    
    # heavyCalculations <- append2list(myHeavyCalculations, heavyCalculations)
    projectionFunctions <- append2list(myProjections, projectionFunctions)
    diffExpFunctions <- append2list(myDiffExpFunctions, diffExpFunctions)
  }
  # .schnappsEnv$projectionFunctions 
  
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
    # myHeavyCalculations <- NULL
    myProjections <- NULL
    myZippedReportFiles <- c()
    base::source(fp, local = TRUE)
    # heavyCalculations <- append2list(myHeavyCalculations, heavyCalculations)
    projectionFunctions <- append2list(myProjections, projectionFunctions)
    zippedReportFiles <- c(zippedReportFiles, myZippedReportFiles)
  }
  
  # overwrite all reactives not needed or modified
  base::source(paste0(packagePath, "/reactives-lite.R"), local = TRUE)

  
  .schnappsEnv$projectionFunctions <- projectionFunctions
  # browser()
  
  # overwrite reactives that should not be calculatated anymore
  for (idx in 1:length(projectionFunctions)) {
    if (!is.null(.schnappsEnv$.SCHNAPPs_LiteData[[ projectionFunctions[[idx]][2] ]])) {
      assign(.schnappsEnv$projectionFunctions[[idx]][2],  as.function(alist(.schnappsEnv$.SCHNAPPs_LiteData[[projectionFunctions[[idx]][2]]])))
    }
  }
  
}) # END SERVER

# shiny::showReactLog()

# enableBookmarking(store = "server")
if (DEBUG) cat(file = stderr(), "end: server-lite\n")