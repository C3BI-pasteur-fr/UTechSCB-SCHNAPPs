
scShinyServer <- shinyServer(function(input, output, session) {
  library(shinyjqui)
  base::cat(file = stderr(), "------ ShinyServer LITE running\n")
  session$onSessionEnded(stopApp)
  
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("starting application", id = "startSCHNAPPs", duration = NULL)
  }
  gmtData <- reactiveVal()
  gmtUserData <- reactiveVal()
  
  if(exists(".SCHNAPPs_GMTData", envir = .schnappsEnv)){
    gmtData(.schnappsEnv[[".SCHNAPPs_GMTData"]])
  }
  
  # Here, we store projections that are created during the session. These can be selections of cells or other values that
  # are not possible to precalculate.
  sessionProjections <- reactiveValues(
    prjs = data.frame()
  )
  
  
  clusterMethodReact <- reactiveValues(
    clusterMethod = "igraph",
    clusterSource = "counts"
  )
  
  # collect copied/renamed projections
  projectionsTable <- reactiveValues(
    newProjections = data.frame()
  )
  
  # TODO needs to be an option
  seed <- 2
  # localContributionDir <- .SCHNAPPs_locContributionDir
  base::set.seed(seed)
  # check that directory is availabl, otherwise create it
  if (DEBUG) {
    if (!dir.exists(normalizePath("~/SCHNAPPsDebug"))) {
      base::dir.create(normalizePath("~/SCHNAPPsDebug"))
    }
    # TODO ??? clean directory??
  }
  
  if (base::exists("devscShinyApp")) {
    if (devscShinyApp) {
      if (dir.exists(paths = normalizePath("~/Rstudio/UTechSCB-SCHNAPPs/inst/app/"))){
        packagePath <- normalizePath("~/Rstudio/UTechSCB-SCHNAPPs/inst/app/")
      } else {
        if (dir.exists(paths = normalizePath("~/Rstudio/Schnapps/inst/app/"))){
          packagePath <- normalizePath("~/Rstudio/Schnapps/inst/app/")
        } else {
          stop("package path not found\n")
        }
      }
      # setwd("~/Rstudio/UTechSCB-SCHNAPPs/")
    }
  } else {
    packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/") %>% normalizePath()
  }
  
  
  # files to be included in report
  # developers can add in outputs.R a variable called "myZippedReportFiles"
  zippedReportFiles <- c(
    "Readme.txt", "report.html", "sessionData.RData",
    "normalizedCounts.csv", "variables.used.txt"
  )
  
  base::options(shiny.maxRequestSize = 2000 * 1024^2)

  
  # TODO check if file exists
  # TODO as parameter to load user specified information
  # TODO have this as an option to load other files
  if (file.exists(normalizePath(paste0(packagePath, "/geneLists.RData")))) {
    base::load(file = normalizePath(paste0(packagePath, "/geneLists.RData")))
  } else {
    if (!exists("geneLists")) {
      geneLists <- list(emtpy = list())
    }
  }
  
  
  # base projections
  # display name, reactive to calculate projections
  projectionFunctions <- list(
    # c("sampleNames", "sample"),
    # c("Gene count", "geneCount"),
    # c("UMI count", "umiCount"),
    # c("before filter", "beforeFilterPrj")
  )
  .schnappsEnv$projectionFunctions <- projectionFunctions
  
  
  ### history setup
  if (base::exists("historyPath", envir = .schnappsEnv)) {
    if (!is.null(x = .schnappsEnv$historyPath)) {
      .schnappsEnv$historyPath = normalizePath(paste0(.schnappsEnv$historyPath, "/hist_",format(Sys.time(), "%Y-%b-%d.%H.%M")))
      if (!dir.exists(.schnappsEnv$historyPath)){
        dir.create(.schnappsEnv$historyPath, recursive = T)
      }  
      if (!base::exists("historyFile", envir = .schnappsEnv)) {
        .schnappsEnv$historyFile = paste0("history.",format(Sys.time(), "%Y-%b-%d.%H.%M"),".Rmd")
      }
      if (is.null(.schnappsEnv$historyFile)) {
        .schnappsEnv$historyFile = "history2.Rmd"
      }
      .schnappsEnv$historyFile <- paste0(.schnappsEnv$historyPath, .Platform$file.sep, basename(.schnappsEnv$historyFile))
      line=paste0("---\ntitle: \"history\"\noutput: html_document\n---\n\n```{r setup, include=FALSE}\nknitr::opts_chunk$set(echo = TRUE)
      suppressMessages(require(shiny))
      suppressMessages(require(shinyTree))
      suppressMessages(require(tibble))
      suppressMessages(require(plotly))
      suppressMessages(require(shinythemes))
      suppressMessages(require(ggplot2))
      suppressMessages(require(DT))
      suppressMessages(require(pheatmap))
      # suppressMessages(require(threejs))
      suppressMessages(require(RColorBrewer))
      # suppressMessages(require(mclust))
      suppressMessages(require(reshape2))
      suppressMessages(require(ggplot2))
      suppressMessages(require(knitr))
      suppressMessages(require(shinyWidgets))
      suppressMessages(require(scater))
      suppressMessages(require(kohonen))
      suppressMessages(require(SingleCellExperiment))
      suppressMessages(require(Matrix))
      suppressMessages(require(colourpicker))
      # suppressMessages(require(shinytest))
      suppressMessages(require(scran))
      suppressMessages(require(ggalluvial))
      suppressMessages(require(BiocSingular))
      suppressMessages(require(dplyr))

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
      .schnappsEnv = list()
      .schnappsEnv$DEBUGSAVE = FALSE
      source(system.file(\"app\", \"serverFunctions.R\", package = \"SCHNAPPs\"))
      DEBUG=FALSE
      \n```\n" )
      write(line,file=.schnappsEnv$historyFile,append=FALSE)
      
    } else {
      rm("historyPath", envir = .schnappsEnv)
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
  
  
  # differential expression functions
  # used in subcluster analysis
  .schnappsEnv$diffExpFunctions <- list()
  diffExpFunctions <- list()
  
  
  projectionColors <- reactiveValues()
  
  
  # load global reactives, modules, etc ----
  base::source(normalizePath(paste0(packagePath, "/reactives.R")), local = TRUE)
  base::source(normalizePath(paste0(packagePath, "/outputs.R")), local = TRUE)
  base::source(normalizePath(paste0(packagePath, "/modulesUI.R")), local = TRUE)
  base::source(normalizePath(paste0(packagePath, "/moduleServer.R")), local = TRUE)
  
  
  
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
    path = c(normalizePath(paste0(packagePath, "/contributions")), localContributionDir), pattern = "reactives.R",
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
  updateSelectizeInput(
    session = session, inputId = "sCA_dgeRadioButton",
    choices = dgeChoices
  )
  .schnappsEnv$diffExpFunctions <- diffExpFunctions
  
  # load contribution outputs ----
  # parse all outputs.R files under contributions to include in application
  uiFiles <- base::dir(
    path = c(normalizePath(paste0(packagePath, "/contributions")), localContributionDir), pattern = "outputs.R",
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
  base::source(normalizePath(paste0(packagePath, "/reactives-lite.R")), local = TRUE)
  cat(file = stderr(), "HALL============================\n")
  # deepDebug()
  .schnappsEnv$projectionFunctions <- projectionFunctions
  
  # overwrite reactives that should not be calculatated anymore
  for (idx in 1:length(projectionFunctions)) {
    if (!is.null(.schnappsEnv$.SCHNAPPs_LiteData[[ projectionFunctions[[idx]][2] ]])) {
      assign(.schnappsEnv$projectionFunctions[[idx]][2],  as.function(alist(.schnappsEnv$.SCHNAPPs_LiteData[[projectionFunctions[[idx]][2]]])))
    }
  }
  
  # Scater QC ----
  output$DE_scaterQC <- renderImage(deleteFile = F, {
    start.time <- base::Sys.time()
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "DE_scaterQC")
      }
    )
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("DE_scaterQC", id = "DE_scaterQC", duration = NULL)
    }
    if (DEBUG) cat(file = stderr(), "output$DE_scaterQC\n")
    scaterReads <- scaterReads()
    if (is.null(scaterReads)) {
      return(list(
        src = "",
        contentType = "image/png",
        width = 10,
        height = 10,
        alt = "Scater plot will be here when 'run scater' is checked"
      ))
    }
    
    DE_scaterPNG()
  })
  
  # introRMD ----
  output$introRMD <- renderUI({
    cat(file = stderr(), paste("wd:", getwd(), "\n"))
    introFile = 'intro.Rmd'
    if (!file.exists(introFile)) return(HTML(paste("introduction file intro.Rmd is missing", getwd(),"\n")))
    HTML(markdown::markdownToHTML(knit(introFile, quiet = TRUE), fragment.only = T))
    # includeHTML("intro.html")
  })
  
   
  
  # colors for samples ----
  sampleCols <- reactiveValues(colPal = get(".SCHNAPPs_LiteData",envir = .schnappsEnv)$sampleCol)
  
  # colors for clusters ----
  clusterCols <- reactiveValues(colPal = get(".SCHNAPPs_LiteData",envir = .schnappsEnv)$clusterCol)
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "startSCHNAPPs")
  }
  
  # add2history
  
}) # END SERVER

# shiny::showReactLog()

# enableBookmarking(store = "server")
if (DEBUG) cat(file = stderr(), "end: server-lite\n")
