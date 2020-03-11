suppressMessages(require(uwot))
suppressMessages(require(tidyverse))
suppressMessages(require(SingleCellExperiment))

# here we define reactive values/variables

# save to history violoin observer ----
observe({
  clicked  = input$save2Histumi
  if (DEBUG) cat(file = stderr(), "observe input$save2HistVio \n")
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "save2Hist")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("save2Hist", id = "save2Hist", duration = NULL)
  }
  
  add2history(type = "renderPlot", input = input, comment = "UMI histogram",  
              plotData = .schnappsEnv[["gQC_plotUmiHist"]])
  
})


# save to history save2HistSample observer ----
observe({
  clicked  = input$save2HistSample
  if (DEBUG) cat(file = stderr(), "observe input$save2HistVio \n")
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "save2Hist")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("save2Hist", id = "save2Hist", duration = NULL)
  }
  
  add2history(type = "renderPlot", input = input, comment = "Sample histogram",  
              plotData = .schnappsEnv[["gQC_plotSampleHist"]])
  
})

# save to history save2HistSample observer ----
observe({
  clicked  = input$save2Histvar
  if (DEBUG) cat(file = stderr(), "observe input$save2HistVio \n")
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "save2Hist")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("save2Hist", id = "save2Hist", duration = NULL)
  }
  
  add2history(type = "renderPlot", input = input, comment = "PC variance",  
              plotData = .schnappsEnv[["gQC_variancePCA"]])
  
})


# gQC_scaterReadsFunc ----
#' gQC_scaterReadsFunc
#' calculate the QC metrix and return updated singleCellExperiment object
#' TODO make the 200 a parameter
gQC_scaterReadsFunc <- function(scEx) {
  if (DEBUG) cat(file = stderr(), "gQC_scaterReadsFunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_scaterReadsFunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_scaterReadsFunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_scaterReadsFunc", id = "gQC_scaterReadsFunc", duration = NULL)
  }
  
  if (is(assays(scEx)[[1]], "dgTMatrix")) {
    assays(scEx)[["counts"]] <- as(assays(scEx)[["counts"]], "dgCMatrix")
  }
  
  # ercc <- rownames(scEx)[grepl("ERCC-", rownames(scEx), ignore.case = TRUE)]
  #
  # mt <- rownames(scEx)[grepl("^MT", rownames(scEx), ignore.case = TRUE)]
  
  cellMet <- scater::perCellQCMetrics(
    scEx
  )
  featureMet <- scater::perFeatureQCMetrics(
    scEx
  )
  filter_by_expr_features <- (cellMet$detected > 200)
  scEx$use <- (
    # sufficient features (genes)
    filter_by_expr_features
    # sufficient molecules counted
    # filter_by_total_counts &
    # sufficient endogenous RNA
    # filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    # filter_by_MT
  )
  
  return(scEx)
}

# scaterReads ----
#' scaterReads
#' singleCellExperiment object/reactive with QC metrix
scaterReads <- reactive({
  if (DEBUG) cat(file = stderr(), "scaterReads started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scaterReads")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scaterReads")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scaterReads", id = "scaterReads", duration = NULL)
  }
  
  scEx <- scEx()
  # scEx_log = scEx_log()
  if (is.null(scEx)) {
    return(NULL)
  }
  retVal <- gQC_scaterReadsFunc(scEx)
  
  exportTestValues(scaterReads = {
    str(retVal)
  })
  return(retVal)
})


# gQC_sampleHistFunc ----
#' gQC_sampleHistFunc
#' create a histogram from samples
gQC_sampleHistFunc <- function(samples, scols) {
  if (DEBUG) cat(file = stderr(), "gQC_sampleHistFunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_sampleHistFunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_sampleHistFunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_sampleHistFunc", id = "gQC_sampleHistFunc", duration = NULL)
  }
  
  counts <- table(samples)
  df <- as.data.frame(counts)
  ggplot(data = df,aes(x=samples, y=Freq, fill=samples)) + geom_bar(stat = "identity")  + 
    scale_color_manual(values=scols) 
  # barplot(counts,
  #         main = "histogram of number of cell per sample",
  #         xlab = "Samples",
  #         col=scols
  # )
}


# projectionTable -----
#' projectionTable
#' input for tableSelectionServer
#' presents all projections
projectionTable <- reactive({
  if (DEBUG) cat(file = stderr(), "projectionTable started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "projectionTable")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "projectionTable")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("projectionTable", id = "projectionTable", duration = NULL)
  }
  
  projections <- projections()
  
  if (is.null(projections)) {
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/projectionTable.RData", list = c(ls()))
  }
  # load(file = "~/SCHNAPPsDebug/projectionTable.RData")
  
  printTimeEnd(start.time, "projectionTable")
  exportTestValues(projectionTable = {
    projections
  })
  return(projections)
})

# tsne ----
#' tsne
#' reactive calculating the tSNE projections

# should be a module :
# https://shiny.rstudio.com/articles/modules.html
# https://stackoverflow.com/questions/43976128/create-a-reactive-function-outside-the-shiny-app
#
# I guess that modules are self-contained and one cannot retrieve information from the parent session
# this means that I cannot create rectivity from within the module to inputs that are outside.

# output$updatetsneParametersButton <- updateButtonUI(name = "updatetsneParameters",
#                                                     variables = c("gQC_tsneDim", "gQC_tsnePerplexity", "gQC_tsneTheta", "gQC_tsneSeed"  ) )
# updateButton(name = "updatetsneParameters",
#                                                   )

tsne <- reactive({
  if (DEBUG) cat(file = stderr(), "tsne started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "tsne")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "tsne")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tsne", id = "tsne", duration = NULL)
  }
  
  pca <- pca()
  # only recalculate when button is pressed.
  input$updatetsneParameters
  gQC_tsneDim <- isolate(input$gQC_tsneDim)
  gQC_tsnePerplexity <- isolate(input$gQC_tsnePerplexity)
  gQC_tsneTheta <- isolate(input$gQC_tsneTheta)
  gQC_tsneSeed <- isolate(input$gQC_tsneSeed)
  
  if (is.null(pca)) {
    if (DEBUG) cat(file = stderr(), "tsne: NULL\n")
    return(NULL)
  }
  
  retVal <- tsneFunc(pca, gQC_tsneDim, gQC_tsnePerplexity, gQC_tsneTheta, gQC_tsneSeed)
  if (is.null(tsne)) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(paste("Problem with tsne:", e),
                       id = "gQC_tsneWarning",
                       type = "error",
                       duration = NULL
      )
    }
    return(NULL)
  }
  
  setRedGreenButton(
    vars = list(
      c("gQC_tsneDim", gQC_tsneDim),
      c("gQC_tsnePerplexity", gQC_tsnePerplexity),
      c("gQC_tsneTheta", gQC_tsneTheta),
      c("gQC_tsneSeed", gQC_tsneSeed)
    ),
    button = "updatetsneParameters"
  )
  
  exportTestValues(tsne = {
    retVal
  })
  return(retVal)
})

tsneFunc <- function(pca, gQC_tsneDim, gQC_tsnePerplexity, gQC_tsneTheta, gQC_tsneSeed) {
  if (DEBUG) cat(file = stderr(), "tsneFunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "tsneFunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "tsneFunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tsneFunc", id = "tsneFunc", duration = NULL)
  }
  
  set.seed(seed = gQC_tsneSeed)
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/tsne.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/tsne.RData')
  suppressMessages(require(parallel))
  suppressMessages(require(Rtsne))
  np <- dim(pca$x)[2]
  tsne <- tryCatch(
    {
      Rtsne::Rtsne(
        pca$x[, 1:np],
        pca = FALSE, dims = gQC_tsneDim,
        perplexity = gQC_tsnePerplexity,
        theta = gQC_tsneTheta,
        check_duplicates = FALSE, num_threads = detectCores()
      )
    },
    error = function(e) {
      return(NULL)
    }
  )
  if (is.null(tsne)) {
    return(NULL)
  }
  retVal <- data.frame(tsne$Y)
  colnames(retVal) <- paste0("tsne", c(1:ncol(retVal)))
  return(retVal)
}


# umapReact ----
#' umapReact
#' reactive for calculating UMAP projection
umapReact <- reactive({
  if (DEBUG) cat(file = stderr(), "umapReact started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "umapReact")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "umapReact")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("umapReact", id = "umapReact", duration = NULL)
  }
  
  # xaxis <- input$um_xaxis
  # yaxis <- input$um_yaxis
  # cellT <- input$um_ct
  # inputCT <- input$um_inputCT
  # sampleRatio <- as.numeric(input$um_sampleRatio)
  # sampleIds <- input$um_sampleIds
  # InlevelOrd <- input$um_levelOrd
  # UMAP1 <- input$um_umap1
  # UMAP2 <- input$um_umap2
  runUMAP <- input$activateUMAP
  scEx_log <- scEx_log()
  pca <- pca()
  
  myseed <- isolate(input$gQC_um_randSeed)
  n_neighbors <- isolate(as.numeric(input$gQC_um_n_neighbors))
  n_components <- isolate(as.numeric(input$gQC_um_n_components))
  n_epochs <- isolate(as.numeric(input$gQC_um_n_epochs))
  # alpha <- isolate(as.numeric(input$um_alpha))
  init <- isolate(input$gQC_um_init)
  min_dist <- isolate(as.numeric(input$gQC_um_min_dist))
  set_op_mix_ratio <- isolate(as.numeric(input$gQC_um_set_op_mix_ratio))
  local_connectivity <- isolate(as.numeric(input$gQC_um_local_connectivity))
  bandwidth <- isolate(as.numeric(input$gQC_um_bandwidth))
  gamma <- isolate(as.numeric(input$um_gamma))
  negative_sample_rate <- isolate(as.numeric(input$gQC_um_negative_sample_rate))
  metric <- isolate(input$gQC_um_metric)
  spread <- isolate(as.numeric(input$gQC_um_spread))
  
  if (is.null(scEx_log)) {
    if (DEBUG) cat(file = stderr(), "output$umap_react:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/umap_react.RData", list = c(ls()))
  }
  # load("~/SCHNAPPsDebug/umap_react.RData")
  if (!runUMAP) {
    if (DEBUG) cat(file = stderr(), "output$umap_react:NULL\n")
    return(NULL)
  }
  umapData <- as.matrix(assays(scEx_log)[[1]])
  compCases <- complete.cases(umapData)
  
  # TODO it might be possible to reuse nearest neighbor information to speeed up recomputations
  # with eg. new seed
  
  set.seed(myseed)
  # embedding <- uwot::umap(t(as.matrix(assays(scEx_log)[[1]])),
  embedding <- uwot::umap(pca$x,
                          n_neighbors = n_neighbors,
                          n_components = n_components,
                          n_epochs = n_epochs,
                          # alpha = alpha,
                          init = init,
                          spread = spread,
                          min_dist = min_dist,
                          set_op_mix_ratio = set_op_mix_ratio,
                          local_connectivity = local_connectivity,
                          bandwidth = bandwidth,
                          # gamma = gamma,
                          negative_sample_rate = negative_sample_rate,
                          metric = metric,
                          n_threads = detectCores()
  )
  embedding <- as.data.frame(embedding)
  colnames(embedding) <- paste0("UMAP", 1:n_components)
  rownames(embedding) <- colnames(scEx_log)
  
  setRedGreenButton(
    vars = list(
      c("gQC_um_randSeed", myseed),
      c("gQC_um_n_neighbors", n_neighbors),
      c("gQC_um_n_components", n_components),
      c("gQC_um_n_epochs", n_epochs),
      # c("um_alpha", alpha),
      c("gQC_um_init", init),
      c("gQC_um_min_dist", min_dist),
      c("gQC_um_set_op_mix_ratio", set_op_mix_ratio),
      c("gQC_um_local_connectivity", local_connectivity),
      c("gQC_um_bandwidth", bandwidth),
      c("um_gamma", gamma),
      c("gQC_um_negative_sample_rate", negative_sample_rate),
      c("gQC_um_metric", metric),
      c("gQC_um_spread", spread)
    ),
    button = "activateUMAP"
  )
  
  
  exportTestValues(umapReact = {
    embedding
  })
  return(embedding)
})


# myProjections ----
myProjections <- list(
  c("tsne", "tsne"),
  c("dbCluster", "dbCluster"),
  c("umap", "umapReact")
)

# myHeavyCalculations ----
# declare function as heavy
# myHeavyCalculations <- list(
#   c("scaterReads", "scaterReads"),
#   c("tsne", "tsne")
# )


#' tsnePlot
#' function that plots in 3D the tsne projection
tsnePlot <- function(projections, dimX, dimY, dimZ, dimCol, scols, ccols) {
  if (DEBUG) cat(file = stderr(), "tsnePlot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "tsnePlot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "tsnePlot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tsnePlot", id = "tsnePlot", duration = NULL)
  }
  
  projections <- as.data.frame(projections)
  if (!all(c(dimX, dimY, dimZ) %in% colnames(projections))) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("Selected projections not available. Did you run normalization?", id = "tsnePlotERROR", type = "error", duration = NULL)
    }
    return(NULL)
  }
  
  projections$dbCluster <- as.factor(projections$dbCluster)
  
  if (dimCol == "sampleNames") {
    myColors <- scols
  } else {
    myColors <- NULL
  }
  if (dimCol == "dbCluster") {
    myColors <- ccols
  }
  
  p <-
    plotly::plot_ly(
      projections,
      x = formula(paste("~ ", dimX)),
      y = formula(paste("~ ", dimY)),
      z = formula(paste("~ ", dimZ)),
      type = "scatter3d",
      color = formula(paste("~ ", dimCol)),
      colors = myColors,
      hoverinfo = "text",
      text = paste("Cluster:", as.numeric(as.character(projections$dbCluster))),
      mode = "markers",
      source = "B",
      marker =
        list(
          line = list(width = 0),
          size = rep(10, nrow(projections)),
          sizeref = 3
        )
    )
  return(p)
}
