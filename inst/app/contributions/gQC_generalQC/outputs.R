# source("moduleServer.R", local = TRUE)
# source("reactives.R", local = TRUE)

# TODO: verify that this anything and then integrate in DUMMY
myZippedReportFiles <- c("gqcProjections.csv")



gQC_X1 <<- "tsne1"
gQC_X2 <<- "tsne2"
gQC_X3 <<- "tsne3"
gQC_col <<- "sampleNames"
observe({
  if (DEBUG) cat(file = stderr(), "observe: gQC_dim3D_x\n")
  gQC_X1 <<- input$gQC_dim3D_x
})
observe({
  if (DEBUG) cat(file = stderr(), "observe: gQC_dim3D_y\n")
  gQC_X2 <<- input$gQC_dim3D_y
})
observe({
  if (DEBUG) cat(file = stderr(), "observe: gQC_dim3D_z\n")
  gQC_X3 <<- input$gQC_dim3D_z
})
observe({
  if (DEBUG) cat(file = stderr(), "observe: gQC_col3D\n")
  gQC_col <<- input$gQC_col3D
})

# gQC_update3DInput ----
#' gQC_update3DInput
#' update axes for tsne display
gQC_update3DInput <- reactive({
  if (DEBUG) cat(file = stderr(), "gQC_update3DInput started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_update3DInput")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_update3DInput")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_update3DInput", id = "gQC_update3DInput", duration = NULL)
  }

  projections <- projections()

  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  # choices = colnames(projections)[unlist(lapply(colnames(projections), function(x) !is.factor(projections[,x])))]
  choices <- colnames(projections)
  # Can also set the label and select items
  updateSelectInput(session, "gQC_dim3D_x",
    choices = choices,
    selected = gQC_X1
  )

  updateSelectInput(session, "gQC_dim3D_y",
    choices = choices,
    selected = gQC_X2
  )
  updateSelectInput(session, "gQC_dim3D_z",
    choices = choices,
    selected = gQC_X3
  )
  updateSelectInput(session, "gQC_col3D",
    choices = colnames(projections),
    selected = gQC_col
  )
})

# gQC_tsne_main ----
output$gQC_tsne_main <- plotly::renderPlotly({
  if (DEBUG) cat(file = stderr(), "gQC_tsne_main started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_tsne_main")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_tsne_main")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_tsne_main", id = "gQC_tsne_main", duration = NULL)
  }

  upI <- gQC_update3DInput()
  projections <- projections()
  dimX <- input$gQC_dim3D_x
  dimY <- input$gQC_dim3D_y
  dimZ <- input$gQC_dim3D_z
  dimCol <- input$gQC_col3D
  scols <- sampleCols$colPal
  ccols <- clusterCols$colPal

  if (is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$gQC_tsne_main:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/gQC_tsne_main.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/gQC_tsne_main.RData")

  retVal <- tsnePlot(projections, dimX, dimY, dimZ, dimCol, scols, ccols)

  exportTestValues(tsnePlot = {
    str(retVal)
  })
  layout(retVal)
})

# gQC_umap_main 2D plot ----
callModule(
  clusterServer,
  "gQC_umap_main",
  projections
)

# gQC_projectionTableMod ----
callModule(
  tableSelectionServer,
  "gQC_projectionTableMod",
  projectionTable
)

# gQC_plotUmiHist ----
output$gQC_plotUmiHist <- renderPlot({
  if (DEBUG) cat(file = stderr(), "gQC_plotUmiHist started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_plotUmiHist")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_plotUmiHist")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_plotUmiHist", id = "gQC_plotUmiHist", duration = NULL)
  }

  scEx <- scEx()
  scols <- sampleCols$colPal

  if (is.null(scEx)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/gQC_plotUmiHist.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/gQC_plotUmiHist.RData")

  dat <- data.frame(counts = Matrix::colSums(assays(scEx)[["counts"]]))
  dat$sample <- colData(scEx)$sampleNames
  ggplot(data = dat, aes(counts, fill = sample)) +
    geom_histogram(bins = 50) +
    labs(title = "Histogram for raw counts", x = "count", y = "Frequency") +
    scale_fill_manual(values = scols, aesthetics = "fill")
})

output$gQC_plotSampleHist <- renderPlot({
  if (DEBUG) cat(file = stderr(), "gQC_plotSampleHist started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_plotSampleHist")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_plotSampleHist")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_plotSampleHist", id = "gQC_plotSampleHist", duration = NULL)
  }

  sampleInf <- sampleInfo()
  scols <- sampleCols$colPal

  if (is.null(sampleInf)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sampleHist.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/sampleHist.RData")
  gQC_sampleHistFunc(sampleInf, scols)
})

output$gQC_variancePCA <- renderPlot({
  if (DEBUG) cat(file = stderr(), "gQC_variancePCA started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_variancePCA")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_variancePCA")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_variancePCA", id = "gQC_variancePCA", duration = NULL)
  }

  if (DEBUG) cat(file = stderr(), "output$gQC_variancePCA\n")
  h2("Variances of PCs")
  pca <- pca()
  if (is.null(pca)) {
    return(NULL)
  }
  barplot(pca$var_pcs, main = "Variance captured by first PCs")
})
