# TODO: make sure this is working
myZippedReportFiles <- c("DGE.csv")




# sCA_dge_plot1 ----
#' sCA_dge_plot1
#' left plot for selection of cells
output$sCA_dge_plot1 <- subCluster2Dplot()
# SUBCLUSTER DGE PLOT2 -----
#' sCA_dge_plot2
#' right plot
output$sCA_dge_plot2 <- subCluster2Dplot()

# sCA_dgeTable ----
#' sCA_dgeTable
#' Table with differential expressed genes

sCA_dgeTableReac <- reactive({
  if (DEBUG) cat(file = stderr(), "sCA_dgeTableReac started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dgeTableReac")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "sCA_dgeTableReac")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dgeTableReac", id = "sCA_dgeTableReac", duration = NULL)
  }

  scEx <- scEx()
  top.genes <- sCA_dge()

  if (is.null(scEx)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/output_dge.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # cp = load(file="~/SCHNAPPsDebug/output_dge.RData")
  featureData <- rowData(scEx)

  top.genes$symbol <-
    featureData[rownames(top.genes), "symbol"]
  if ("Description" %in% colnames(featureData)) {
    top.genes$Description <- featureData[rownames(top.genes), "Description"]
  }
  rownames(top.genes) <- make.unique(top.genes$symbol, sep="___")
  if (dim(top.genes)[1] > 0) {
    return(top.genes)
  } else {
    return(NULL)
  }

})

# dge table ----
callModule(
  tableSelectionServer,
  "sCA_dgeTable",
  sCA_dgeTableReac)


# sub cluster analysis ( used for 2 panels )
output$sCA_dgeClustersSelection <- renderUI({
  if (DEBUG) cat(file = stderr(), "sCA_dgeClustersSelection started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dgeClustersSelection")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "sCA_dgeClustersSelection")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dgeClustersSelection", id = "sCA_dgeClustersSelection", duration = NULL)
  }

  projections <- projections()
  up1 <- updateInputSubclusterAxes()

  if (DEBUG) cat(file = stderr(), "output$sCA_dgeClustersSelection\n")
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dgeClustersSelection.RData", list = c(ls(envir = globalenv(), ls(), "subClusterClusters")))
  }
  # load(file="~/SCHNAPPsDebug/sCA_dgeClustersSelection.RData")


  if (is.null(projections)) {
    tags$span(style="color:red", "Please load data first")
  } else {
    noOfClusters <- levels(as.factor(projections$dbCluster))
    # noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
    selectizeInput(
      "sCA_dgeClustersSelection",
      label = "Cluster",
      choices = noOfClusters,
      selected = subClusterClusters,
      multiple = TRUE
    )
  }
})




