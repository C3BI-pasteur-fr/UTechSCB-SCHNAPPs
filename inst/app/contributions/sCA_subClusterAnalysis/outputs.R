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
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_dgeTableReac")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dgeTableReac", id = "sCA_dgeTableReac", duration = NULL)
  }

  scEx <- scEx()
  top.genes <- sCA_dge()

  if (is.null(scEx) || is.null(top.genes)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/output_dge.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # cp = load(file="~/SCHNAPPsDebug/output_dge.RData")
  featureData <- rowData(scEx)

  top.genes$symbol <-
    featureData[rownames(top.genes), "symbol"]
  if ("Description" %in% colnames(featureData)) {
    top.genes$Description <- featureData[rownames(top.genes), "Description"]
  }
  rownames(top.genes) <- make.unique(as.character(top.genes$symbol), sep = "___")
  if (dim(top.genes)[1] > 0) {
    # change inf to high/low number
    infIdx <- which(is.infinite(top.genes$avg_diff))
    top.genes$avg_diff[infIdx[top.genes$avg_diff[infIdx] > 0]] <- 9999999
    top.genes$avg_diff[infIdx[top.genes$avg_diff[infIdx] < 0]] <- -9999999
    top.genes$Description[is.na(top.genes$Description)] = ""
    return(top.genes)
  } else {
    return(NULL)
  }
})

# dge table ----
callModule(
  tableSelectionServer,
  "sCA_dgeTable",
  sCA_dgeTableReac
)


# sub cluster analysis ( used for 2 panels )
output$sCA_dgeClustersSelection <- renderUI({
  if (DEBUG) cat(file = stderr(), "sCA_dgeClustersSelection started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dgeClustersSelection")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_dgeClustersSelection")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dgeClustersSelection", id = "sCA_dgeClustersSelection", duration = NULL)
  }

  projections <- projections()
  up1 <- updateInputSubclusterAxes()

  if (DEBUG) cat(file = stderr(), "output$sCA_dgeClustersSelection\n")
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dgeClustersSelection.RData", list = c(ls(envir = globalenv(), ls(), ls(envir = .schnappsEnv))))
  }
  # load(file="~/SCHNAPPsDebug/sCA_dgeClustersSelection.RData")


  if (is.null(projections)) {
    tags$span(style = "color:red", "Please load data first")
  } else {
    noOfClusters <- levels(as.factor(projections$dbCluster))
    # noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
    selectizeInput(
      "sCA_dgeClustersSelection",
      label = "Cluster",
      choices = noOfClusters,
      selected = .schnappsEnv$subClusterClusters,
      multiple = TRUE
    )
  }
})

output$sCA_volcanoPlot <- plotly::renderPlotly({
  if (DEBUG) cat(file = stderr(), "sCA_volcanoPlot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_volcanoPlot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_volcanoPlot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_volcanoPlot", id = "sCA_volcanoPlot", duration = NULL)
  }
  
  DGEdata <- sCA_dgeTableReac()
  
  if (is.null(DGEdata)) {
    if (DEBUG) cat(file = stderr(), "output$sCA_volcanoPlot:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_volcanoPlot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/sCA_volcanoPlot.RData")
  
  if ("p_val" %in% colnames(DGEdata)){
    pval = "p_val"
  }
  if ("avg_diff" %in% colnames(DGEdata)){
    effect_size = "avg_diff"
  }
  DGEdata[is.na(DGEdata[,pval]),pval] = 1
  DGEdata[is.na(DGEdata[,effect_size]),effect_size] = 0
  require(manhattanly)
  DGEdata$EFFECTSIZE = DGEdata[,effect_size]
  DGEdata$P = DGEdata[, pval]
  TEXT <- paste(paste("symbol: ",DGEdata$symbol), sep = "<br>")
  
  retVal = volcanoly(DGEdata, snp="symbol")

  # magrittr::%<>%
  retVal %<>% plotly::add_trace(x = DGEdata$EFFECTSIZE, y = -log10(DGEdata$p_val),
                         type = "scatter",
                         mode = "markers",
                         text = TEXT,
                         marker = list(color = "grey",
                                       size = 5),
                         name = "")
# retVal
  
  exportTestValues(dgeVolcanoPlot = {
    str(retVal)
  })
  layout(retVal)
  
})
