#' sCA_selectedDge
#' stores table with differentially expressed genes
#' is used to export file and write csv file
sCA_selectedDge <- reactiveValues(
  sCA_dgeTable = data.frame()
)

#' sCA_getCells
#' get list of cells from selctions
sCA_getCells <- function(projections, cl1, db1, db2) {
  if (DEBUG) cat(file = stderr(), "sCA_getCells started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_getCells")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "sCA_getCells")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_getCells", id = "sCA_getCells", duration = NULL)
  }

    dbCluster = projections$dbCluster
  subsetData <- subset(projections, dbCluster %in% cl1)
  if (class(subsetData[,db1$mapping$x]) == "logical") {
    subsetData[,db1$mapping$x] = as.numeric(subsetData[,db1$mapping$x]) + 1
  }
  if (class(subsetData[,db1$mapping$y]) == "logical") {
    subsetData[,db1$mapping$y] = as.numeric(subsetData[,db1$mapping$y]) + 1
  }
  cells.1 <- rownames(shiny::brushedPoints(subsetData, db1))
  cells.2 <- rownames(shiny::brushedPoints(subsetData, db2))
  retVal <- list(c1 = cells.1, c2 = cells.2)
  return(retVal)
}

#' define different methods for calculating diff. expressed genes
#' first entry is displayed in Radio box, second is function to be called.
myDiffExpFunctions = list(
  c("Chi-square test of an estimated binomial distribution", "sCA_dge_CellViewfunc"),
  c("t-test", "sCA_dge_ttest")
)

#' sCA_dge_CellViewfunc
#' calculate differentically expressed genes given 2 sets of cells
sCA_dge_CellViewfunc <- function(scEx_log, cells.1, cells.2) {
  if (DEBUG) cat(file = stderr(), "sCA_dge_CellViewfunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dge_CellViewfunc")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "sCA_dge_CellViewfunc")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dge_CellViewfunc", id = "sCA_dge_CellViewfunc", duration = NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dge_CellViewfunc.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/SCHNAPPsDebug/sCA_dge_CellViewfunc.RData')
  
  featureData <- rowData(scEx_log)
  scEx_log <- as.matrix(assays(scEx_log)[[1]])
  subsetExpression <- scEx_log[complete.cases(scEx_log[, union(cells.1, cells.2)]),]
  genes.use <- rownames(subsetExpression)
  # expMean exponential mean
  data.1 <- apply(subsetExpression[genes.use, cells.1], 1, expMean)
  data.2 <- apply(subsetExpression[genes.use, cells.2], 1, expMean)
  total.diff <- (data.1 - data.2)

  genes.diff <- names(which(abs(total.diff) > .2))
  genes.use <- ainb(genes.use, genes.diff)

  retVal <-
    DiffExpTest(subsetExpression, cells.1, cells.2, genes.use = genes.use)
  retVal[, "avg_diff"] <- total.diff[rownames(retVal)]
  retVal$symbol <-
    featureData[rownames(retVal), "symbol"]
  return(retVal)
}

#' sCA_dge_ttest
#' t-test on selected cells
sCA_dge_ttest <- function(scEx_log, cells.1, cells.2) {
  if (DEBUG) cat(file = stderr(), "sCA_dge_ttest started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dge_ttest")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "sCA_dge_ttest")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dge_ttest", id = "sCA_dge_ttest", duration = NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dge_ttest.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/SCHNAPPsDebug/sCA_dge_ttest.RData')
  
  featureData <- rowData(scEx_log)
  scEx_log <- as.matrix(assays(scEx_log)[[1]])
  subsetExpression <- scEx_log[complete.cases(scEx_log[, union(cells.1, cells.2)]),]
  genes.use <- rownames(subsetExpression)

  p_val <- apply(subsetExpression, 1, function(x) t.test(x[cells.1], x[cells.2])$p.value)
  p_val <- p_val[!is.na(p_val)]

  retVal = data.frame(p_val = p_val, symbol = featureData[names(p_val), "symbol"])
  return(retVal)
}

#' sCA_dge
#' manage calculation for differential expression analysis
sCA_dge <- reactive({
  if (DEBUG) cat(file = stderr(), "sCA_dge started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dge")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "sCA_dge")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dge", id = "sCA_dge", duration = NULL)
    removeNotification(id = "dgewarning")
  }

  scEx_log <- scEx_log()
  projections <- projections()
  cl1 <- input$sCA_dgeClustersSelection
  db1 <- input$db1
  db2 <- input$db2
  method <- input$sCA_dgeRadioButton

  if (is.null(scEx_log) | is.null(projections)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dge.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/SCHNAPPsDebug/sCA_dge.RData')

  methodIdx <- ceiling(which(unlist(diffExpFunctions)== method)/2)
  dgeFunc <- diffExpFunctions[[methodIdx]][2]
  gCells <- sCA_getCells(projections, cl1, db1, db2)
  retVal <- do.call(dgeFunc, args = list(scEx_log = scEx_log,
                                         cells.1 = gCells$c1, cells.2 = gCells$c2))

  if (nrow(retVal) == 0) {
    if (DEBUG) cat(file = stderr(), "dge: nothing found\n")
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("dge: nothing found", id = "dgewarning", duration = 10, type = "warning")
    }
  }
  # update reactiveValue
  sCA_selectedDge$sCA_dgeTable <- retVal

  exportTestValues(sCA_dge = {retVal})
  return(retVal)
})

# using these global variables allows us to store a set values even when the projections are changing
subClusterDim1 <- "PC1"
subClusterDim2 <- "PC2"
subClusterClusters <<- NULL

observe({
  if (DEBUG) cat(file = stderr(), "observe: sCA_subscluster_x1\n")
  subClusterDim1 <<- input$sCA_subscluster_x1
})

observe({
  if (DEBUG) cat(file = stderr(), "observe: sCA_subscluster_y1\n")
  subClusterDim2 <<- input$sCA_subscluster_y1
})

#' TODO
#' if this observer is really needed we need to get rid of projections
observe({
  if (DEBUG) cat(file = stderr(), "observe: projections\n")
  projections <- projections()
  if (!is.null(projections)) {
    noOfClusters <- levels(as.factor(projections$dbCluster))
    # noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
    if (is.null(subClusterClusters)){
       subClusterClusters <<- noOfClusters
    }
  }
})

observe({
  if (DEBUG) cat(file = stderr(), "observe: sCA_dgeClustersSelection\n")
  subClusterClusters <<- input$sCA_dgeClustersSelection
})


# subcluster axes ----
# update axes in subcluster analysis
updateInputSubclusterAxes <- reactive({
  if (DEBUG) cat(file = stderr(), "updateInputSubclusterAxes started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "updateInputSubclusterAxes")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "updateInputSubclusterAxes")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("updateInputSubclusterAxes", id = "updateInputSubclusterAxes", duration = NULL)
  }

    projections <- projections()
  # we combine the group names with the projections to add ability to select groups
  # gn <- groupNames$namesDF
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/updateInputSubclusterAxes.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/updateInputSubclusterAxes.RData")
  # if (length(gn) > 0) {
  #   projections <- cbind(projections, gn[rownames(projections), ] * 1)
  # }
  # Can also set the label and select items
  updateSelectInput(session, "sCA_subscluster_x1",
                    choices = colnames(projections),
                    selected = subClusterDim1
  )

  updateSelectInput(session, "sCA_subscluster_y1",
                    choices = colnames(projections),
                    selected = subClusterDim2
  )
})



#' subCluster2Dplot
#' plots cells in 2D for the subcluster anlaysis. The handling of the selection is done
#' outside this function
subCluster2Dplot <- function() {
  if (DEBUG) cat(file = stderr(), "subCluster2Dplot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "subCluster2Dplot")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "subCluster2Dplot")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("subCluster2Dplot", id = "subCluster2Dplot", duration = NULL)
  }

  renderPlot({
    if (DEBUG) cat(file = stderr(), "output$sCA_dge_plot2\n")

    projections <- projections()
    x1 <- input$sCA_subscluster_x1
    y1 <- input$sCA_subscluster_y1
    c1 <- input$sCA_dgeClustersSelection

    if (is.null(projections)) {
      return(NULL)
    }
    if (DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/sCA_dge_plot2.RData", list = c(ls(envir = globalenv(), ls())))
    }
    # load(file="~/SCHNAPPsDebug/sCA_dge_plot2.RData")


    subsetData <- subset(projections, dbCluster %in% c1)
    p1 <-
      ggplot(subsetData,
             aes_string(x = x1, y = y1),
             color = "dbCluster"
      ) +
      geom_point(aes(colour = dbCluster)) +
      geom_point(
        shape = 1,
        size = 4,
        aes(colour = dbCluster)
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          size = 12,
          vjust = 0.5
        ),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16),
        legend.position = "none"
      ) +
      ggtitle(c1)
    p1
  })
}
