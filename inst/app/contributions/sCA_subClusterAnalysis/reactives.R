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
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_getCells")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_getCells", id = "sCA_getCells", duration = NULL)
  }

  # save(file = "~/SCHNAPPsDebug/sCA_getCells.RData", list = c( ls()))
  # cp =load("~/SCHNAPPsDebug/sCA_getCells.RData")
  # dbCluster = projections$dbCluster
  subsetData <- projections[cl1,]
  
    if (is(subsetData[,db1$mapping$x], "logical")) {
      subsetData[,db1$mapping$x] = as.numeric(subsetData[,db1$mapping$x]) + 1
    }
    if (is(subsetData[, db1$mapping$y], "logical")) {
      subsetData[, db1$mapping$y] <- as.numeric(subsetData[, db1$mapping$y]) + 1
    }
  

  # factors and brushedPoints don't work together.
  # so we change a factor into a numeric
  # TODO WHY is discrete_limits set/misused?????
  # ggplot is not displaying levels for which there are no values
  # thus, the numbering might be off 
  if (is(subsetData[, db1$mapping$x], "factor")) {
    subsetData[, db1$mapping$x] <- as.numeric(subsetData[, db1$mapping$x])
    db1$domain$discrete_limits <- NULL
  }
  if (is(subsetData[, db1$mapping$y], "factor")) {
    subsetData[, db1$mapping$y] <- as.numeric(subsetData[, db1$mapping$y])
    db1$domain$discrete_limits <- NULL
  }
  db1$domain$discrete_limits <- NULL
  cells.1 <- rownames(shiny::brushedPoints(df = subsetData, brush = db1))
  
  cells.2 = c()
  if (!is.null(db2)) {
    db2$domain$discrete_limits <- NULL

  # factors and brushedPoints don't work together.
  # so we change a factor into a numeric
  if (is(subsetData[, db2$mapping$x], "factor")) {
    subsetData[, db2$mapping$x] <- as.numeric(subsetData[, db2$mapping$x])
    db2$domain$discrete_limits <- NULL
  }
  if (is(subsetData[, db2$mapping$y], "factor")) {
    subsetData[, db2$mapping$y] <- as.numeric(subsetData[, db2$mapping$y])
    db2$domain$discrete_limits <- NULL
  }
    cells.2 <- rownames(shiny::brushedPoints(df = subsetData, brush = db2))
  } else {
    cells.2 <- rownames(subsetData)[!rownames(subsetData) %in% cells.1]
  }
  
  # cells.2 <- rownames(shiny::brushedPoints(df = subsetData, brush = db2))
  retVal <- list(c1 = cells.1, c2 = cells.2)
  # save(file = "~/SCHNAPPsDebug/sCA_getCells.RData", list = c( ls()))
  # cp = load("~/SCHNAPPsDebug/sCA_getCells.RData")
  return(retVal)
}

#' define different methods for calculating diff. expressed genes
#' first entry is displayed in Radio box, second is function to be called.
myDiffExpFunctions <- list(
  c("Chi-square test of an estimated binomial distribution", "sCA_dge_CellViewfunc"),
  c("t-test", "sCA_dge_ttest"),
  c("DESeq2", "sCA_dge_deseq2"),
  c("seurat:wilcox", "sCA_dge_s_wilcox"),
  c("seurat:bimod", "sCA_dge_s_bimod"),
  c("seurat:t-test", "sCA_dge_s_t"),
  c("seurat:LR", "sCA_dge_s_LR"),
  c("seurat:neg-binomial", "sCA_dge_s_negbinom"),
  c("seurat:Poisson", "sCA_dge_s_poisson")
)






# ! TODO
# some dge functions are based on counts most on log-counts
# this is hard-coded here, but should be a parameter somehow like myDiffExpFunctions
# such it can be used by contributed functions


#' Seurat FindMarkers
#'
#' cellMeta = colData(scEx)
sCA_seuratFindMarkers <- function(scEx, cells.1, cells.2, test="wilcox", normFact = 1){
  if (DEBUG) cat(file = stderr(), "sCA_seuratFindMarkers started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_seuratFindMarkers")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_seuratFindMarkers")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_seuratFindMarkers", id = "sCA_seuratFindMarkers", duration = NULL)
  }
  if (!"DESeq2" %in% rownames(installed.packages())) {
    warning("Please install DESeq2 - learn more at https://bioconductor.org/packages/release/bioc/html/DESeq2.html")
    showNotification("Please install DESeq2", id = "sCA_dge_deseq2NOTFOUND", duration = NULL, type = "error")
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_seuratFindMarkers.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/sCA_seuratFindMarkers.RData')

  cellMeta <- colData(scEx)
  rData <- rowData(scEx)
  meta.data <- cellMeta[, "sampleNames", drop = FALSE]
  # creates object @assays$RNA@data and @assays$RNA@counts
  seurDat <- CreateSeuratObject(
    counts = assays(scEx)[[1]],
    meta.data = meta.data
  )
  # we remove e.g. "genes" from total seq (CD3-TotalSeqB)
  useGenes = which(rownames(seurDat@assays$RNA@data) %in% rownames(as(assays(scEx)[[1]], "dgCMatrix")))
  seurDat@assays$RNA@data = as(assays(scEx)[[1]], "dgCMatrix")[useGenes,]
  
  markers <- Seurat::FindMarkers(seurDat@assays$RNA@data/normFact, 
                                 cells.1 = cells.1,
                                 cells.2 = cells.2,
                                 min.pct = 0,
                                 test.use = test,
                                 logfc.threshold = 0.001
                                 # test.use = "wilcox" # p_val  avg_logFC pct.1 pct.2    p_val_adj
                                 # test.use = "bimod"  # p_val  avg_logFC pct.1 pct.2    p_val_adj
                                 # test.use = "roc"    # myAUC   avg_diff power pct.1 pct.2
                                 # test.use = "t"      # p_val avg_logFC pct.1 pct.2 p_val_adj
                                 # test.use = "negbinom" # needs UMI; p_val  avg_logFC pct.1 pct.2    p_val_adj
                                 # test.use = "poisson"  # needs UMI; p_val  avg_logFC pct.1 pct.2    p_val_adj
                                 # test.use = "LR"     # p_val  avg_logFC pct.1 pct.2 p_val_adj
                                 # test.use = "MAST" # not working: Assay in position 1, with name et is unlogged. Set `check_sanity = FALSE` to override and then proceed with caution.
                                 # test.use = "DESeq2" # needs UMI # done separately because the estimating process isn't working with 0s
  )
  if (nrow(markers) > 0) {
    markers$symbol <- rData[rownames(markers), "symbol"]
  }
  if ("Description" %in% colnames(rData)) {
    markers$Description <- rData[rownames(markers), "Description"]
  }
  return(markers)
}


sCA_dge_s_wilcox <- function(scEx_log, cells.1, cells.2){
  normFact = .schnappsEnv$normalizationFactor
  sCA_seuratFindMarkers(scEx_log, cells.1, cells.2, test="wilcox", normFact)
}
sCA_dge_s_bimod <- function(scEx_log, cells.1, cells.2){
  normFact = .schnappsEnv$normalizationFactor
  sCA_seuratFindMarkers(scEx_log, cells.1, cells.2, test="bimod")
}
sCA_dge_s_t <- function(scEx_log, cells.1, cells.2){
  normFact = .schnappsEnv$normalizationFactor
  sCA_seuratFindMarkers(scEx_log, cells.1, cells.2, test="t", normFact)
}
sCA_dge_s_LR<- function(scEx_log, cells.1, cells.2){
  normFact = .schnappsEnv$normalizationFactor
  sCA_seuratFindMarkers(scEx_log, cells.1, cells.2, test="LR", normFact)
}
sCA_dge_s_negbinom <- function(scEx_log, cells.1, cells.2){
  sCA_seuratFindMarkers(scEx_log, cells.1, cells.2, test="negbinom", normFact = 1)
}
sCA_dge_s_poisson <- function(scEx_log, cells.1, cells.2){
  sCA_seuratFindMarkers(scEx_log, cells.1, cells.2, test="poisson", normFact = 1)
}



#' sCA_dge_deseq2
#' calculate dge using DESeq2
sCA_dge_deseq2 <- function(scEx_log, cells.1, cells.2) {
  if (DEBUG) cat(file = stderr(), "sCA_dge_deseq2 started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dge_deseq2")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_dge_deseq2")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dge_deseq2", id = "sCA_dge_deseq2", duration = NULL)
  }
  if (!"DESeq2" %in% rownames(installed.packages())) {
    warning("Please install DESeq2 - learn more at https://bioconductor.org/packages/release/bioc/html/DESeq2.html")
    showNotification("Please install DESeq2", id = "sCA_dge_deseq2NOTFOUND", duration = NULL, type = "error")
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dge_deseq2.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/sCA_dge_deseq2.RData')


  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  group.info$wellKey <- rownames(x = group.info)
  # TODO how to handle data / transformation ?
  # it uses non-transformed org data. What happens if we loaded normalized data?
  # can back tronsform the data in counts?
  if (is.null(.schnappsEnv$normalizationFactor)) {
    .schnappsEnv$normalizationFactor = 1
  }
  data.use = assays(scEx_log)[[1]][,rownames(group.info)]
  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data.use,
    colData = group.info,
    design = ~group
  )
  dds1 <- DESeq2::estimateSizeFactors(object = dds1, type = "poscounts")
  dds1 <- DESeq2::estimateDispersions(object = dds1, fitType = "local")
  dds1 <- DESeq2::nbinomWaldTest(object = dds1)
  res <- DESeq2::results(
    object = dds1,
    contrast = c("group", "Group1", "Group2"),
    alpha = 0.05
  )
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
  res$padj[is.na(res$padj)] <- 1
  res$pvalue[is.na(res$pvalue)] <- 1
  featureData <- rowData(scEx_log)

  to.return <- data.frame(
    p_val = res$padj,
    avg_diff = res$log2FoldChange, row.names = rownames(res), symbol = rownames(res)
  )


  return(to.return)
}



#' sCA_dge_CellViewfunc
#' calculate differentically expressed genes given 2 sets of cells
sCA_dge_CellViewfunc <- function(scEx_log, cells.1, cells.2) {
  if (DEBUG) cat(file = stderr(), "sCA_dge_CellViewfunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dge_CellViewfunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_dge_CellViewfunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dge_CellViewfunc", id = "sCA_dge_CellViewfunc", duration = NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dge_CellViewfunc.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/sCA_dge_CellViewfunc.RData')

  featureData <- rowData(scEx_log)
  scEx_log <- as.matrix(assays(scEx_log)[[1]])
  subsetExpression <- scEx_log[complete.cases(scEx_log[, union(cells.1, cells.2)]), ]
  genes.use <- rownames(subsetExpression)
  # expMean exponential mean
  dat <- subsetExpression[genes.use, cells.1]
  data.1 <- apply(dat, 1, function(x) expMean(x, .schnappsEnv$normalizationFactor))
  dat <- subsetExpression[genes.use, cells.2]
  data.2 <- apply(dat, 1, function(x) expMean(x, .schnappsEnv$normalizationFactor))
  total.diff <- (data.1 - data.2)

  genes.diff <- names(which(abs(total.diff) > .2))
  genes.use <- ainb(genes.use, genes.diff)

  retVal <-
    DiffExpTest(subsetExpression, cells.1, cells.2, genes.use = genes.use)
  retVal[is.na(retVal[, "p_val"]), ] <- 1
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
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_dge_ttest")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dge_ttest", id = "sCA_dge_ttest", duration = NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dge_ttest.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/sCA_dge_ttest.RData')

  featureData <- rowData(scEx_log)
  scEx_log <- as.matrix(assays(scEx_log)[[1]])
  subsetExpression <- scEx_log[complete.cases(scEx_log[, union(cells.1, cells.2)]), ]
  genes.use <- rownames(subsetExpression)

  p_val <- apply(subsetExpression, 1, function(x) t.test(x[cells.1], x[cells.2])$p.value)
  p_val[is.na(p_val)] <- 1
  dat <- subsetExpression[genes.use, cells.1]
  data.1 <- apply(dat, 1, function(x) expMean(x, normFactor = .schnappsEnv$normalizationFactor))
  dat <- subsetExpression[genes.use, cells.2]
  data.2 <- apply(dat, 1, function(x) expMean(x, normFactor = .schnappsEnv$normalizationFactor))
  avg_diff <- (data.1 - data.2)

  retVal <- data.frame(p_val = p_val, avg_diff = avg_diff, symbol = featureData[names(p_val), "symbol"])
  return(retVal)
}

# sCA_dge reactive ----
#' manage calculation for differential expression analysis
sCA_dge <- reactive({
  if (DEBUG) cat(file = stderr(), "sCA_dge started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dge")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_dge")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dge", id = "sCA_dge", duration = NULL)
    removeNotification(id = "dgewarning")
  }

  scEx_log <- scEx_log()
  scEx <- scEx()
  projections <- projections()
  # cl1 <- input$sCA_dgeClustersSelection
  # selectedCells <- isolate(dePanelCellSelection())
  # cellNs <- isolate(selectedCells$cellNames())
  # sampdesc <- isolate(selectedCells$selectionDescription())
  
  clicked <- input$updateDGEParameters
  selectedCells <- isolate(sCA_dataInp())
  cellNs <- isolate(selectedCells$cellNames())
  sampdesc <- isolate(selectedCells$selectionDescription())
  prj <- isolate(selectedCells$ProjectionUsed())
  prjVals <- isolate(selectedCells$ProjectionValsUsed())
  
  db1 <- isolate(input$db1)
  db2 <- isolate(input$db2)
  method <- isolate(input$sCA_dgeRadioButton)

  if (is.null(scEx_log) | is.null(projections)  || is.null(db1)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dge.RData", list = c(ls(), ".schnappsEnv"))
  }
  # load(file='~/SCHNAPPsDebug/sCA_dge.RData')
  # require(archivist)
  # library(tools)
  # lazyLoad = local({load("~/SCHNAPPsDebug/sCA_dge.RData"); environment()})
  # tools:::makeLazyLoadDB(lazyLoad, "Huge")
  # lazyLoad("Huge")
  # objNames <- ls()
  # browser()
  methodIdx <- ceiling(which(unlist(.schnappsEnv$diffExpFunctions) == method) / 2)
  dgeFunc <- .schnappsEnv$diffExpFunctions[[methodIdx]][2]
  gCells <- sCA_getCells(projections, cl1 = cellNs, db1, db2)
  
  # in case we need counts and not normalized counts
  if (dgeFunc %in% c("sCA_dge_deseq2", "sCA_dge_s_negbinom", "sCA_dge_s_poisson")) {
    scEx_log = scEx
  }
  retVal <- do.call(dgeFunc, args = list(
    scEx_log = scEx_log,
    cells.1 = gCells$c1, cells.2 = gCells$c2
  ))

  if (nrow(retVal) == 0) {
    if (DEBUG) cat(file = stderr(), "dge: nothing found\n")
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("dge: nothing found", id = "dgewarning", duration = 10, type = "warning")
    }
  }
  # update reactiveValue
  sCA_selectedDge$sCA_dgeTable <- retVal

  setRedGreenButton(
    vars = list(
      # c("sCA_dataInpSelected_cells", isolate(sCA_dataInp()$selectedCells())),
      c("db1", isolate(input$db1)),
      c("db2", isolate(input$db2)),
      c("sCA_dgeRadioButton", isolate(input$sCA_dgeRadioButton)),
      c("sCA_dataInput-Mod_PPGrp", prjVals),
      c("sCA_dataInput-Mod_clusterPP", prj)
      
    ),
    button = "updateDGEParameters"
  )
  
  exportTestValues(sCA_dge = {
    retVal
  })
  return(retVal)
})

# using these global variables allows us to store a set values even when the projections are changing
.schnappsEnv$subClusterDim1 <- "PC1"
.schnappsEnv$subClusterDim2 <- "PC2"
.schnappsEnv$subClusterClusters <- NULL

observe(label = "ob7", {
  if (DEBUG) cat(file = stderr(), "observe: sCA_subscluster_x1\n")
  .schnappsEnv$subClusterDim1 <- input$sCA_subscluster_x1
})

observe(label = "ob8", {
  if (DEBUG) cat(file = stderr(), "observe: sCA_subscluster_y1\n")
  .schnappsEnv$subClusterDim2 <- input$sCA_subscluster_y1
})

#' TODO
#' if this observer is really needed we need to get rid of projections
observe(label = "ob9", {
  if (DEBUG) cat(file = stderr(), "observe: projections\n")
  projections <- projections()
  if (!is.null(projections)) {
    noOfClusters <- levels(as.factor(projections$dbCluster))
    # noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
    if (is.null(.schnappsEnv$subClusterClusters)) {
      .schnappsEnv$subClusterClusters <- noOfClusters
    }
  }
})

observe(label = "ob10", {
  if (DEBUG) cat(file = stderr(), "observe: sCA_dgeClustersSelection\n")
  .schnappsEnv$subClusterClusters <- input$sCA_dgeClustersSelection
})


# subcluster axes ----
# update axes in subcluster analysis
observeEvent(projections(), {
  if (DEBUG) cat(file = stderr(), "updateInputSubclusterAxes started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "updateInputSubclusterAxes")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "updateInputSubclusterAxes")
    }
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
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/updateInputSubclusterAxes.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/updateInputSubclusterAxes.RData")
  # if (length(gn) > 0) {
  #   projections <- cbind(projections, gn[rownames(projections), ] * 1)
  # }
  # Can also set the label and select items
  updateSelectInput(session, "sCA_subscluster_x1",
    choices = colnames(projections),
    selected = .schnappsEnv$subClusterDim1
  )

  updateSelectInput(session, "sCA_subscluster_y1",
    choices = colnames(projections),
    selected = .schnappsEnv$subClusterDim2
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
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "subCluster2Dplot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("subCluster2Dplot", id = "subCluster2Dplot", duration = NULL)
  }

  renderPlot({
    if (DEBUG) cat(file = stderr(), "output$sCA_dge_plot2\n")

    projections <- projections()
    x1 <- input$sCA_subscluster_x1
    y1 <- input$sCA_subscluster_y1
    # c1 <- input$sCA_dgeClustersSelection
    # selectedCells <- isolate(dePanelCellSelection())
    # cellNs <- isolate(selectedCells$cellNames())
    # sampdesc <- isolate(selectedCells$selectionDescription())
    selectedCells <- sCA_dataInp()
    cellNs <- selectedCells$cellNames()
    sampdesc <- selectedCells$selectionDescription()
    prjs <- selectedCells$ProjectionUsed()
    if (is.null(projections)) {
      return(NULL)
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/sCA_dge_plot2.RData", list = c(ls()))
    }
    # cp = load(file="~/SCHNAPPsDebug/sCA_dge_plot2.RData")

    subsetData <- projections[cellNs,]
#     xAxis <- list(
#       title = x1,
#       titlefont = f
#     )
#     yAxis <- list(
#       title = y1,
#       titlefont = f
#     )
#     
#     p1 <- plotly::plot_ly(
#       data = subsetData, source = "subset",
#       key = rownames(subsetData)
#     ) %>%
#       add_trace(
#         x = ~ get(x1),
#         # x = ~ get(dimX),
#         # y = ~ get(dimY),
#         y = ~ get(y1),
#         type = "scatter", mode = "markers",
#         text = ~ paste(1:nrow(subsetData), " ", rownames(subsetData), "<br />", subsetData$exprs),
#         # color = ~ get(dimCol),
#         colors = colors,
#         showlegend = TRUE
#       ) %>%
#       layout(
#         xaxis = xAxis,
#         yaxis = yAxis,
#         # title = "gtitle",
#         dragmode = "select"
#       )
# p1
    p1 <-
      ggplot(subsetData,
             aes_string(x = x1, y = y1),
             color = prjs
      ) +
      geom_point(aes(colour = get(prjs))) +
      geom_point(
        shape = 1,
        size = 4,
        aes(colour = get(prjs))
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
      ) + ggtitle(sampdesc) 
    if (is.factor(subsetData[,x1])) {
      p1 <- p1 + scale_x_discrete(drop=FALSE) 
    }
    if (is.factor(subsetData[,y1])) {
      p1 <- p1 + scale_y_discrete(drop=FALSE) 
    }
    
    # + scale_y_discrete(drop=FALSE)
    p1
  })
}

# save to history violoin observer ----
observe({
  clicked  = input$save2HistVolc
  if (DEBUG) cat(file = stderr(), "observe input$save2HistVolc \n")
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
  
  if (is.null(clicked)) return()
  if (clicked < 1) return()
  add2history(type = "renderPlotly", input = input, comment = "volcano plot",  
              plotData = .schnappsEnv[["sCA_volcanoPlot"]])
  
})


