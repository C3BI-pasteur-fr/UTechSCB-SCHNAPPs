# TODO: make sure this is working
myZippedReportFiles <- c("DGE.csv")
require(Seurat)


observe(label ="obs_sCA_volc_pval", x = {
  .schnappsEnv$defaultValues[["sCA_volc_pval"]] = input$sCA_volc_pval
})
observe(label ="obs_sCA_volc_effectLimit", x = {
  .schnappsEnv$defaultValues[["sCA_volc_effectLimit"]] = input$sCA_volc_effectLimit
})
observe(label ="obs_scDEA_zingeR.edgeR", x = {
  .schnappsEnv$defaultValues[["scDEA_zingeR.edgeR"]] = input$scDEA_zingeR.edgeR
})
observe(label ="obs_scDEA_Seurat", x = {
  .schnappsEnv$defaultValues[["scDEA_Seurat"]] = input$scDEA_Seurat
})
observe(label ="obs_scDEA_limma", x = {
  .schnappsEnv$defaultValues[["scDEA_limma"]] = input$scDEA_limma
})
observe(label ="obs_scDEA_Wilcoxon", x = {
  .schnappsEnv$defaultValues[["scDEA_Wilcoxon"]] = input$scDEA_Wilcoxon
})
observe(label ="obs_scDEA_Ttest", x = {
  .schnappsEnv$defaultValues[["scDEA_Ttest"]] = input$scDEA_Ttest
})
observe(label ="obs_scDEA_scDD", x = {
  .schnappsEnv$defaultValues[["scDEA_scDD"]] = input$scDEA_scDD
})
observe(label ="obs_scDEA_monocle", x = {
  .schnappsEnv$defaultValues[["scDEA_monocle"]] = input$scDEA_monocle
})
observe(label ="obs_scDEA_MAST", x = {
  .schnappsEnv$defaultValues[["scDEA_MAST"]] = input$scDEA_MAST
})
observe(label ="obs_scDEA_DESeq2", x = {
  .schnappsEnv$defaultValues[["scDEA_DESeq2"]] = input$scDEA_DESeq2
})
observe(label ="obs_scDEA_DEsingle", x = {
  .schnappsEnv$defaultValues[["scDEA_DEsingle"]] = input$scDEA_DEsingle
})
observe(label ="obs_scDEA_BPSC", x = {
  .schnappsEnv$defaultValues[["scDEA_BPSC"]] = input$scDEA_BPSC
})
observe(label ="obs_scDEA_parallel", x = {
  .schnappsEnv$defaultValues[["scDEA_parallel"]] = input$scDEA_parallel
})
observe(label ="obs_sCA_dgeRadioButton", x = {
  .schnappsEnv$defaultValues[["sCA_dgeRadioButton"]] = input$sCA_dgeRadioButton
})
observe(label ="obs_sCA_subscluster_y1", x = {
  .schnappsEnv$defaultValues[["sCA_subscluster_y1"]] = input$sCA_subscluster_y1
})
observe(label ="obs_sCA_subscluster_x1", x = {
  .schnappsEnv$defaultValues[["sCA_subscluster_x1"]] = input$sCA_subscluster_x1
})


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
  
  if (is.null(scEx) | is.null(top.genes)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/output_dge.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/output_dge.RData")
  featureData <- rowData(scEx)
  
  top.genes$symbol <-
    featureData[rownames(top.genes), "symbol"]
  if ("Description" %in% colnames(featureData)) {
    top.genes$Description <- featureData[rownames(top.genes), "Description"]
  }
  rownames(top.genes) <- make.unique(as.character(top.genes$symbol), sep = "_#_")
  if (dim(top.genes)[1] > 0) {
    # change inf to high/low number
    infIdx <- which(is.infinite(top.genes$avg_diff))
    if (length(infIdx) > 0) {
      top.genes$avg_diff[infIdx[top.genes$avg_diff[infIdx] > 0]] <- 9999999
      top.genes$avg_diff[infIdx[top.genes$avg_diff[infIdx] < 0]] <- -9999999
      top.genes$Description[is.na(top.genes$Description)] <- ""
    }
    return(top.genes)
  } else {
    return(NULL)
  }
})

# dge table ----
callModule(
  tableSelectionServer,
  "sCA_dgeTable",
  sCA_dgeTableReac, caption = "Tables with differentially expressed genes"
)

sCA_dataInp <- callModule(
  cellSelectionModule,
  "sCA_dataInput"
)

# observe button ----
observeEvent(input$updateDGEParameters, handlerExpr = {
  if (DEBUG) cat(file = stderr(), "observe updateDGEParameters\n")
  sCA_dge()
  if (DEBUG) cat(file = stderr(), "observe updateDGEParameters end\n")
})

# observer for button ----
observe(label = "ob_DGEParams", {
  if (DEBUG) cat(file = stderr(), "observe DGEVars\n")
  input$updateDGEParameters
  selectedCells <- sCA_dataInp()
  prj <- selectedCells$ProjectionUsed()
  prjVals <- selectedCells$ProjectionValsUsed()
  
  setRedGreenButtonCurrent(
    vars = list(
      c("db1", input$db1),
      c("db2", input$db2),
      c("sCA_dgeRadioButton", input$sCA_dgeRadioButton),
      c("sCA_dataInput-Mod_PPGrp", prjVals),
      c("sCA_dataInput-Mod_clusterPP", prj)
    )
  )
  
  updateButtonColor(buttonName = "updateDGEParameters", parameters = c(
    "db1", "db2", "sCA_dgeRadioButton", 
    "sCA_dataInput-Mod_PPGrp", "sCA_dataInput-Mod_clusterPP"
  ))
  
})


# sub cluster analysis ( used for 2 panels )

# output$sCA_dgeClustersSelection <- renderUI({
#   if (DEBUG) cat(file = stderr(), "sCA_dgeClustersSelection started.\n")
#   start.time <- base::Sys.time()
#   on.exit({
#     printTimeEnd(start.time, "sCA_dgeClustersSelection")
#     if (!is.null(getDefaultReactiveDomain())) {
#       removeNotification(id = "sCA_dgeClustersSelection")
#     }
#   })
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("sCA_dgeClustersSelection", id = "sCA_dgeClustersSelection", duration = NULL)
#   }
# 
#   projections <- projections()
#   up1 <- updateInputSubclusterAxes()
# 
#   if (DEBUG) cat(file = stderr(), "output$sCA_dgeClustersSelection\n")
#   if (.schnappsEnv$DEBUGSAVE) {
#     save(file = "~/SCHNAPPsDebug/sCA_dgeClustersSelection.RData", list = c(ls(envir = globalenv(), ls(), ls(envir = .schnappsEnv))))
#   }
#   # load(file="~/SCHNAPPsDebug/sCA_dgeClustersSelection.RData")
# 
# 
#   if (is.null(projections)) {
#     tags$span(style = "color:red", "Please load data first")
#   } else {
#     noOfClusters <- levels(as.factor(projections$dbCluster))
#     # noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
#     selectizeInput(
#       "sCA_dgeClustersSelection",
#       label = "Cluster",
#       choices = noOfClusters,
#       selected = .schnappsEnv$subClusterClusters,
#       multiple = TRUE
#     )
#   }
# })


output$sCA_volc_selected <- renderText({
  if (DEBUG)
    cat(file = stderr(), "sCA_volc_selected started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_volc_selected")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_volc_selected")
    }
  })
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_volc_selected", id = "sCA_volc_selected", duration = NULL)
  }
  DGEdata <- sCA_dgeTableReac()
  req(DGEdata)
  brushedPs <- plotly::event_data("plotly_selected")
  
  if (is.null(brushedPs)) {
    if (DEBUG) cat(file = stderr(), "cluster: selectedCellNames: brush null\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_volc_selected.RData", list = ls())
  }
  # cp =load(file=paste0("~/SCHNAPPsDebug/sCA_volc_selected.RData"))
  # cells.names <- brushedPs$key
  # cells.names <- unique(cells.names[!is.na(cells.names)])
  # if (DEBUG) {
  #   cat(file = stderr(), paste("curveNumbers:", unique(brushedPs$curveNumber), "\n"))
  # }
  retVal <- paste(rownames(DGEdata)[brushedPs[brushedPs$curveNumber == 0, "pointNumber"] + 1], collapse = ", ")
  return(retVal)
})


### sCA_cells used started ----
observe({
  if (DEBUG) cat(file = stderr(), "sCA_cells used started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_cells_used")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_cells_used")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_cells used", id = "sCA_cells_used", duration = NULL)
    removeNotification(id = "dgewarning")
  }
  
  

  db1 <- isolate(input$db1)
  db2 <- isolate(input$db2)
  db1x <- isolate(input$sCA_subscluster_x1)
  db1y <- isolate(input$sCA_subscluster_y1)
  projections <- projections()
  if ( is.null(projections)  || is.null(db1)) {
    return(NULL)
  }
  
  selectedCells <- sCA_dataInp()
  sampdesc <- isolate(selectedCells$selectionDescription())
  prj <- isolate(selectedCells$ProjectionUsed())
  cellNs <- isolate(selectedCells$cellNames())
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_cells_used.RData", list = c(ls(), ".schnappsEnv"))
  }
  # cp = load(file='~/SCHNAPPsDebug/sCA_cells_used.RData')
  # browser()
  gCells <- sCA_getCells(projections, cl1 = cellNs, db1, db2, db1x, db1y)
  
  sCA_dge_Ncells(c(cells.1 = gCells$c1, cells.2 = gCells$c2))
})

output$sCA_dge_plot1_Ncells = renderText(
  sCA_dge_Ncells()[[1]]
)
output$sCA_dge_plot2_Ncells = renderText(
  sCA_dge_Ncells()[[2]]
)

output$sCA_volcanoPlot <- plotly::renderPlotly({
  require(manhattanly)
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
  pvalT <- input$sCA_volc_pval
  effT <- input$sCA_volc_effectLimit
  
  if (is.null(DGEdata)) {
    if (DEBUG) cat(file = stderr(), "output$sCA_volcanoPlot:NULL\n")
    return(NULL)
  }
  # reset selection
  # js$sCA_volcanoPlot_resetClick()
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_volcanoPlot.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/sCA_volcanoPlot.RData")
  if ("p_val_adj" %in% colnames(DGEdata)) {
    pval <- "p_val_adj"
  } else
    if ("p_val" %in% colnames(DGEdata)) {
      pval <- "p_val"
    }
  effect_size = NULL
  if ("avg_log2FC" %in% colnames(DGEdata)) {
    effect_size <- "avg_log2FC"
  }
  if ("avg_diff" %in% colnames(DGEdata)) {
    effect_size <- "avg_diff"
  }
  DGEdata[is.na(DGEdata[, pval]), pval] <- 1
  if (!is.null(effect_size)){
    DGEdata[is.na(DGEdata[, effect_size]), effect_size] <- 0
    DGEdata$EFFECTSIZE <- DGEdata[, effect_size]
  }
  DGEdata$P <- DGEdata[, pval]
  TEXT <- paste(paste("symbol: ", DGEdata$symbol), sep = "<br>")
  
  retVal <- volcanoly(DGEdata, snp = "symbol", genomewideline = pvalT, effect_size_line = c(-effT, effT)) %>%
    layout(
      dragmode = "select"
    )
  
  if (!is.null(effect_size)){
    # magrittr::%<>%
    retVal %<>% plotly::add_trace(
      x = DGEdata$EFFECTSIZE, y = -log10(DGEdata[, pval]),
      type = "scatter",
      mode = "markers",
      text = TEXT,
      opacity = 0,
      marker = list(size = 5),
      name = ""
    )
    # retVal
  }
  
  .schnappsEnv[["sCA_volcanoPlot"]] <- retVal
  
  exportTestValues(dgeVolcanoPlot = {
    str(retVal)
  })
  retVal
})
