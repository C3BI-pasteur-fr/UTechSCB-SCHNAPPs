source(paste0(packagePath,  "/reactives.R"), local = TRUE)

# since DE_scaterPNG is not used frequently it is not included in the heavyCalculations
# list
# myHeavyCalculations = list(c("DE_scaterPNG", "DE_scaterPNG"))

# Expression ------------------------------------------------------------------
callModule(
  clusterServer,
  "DE_expclusters",
  projections,
  reactive(input$DE_gene_id)
)

# DE_updateInputExpPanel ----
#' DE_updateInputExpPanel
#' update x/y coordinates that can be chosen based on available
#' projections

.schnappsEnv$DE_X1 <- "tsne1"
.schnappsEnv$DE_Y1 <- "tsne1"
observe({
  if (DEBUG) cat(file = stderr(), "observe: DE_dim_x\n")
  .schnappsEnv$DE_X1 <- input$DE_dim_x
})
observe({
  if (DEBUG) cat(file = stderr(), "observe: DE_dim_y\n")
  .schnappsEnv$DE_Y1 <- input$DE_dim_y
})

DE_updateInputExpPanel <- reactive({
  if (DEBUG) cat(file = stderr(), "DE_updateInputExpPanel started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_updateInputExpPanel")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_updateInputExpPanel")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_updateInputExpPanel", id = "DE_updateInputExpPanel", duration = NULL)
  }
  
  projections <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  
  # Can also set the label and select items
  updateSelectInput(session, "DE_dim_x",
                    choices = colnames(projections),
                    selected = .schnappsEnv$DE_X1
  )
  
  # Can also set the label and select items
  updateSelectInput(session, "DE_dim_y",
                    choices = colnames(projections),
                    selected = .schnappsEnv$DE_Y1
  )
  return(TRUE)
})




# EXPLORE TAB VIOLIN PLOT ----
# TODO module for violin plot  ??
output$DE_gene_vio_plot <- renderPlot({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_gene_vio_plot")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_gene_vio_plot", id = "DE_gene_vio_plot", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$DE_gene_vio_plot\n")
  
  scEx_log <- scEx_log()
  projections <- projections()
  g_id <- input$DE_gene_id
  ccols <- clusterCols$colPal
  
  if (is.null(scEx_log) | is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$DE_gene_vio_plot:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_gene_vio_plot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/DE_gene_vio_plot.RData")
  
  
  p1 <- DE_geneViolinFunc(scEx_log, g_id, projections, ccols)
  
  printTimeEnd(start.time, "DE_gene_vio_plot")
  exportTestValues(DE_gene_vio_plot = {
    p1
  })
  return(p1)
})


### Panel Plot ----
#' DE_clusterSelectionPanelPlot
#' update selection options for clusters
#' since we allow also "all" we have to have a different strategy for the input
#' it is debateable whether this is usefull to have a different strategy, but for now
#' we leave it as it.
.schnappsEnv$DE_cl1 <- "All"
observe({
  if (DEBUG) cat(file = stderr(), "observe: DE_clusterSelectionPanelPlot\n")
  .schnappsEnv$DE_cl1 <- input$DE_clusterSelectionPanelPlot
})
output$DE_clusterSelectionPanelPlot <- renderUI({
  if (DEBUG) cat(file = stderr(), "output$DE_clusterSelectionPanelPlot\n")
  projections <- projections()
  upI <- DE_updateInputExpPanel()
  if (is.null(projections)) {
    HTML("Please load data")
  } else {
    noOfClusters <- levels(as.factor(projections$dbCluster))
    selectInput(
      "DE_clusterSelectionPanelPlot",
      label = "Cluster",
      choices = c(c("All"), noOfClusters),
      selected = .schnappsEnv$DE_cl1
    )
  }
})

# DE_panelPlot ----
#' DE_panelPlot
#' plot multiple panels for a given list of genes
#' If the x-axis is a categorical value and the y-axis is UMI.counts the y-axis related to
#' the count for that gene. Otherwise, all genes are used.
#' normalized counts are used for plotting
output$DE_panelPlot <- renderPlot({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_panelPlot")
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_panelPlot", id = "DE_panelPlot", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$DE_panelPlot\n")
  
  scEx_log <- scEx_log()
  projections <- projections()
  genesin <- input$DE_panelplotids
  cl4 <- input$DE_clusterSelectionPanelPlot
  dimx4 <- input$DE_dim_x
  dimy4 <- input$DE_dim_y
  sameScale <- input$DE_panelplotSameScale
  nCol <- as.numeric(input$DE_nCol)
  
  if (is.null(scEx_log) | is.null(projections) | is.null(cl4)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_panelPlot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/DE_panelPlot.RData")
  
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")
  genesin <- genesin[[1]]
  
  featureData <- rowData(scEx_log)
  # featureData$symbol = toupper(featureData$symbol)
  genesin <- genesin[which(genesin %in% toupper(featureData$symbol))]
  if (length(genesin)<1) {return (NULL)}
  par(mfrow = c(ceiling(length(genesin) / 4), 4), mai = c(0., .3, .3, .3))
  rbPal <- colorRampPalette(c("#f0f0f0", "red"))
  ylim <- c(min(projections[, dimy4]), max(projections[, dimy4]))
  if (is(projections[, dimx4], "factor") & dimy4 == "UMI.count") {
    ymax <- 0
    for (i in 1:length(genesin)) {
      geneIdx <- which(toupper(featureData$symbol) == genesin[i])
      ymax <- max(ymax, max(Matrix::colSums(assays(scEx_log)[["logcounts"]][geneIdx, , drop = FALSE])))
    }
    ylim <- c(0, ymax)
    if(!sameScale){
      ylim <- NULL
    }
  }
  plotList = list()
  plotIdx = 0
  if (cl4 == "All") {
    for (i in 1:length(genesin)) {
      geneIdx <- which(toupper(featureData$symbol) == genesin[i])
      Col <- rbPal(10)[
        as.numeric(
          cut(
            as.numeric(
              assays(scEx_log)[[1]][
                rownames(featureData[geneIdx, ]),
                ]
            ),
            breaks = 10
          )
        )
        ]
      plotIdx = plotIdx +1
      
      plotList[[plotIdx]] = ggplot(projections, aes_string(x=dimx4, y=dimy4)) 
      if (is(projections[, dimx4], "factor") & dimy4 == "UMI.count") {
        projections[, dimy4] <- Matrix::colSums(assays(scEx_log)[["logcounts"]][geneIdx, , drop = FALSE])
        plotList[[plotIdx]] = plotList[[plotIdx]] + geom_boxplot(show.legend = FALSE) + ggtitle(genesin[i])
      } else{
        plotList[[plotIdx]] = plotList[[plotIdx]] + geom_point(color = Col, show.legend = FALSE) + ggtitle(genesin[i])
        
      }
      if (!is.null(ylim)) {
        plotList[[plotIdx]] = plotList[[plotIdx]]  + ylim(ylim)
      }
      # plot(projections[, dimx4], projections[, dimy4],
      #      col = Col, pch = 16, frame.plot = TRUE, ann = FALSE, ylim = ylim
      # )
      # title(genesin[i], line = -1.2, adj = 0.05, cex.main = 2)
      if (DEBUG) cat(file = stderr(), genesin[i])
    }
  } else {
    for (i in 1:length(genesin)) {
      geneIdx <- which(toupper(featureData$symbol) == genesin[i])
      subsetTSNE <- subset(projections, dbCluster == cl4)
      
      Col <- rbPal(10)[
        as.numeric(
          cut(
            as.numeric(
              assays(scEx_log)[[1]][
                rownames(featureData[geneIdx, ]),
                ]
            ),
            breaks = 10
          )
        )
        ]
      
      names(Col) <- rownames(projections)
      plotCol <- Col[rownames(subsetTSNE)]
      if (is(projections[, dimx4], "factor") & dimy4 == "UMI.count") {
        projections[, dimy4] <- Matrix::colSums(assays(scEx_log)[["logcounts"]][geneIdx, , drop = FALSE])
        subsetTSNE <- subset(projections, dbCluster == cl4)
      }
      
      plotIdx = plotIdx + 1
      
      plotList[[plotIdx]] = ggplot(subsetTSNE, aes_string(x=dimx4, y=dimy4)) 
      if (is(subsetTSNE[, dimx4], "factor") & dimy4 == "UMI.count") {
        subsetTSNE[, dimy4] <- Matrix::colSums(assays(scEx_log)[["logcounts"]][geneIdx, rownames(subsetTSNE), drop = FALSE])
        plotList[[plotIdx]] = plotList[[plotIdx]] + geom_boxplot(show.legend = FALSE) + ggtitle(genesin[i])
      } else{
        plotList[[plotIdx]] = plotList[[plotIdx]] + geom_point(color = plotCol, show.legend = FALSE) + ggtitle(genesin[i])
        
      }
      if (!is.null(ylim)) {
        plotList[[plotIdx]] = plotList[[plotIdx]]  + ylim(ylim)
      }
      # plot(subsetTSNE[, dimx4], subsetTSNE[, dimy4],
      #      col = plotCol, pch = 16, frame.plot = TRUE,
      #      ann = FALSE, ylim = ylim
      # )
      # title(genesin[i], line = -1.2, adj = 0.05, cex.main = 2)
      if (DEBUG) cat(file = stderr(), cl4)
    }
  }
  retVal = ggarrange(plotlist = plotList, ncol = nCol, nrow = ceil(length(plotList)/nCol))
  printTimeEnd(start.time, "DE_panelPlot")
  exportTestValues(DE_panelPlot = {ls()})
  .schnappsEnv[["DE_panelPlot"]] <- retVal
  retVal
  
})


# Scater QC ----
output$DE_scaterQC <- renderImage({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_scaterQC")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_scaterQC", id = "DE_scaterQC", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$DE_scaterQC\n")
  scaterReads <- scaterReads()
  if (is.null(scaterReads)) {
    return(NULL)
  }
  
  DE_scaterPNG()
})

# DE_tsne_plt ----
# tSNE plot within Data exploration - Expressoin
output$DE_tsne_plt <- plotly::renderPlotly({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_tsne_plt")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_tsne_plt", id = "DE_tsne_plt", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$DE_tsne_plt\n")
  
  scEx_log <- scEx_log()
  g_id <- input$DE_gene_id
  projections <- projections()
  
  if (is.null(scEx_log) | is.null(projections)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_tsne_plt.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/DE_tsne_plt.RData")
  
  retVal <- DE_dataExpltSNEPlot(scEx_log, g_id, projections)
  
  printTimeEnd(start.time, "DE_dataExpltSNEPlot")
  exportTestValues(DE_dataExpltSNEPlot = {str(retVal)})
  retVal
})

# download RDS ----
output$DE_downloadPanel <- downloadHandler(
  filename = paste0("panelPlot.", Sys.Date(), ".Zip"),
  content = function(file) {
    if (DEBUG) cat(file = stderr(), paste("DE_downloadPanel: \n"))
    
    scEx <- scEx()
    projections <- projections()
    scEx_log <- scEx_log()
    pca <- pca()
    tsne <- tsne()
    
    if (is.null(scEx)) {
      return(NULL)
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/RDSsave.RData", list = c(ls(), ls(envir = globalenv())))
    }
    # load(file='~/SCHNAPPsDebug/RDSsave.RData')
    
    reducedDims(scEx) <- SimpleList(PCA = pca$x, TSNE = tsne)
    assays(scEx)[["logcounts"]] = assays(scEx_log)[[1]]
    colData(scEx)[["before.Filter"]] = projections$before.filter
    colData(scEx)[["dbCluster"]] = projections$dbCluster
    colData(scEx)[["UmiCountPerGenes"]] = projections$UmiCountPerGenes
    colData(scEx)[["UmiCountPerGenes2"]] = projections$UmiCountPerGenes2
    
    save(file = file, list = c("scEx"))
    if (DEBUG) cat(file = stderr(), paste("RDSsave:done \n"))
    
    # write.csv(as.matrix(exprs(scEx)), file)
  }
)

# observe({
#   if (DEBUG) cat(file = stderr(), paste0("observe: dimension_x\n"))
#   .schnappsEnv$DE_dimension_y <- input$DE_dimension_y
# })
# observe({
#   updateSelectInput(session, "DE_dimension_y",
#                     choices = c(colnames(projections), "UmiCountPerGenes"),
#                     selected = .schnappsEnv$DE_dimension_y
#   )
#   
# })
# 
# output$DE_sortedPlot <- plotly::renderPlotly({
#   if (DEBUG) cat(file = stderr(), paste("Module: output$DE_sortedPlot\n"))
#   start.time <- base::Sys.time()
#   
#   # remove any notification on exit that we don't want
#   on.exit(
#     if (!is.null(getDefaultReactiveDomain())) {
#       removeNotification(id = "DE_sortedPlot")
#     }
#   )
#   # show in the app that this is running
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("DE_sortedPlot", id = "DE_sortedPlot", duration = NULL)
#   }
#   
#   scEx_log <- scEx_log()
#   tdata <- tData()
#   projections <- projections()
#   DE_dimension_y <- input$DE_dimension_y
#   
#   if (is.null(scEx_log) | is.null(scEx_log) | is.null(tdata)) {
#     if (DEBUG) cat(file = stderr(), paste("output$clusterPlot:NULL\n"))
#     return(NULL)
#   }
#   
#   
#   
#   # event_register(p1, 'plotly_selected')
#   printTimeEnd(start.time, "DE_sortedPlot")
#   exportTestValues(DE_sortedPlot = {
#     p1
#   })
#   suppressMessages(p1)
#   
# })
# 
# output$DE_SelectionText <- renderText({
#   if (DEBUG) cat(file = stderr(), "DE_SelectionText\n")
#   start.time <- base::Sys.time()
#   on.exit({
#     printTimeEnd(start.time, "DE_SelectionText")
#     if (!is.null(getDefaultReactiveDomain())) {
#       removeNotification(id = "DE_SelectionText")
#     }
#   })
#   # show in the app that this is running
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("DE_SelectionText", id = "DE_SelectionText", duration = NULL)
#   }
#   
#   selectedCells <- DE_selectedCells()$selectedCells
#   
#   if (is.null(selectedCells)) {
#     return(NULL)
#   }
#   if (.schnappsEnv$DEBUGSAVE) {
#     save(file = "~/SCHNAPPsDebug/DE_SelectionText.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file="~/SCHNAPPsDebug/DE_SelectionText.RData")
#   inpClusters <- levels(projections$dbCluster)
#   
#   retVal <- paste(selectedCells)
#   
#   exportTestValues(DummyReactive = {
#     retVal
#   })
#   return(retVal)
# })
# 
#   
#   