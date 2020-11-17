# source(paste0(packagePath, "/reactives.R"), local = TRUE)

# since DE_scaterPNG is not used frequently it is not included in the heavyCalculations
# list
# myHeavyCalculations = list(c("DE_scaterPNG", "DE_scaterPNG"))

# Expression ------------------------------------------------------------------
callModule(
  clusterServer,
  "DE_expclusters",
  deProjTable,
  reactive(input$DE_gene_id)
)

# DE_updateInputExpPanel ----
#' DE_updateInputExpPanel
#' update x/y coordinates that can be chosen based on available
#' projections

.schnappsEnv$DE_X1 <- "tsne1"
.schnappsEnv$DE_Y1 <- "tsne1"
observe(label = "ob17", {
  if (DEBUG) cat(file = stderr(), "observe: DE_dim_x\n")
  .schnappsEnv$DE_X1 <- input$DE_dim_x
})
observe(label = "ob18", {
  if (DEBUG) cat(file = stderr(), "observe: DE_dim_y\n")
  .schnappsEnv$DE_Y1 <- input$DE_dim_y
})

observe({
  if (DEBUG) cat(file = stderr(), "DE_updateInputExpPanel started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_updateInputExpPanel")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_updateInputExpPanel")
    }
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
  # return(TRUE)
})

observe({
  scEx <- scEx()
  setRedGreenButtonCurrent(
    vars = list(
      c("scaterRan", 0)
    )
  )
  
  updateButtonColor(buttonName = "runScater", parameters = c(
    "scaterRan"
  ))
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

  selectedCells <- DE_Exp_dataInput()
  cellNs <- selectedCells$cellNames()
  sampdesc <- selectedCells$selectionDescription()
  prj <- isolate(selectedCells$ProjectionUsed())
  prjVals <- isolate(selectedCells$ProjectionValsUsed())
  
  if (is.null(scEx_log) | is.null(projections) | is.null(cellNs)) {
    if (DEBUG) cat(file = stderr(), "output$DE_gene_vio_plot:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_gene_vio_plot.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/DE_gene_vio_plot.RData")


  p1 <- DE_geneViolinFunc(scEx_log[, cellNs], g_id, projections[cellNs, ], ccols)

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
observe(label = "ob19", {
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



dePanelCellSelection <- callModule(
  cellSelectionModule,
  "DE_PanelPlotCellSelection"
)

DE_Exp_dataInput <- callModule(
  cellSelectionModule,
  "DE_Exp_dataInput"
)

# Panel plot button observer ----
#observe: cellNameTable_rows_selected ----
  observe(label = "ob_panelPlotParams", {
    if (DEBUG) cat(file = stderr(), "observe ob_panelPlotParams\n")
    
    input$updatePanelPlot
    setRedGreenButtonCurrent(
      vars = list(
        c("DE_panelplotids", (input$DE_panelplotids)),
        c("DE_dim_x", (input$DE_dim_x)),
        c("DE_dim_y", (input$DE_dim_y)),
        c("DE_panelplotSameScale", (input$DE_panelplotSameScale)),
        c("DE_nCol", (input$DE_nCol)),
        c("DE_PanelPlotCellSelection-Mod_clusterPP", (input$`DE_PanelPlotCellSelection-Mod_clusterPP`)),
        c("DE_PanelPlotCellSelection-Mod_PPGrp", (input$`DE_PanelPlotCellSelection-Mod_PPGrp`))
      )
    )
    
    updateButtonColor(buttonName = "updatePanelPlot", parameters = c(
      "DE_panelplotids", "DE_dim_x", "DE_dim_y", "DE_panelplotSameScale", 
      "DE_nCol", "DE_PanelPlotCellSelection-Mod_clusterPP", "DE_PanelPlotCellSelection-Mod_PPGrp"
    ))
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
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_panelPlot")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_panelPlot", id = "DE_panelPlot", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$DE_panelPlot\n")

  clicked <- input$updatePanelPlot
  scEx_log <- scEx_log()
  projections <- projections()
  DE_updateInputPPt()
  genesin <- isolate(input$DE_panelplotids)
  # cl4 <- input$DE_clusterSelectionPanelPlot
  # ppgrp <- isolate(input$DE_PPGrp)
  # ppCluster <- isolate(input$DE_clusterPP)

  selectedCells <- isolate(dePanelCellSelection())
  cellNs <- isolate(selectedCells$cellNames())
  sampdesc <- isolate(selectedCells$selectionDescription())

  dimx4 <- isolate(input$DE_dim_x)
  dimy4 <- isolate(input$DE_dim_y)
  sameScale <- isolate(input$DE_panelplotSameScale)
  nCol <- isolate(as.numeric(input$DE_nCol))

  if (is.null(scEx_log) | is.null(projections) | is.null(cellNs)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_panelPlot.RData", list = c("scEx_log", "projections", "genesin",
                                                               "dimx4", "dimy4", "sameScale", "nCol", "sampdesc" , "cellNs"
                                                               ), compress = F)
  }
  # cp = load(file="~/SCHNAPPsDebug/DE_panelPlot.RData")

  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")
  genesin <- genesin[[1]]

  if (DEBUG) cat(file = stderr(), paste("output:sampdesc",sampdesc,"\n"))
  retVal <- panelPlotFunc(scEx_log, projections, genesin, dimx4, dimy4, sameScale, nCol, sampdesc, cellNs )

  setRedGreenButton(
    vars = list(
      c("DE_panelplotids", isolate(input$DE_panelplotids)),
      c("DE_dim_x", isolate(input$DE_dim_x)),
      c("DE_dim_y", isolate(input$DE_dim_y)),
      c("DE_panelplotSameScale", isolate(input$DE_panelplotSameScale)),
      c("DE_nCol", isolate(input$DE_nCol)),
      c("DE_PanelPlotCellSelection-Mod_clusterPP", isolate(input$`DE_PanelPlotCellSelection-Mod_clusterPP`)),
      c("DE_PanelPlotCellSelection-Mod_PPGrp", isolate(input$`DE_PanelPlotCellSelection-Mod_PPGrp`))
    ),
    button = "updatePanelPlot"
  )
  
  
  printTimeEnd(start.time, "DE_panelPlot")
  exportTestValues(DE_panelPlot = {
    ls()
  })
  af = panelPlotFunc
  # remove env because it is too big
  environment(af) = new.env(parent = emptyenv())
  
  .schnappsEnv[["DE_panelPlot"]] <- list(panelPlotFunc = af,
                                         scEx_log = scEx_log, 
                                         projections=projections, 
                                         genesin=genesin, dimx4=dimx4, 
                                         dimy4=dimy4, sameScale=sameScale, 
                                         nCol=nCol, sampdesc=sampdesc,
                                         cellNs=cellNs
                                         )
  retVal
})


# Scater QC ----
output$DE_scaterQC <- renderImage(deleteFile = T, {
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

# DE_tsne_plt ----
# tSNE plot within Data exploration - Expressoin
output$DE_tsne_plt <- plotly::renderPlotly({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_tsne_plt")
    }
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_tsne_plt", id = "DE_tsne_plt", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$DE_tsne_plt\n")

  scEx_log <- scEx_log()
  g_id <- input$DE_gene_id
  projections <- projections()
  selectedCells <- DE_Exp_dataInput()
  cellNs <- selectedCells$cellNames()
  sampdesc <- selectedCells$selectionDescription()
  prj <- isolate(selectedCells$ProjectionUsed())
  prjVals <- isolate(selectedCells$ProjectionValsUsed())
  
  if (is.null(scEx_log) | is.null(projections) | is.null(cellNs) ) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_tsne_plt.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/DE_tsne_plt.RData")

  retVal <- DE_dataExpltSNEPlot(scEx_log[,cellNs], g_id, projections[cellNs, ])

  printTimeEnd(start.time, "DE_dataExpltSNEPlot")
  exportTestValues(DE_dataExpltSNEPlot = {
    str(retVal)
  })
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
    pca <- pcaReact()
    # TODO should be taken from projections.
    tsne <- tsne()

    if (is.null(scEx)) {
      return(NULL)
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/DE_downloadPanel.RData", list = c(ls()))
    }
    # load(file='~/SCHNAPPsDebug/DE_downloadPanel.RData')

    reducedDims(scEx) <- SimpleList(PCA = pca$x, TSNE = tsne)
    assays(scEx)[["logcounts"]] <- assays(scEx_log)[[1]]
    colData(scEx)[["before.Filter"]] <- projections$before.filter
    colData(scEx)[["dbCluster"]] <- projections$dbCluster
    colData(scEx)[["UmiCountPerGenes"]] <- projections$UmiCountPerGenes
    colData(scEx)[["UmiCountPerGenes2"]] <- projections$UmiCountPerGenes2

    save(file = file, list = c("scEx"))
    if (DEBUG) cat(file = stderr(), paste("DE_downloadPanel:done \n"))

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
#     save(file = "~/SCHNAPPsDebug/DE_SelectionText.RData", list = c(ls()))
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

# observers for input parameters and

observe(label = "observe DE_seuratRefBased", {
  if (DEBUG) cat(file = stderr(), "observe DE_seuratRefBased\n")
  if (is.null(input$updateNormalization)) {
    return(NULL)
  }
  if (!input$normalizationRadioButton == "DE_seuratRefBased") {
    return(NULL)
  }
  out <- DE_seuratRefBased()
  if (is.null(out)) {
    # set one parameter to something not possible to deactivate button/or choose different
    .schnappsEnv$DE_seuratRefBased_nfeatures <- "NA"
  }

  # currentValues
 setRedGreenButtonCurrent(
    vars = list(
      c("DE_seuratRefBased_nfeatures", input$DE_seuratRefBased_nfeatures),
      c("DE_seuratRefBased_k.filter", input$DE_seuratRefBased_k.filter),
      c("DE_seuratRefBased_scaleFactor", input$DE_seuratRefBased_scaleFactor)
    )
  )
  
  updateButtonColor(buttonName = "updateNormalization", parameters = c(
    "DE_seuratRefBased_nfeatures", "DE_seuratRefBased_k.filter",
    "DE_seuratRefBased_scaleFactor"
  ))


  # Here, we create the actual button
  # output$updateNormalizationButton <- updateButtonUI(input = input, name = "updateNormalization",
  #                                                    variables = c("nfeatures", "k.filter", "scalingFactor"  ) )
})


# obs.updateNormalization DE_seuratRefBased ----
observe(label = "obs.updateNormalization", {
  buttonPressed <- input$updateNormalization
  radioButtonVal <- isolate(input$normalizationRadioButton)
  if (!exists("DE_seuratRefBasedButtonOldVal", envir = .schnappsEnv)) {
    .schnappsEnv$DE_seuratRefBasedButtonOldVal <- 0
  }
  if (is.null(radioButtonVal)) {
    radioButtonVal <- ""
  }
  if (is.null(buttonPressed)) {
    buttonPressed <- 0
  }

  # changing the reactive DE_logGeneNormalizationButton will trigger the recalculation
  if (radioButtonVal == "DE_seuratRefBased" &
    !.schnappsEnv$DE_seuratRefBasedButtonOldVal == buttonPressed) {
    cat(file = stderr(), green(paste("\n=====changing value\n")))
    DE_seuratRefBasedButton(buttonPressed)
    .schnappsEnv$DE_seuratRefBasedButtonOldVal <- buttonPressed
  }
})

# obs.updateNormalization DE_seuratSCtransformButton ----
observe(label = "ob DE_seuratSCtransformButtonOldVal", {
  buttonPressed <- input$updateNormalization
  radioButtonVal <- isolate(input$normalizationRadioButton)
  if (!exists("DE_seuratSCtransformButtonOldVal", envir = .schnappsEnv)) {
    .schnappsEnv$DE_seuratSCtransformButtonOldVal <- 0
  }
  if (is.null(radioButtonVal)) {
    radioButtonVal <- ""
  }
  if (is.null(buttonPressed)) {
    buttonPressed <- 0
  }

  # changing the reactive DE_logGeneNormalizationButton will trigger the recalculation
  if (radioButtonVal == "DE_seuratSCtransform" &
    !.schnappsEnv$DE_seuratSCtransformButtonOldVal == buttonPressed) {
    cat(file = stderr(), green(paste("\n=====changing value\n")))
    DE_seuratSCtransformButton(buttonPressed)
    .schnappsEnv$DE_seuratSCtransformButtonOldVal <- buttonPressed
  }
})

observe(label = "observe DE_seuratSCtransform", {
  if (DEBUG) cat(file = stderr(), "observe DE_seuratSCtransform\n")
  if (is.null(input$updateNormalization)) {
    return(NULL)
  }
  if (!input$normalizationRadioButton == "DE_seuratSCtransform") {
    return(NULL)
  }
  out <- DE_seuratSCtransform()
  if (is.null(out)) {
    # set one parameter to something not possible to deactivate button/or choose different
    .schnappsEnv$DE_seuratRefBased_nfeatures <- "NA"
  }

  setRedGreenButtonCurrent(
    vars = list(
      c("DE_seuratSCtransform_nfeatures", input$DE_seuratSCtransform_nfeatures),
      c("DE_seuratSCtransform_k.filter", input$DE_seuratSCtransform_k.filter),
      c("DE_seuratSCtransform_scaleFactor", input$DE_seuratSCtransform_scaleFactor)
    )
  )

  updateButtonColor(buttonName = "updateNormalization", parameters = c(
    "DE_seuratSCtransform_nfeatures",
    "DE_seuratSCtransform_k.filter",
    "DE_seuratSCtransform_scaleFactor"
  ))
})


# obs.updateNormalization DE_seuratSCtransformButton ----
observe(label = "ob DE_seuratStandardButton", {
  buttonPressed <- input$updateNormalization
  radioButtonVal <- isolate(input$normalizationRadioButton)
  if (!exists("DE_seuratStandardButtonOldVal", envir = .schnappsEnv)) {
    .schnappsEnv$DE_seuratStandardButtonOldVal <- 0
  }
  if (is.null(radioButtonVal)) {
    radioButtonVal <- ""
  }
  if (is.null(buttonPressed)) {
    buttonPressed <- 0
  }
  
  # changing the reactive DE_logGeneNormalizationButton will trigger the recalculation
  if (radioButtonVal == "DE_seuratStandard" &
      !.schnappsEnv$DE_seuratStandardButtonOldVal == buttonPressed) {
    cat(file = stderr(), green(paste("\n=====changing value\n")))
    DE_seuratStandardButton(buttonPressed)
    .schnappsEnv$DE_seuratStandardButtonOldVal <- buttonPressed
  }
})

# obs.updateNormalization DE_seuratLogNormButton ----
observe(label = "ob DE_seuratLogNormButton", {
  buttonPressed <- input$updateNormalization
  radioButtonVal <- isolate(input$normalizationRadioButton)
  if (!exists("DE_seuratLogNormButtonOldVal", envir = .schnappsEnv)) {
    .schnappsEnv$DE_seuratLogNormButtonOldVal <- 0
  }
  if (is.null(radioButtonVal)) {
    radioButtonVal <- ""
  }
  if (is.null(buttonPressed)) {
    buttonPressed <- 0
  }
  
  # changing the reactive DE_logGeneNormalizationButton will trigger the recalculation
  if (radioButtonVal == "DE_seuratLogNorm" &
      !.schnappsEnv$DE_seuratLogNormButtonOldVal == buttonPressed) {
    cat(file = stderr(), green(paste("\n=====changing value\n")))
    DE_seuratLogNormButton(buttonPressed)
    .schnappsEnv$DE_seuratLogNormButtonOldVal <- buttonPressed
  }
})

observe(label = "observe DE_seuratStandard", {
  if (DEBUG) cat(file = stderr(), "observe DE_seuratStandard\n")
  if (is.null(input$updateNormalization)) {
    return(NULL)
  }
  if (!input$normalizationRadioButton == "DE_seuratStandard") {
    return(NULL)
  }
  out <- DE_seuratStandard()
  if (is.null(out)) {
    # set one parameter to something not possible to deactivate button/or choose different
    .schnappsEnv$DE_seuratRefBased_nfeatures <- "NA"
  }

  setRedGreenButtonCurrent(
    vars = list(
      c("DE_seuratStandard_dims", input$DE_seuratStandard_dims),
      c("DE_seuratStandard_anchorF", input$DE_seuratStandard_anchorF),
      c("DE_seuratStandard_kF", input$DE_seuratStandard_kF),
      c("DE_seuratStandard_k.weight", input$DE_seuratStandard_k.weight)
    )
  )

  updateButtonColor(buttonName = "updateNormalization", parameters = c(
    "DE_seuratStandard_dims",
    "DE_seuratStandard_anchorF",
    "DE_seuratStandard_kF",
    "DE_seuratStandard_k.weight"
  ))
})

# obs.updateNormalization DE_logNormalization ----
observe(label = "obs.updateNormalization", {
  buttonPressed <- input$updateNormalization
  radioButtonVal <- isolate(input$normalizationRadioButton)
  if (!exists("DE_logNormalizationButtonOldVal", envir = .schnappsEnv)) {
    .schnappsEnv$DE_logNormalizationButtonOldVal <- 0
  }
  if (is.null(radioButtonVal)) {
    radioButtonVal <- ""
  }
  if (is.null(buttonPressed)) {
    buttonPressed <- 0
  }
  
  # changing the reactive DE_logGeneNormalizationButton will trigger the recalculation
  if (radioButtonVal == "DE_logNormalization" &
      !.schnappsEnv$DE_logNormalizationButtonOldVal == buttonPressed) {
    cat(file = stderr(), green(paste("\n=====changing value\n")))
    DE_logNormalizationButton(buttonPressed)
    .schnappsEnv$DE_logNormalizationButtonOldVal <- buttonPressed
  }
})

# if no parameters the second observer is not needed


# obs.updateNormalization DE_logGeneNormalization ----
observe(label = "obs.updateNormalization", {
  buttonPressed <- input$updateNormalization
  radioButtonVal <- isolate(input$normalizationRadioButton)
  if (!exists("DE_logGeneNormalizationButtonOldVal", envir = .schnappsEnv)) {
    .schnappsEnv$DE_logGeneNormalizationButtonOldVal <- 0
  }
  if (is.null(radioButtonVal)) {
    radioButtonVal <- ""
  }
  if (is.null(buttonPressed)) {
    buttonPressed <- 0
  }

  # changing the reactive DE_logGeneNormalizationButton will trigger the recalculation
  if (radioButtonVal == "DE_logGeneNormalization" &
    !.schnappsEnv$DE_logGeneNormalizationButtonOldVal == buttonPressed) {
    cat(file = stderr(), green(paste("\n=====changing value\n")))
    DE_logGeneNormalizationButton(buttonPressed)
    .schnappsEnv$DE_logGeneNormalizationButtonOldVal <- buttonPressed
  }
})

# observe parameters for logGene ----
observe(label = "oblogGene", {
  temp <- input$DE_geneIds_norm
  if (DEBUG) green(cat(file = stderr(), "observe DE_logGeneNormalization\n"))
  if (is.null(input$updateNormalization)) {
    if (DEBUG) cat(file = stderr(), "observe DE_logGeneNormalization, input$updateNormalization NULL\n")
    return(NULL)
  }
  if (!input$normalizationRadioButton == "DE_logGeneNormalization") {
    if (DEBUG) {
      cat(file = stderr(), paste(
        "observe DE_logGeneNormalization, input$normalizationRadioButton good",
        input$updateNormalization, "\n"
      ))
    }
    return(NULL)
  }
  out <- isolate(DE_logGeneNormalization())

  if (is.null(out)) {
    # set one parameter to something not possible to deactivate button/or choose different
    .schnappsEnv$calculated_DE_geneIds_norm <- "NOT AVAILABLE"
  }

  setRedGreenButtonCurrent(
    vars = list(
      c("DE_geneIds_norm", input$DE_geneIds_norm)
    )
  )
 # Here, we create the actual button
  updateButtonColor(buttonName = "updateNormalization", parameters = c("DE_geneIds_norm"))
  # output$updateNormalizationButton <- updateButtonUI(input = input, name = "updateNormalization",
  #                                                    variables = c("DE_geneIds_norm"))
})

# observer for normalization button ----
# checks the radio button for changes.
observe(label = "ob12", {
  if (DEBUG) cat(file = stderr(), "observe normalizationRadioButton\n")
  out <- scEx_log()
  radioButtonValue <- input$normalizationRadioButton

  if (is.null(out)) {
    # set one parameter to something not possible to deactivate button/or choose different
    .schnappsEnv$calculated_normalizationRadioButton <- "NA"
  }
  if (DEBUG) {
    cat(file = stderr(), paste(
      "observe normalizationRadioButton: ",
      radioButtonValue, "\n"
    ))
  }
  assign("normalizationRadioButton", radioButtonValue, envir = .schnappsEnv)
  # Here, we create the actual button
  updateButtonColor(buttonName = "updateNormalization", parameters = c("normalizationRadioButton"))

  # output$updateNormalizationButton <- updateButtonUI(input = input, name = "updateNormalization",
  #                                                    variables = c("normalizationRadioButton"))
})
