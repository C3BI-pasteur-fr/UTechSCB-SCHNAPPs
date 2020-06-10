
myZippedReportFiles <- c("output_coE_topExpGenes.csv")


# All clusters heat map ------
callModule(
  pHeatMapModule,
  "coExpHeatmapModule",
  coE_heatmapReactive
)


# 2D plot with selection of cells ------
# assigning it to a variable allows us to interact with the plot and collect selection events
coE_selctedCluster <-
  callModule(
    clusterServer,
    "coE_selected",
    projections # ,
    # reactive(input$coE_gene_id_sch)
  )

# selected clusters heatmap module -----
callModule(
  pHeatMapModule,
  "coE_heatmapSelectedModule",
  coE_heatmapSelectedReactive
)

# max expressed genes ----
callModule(
  tableSelectionServer,
  "coE_topExpGenes",
  coE_topExpGenesTable
)

# max expressed genes ----
callModule(
  tableSelectionServer,
  "coE_topCCGenes",
  coE_topExpCCTable
)
# observe: updateHeatMapSelectedParameters ----
observe(label = "ob_HeatMapSelectedParams", {
  if (DEBUG) cat(file = stderr(), "observe HeatMapSelected\n")
  input$updateHeatMapSelectedParameters
  setRedGreenButtonCurrent(
    vars = list(
      c("coE_heatmapselected_geneids", input$coE_heatmapselected_geneids),
      c("coE_heatmapselected_cells", coE_selctedCluster()$selectedCells()),
      c("coE_heatmapselected_sampcolPal", sampleCols$colPal),
      c("coE_heatmapselected_cluscolPal", clusterCols$colPal)
    )
  )
  
  updateButtonColor(buttonName = "updateHeatMapSelectedParameters", parameters = c(
    "coE_heatmapselected_geneids", "coE_heatmapselected_cells",
    "coE_heatmapselected_sampcolPal", "coE_heatmapselected_cluscolPal"
  ))
})

# observe: updatetopCCGenesSelectedParameters ----
observe(label = "ob_updatetopCCGenesSelectedParams", {
  if (DEBUG) cat(file = stderr(), "observe ob_updatetopCCGenesSelectedParams\n")
  input$updatetopCCGenesSelectedParameters
  setRedGreenButtonCurrent(
    vars = list(
      c("coE_heatmapselected_geneids", input$coE_heatmapselected_geneids),
      c("coE_heatmapselected_cells", coE_selctedCluster()$selectedCells())
    )
  )
  
  updateButtonColor(buttonName = "updatetopCCGenesSelectedParameters", parameters = c(
    "coE_heatmapselected_geneids", "coE_heatmapselected_cells"
    ))
})

# observe: updateMinExprSelectedParameters ----
observe(label = "ob_MinExprParams", {
  if (DEBUG) cat(file = stderr(), "observe MinExprVars\n")
  input$updateMinExprSelectedParameters
  setRedGreenButtonCurrent(
    vars = list(
      c("coEtgPerc",(input$coEtgPerc)),
      c("coEtgMinExpr",(input$coEtgMinExpr)),
      c("coE_heatmapselected_cells", coE_selctedCluster()$selectedCells())
    )
  )
  
  updateButtonColor(buttonName = "updateMinExprSelectedParameters", parameters = c(
    "coEtgPerc", "coEtgMinExpr", "coE_heatmapselected_cells"
  ))
})


# EXPLORE TAB VIOLIN PLOT ------------------------------------------------------------------
# output$coE_geneGrp_vio_plot <- plotly::renderPlotly({
output$coE_geneGrp_vio_plot <- renderPlot({
  if (DEBUG) cat(file = stderr(), "coE_geneGrp_vio_plot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_geneGrp_vio_plot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_geneGrp_vio_plot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_geneGrp_vio_plot", id = "coE_geneGrp_vio_plot", duration = NULL)
  }

  projections <- projections()
  scEx_log <- scEx_log()
  geneListStr <- input$coE_geneGrpVioIds
  projectionVar <- input$coE_dimension_xVioiGrp
  minExpr <- input$coEminExpr
  coE_showPermutations <- input$coE_showPermutations
  # colPal = coE_geneGrp_vioFunc # TODO must be wrong
  sampCol <- sampleCols$colPal
  ccols <- clusterCols$colPal

  # upI <- coE_updateInputXviolinPlot() # no need to check because this is done in projections
  if (is.null(projections) | is.null(scEx_log)) {
    if (DEBUG) cat(file = stderr(), "output$coE_geneGrp_vio_plot:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_geneGrp_vio_plot.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/coE_geneGrp_vio_plot.RData")

  featureData <- rowData(scEx_log)
  retVal <- coE_geneGrp_vioFunc(
    genesin = geneListStr,
    projections = projections,
    scEx = scEx_log,
    featureData = featureData,
    minExpr = minExpr,
    dbCluster = projectionVar,
    coE_showPermutations = coE_showPermutations,
    sampCol = sampCol,
    ccols = ccols
  )
  .schnappsEnv[["coE_geneGrp_vio_plot"]] <- retVal
   exportTestValues(coE_geneGrp_vio_plot = {
    retVal
  })
  return(retVal)
})

# EXPLORE GROUPEDB VIOLIN PLOT ------------------------------------------------------------------
output$coE_geneGrp_vio_plot2 <- plotly::renderPlotly({
# output$coE_geneGrp_vio_plot2 <- renderPlot({
  if (DEBUG) cat(file = stderr(), "coE_geneGrp_vio_plot2 started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_geneGrp_vio_plot2")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_geneGrp_vio_plot2")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_geneGrp_vio_plot2", id = "coE_geneGrp_vio_plot2", duration = NULL)
  }
  
  projections <- projections()
  scEx_log <- scEx_log()
  geneListStr <- input$coE_geneGrpVioIds2
  projectionVar <- input$coE_dimension_xVioiGrp2
  minExpr <- input$coEminExpr2
  # colPal = coE_geneGrp_vioFunc # TODO must be wrong
  sampCol <- sampleCols$colPal
  ccols <- clusterCols$colPal
  
  # upI <- coE_updateInputXviolinPlot() # no need to check because this is done in projections
  if (is.null(projections) | is.null(scEx_log)) {
    if (DEBUG) cat(file = stderr(), "output$coE_geneGrp_vio_plot2:NULL\n")
    return(NULL)
  }
  if (!projectionVar %in% colnames(projections)) {
    if (DEBUG) cat(file = stderr(), "coE_geneGrp_vio_plot2: projectionVar not known: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_geneGrp_vio_plot2.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/coE_geneGrp_vio_plot2.RData")
  
  featureData <- rowData(scEx_log)
  retVal <- coE_geneGrp_vioFunc2(
    genesin = geneListStr,
    projections = projections,
    scEx = scEx_log,
    featureData = featureData,
    minExpr = minExpr,
    dbCluster = projectionVar,
    sampCol = sampCol,
    ccols = ccols
  )
  .schnappsEnv[["coE_geneGrp_vio_plot"]] <- retVal
  exportTestValues(coE_geneGrp_vio_plot = {
    retVal
  })
  return(retVal)
})



require(ggplot2)
require(ggalluvial)
# Alluvial plot of two factors
output$alluvial_plot <- renderPlot({
  if (DEBUG) {
    cat(file = stderr(), "alluvial_plot started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "alluvial_plot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "alluvial_plot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("alluvial_plot", id = "alluvial_plot", duration = NULL)
  }
  
  # load reactive data
  if (DEBUG) {
    cat(file = stderr(), "alluvial_plot started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "alluvial_plot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "alluvial_plot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("alluvial_plot", id = "alluvial_plot", duration = NULL)
  }
  
  projections <- projections()
  alluiv1 <- input$alluiv1
  alluiv2 <- input$alluiv2
  
  # return if nothing to be computed
  if (is.null(projections) ) {
    return(NULL)
  }
  if (alluiv1 == alluiv2 ) {
    return(NULL)
  }
  if(!alluiv1 %in% colnames(projections)) {
    return(NULL)
  }
  if(!alluiv2 %in% colnames(projections)) {
    return(NULL)
  }
  # some debugging messages
  if (DEBUG) cat(file = stderr(), paste("alluvial_plot:\n"))
  # for development and debugging purposes
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/alluvial_plot.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/alluvial_plot.RData")
  
  dat = projections[,c(alluiv1, alluiv2)]
  # dat$cells = rownames(projections)
  gg = ggplot(as.data.frame(dat),
              aes_string(  axis1 = alluiv1, axis2 = alluiv2)) +
    geom_alluvium(aes_string(fill = alluiv1), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", infer.label = TRUE) +
    scale_x_discrete(limits = c(alluiv1, alluiv2), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    ggtitle(paste("Alluvial plot of ", alluiv1, "and", alluiv2))
  
 
  # create and return the plot
  # ggalluvial::(dummyNRow)
  return(gg)
})

# observe alluiv ----
observe({
  projections <- projections()
  
  req(projections)
  # save(file = "~/SCHNAPPsDebug/alluvial_plot2.RData", list = c(ls(), ls(envir = globalenv())))
  # load(file="~/SCHNAPPsDebug/alluvial_plot2.RData")
  facs = which(lapply(projections, class) == "factor")
  idx=1
  for (ff in facs){
    if(length(levels(projections[,ff]))>100){
      facs = facs[-idx]
    }
    idx = idx + 1
  }
  
  updateSelectInput(session, "alluiv1",
                    choices = c(colnames(projections)[facs]),
                    selected = .schnappsEnv$alluiv1
  )
  updateSelectInput(session, "alluiv2",
                    choices = c(colnames(projections)[facs]),
                    selected = .schnappsEnv$alluiv2
  )
})


