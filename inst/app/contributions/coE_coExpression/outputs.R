
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

# SOM heatmap module -----
callModule(
  pHeatMapModule,
  "coE_heatmapSOM",
  coE_heatmapSOMReactive
)

# EXPLORE TAB VIOLIN PLOT ------------------------------------------------------------------
# output$coE_geneGrp_vio_plot <- plotly::renderPlotly({
output$coE_geneGrp_vio_plot <- renderPlot({
  if (DEBUG) cat(file = stderr(), "coE_geneGrp_vio_plot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_geneGrp_vio_plot")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "coE_geneGrp_vio_plot")
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

  upI <- coE_updateInputXviolinPlot() # no need to check because this is done in projections
  if (is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$coE_geneGrp_vio_plot:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_geneGrp_vio_plot.RData", list = c(ls(), ls(envir = globalenv()), ls(.schnappsEnv)))
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

  exportTestValues(coE_geneGrp_vio_plot = {
    retVal
  })
  return(retVal)
})
