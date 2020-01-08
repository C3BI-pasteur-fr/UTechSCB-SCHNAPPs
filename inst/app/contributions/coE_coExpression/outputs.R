
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

coE_SOM_dataInput <- callModule(
  cellSelectionModule,
  "coE_SOM_dataInput"
)


# SOM heatmap module -----
callModule(
  pHeatMapModule,
  "coE_heatmapSOM",
  coE_heatmapSOMReactive
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

  upI <- coE_updateInputXviolinPlot() # no need to check because this is done in projections
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

# # observer for coE_clusterValSOM ----
# observeEvent(
#   label = "ob15",
#   eventExpr = input$coE_clusterSOM,
#   handlerExpr = {
#     projections <- projections()
#     if (DEBUG) cat(file = stderr(), "observeEvent: input$coE_clusterSOM.\n")
#     # Can use character(0) to remove all choices
#     if (is.null(projections)) {
#       return(NULL)
#     }
#     if (!input$coE_clusterSOM %in% colnames(projections)) {
#       return(NULL)
#     }
#     choicesVal <- levels(projections[, input$coE_clusterSOM])
#     updateSelectInput(
#       session,
#       "coE_clusterValSOM",
#       choices = choicesVal,
#       selected = .schnappsEnv$coE_SOMSelection
#     )
#   }
# )

# observer of button Color SOM ----
observe(label = "ob_somParameter", 
        {
          if (DEBUG) cat(file = stderr(), "ob_somParameter\n")
          # browser()
          input$updateSOMParameters
          setRedGreenButtonCurrent(
            vars = list(
              c("coE_geneSOM", input$coE_geneSOM),
              c("coE_dimSOM", input$coE_dimSOM),
              c("coE_SOM_dataInput-Mod_PPGrp", input$'coE_SOM_dataInput-Mod_PPGrp'),
              c("coE_SOM_dataInput-Mod_clusterPP", input$'coE_SOM_dataInput-Mod_clusterPP')
            )
          )
          updateButtonColor(buttonName = "updateSOMParameters", parameters = c(
            "coE_geneSOM", "coE_dimSOM",
            "coE_SOM_dataInput-Mod_PPGrp", "coE_SOM_dataInput-Mod_clusterPP"
          ))
          
        })


