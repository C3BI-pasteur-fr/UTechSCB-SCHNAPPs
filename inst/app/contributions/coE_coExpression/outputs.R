
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
  coE_topExpGenesTable, caption = "Tables with highest expressed genes"
)

# max expressed genes ----
callModule(
  tableSelectionServer,
  "coE_topCCGenes",
  coE_topExpCCTable, caption = "Tables with highest expressed correlated genes"
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

observe({
  if (DEBUG) cat(file = stderr(), "coE_geneGrp_vioOBS started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_geneGrp_vioOBS")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_geneGrp_vioOBS")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_geneGrp_vioOBS", id = "coE_geneGrp_vioOBS", duration = NULL)
  }
  scEx_log <- scEx_log()
  # genesin <- input$coE_geneGrpVioIds
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_geneGrp_vioOBS.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/coE_geneGrp_vioOBS.RData")
  if(is.null(scEx_log)) return(NULL)

  # genesin <- toupper(genesin)
  # genesin <- gsub(" ", "", genesin, fixed = TRUE)
  # genesin <- strsplit(genesin, ",")[[1]]
  featureData <- rowData(scEx_log)
  
  # map <- rownames(featureData[which(toupper(featureData$symbol) %in% genesin), ])
  # minExp <- min(assays(scEx_log)[[1]][map, , drop = F])
  minExp <- min(assays(scEx_log)[[1]])
  maxExp <- max(assays(scEx_log)[[1]])
  step = (maxExp - minExp) / 1000000
  updateSliderInput(session,
                    inputId = "coEminMaxExpr",
                    step = step,
                    min = signif(minExp + step, digits = 4),
                    max = signif(maxExp, digits = 4)
  )
  updateSliderInput(session,
                    inputId = "coEminMaxExpr2",
                    step = (maxExp - minExp) / 1000000,
                    min = signif(minExp + step, digits = 4),
                    max = signif(maxExp, digits = 4)
  )
  
  
})



coeMinMax = reactive({
  if (is.null(input$coEminMaxExpr))
    list(x = NA, y = NA)
  else
    input$coEminMaxExpr
}) %>% debounce(1000)

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
  showExpression <- input$coE_showExpression
  minMaxExpr <- coeMinMax()
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
    minMaxExpr = minMaxExpr,
    dbCluster = projectionVar,
    coE_showPermutations = coE_showPermutations,
    sampCol = sampCol,
    ccols = ccols,
    showExpression = showExpression
  )
  
  if(is.null(retVal)) return(NULL)
  
  # serializePlots <- function(data) {
  #   lapply(data, function(plot) {
  #     rawToChar(serialize(plot, NULL, ascii = TRUE))
  #   })
  # }
  
  # removed because it was too long and took too much memory
  # .schnappsEnv[["coE_geneGrp_vio_plot"]] <- serializePlots(retVal)
  af = coE_geneGrp_vioFunc
  # remove env because it is too big
  specEnv = emptyenv()
  eval(envir = specEnv, {
    finner = finner;
    combinePermutations = combinePermutations
  })
  environment(af) = new.env(parent = specEnv)
  .schnappsEnv[["coE_geneGrp_vio_plot"]] <- list(plotFunc = af,
                                                 genesin = geneListStr,
                                                 projections = projections,
                                                 scEx = scEx_log,
                                                 featureData = featureData,
                                                 minMaxExpr = minMaxExpr,
                                                 dbCluster = projectionVar,
                                                 coE_showPermutations = coE_showPermutations,
                                                 sampCol = sampCol,
                                                 ccols = ccols,
                                                 showExpression = showExpression
  )
  # .schnappsEnv[["coE_geneGrp_vio_plot"]] <- retVal
  #  exportTestValues(coE_geneGrp_vio_plot = {
  #    serializePlots(retVal)
  # })
  return(retVal)
})

# EXPLORE GROUPEDB VIOLIN PLOT ------------------------------------------------------------------
coeMinMax2 = reactive({
  if (is.null(input$coEminMaxExpr2))
    list(x = NA, y = NA)
  else
    input$coEminMaxExpr2 
})%>% debounce(1000)

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
  minMaxExpr <- coeMinMax2()
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
    minMaxExpr = minMaxExpr,
    dbCluster = projectionVar,
    sampCol = sampCol,
    ccols = ccols
  )
  
  if(is.null(retVal)) return(NULL)
  
  af = coE_geneGrp_vioFunc2
  # remove env because it is too big
  environment(af) = new.env(parent = emptyenv())
  .schnappsEnv[["coE_geneGrp_vio_plot"]] <- list(plotFunc = af,
                                                 genesin = geneListStr,
                                                 projections = projections,
                                                 scEx = scEx_log,
                                                 featureData = featureData,
                                                 minMaxExpr = minMaxExpr,
                                                 dbCluster = projectionVar,
                                                 sampCol = sampCol,
                                                 ccols = ccols)
  
  # .schnappsEnv[["coE_geneGrp_vio_plot"]] <- retVal
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
  # for development and debugging purposes
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/alluvial_plot.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/alluvial_plot.RData")
  
  dat = projections[,c(alluiv1, alluiv2)]
  # dat$cells = rownames(projections)
  gg = alluvialPlotFunc(dat, alluiv1, alluiv2) 
  af = alluvialPlotFunc
  environment(af) = new.env(parent = emptyenv())
  # for saving to history
  .schnappsEnv[["coE_alluvialPlot"]] <- list(plotFunc = af,
                                             dat = dat, 
                                             alluiv1=alluiv1, 
                                             alluiv2=alluiv2)
  
  return(gg)
})

# observe alluiv ----
observe({
  # projections <- projections()
  # 
  # req(projections)
  projF = projFactors()
  # save(file = "~/SCHNAPPsDebug/alluvial_plot2.RData", list = c(ls(), ls(envir = globalenv())))
  # load(file="~/SCHNAPPsDebug/alluvial_plot2.RData")
  # facs = which(lapply(projections, class) == "factor")
  # idx=1
  # for (ff in facs){
  #   if(length(levels(projections[,ff]))>100){
  #     facs = facs[-idx]
  #   }
  #   idx = idx + 1
  # }
  
  updateSelectInput(session, "alluiv1",
                    choices = projF,
                    selected = .schnappsEnv$alluiv1
  )
  updateSelectInput(session, "alluiv2",
                    choices = projF,
                    selected = .schnappsEnv$alluiv2
  )
})


