# coExpression outputs.R

myZippedReportFiles <- c("output_coE_topExpGenes.csv")


# All clusters heatmap ------
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

callModule(
  tableSelectionServer,
  "coE_scranFindMarkerTable",
  coE_scranFindMarkerTableReact, caption = "Tables with marker genes per cluster"
  
)

# observe: updateHeatMapSelectedParameters ----
observe(label = "ob_HeatMapSelectedParams", {
  if (DEBUG) cat(file = stderr(), "observe HeatMapSelected\n")
  input$updateHeatMapSelectedParameters
  setRedGreenButtonCurrent(
    vars = list(
      c("coE_heatmapselected_geneids", input$coE_heatmapselected_geneids),
      c("coE_heatmapselected_cells", coE_selctedCluster()$selectedCells()),
      c("coE_heatmapselected_sampcolPal", projectionColors$sampleNames),
      c("coE_heatmapselected_cluscolPal", projectionColors$dbCluster)
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

observe(label = "ob16", {
  if (DEBUG) cat(file = stderr(), paste0("observe: coE_dimension_xVioiGrp\n"))
  .schnappsEnv$coE_dimension_xVioiGrp <- input$coE_dimension_xVioiGrp
  .schnappsEnv$defaultValues$coE_dimension_xVioiGrp <- input$coE_dimension_xVioiGrp
  .schnappsEnv$defaultValues$coE_scranFactor <- input$coE_scranFactor
})

observe(label = "ob16a", {
  if (DEBUG) cat(file = stderr(), paste0("observe: coE_dimension_xVioiGrp2\n"))
  .schnappsEnv$coE_dimension_xVioiGrp2 <- input$coE_dimension_xVioiGrp2
  .schnappsEnv$defaultValues$coE_dimension_xVioiGrp2 <- input$coE_dimension_xVioiGrp2
})
observe(label = "ob16b", {
  if (DEBUG) cat(file = stderr(), paste0("observe: alluiv1\n"))
  .schnappsEnv$alluiv1 <- input$alluiv1
  .schnappsEnv$defaultValues$alluiv1 <- input$alluiv1
})
observe(label = "ob16c", {
  if (DEBUG) cat(file = stderr(), paste0("observe: alluiv2\n"))
  .schnappsEnv$alluiv2 <- input$alluiv2
  .schnappsEnv$defaultValues$alluiv2 <- input$alluiv2
})

observe(label ="obs_coE_geneGrpVioIds2", x = {
  .schnappsEnv$defaultValues[["coE_geneGrpVioIds2"]] = input$coE_geneGrpVioIds2
})
observe(label ="obs_coEminMaxExprValue", x = {
  .schnappsEnv$defaultValues[["coEminMaxExprValue"]] = input$coEminMaxExprValue
})
observe(label ="obs_coE_geneGrpVioIds", x = {
  .schnappsEnv$defaultValues[["coE_geneGrpVioIds"]] = input$coE_geneGrpVioIds
})
observe(label ="obs_coE_scale", x = {
  .schnappsEnv$defaultValues[["coE_scale"]] = input$coE_scale
})
observe(label ="obs_coE_showExpression", x = {
  .schnappsEnv$defaultValues[["coE_showExpression"]] = input$coE_showExpression
})
observe(label ="obs_coE_showPermutations", x = {
  .schnappsEnv$defaultValues[["coE_showPermutations"]] = input$coE_showPermutations
})
observe(label ="obs_coEtgPerc", x = {
  .schnappsEnv$defaultValues[["coEtgPerc"]] = input$coEtgPerc
})
observe(label ="obs_coEtgMinExpr", x = {
  .schnappsEnv$defaultValues[["coEtgMinExpr"]] = input$coEtgMinExpr
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
                    min = signif(minExp, digits = 4),
                    max = signif(maxExp, digits = 4),
                    value = c(minExp + step, maxExp)
  )
  updateSliderInput(session,
                    inputId = "coEminMaxExpr2",
                    step = (maxExp - minExp) / 1000000,
                    min = signif(minExp , digits = 4),
                    max = signif(maxExp, digits = 4), 
                    value = c(minExp + step, maxExp)
  )
  
  
})



coeMinMax = reactive({
  if (is.null(input$coEminMaxExpr))
    list(x = NA, y = NA)
  else
    input$coEminMaxExpr
}) %>% debounce(1000)


output$coE_objSize <- renderText({
  if (DEBUG) cat(file = stderr(), "coE_objSize started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_objSize")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_objSize")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_objSize", id = "coE_objSize", duration = NULL)
  }
  # deepDebug()
  scEx_log <- scEx_log()
  projections <- projections()
  direction <- isolate(input$coE_direction)
  prjFact <- input$coE_scranFactor
  
  lfc <- isolate(input$coE_lfc)
  if (is.null(scEx_log)) {
    return(NULL)
  }
  if(!prjFact %in% colnames(projections)) return("please select factor")
  objSize = pryr::object_size(scEx_log, 
                        projections[,prjFact],
                        direction = direction,
                        lfc = lfc)
  
  paste("size of input object for one factor: ", format(objSize))
})

# coE_dotPlot_GeneSets ----
output$coE_dotPlot_GeneSets <- renderPlotly({
  # output$coE_dotPlot_GeneSets <- renderPlotly({
    if (DEBUG) cat(file = stderr(), "coE_dotPlot_GeneSets started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_dotPlot_GeneSets")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_dotPlot_GeneSets")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_dotPlot_GeneSets", id = "coE_dotPlot_GeneSets", duration = NULL)
  }
  projections <- projections()
  scEx_log <- scEx_log()
  clusters = input$coE_dimension_ydotPlotClusters
  geneSets = input$coE_dotPlot_geneSets
  col = input$coE_dotPlot_col
  col.min = input$coE_dotPlot_col.min
  col.max = input$coE_dotPlot_col.max
  dot.min = input$coE_dotPlot_dot.min
  dot.scale = input$coE_dotPlot_dot.scale
  scale.by = input$coE_dotPlot_scale.by
  
  gmtData = gmtData()
  
  if (is.null(projections) | is.null(scEx_log) | is.null(gmtData) | isEmpty(gmtData) | isEmpty(geneSets) | !all(geneSets %in% names(gmtData))) {
    if (DEBUG) cat(file = stderr(), "output$coE_dotPlot_GeneSets:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_dotPlot_GeneSets.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/coE_dotPlot_GeneSets.RData")
  featureData <- rowData(scEx_log)
  retVal <- coE_dotPlot_GeneSets(
    projections = projections,
    scEx_log <- scEx_log,
    clusters = clusters,
    geneSets = geneSets,
    gmtData = gmtData,
    col = col,
    col.min = col.min,
    col.max = col.max,
    dot.min = dot.min,
    dot.scale = dot.scale,
    scale.by = scale.by
  )
  

  if(is.null(retVal)) return(NULL)
  if(length(levels(projections[,clusters]))<3){
    retVal + scale_color_gradient(low = "#808080", high = "white", limits=c(-2, 2))
  }
  af = coE_dotPlot_GeneSets
  # remove env because it is too big
  specEnv = emptyenv()
  environment(af) = new.env(parent = specEnv)
  .schnappsEnv[["coE_dotPlot_GeneSets"]] <- list(plotFunc = af,
                                                 projections = projections,
                                                 scEx_log = scEx_log,
                                                 clusters = clusters,
                                                 geneSets = geneSets,
                                                 gmtData = gmtData,
                                                 col = col,
                                                 col.min = col.min,
                                                 col.max = col.max,
                                                 dot.min = dot.min,
                                                 dot.scale = dot.scale,
                                                 scale.by = scale.by
                                                 
  )
  return(retVal %>% ggplotly())
})

output$scranFindMarkersSelected <- renderText({
  scEx_log = scEx_log()
  req(scEx_log)
  featureData <- rowData(scEx_log)
  wmarkers = scranFindMarkerFullReactiveTable()
  nFindCluster <- isolate(input$coE_nFindMarker)
  
  req(wmarkers)
  
  markerlist = lapply(wmarkers,FUN = function(x){
    rownames(x)[order(x$p.value)[1:nFindCluster]]}) %>% unlist %>% unique
  return(paste(featureData[markerlist, "symbol"], collapse  = ", "))
})

# coE_dotPlot_GeneSetsModuleScore ----
output$coE_dotPlot_GeneSetsModuleScore <- renderPlotly({
  # output$coE_dotPlot_GeneSets <- renderPlotly({
  if (DEBUG) cat(file = stderr(), "coE_dotPlot_GeneSetsModuleScore started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_dotPlotModuleScore_GeneSets")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_dotPlotModuleScore_GeneSets")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_dotPlotModuleScore_GeneSets", id = "coE_dotPlotModuleScore_GeneSets", duration = NULL)
  }
  projections <- projections()
  scEx_log <- scEx_log()
  clusters = input$coE_dimension_ydotPlotModuleScoreClusters
  geneSets = input$coE_dotPlotModuleScore_geneSets
  col = input$coE_dotPlotModuleScore_col
  col.min = input$coE_dotPlotModuleScore_col.min
  col.max = input$coE_dotPlotModuleScore_col.max
  dot.min = input$coE_dotPlotModuleScore_dot.min
  dot.scale = input$coE_dotPlotModuleScore_dot.scale
  scale.by = input$coE_dotPlotModuleScore_scale.by
  
  gmtData = gmtData()
  
  if (is.null(projections) | is.null(scEx_log) | is.null(gmtData) | isEmpty(gmtData) | isEmpty(geneSets) | !all(geneSets %in% names(gmtData))) {
    if (DEBUG) cat(file = stderr(), "output$coE_dotPlotModuleScore_GeneSets:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_dotPlotModuleScore_GeneSets.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/coE_dotPlotModuleScore_GeneSets.RData")
  featureData <- rowData(scEx_log)
  retVal <- coE_dotPlot_GeneSetsModuleScore(
    projections = projections,
    scEx_log <- scEx_log,
    clusters = clusters,
    geneSets = geneSets,
    gmtData = gmtData,
    col = col,
    col.min = col.min,
    col.max = col.max,
    dot.min = dot.min,
    dot.scale = dot.scale,
    scale.by = scale.by
  )
  
  
  if(is.null(retVal)) return(NULL)
  af = coE_dotPlot_GeneSetsModuleScore
  # remove env because it is too big
  specEnv = emptyenv()
  environment(af) = new.env(parent = specEnv)
  .schnappsEnv[["coE_dotPlotModuleScore_GeneSets"]] <- list(plotFunc = af,
                                                 projections = projections,
                                                 scEx_log = scEx_log,
                                                 clusters = clusters,
                                                 geneSets = geneSets,
                                                 gmtData = gmtData,
                                                 col = col,
                                                 col.min = col.min,
                                                 col.max = col.max,
                                                 dot.min = dot.min,
                                                 dot.scale = dot.scale,
                                                 scale.by = scale.by
                                                 
  )
  return(retVal %>% ggplotly())
})

# save to history dotplot ---d-
observe(label = "save2histDotPlot", {
  clicked  = input$save2histDotPlot
  if (DEBUG) cat(file = stderr(), "observe save2histDotPlot \n")
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
  add2history(type = "save", input = isolate( reactiveValuesToList(input)), 
              comment = paste("# DotPlot genes \n",
                              "fun = plotData$plotData$plotFunc\n", 
                              "environment(fun) = environment()\n",
                              "plotData$plotData$outfile=NULL\n",
                              "print(do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)]))\n"
              ),
              plotData = .schnappsEnv[["coE_dotPlot_GeneSets"]])
  
})

# save to history dotplot ---d-
observe(label = "save2histDotPlotModuleScore", {
  clicked  = input$save2histDotPlotModuleScore
  if (DEBUG) cat(file = stderr(), "observe save2histDotPlotModuleScore \n")
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
  add2history(type = "save", input = isolate( reactiveValuesToList(input)), 
              comment = paste("# DotPlotModuleScore\n",
                              "fun = plotData$plotData$plotFunc\n", 
                              "environment(fun) = environment()\n",
                              "plotData$plotData$outfile=NULL\n",
                              "print(do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)]))\n"
              ),
              plotData = .schnappsEnv[["coE_dotPlotModuleScore_GeneSets"]])
  
})

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
  coE_scale <- input$coE_scale
  minMaxExpr <- coeMinMax()
  coE_showPermutations <- input$coE_showPermutations
  # colPal = coE_geneGrp_vioFunc # TODO must be wrong
  # sampCol <- projectionColors$sampleNames
  # ccols <- projectionColors$dbCluster
  pc = projectionColors %>% reactiveValuesToList()
  
  # upI <- coE_updateInputXviolinPlot() # no need to check because this is done in projections
  if (is.null(projections) | is.null(scEx_log)) {
    if (DEBUG) cat(file = stderr(), "output$coE_geneGrp_vio_plot:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_geneGrp_vio_plot.RData", list = c(ls()))
  }
  # cp=load(file="~/SCHNAPPsDebug/coE_geneGrp_vio_plot.RData")
  
  featureData <- rowData(scEx_log)
  retVal <- coE_geneGrp_vioFunc(
    genesin = geneListStr,
    projections = projections,
    scEx = scEx_log,
    featureData = featureData,
    minMaxExpr = minMaxExpr,
    dbCluster = projectionVar,
    coE_showPermutations = coE_showPermutations,
    projectionColors = pc,
    showExpression = showExpression,
    coE_scale = coE_scale
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
                                                 projectionColors=pc,
                                                 showExpression = showExpression,
                                                 coE_scale = coE_scale
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
}) %>% debounce(1000)

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
  # sampCol <- projectionColors$sampleNames
  # ccols <- projectionColors$dbCluster
  pc = projectionColors %>% reactiveValuesToList()
  
  # upI <- coE_updateInputXviolinPlot() # no need to check because this is done in projections
  if (is.null(projections) | is.null(scEx_log)) {
    if (DEBUG) cat(file = stderr(), "output$coE_geneGrp_vio_plot2:NULL\n")
    return(NULL)
  }
  if (!all(projectionVar %in% colnames(projections))) {
    if (DEBUG) cat(file = stderr(), "coE_geneGrp_vio_plot2: projectionVar not known: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_geneGrp_vio_plot2.RData", list = c(ls()))
  }
  #cp = load(file="~/SCHNAPPsDebug/coE_geneGrp_vio_plot2.RData")

  featureData <- rowData(scEx_log)
  retVal <- coE_geneGrp_vioFunc2(
    genesin = geneListStr,
    projections = projections,
    scEx = scEx_log,
    featureData = featureData,
    minMaxExpr = minMaxExpr,
    dbCluster = projectionVar,
    projectionColors = pc
    
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
                                                 projectionColors = pc
  )
  
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
  # cp =load(file="~/SCHNAPPsDebug/alluvial_plot.RData")
  
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

observe({
  if (DEBUG) cat(file = stderr(), "observe gmtData: coE_dotPlot_geneSets\n")
  
  gmtList = gmtData()
  
  updateSelectizeInput(session, "coE_dotPlot_geneSets",
                    choices = names(gmtList),
                    selected = .schnappsEnv$coE_dotPlot_geneSets, server = TRUE
  )
  updateSelectizeInput(session, "coE_dotPlotModuleScore_geneSets",
                    choices = names(gmtList),
                    selected = .schnappsEnv$coE_dotPlotModuleScore_geneSets, server = TRUE
  )
})


# observe obs_coE_heatmap_geneids ----
observe(label = "obs_coE_heatmap_geneids", x= {
  .schnappsEnv$defaultValues[["coE_heatmap_geneids"]] = input$coE_heatmap_geneids
  .schnappsEnv$defaultValues[["coE_subSampleFactor"]] = input$coE_subSampleFactor
  .schnappsEnv$defaultValues[["coE_nSubsample"]] = input$coE_nSubsample
  .schnappsEnv[["coE_heatmap_geneids"]] = input$coE_heatmap_geneids
  .schnappsEnv[["coE_subSampleFactor"]] = input$coE_subSampleFactor
  .schnappsEnv[["coE_nSubsample"]] = input$coE_nSubsample
  
})

# observe alluiv ----
observe({
  projF = projFactors()
  # browser()
  updateSelectInput(session, "alluiv1",
                    choices = projF,
                    selected = .schnappsEnv$alluiv1
  )
  updateSelectInput(session, "coE_subSampleFactor",
                    choices = projF
                    ,
                    selected = .schnappsEnv$coE_subSampleFactor
  )
  updateSelectInput(session, "alluiv2",
                    choices = projF,
                    selected = .schnappsEnv$alluiv2
  )
  updateSelectInput(session, "coE_dimension_ydotPlotClusters",
                    choices = projF,
                    selected = .schnappsEnv$coE_dimension_ydotPlotClusters)
  updateSelectInput(session, "coE_dimension_ydotPlotModuleScoreClusters",
                    choices = projF,
                    selected = .schnappsEnv$coE_dimension_ydotPlotModuleScoreClusters)
})
