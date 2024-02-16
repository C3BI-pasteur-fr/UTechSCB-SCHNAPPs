require(parallel)


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

observe(label ="obs_DE_logNormalization_sf", x = {
  .schnappsEnv$defaultValues[["DE_logNormalization_sf"]] = input$DE_logNormalization_sf
})
observe(label ="obs_DE_seuratRefBased_splitby", x = {
  .schnappsEnv$defaultValues[["DE_seuratRefBased_splitby"]] = input$DE_seuratRefBased_splitby
})
observe(label ="obs_DE_seuratSCtransform_split.by", x = {
  .schnappsEnv$defaultValues[["DE_seuratSCtransform_split.by"]] = input$DE_seuratSCtransform_split.by
})
observe(label ="obs_DE_seuratSCtransform_vars2regress", x = {
  .schnappsEnv$defaultValues[["DE_seuratSCtransform_vars2regress"]] = input$DE_seuratSCtransform_vars2regress
})
observe(label ="obs_DE_seuratStandard_splitby", x = {
  .schnappsEnv$defaultValues[["DE_seuratStandard_splitby"]] = input$DE_seuratStandard_splitby
})
observe(label ="obs_DE_seuratLogNorm_var2reg", x = {
  .schnappsEnv$defaultValues[["DE_seuratLogNorm_var2reg"]] = input$DE_seuratLogNorm_var2reg
})
observe(label ="obs_DE_logNormalization_sf", x = {
  .schnappsEnv$defaultValues[["DE_logNormalization_sf"]] = input$DE_logNormalization_sf
})
observe(label ="obs_DE_panelplotids", x = {
  .schnappsEnv$defaultValues[["DE_panelplotids"]] = input$DE_panelplotids
})
observe(label ="obs_DE_seuratSCTnorm_var2reg", x = {
  .schnappsEnv$defaultValues[["DE_seuratSCTnorm_var2reg"]] = input$DE_seuratSCTnorm_var2reg
})
observe(label ="obs_DE_nCol", x = {
  .schnappsEnv$defaultValues[["DE_nCol"]] = input$DE_nCol
})
observe(label ="obs_DE_panelplotSameScale", x = {
  .schnappsEnv$defaultValues[["DE_panelplotSameScale"]] = input$DE_panelplotSameScale
})
observe(label ="obs_DE_expclusters_col", x = {
  .schnappsEnv$defaultValues[["DE_expclusters_col"]] = input$DE_expclusters_col
})
observe(label ="obs_DE_gene_id", x = {
  .schnappsEnv$defaultValues[["DE_gene_id"]] = input$DE_gene_id
})

observe(label ="obs_DE_pFact_dim_x", x = {
  .schnappsEnv$obs_DE_pFact_dim_x = input$DE_pFact_dim_x
  .schnappsEnv$defaultValues[["DE_pFact_dim_x"]] = input$DE_pFact_dim_x
})
observe(label ="obs_DE_pFact_dim_y", x = {
  .schnappsEnv$DE_pFact_dim_y = input$DE_pFact_dim_y
  .schnappsEnv$defaultValues[["DE_pFact_dim_y"]] = input$DE_pFact_dim_y
})
observe(label ="obs_DE_panelplotFactSameScale", x = {
  .schnappsEnv$DE_panelplotFactSameScale = input$DE_panelplotFactSameScale
  .schnappsEnv$defaultValues[["DE_panelplotFactSameScale"]] = input$DE_panelplotFactSameScale
})
observe(label ="obs_DE_panelplotFactPvalue", x = {
  .schnappsEnv$DE_panelplotFactPvalue = input$DE_panelplotFactPvalue
  .schnappsEnv$defaultValues[["DE_panelplotFactPvalue"]] = input$DE_panelplotFactPvalue
})
observe(label ="obs_DE_pFactnCol", x = {
  .schnappsEnv$DE_pFactnCol = input$DE_pFactnCol
  .schnappsEnv$defaultValues[["DE_pFactnCol"]] = input$DE_pFactnCol
})
observe(label ="obs_DE_pFactIds", x = {
  .schnappsEnv$DE_pFactIds = input$DE_pFactIds
  .schnappsEnv$defaultValues[["DE_pFactIds"]] = input$DE_pFactIds
})
observe(label = "ob19", {
  if (DEBUG) cat(file = stderr(), "observe: DE_clusterSelectionPanelPlot\n")
  .schnappsEnv$DE_cl1 <- input$DE_clusterSelectionPanelPlot
})


observe(label = "DE_seuratLogNorm_var2regOBSinp", {
  if (DEBUG) cat(file = stderr(), paste0("observe: DE_seuratLogNorm_var2regOBSinp\n"))
  .schnappsEnv$DE_seuratLogNorm_var2reg <- input$DE_seuratLogNorm_var2reg
})
observe(label = "DE_seuratSCtransform_vars2regressOBSinp", {
  if (DEBUG) cat(file = stderr(), paste0("observe: DE_seuratSCtransform_vars2regress\n"))
  .schnappsEnv$DE_seuratSCtransform_vars2regress <- input$DE_seuratSCtransform_vars2regress
})
observe(label = "DE_seuratSCtransform_split.byOBSinp", {
  if (DEBUG) cat(file = stderr(), paste0("observe: DE_seuratSCtransform_split.by\n"))
  .schnappsEnv$DE_seuratSCtransform_split.by <- input$DE_seuratSCtransform_split.by
})
observe(label = "DE_seuratStandard_splitbyOBSinp", {
  if (DEBUG) cat(file = stderr(), paste0("observe: DE_seuratStandard_splitby\n"))
  .schnappsEnv$DE_seuratStandard_splitby <- input$DE_seuratStandard_splitby
})




observe(label = "ob17x", {
  if (DEBUG) cat(file = stderr(), "observe: DE_expclusters_x\n")
  .schnappsEnv$DE_expclusters_x <- input$DE_expclusters_x
  .schnappsEnv$defaultValues$DE_expclusters_x <- input$DE_expclusters_x
})
observe(label = "ob17y", {
  if (DEBUG) cat(file = stderr(), "observe: DE_expclusters_y\n")
  .schnappsEnv$DE_expclusters_y <- input$DE_expclusters_y
  .schnappsEnv$defaultValues$DE_expclusters_y <- input$DE_expclusters_y
})
observe(label = "ob17z", {
  if (DEBUG) cat(file = stderr(), "observe: DE_expclusters_z\n")
  .schnappsEnv$DE_expclusters_z <- input$DE_expclusters_z
  .schnappsEnv$defaultValues$DE_expclusters_z <- input$DE_expclusters_z
})
# observe(label = "ob17c", {
#   if (DEBUG) cat(file = stderr(), "observe: DE_expclusters_col\n")
#   .schnappsEnv$DE_expclusters_col <- input$DE_expclusters_col
# })
observe(label = "ob17d", {
  if (DEBUG) cat(file = stderr(), "observe: DE_gene_vio_x\n")
  .schnappsEnv$DE_gene_vio_x <- input$DE_gene_vio_x
  .schnappsEnv$defaultValues$DE_gene_vio_x <- input$DE_gene_vio_x
})
.schnappsEnv$DE_dim_x <- "tsne1"
.schnappsEnv$DE_dim_y <- "tsne1"
observe(label = "ob17", {
  if (DEBUG) cat(file = stderr(), "observe: DE_dim_x\n")
  .schnappsEnv$DE_dim_x <- input$DE_dim_x
  .schnappsEnv$defaultValues$DE_dim_x <- input$DE_dim_x
})
observe(label = "ob18", {
  if (DEBUG) cat(file = stderr(), "observe: DE_dim_y\n")
  .schnappsEnv$DE_dim_y <- input$DE_dim_y
  .schnappsEnv$defaultValues$DE_dim_y <- input$DE_dim_y
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
  projFactors <- projFactors()
  
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  
  # Can also set the label and select items
  updateSelectInput(session, "DE_expclusters_x",
                    choices = colnames(projections),
                    selected = .schnappsEnv$DE_expclusters_x
  )
  updateSelectInput(session, "DE_expclusters_y",
                    choices = colnames(projections),
                    selected = .schnappsEnv$DE_expclusters_y
  )
  updateSelectInput(session, "DE_expclusters_z",
                    choices = colnames(projections),
                    selected = .schnappsEnv$DE_expclusters_z
  )
  # updateSelectInput(session, "DE_expclusters_col",
  #                   choices = colnames(projections),
  #                   selected = .schnappsEnv$DE_expclusters_col
  # )
  updateSelectInput(session, "DE_gene_vio_x",
                    choices = projFactors,
                    selected = .schnappsEnv$DE_gene_vio_x
  )
  
  updateSelectInput(session, "DE_dim_x",
                    choices = colnames(projections),
                    selected = .schnappsEnv$DE_dim_x
  )
  
  # Can also set the label and select items
  updateSelectInput(session, "DE_dim_y",
                    choices = colnames(projections),
                    selected = .schnappsEnv$DE_dim_y
  )
  updateSelectInput(session,"DE_pFact_dim_x",
                    choices = colnames(projections),
                    selected = .schnappsEnv$DE_pFact_dim_x)
  updateSelectInput(session,"DE_pFact_dim_y",
                    choices = colnames(projections),
                    selected = .schnappsEnv$DE_pFact_dim_y)
  updateSelectInput(session, "DE_pFactIds",
                    choices = projFactors,
                    selected = .schnappsEnv$DE_pFactIds)
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


## observe DE_seuratLogNorm_var2regOBS ----
observe(label = "DE_seuratLogNorm_var2regOBS", {
  scEx <- scEx()
  tmp <- input$normalizationRadioButton
  if (DEBUG) cat(file = stderr(), "observe: DE_seuratLogNorm_var2regOBS\n")
  # Can use character(0) to remove all choices
  if (is.null(scEx)) {
    return(NULL)
  }
  choicesVal = names(Filter(is.factor, colData(scEx)))
  cdat =  colData(scEx)
  choicesVal = choicesVal[unlist(lapply(choicesVal, FUN = function(x) {length(levels(cdat[,x]))>1}))]
  choicesVal = c("", choicesVal)
  # save(file = "~/SCHNAPPsDebug/DE_seuratLogNorm_var2regOBS.RData", list = c(ls(), ".schnappsEnv"))
  # cp = load(file="~/SCHNAPPsDebug/DE_seuratLogNorm_var2regOBS.RData")
  # deepDebug()
  updateSelectInput(
    session,
    "DE_seuratLogNorm_var2reg",
    choices = choicesVal
    ,
    selected = .schnappsEnv$DE_seuratLogNorm_var2reg
  )
  updateSelectInput(
    session,
    "DE_seuratSCtransform_vars2regress",
    choices = colData(scEx) %>% names()
    ,
    selected = .schnappsEnv$DE_seuratSCtransform_vars2regress
  )
  updateSelectInput(
    session,
    "DE_seuratSCtransform_split.by",
    choices = choicesVal
    ,
    selected = .schnappsEnv$DE_seuratSCtransform_split.by
  )
  updateSelectInput(
    session,
    "DE_seuratSCTnorm_var2reg",
    choices = choicesVal
    ,
    selected = .schnappsEnv$DE_seuratSCtransform_split.by
  )
  updateSelectInput(
    session,
    "DE_seuratStandard_splitby",
    choices = choicesVal
    ,
    selected = .schnappsEnv$DE_seuratSCtransform_split.by
  )
  updateSelectInput(
    session,
    "DE_seuratRefBased_splitby",
    choices = choicesVal
    ,
    selected = .schnappsEnv$DE_seuratSCtransform_split.by
  )
  # save(file = "~/SCHNAPPsDebug/DE_seuratLogNorm_var2regOBS2.RData", list = c(ls(), ".schnappsEnv"))
  # cp = load(file="~/SCHNAPPsDebug/DE_seuratLogNorm_var2regOBS2.RData")
  
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
  pc <- projectionColors %>% reactiveValuesToList()
  x = input$DE_gene_vio_x
  
  selectedCells <- DE_Exp_dataInput()
  if(is.null(selectedCells)) return(NULL)
  cellNs <- selectedCells$cellNames()
  sampdesc <- selectedCells$selectionDescription()
  prj <- isolate(selectedCells$ProjectionUsed())
  prjVals <- isolate(selectedCells$ProjectionValsUsed())
  
  if (is.null(scEx_log) | is.null(projections) | is.null(cellNs)) {
    if (DEBUG) cat(file = stderr(), "output$DE_gene_vio_plot:NULL\n")
    return(NULL)
  }
  debugControl("DE_gene_vio_plot", list = c(ls()))
  # cp = load(file="~/SCHNAPPsDebug/DE_gene_vio_plot.RData")
  
  
  p1 <- DE_geneViolinFunc(scEx_log = scEx_log[, cellNs], g_id = g_id, projections = projections[cellNs, ], ccols = pc[[x]], x)
  
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
output$DE_clusterSelectionPanelPlot <- renderUI({
  if (DEBUG) cat(file = stderr(), "output$DE_clusterSelectionPanelPlot\n")
  projections <- projections()
  upI <- DE_updateInputExpPanel()
  if (is.null(projections)) {
    HTML("Please load data")
  } else {
    noOfClusters <- levels(as.factor(projections$dbCluster))
    sc_selectInput(
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

dePanelFactCellSelection <- callModule(
  cellSelectionModule,
  "DE_PanelPlotFactCellSelection"
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
      c("DE_PanelPlotCellSelection-Mod_PPGrp", (input$`DE_PanelPlotCellSelection-Mod_PPGrp`)),
      c("DE_panelplotPvalue", (input$DE_panelplotPvalue))
    )
  )
  
  updateButtonColor(buttonName = "updatePanelPlot", parameters = c(
    "DE_panelplotids", "DE_dim_x", "DE_dim_y", "DE_panelplotSameScale", 
    "DE_nCol", "DE_PanelPlotCellSelection-Mod_clusterPP", "DE_PanelPlotCellSelection-Mod_PPGrp", "DE_panelplotPvalue"
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
  applyPvalue <- isolate(input$DE_panelplotPvalue)
  scEx_log <- scEx_log()
  projections <- projections()
  # DE_updateInputPPt()
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
  # debugControl("DE_panelPlot", list = c("scEx_log", "projections", "genesin",
  # "dimx4", "dimy4", "sameScale", "nCol", "sampdesc" , "cellNs"))
  # cp = load(file="~/SCHNAPPsDebug/DE_panelPlot.RData")
  
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")
  genesin <- genesin[[1]]
  
  # if (DEBUG) cat(file = stderr(), paste("output:sampdesc",sampdesc,"\n"))
  retVal <- panelPlotFunc_m(scEx_log, projections, genesin, dimx4, dimy4, sameScale, nCol, sampdesc, cellNs, applyPvalue = applyPvalue) 
  
  setRedGreenButton(
    vars = list(
      c("DE_panelplotids", isolate(input$DE_panelplotids)),
      c("DE_dim_x", isolate(input$DE_dim_x)),
      c("DE_dim_y", isolate(input$DE_dim_y)),
      c("DE_panelplotSameScale", isolate(input$DE_panelplotSameScale)),
      c("DE_nCol", isolate(input$DE_nCol)),
      c("DE_PanelPlotCellSelection-Mod_clusterPP", isolate(input$`DE_PanelPlotCellSelection-Mod_clusterPP`)),
      c("DE_PanelPlotCellSelection-Mod_PPGrp", isolate(input$`DE_PanelPlotCellSelection-Mod_PPGrp`)),
      c("DE_panelplotPvalue", isolate(input$DE_panelplotPvalue))
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
  
  .schnappsEnv[["DE_panelPlot"]] <- list(plotFunc = af,
                                         scEx_log = scEx_log, 
                                         projections=projections, 
                                         genesin=genesin, dimx4=dimx4, 
                                         dimy4=dimy4, sameScale=sameScale, 
                                         nCol=nCol, sampdesc=sampdesc,
                                         cellNs=cellNs,
                                         applyPvalue = applyPvalue
  )
  retVal
})




# DE_panelPlotFact ----
#' DE_panelPlotFact
#' plot multiple panels for a given list of genes
#' If the x-axis is a categorical value and the y-axis is UMI.counts the y-axis related to
#' the count for that gene. Otherwise, all genes are used.
#' normalized counts are used for plotting
output$DE_panelPlotFact <- renderPlot({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_panelPlotFact")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_panelPlotFact", id = "DE_panelPlotFact", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$DE_panelPlotFact\n")
  
  clicked <- input$updatePanelPlotFact
  applyPvalue <- isolate(input$DE_panelplotFactPvalue)
  scEx_log <- scEx_log()
  projections <- projections()
  DE_updateInputPPt()
  factsin <- isolate(input$DE_pFactIds)
  # cl4 <- input$DE_clusterSelectionPanelPlot
  # ppgrp <- isolate(input$DE_PPGrp)
  # ppCluster <- isolate(input$DE_clusterPP)
  
  selectedCells <- isolate(dePanelFactCellSelection())
  cellNs <- isolate(selectedCells$cellNames())
  sampdesc <- isolate(selectedCells$selectionDescription())
  
  dimx4 <- isolate(input$DE_pFact_dim_x)
  dimy4 <- isolate(input$DE_pFact_dim_y)
  sameScale <- isolate(input$DE_panelplotFactSameScale)
  nCol <- isolate(as.numeric(input$DE_pFactnCol))
  
  if (is.null(scEx_log) | is.null(projections) | is.null(cellNs)) {
    return(NULL)
  }
  # debugControl("DE_panelPlotFact", list = c("scEx_log", "projections", "genesin",
  # "dimx4", "dimy4", "sameScale", "nCol", "sampdesc" , "cellNs"))
  # cp = load(file="~/SCHNAPPsDebug/DE_panelPlotFact.RData")
  
  # genesin <- toupper(genesin)
  # genesin <- gsub(" ", "", genesin, fixed = TRUE)
  # genesin <- strsplit(genesin, ",")
  # genesin <- genesin[[1]]
  
  # if (DEBUG) cat(file = stderr(), paste("output:sampdesc",sampdesc,"\n"))
  retVal <- panelPlotFactFunc(scEx_log, projections, factsin, dimx4, dimy4, sameScale, nCol, sampdesc, cellNs, applyPvalue = applyPvalue) 
  
  setRedGreenButton(
    vars = list(
      c("DE_pFactIds", isolate(input$DE_pFactIds)),
      c("DE_pFact_dim_x", isolate(input$DE_pFact_dim_x)),
      c("DE_pFact_dim_y", isolate(input$DE_pFact_dim_y)),
      c("DE_panelplotFactSameScale", isolate(input$DE_panelplotFactSameScale)),
      c("DE_pFactnCol", isolate(input$DE_pFactnCol)),
      c("DE_PanelPlotFactCellSelection-Mod_clusterPP", isolate(input$`DE_PanelPlotFactCellSelection-Mod_clusterPP`)),
      c("DE_PanelPlotFactCellSelection-Mod_PPGrp", isolate(input$`DE_PanelPlotFactCellSelection-Mod_PPGrp`)),
      c("DE_panelplotFactPvalue", isolate(input$DE_panelplotFactPvalue))
    ),
    button = "updatePanelPlotFact"
  )
  
  
  printTimeEnd(start.time, "DE_panelPlotFact")
  exportTestValues(DE_panelPlotFact = {
    ls()
  })
  af = panelPlotFunc
  # remove env because it is too big
  environment(af) = new.env(parent = emptyenv())
  
  .schnappsEnv[["DE_panelPlotFact"]] <- list(plotFunc = af,
                                         scEx_log = scEx_log, 
                                         projections=projections, 
                                         dimx4=dimx4, 
                                         dimy4=dimy4, sameScale=sameScale, 
                                         nCol=nCol, sampdesc=sampdesc,
                                         cellNs=cellNs,
                                         applyPvalue = applyPvalue
  )
  retVal
})

#
# 
# 
# Scater QC ----
#
# 
# 

emptyImage = list(
  src =  normalizePath("www/images/schnappsLogo.png",mustWork = F),
  contentType = "image/png",
  width = 500,
  height = 500,
  alt = "Scater plot will be here when 'apply changes' is checked"
)

createScaterPNG <- function(scaterReads, n, scols, width=NULL, height=NULL, DEBUG, outfile) {
  if (DEBUG) cat(file = stderr(), "function: createScaterPNG\n")
  # save(file = "~/SCHNAPPsDebug/createScaterPNG.RData", list = c(ls()))
  # cp=load(file='~/SCHNAPPsDebug/createScaterPNG.RData')
  if (DEBUG) cat(file = stderr(), "function: createScaterPNG2\n")
  p1 = pltHighExp( scaterReads, n, scols) 
  if (DEBUG) cat(file = stderr(), "function: createScaterPNG3\n")
  # calculations
  if (is.null(width)) {
    width <- 96 * 7
  }
  if (is.null(height)) {
    height <- 96 * 7
  }
  myPNGwidth <- width / 96
  myPNGheight <- height / 96
  
  if (DEBUG) cat(file = stderr(), "function: createScaterPNG4\n")
  tryCatch(
    ggsave(file = normalizePath(outfile, mustWork = FALSE), plot = p1, width = myPNGwidth, height = myPNGheight, units = "in"),
    error = function(e) {
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("Problem saving ggplot", type = "warning", duration = NULL)
      }
      return(emptyImage)
    }
  )
  retVal <- list(
    src = normalizePath(outfile, mustWork = FALSE),
    contentType = "image/png",
    width = width,
    height = height,
    alt = "Scater plot should be here"
  )
  if (DEBUG) cat(file = stderr(), "function: createScaterPNG5\n")
  
  return(retVal)
}

# click on start button
observeEvent(input$runScater,{
  if (.schnappsEnv$DEBUG) cat(file = stderr(), "observeEvent: detachedProc$runScater\n")
  if (!is.null(detachedProc$process)){
    return()
  }
  start.time <- base::Sys.time()
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_scaterPNG", id = "DE_scaterPNG", duration = NULL)
    removeNotification(id="DE_scaterPNG_Error")
  }
  scaterReads <- isolate(scaterReads())
  scols <- isolate(projectionColors$sampleNames)
  maxMemory = isolate(input$maxMemory)
  
  if (is.null(scaterReads)){
    removeNotification(id="DE_scaterPNG")
    detachedProc$result <- emptyImage
    return()
  }
  # width <- session$clientData$output_plot_width
  # height <- session$clientData$output_plot_height
  width <- NULL
  height <- NULL
  # outfile has to be set outside of the future since it will be removed after the session closes.
  outfile <- paste0(tempdir(), "/scaterPlot.png")
  
  n <- min(nrow(scaterReads), 50)
  # browser()
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scater.RData", list = c(ls()))
  }
  # cp=load(file='~/SCHNAPPsDebug/scater.RData')
  detachedProc$result <- emptyImage
  pl=plan()
  # span process
  # detachedProc$process <- mcparallel({
  #   createScaterPNG(scaterReads, n, scols, width=width, height=height)
  # })
  #span the process/function call
  # This should be set globally as it is also not reset here
  # options(future.globals.maxSize= maxMemory * 1024^3)

  detachedProc$process <- tryCatch({
    future({
      # detachedProc$process$pid = Sys.getpid()
      createScaterPNG(scaterReads=scaterReads, n=n, scols=scols, width=NULL, height=NULL, DEBUG = DEBUG, outfile=outfile)
    },seed=NULL,
    packages = "scater",
    globals = list(createScaterPNG=createScaterPNG, emptyImage=emptyImage, DEBUG=.schnappsEnv$DEBUG, pltHighExp=pltHighExp,
                   scaterReads=scaterReads, n=n, scols=scols, width=NULL, height=NULL, outfile=outfile), # we specify all variables with the function call
    lazy = FALSE, #start immediatly
    stdout = structure(TRUE, drop = TRUE)
    )},
    error = function(e) {
      cat(file = stderr(), paste("\n\n!!!Error during detach process:", e, "\n\nDo you need to increase the memory?\n\n"))
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("DE_scaterPNG ERROR", id = "DE_scaterPNG_Error", duration = NULL, type = "error")
      }
      return(NULL)
    }
  )
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scater2.RData", list = c(ls()))
  }
  # cp=load(file='~/SCHNAPPsDebug/scater2.RData')
  # browser()
  if(is.null(detachedProc$process)) return(NULL)
  activateObserver(1)
  if(!is.null(detachedProc$process$process))
    cat(file = stderr(), paste("input$start",detachedProc$process$process$get_pid(),"me:",Sys.getpid(),"\n"))
  if("callr" %in% class(pl)){
    detachedProc$PID = detachedProc$process$process$get_pid()
  }else{
    # if("multisession" %in% class(pl)){
    #   currentWorkerPIDs = getWorkerPIDs()
    # }
    #
    # please use callr otherwise we cannot kill process (for now)
    # 
    detachedProc$PID = NULL
  }
  
  detachedProc$startTime = start.time
  detachedProc$msg <- sprintf("%1$s started", detachedProc$process$pid)
})
#
# Stop the process
#
observeEvent(input$stopScater, {
  if(.schnappsEnv$DEBUG) cat(file = stderr(), "input$stopScater\n")

  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/stopScater.RData", list = c(ls()))
  }
  # cp=load(file='~/SCHNAPPsDebug/stopScater.RData')
  if (!is.null(detachedProc$PID)) {
    if("running" == detachedProc$process$state){
      #For windows
      # system(sprintf("taskkill /F /PID %s", v[[i]]))
      
      #For Linux
      system(sprintf("kill -9 %s", detachedProc$PID))
      activateObserver(0)
      detachedProc$PID = NULL
      detachedProc$process = NULL
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "DE_scaterQC")
        removeNotification(id = "DE_scaterPNG")
      }
      
    }
  }
})

#
# Handle process event
#
observe({
  if (.schnappsEnv$DEBUG) cat(file = stderr(), "observeEvent: detachedProc$process\n")
  # this will re-execute the collection process of mcparallel
  # if(!is.null(detachedProc$process))
  if(activateObserver()>0)
    invalidateLater(500, session)
  # browser()
  if (.schnappsEnv$DEBUG) cat(file = stderr(), "observeEvent: detachedProc$process2\n")
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/processScater.RData", list = c(ls()))
  }
  # cp=load(file='~/SCHNAPPsDebug/processScater.RData')
  
    isolate({
    if(resolved(detachedProc$process))
      if(!is.null(detachedProc$process)){
        # browser()
        detachedProc$result <- value(detachedProc$process)
        result = detachedProc$result
        # save(file = "~/SCHNAPPsDebug/createScaterPNGprocess.RData", list = c("result"))
        # cp=load(file='~/SCHNAPPsDebug/createScaterPNGprocess.RData')
        detachedProc$process <- NULL
        detachedProc$PID = NULL
        activateObserver(0)
        if (.schnappsEnv$DEBUG) cat(file = stderr(), "observeEvent: detachedProc$process3\n")
        printTimeEnd(detachedProc$startTime, "DE_scaterPNG")
        if (.schnappsEnv$DEBUG) cat(file = stderr(), "observeEvent: detachedProc$process4\n")
        if (!is.null(getDefaultReactiveDomain())) {
          removeNotification( id = "DE_scaterPNG")
        }
        
      }
    if (.schnappsEnv$DEBUG) cat(file = stderr(), "observeEvent: detachedProc$process5\n")
    
  })
})

output$DE_scaterQC <- renderImage(deleteFile = F, {
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_scaterQC")
    }
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_scaterQC", id = "DE_scaterQC", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "renderImage. output$DE_scaterQC\n")
  # browser()
  result = detachedProc$result
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/result1Scater.RData", list = c(ls()))
  }
  # cp=load(file='~/SCHNAPPsDebug/result1Scater.RData')
  
  scaterReads <- isolate(scaterReads())
  if (is.null(scaterReads) | is.null(result)) {
    return(emptyImage)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/result2Scater.RData", list = c(ls()))
  }
  # cp=load(file='~/SCHNAPPsDebug/result2Scater.RData')
  af = pltHighExp
  # remove env because it is too big
  environment(af) = new.env(parent = emptyenv())
  n <- min(nrow(scaterReads), 50)
  scols <- isolate(projectionColors$sampleNames)
  .schnappsEnv[["DE_scaterPNG"]] <- list(plotFunc = af,
                                         # plotHighestExprs = plotHighestExprs,
                                         scaterReads = scaterReads, 
                                         n = n,
                                         scols = scols
  )
  setRedGreenButton(
    vars = list(
      c("scaterRan", 1)
    ),
    button = "runScater"
  )
  
  exportTestValues(DE_scaterPNG = {
    result
  })
  result
})
#  Cannot call `bindCache()` on this object because it is marked as not cacheable.
# %>% bindCache(scaterReads())

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
  if(is.null(selectedCells)) return(NULL)
  cellNs <- selectedCells$cellNames()
  sampdesc <- selectedCells$selectionDescription()
  prj <- isolate(selectedCells$ProjectionUsed())
  prjVals <- isolate(selectedCells$ProjectionValsUsed())
  x = input$DE_expclusters_x
  y = input$DE_expclusters_y
  z = input$DE_expclusters_z
  # dimCol <- input$DE_expclusters_col
  # scols <- projectionColors$sampleNames
  # ccols <- projectionColors$dbCluster
  pc = projectionColors %>% reactiveValuesToList()
  if (is.null(scEx_log) | is.null(projections) | is.null(cellNs) ) {
    return(NULL)
  }
  # debugControl("DE_tsne_plt", list = c(ls()))
  # cp = load(file="~/SCHNAPPsDebug/DE_tsne_plt.RData")
  
  # projections$ExpressionColor = 
  dimCol = "ExpressionColor"
  featureData <- rowData(scEx_log)
  geneid <- geneName2Index(g_id, featureData)
  if (length(geneid) == 0) {
    return(NULL)
  }
  
  if (length(geneid) == 1) {
    projections$ExpressionColor <- assays(scEx_log)[[1]][geneid, ]
  } else {
    projections$ExpressionColor <- Matrix::colSums(assays(scEx_log)[[1]][geneid, ])
  }
  
  retVal <- tsnePlot(projections, x,y,z, dimCol, projColors=pc) 
  
  # retVal <- DE_dataExpltSNEPlot(scEx_log[,cellNs], g_id, projections[cellNs, ], x,y,z)
  
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
    
    if (is.null(scEx) | is.null(scEx_log)) {
      return(NULL)
    }
    # debugControl("DE_downloadPanel", list = c(ls()))
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

# observers for input parameters and

observe(label = "observe DE_seuratRefBased", {
  if (DEBUG) cat(file = stderr(), "observe DE_seuratRefBased\n")
  if (is.null(input$updateNormalization)) {
    return(NULL)
  }
  if (!input$normalizationRadioButton == "DE_seuratRefBased") {
    return(NULL)
  }
  if (!input$whichscLog == "calcLog") {
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


# obs.updateNormalization DE_seuratSCTnormButton ----
observe(label = "ob DE_seuratSCTnormButtonOldVal", {
  buttonPressed <- input$updateNormalization
  radioButtonVal <- isolate(input$normalizationRadioButton)
  if (!exists("DE_seuratSCTnormButtonOldVal", envir = .schnappsEnv)) {
    .schnappsEnv$DE_seuratSCTnormButtonOldVal <- 0
  }
  if (is.null(radioButtonVal)) {
    radioButtonVal <- ""
  }
  if (is.null(buttonPressed)) {
    buttonPressed <- 0
  }
  
  # changing the reactive DE_logGeneNormalizationButton will trigger the recalculation
  if (radioButtonVal == "DE_seuratSCTnorm" &
      !.schnappsEnv$DE_seuratSCTnormButtonOldVal == buttonPressed) {
    cat(file = stderr(), green(paste("\n=====changing value\n")))
    DE_seuratSCTnormButton(buttonPressed)
    .schnappsEnv$DE_seuratSCTnormButtonOldVal <- buttonPressed
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
  if (!input$whichscLog == "calcLog") {
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
  if (!input$whichscLog == "calcLog") {
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
  if (!input$whichscLog == "calcLog") {
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
