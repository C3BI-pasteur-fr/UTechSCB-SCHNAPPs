require(ElPiGraph.R)
require(plyr)
library(dplyr)
library(SCORPIUSbj)


# Scorpius ----------------------------------------------------------------
## observe scorpius Proj ---
# update projections: add new projections
observe(label = "ob20sc", {
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "observeScorpiusProj")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "observeScorpiusProj")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("observeScorpiusProj", id = "observeScorpiusProj", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "observeScorpiusProj 20sc started.\n")
  
  traj <- scorpiusTrajectory()
  scEx_log <- Scorpius_scEx_log()
  scExNames <- isolate(colnames(scEx_log()))
  
  isolate({
    prjs <- sessionProjections$prjs
  })
  if (is.null(scEx_log) || is.null(traj)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/observeScorpiusProj.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/observeScorpiusProj.RData")
  cn <- "traj_scorpius"
  # if (cn %in% colnames(prjs)) {
  #   return(NULL)
  # }
  # browser()
  if (ncol(prjs) > 0) {
    # # make sure we are working with the correct cells. This might change when cells were removed.
    # prjs <- prjs[colnames(scEx_log), , drop = FALSE]
    # didn't find a way to easily overwrite columns
    
    # if (cn %in% colnames(prjs)) {
    prjs[,cn] <- -1
    prjs[rownames(traj), cn] <- traj$time
    # } else {
    #   prjs <- base::cbind(prjs, traj$time, deparse.level = 0)
    #   colnames(prjs)[ncol(prjs)] <- cn
    # }
    
  } else {
    prjs <- data.frame(row.names = scExNames)
    prjs[,cn] <- -1
    prjs[colnames(scEx_log),cn] <- traj$time
    # colnames(prjs)[ncol(prjs)] <- cn
  }
  sessionProjections$prjs <- prjs
})

observe(label = "ob_scorpButton",priority = 99,
        {
          if (DEBUG) cat(file = stderr(), "observe ob_scorpButton\n")
          input$updatetScorpiusParameters
          setRedGreenButtonCurrent(
            vars = list(
              c("scorpDimX", input$dimScorpiusX),
              c("scorpDimY", input$dimScorpiusY),
              c("scorpMaxGenes", input$scorpMaxGenes),
              c("scorpRepeat", input$scorpRepeat),
              c("scorpInFile", input$trajInputFile)
            )
          )
          
          updateButtonColor(buttonName = "updatetScorpiusParameters", parameters = c(
            "scorpDimX", "scorpDimY",
            "scorpMaxGenes", "scorpRepeat" , "scorpInFile"
          ))
        })

observe(label = "ob1", {
  if (DEBUG) cat(file = stderr(), "observe ob1 scorpius\n")
  # .schnappsEnv$dimScorpiusX <- input$dimScorpiusX
  # .schnappsEnv$dimScorpiusY <- input$dimScorpiusY
  # .schnappsEnv$dimScorpiusCol <- input$dimScorpiusCol
})
# updateScorpiusInput <- reactive({
observe(label = "ob2", {
  projections <- scorpius_projections()
  
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  
  # Can also set the label and select items
  if (is.null(.schnappsEnv$dimScorpiusX)) {
    .schnappsEnv$dimScorpiusX <- colnames(projections)[1]
  }
  if (is.null(.schnappsEnv$dimScorpiusY)) {
    .schnappsEnv$dimScorpiusY <- colnames(projections)[2]
  }
  if (is.null(.schnappsEnv$dimScorpiusCol)) {
    .schnappsEnv$dimScorpiusCol <- "dbCluster"
  }
  updateSelectInput(session, "dimScorpiusX",
                    choices = colnames(projections),
                    selected = defaultValue("dimScorpiusX", "notyet")
  )
  
  updateSelectInput(session, "dimScorpiusY",
                    choices = colnames(projections),
                    selected = defaultValue("dimScorpiusY", "notyet")
  )
  updateSelectInput(session, "dimScorpiusCol",
                    choices = colnames(projections),
                    selected = defaultValue("dimScorpiusCol", "notyet")
  )
  # updateNumericInput(session, "scorpMaxGenes",
  #                   value = 500
  # )
})
## elpi observers ----
observe(label = "ob3", {
  if (DEBUG) cat(file = stderr(), paste0("observe: dimElpi\n"))
})
# observe(label = "ob4", {
#   if (DEBUG) cat(file = stderr(), paste0("observe: dimElpiX\n"))
#   .schnappsEnv$dimElpi <- input$dimElpi
#   .schnappsEnv$dimElpiX <- input$dimElpiX
#   .schnappsEnv$dimElpiY <- input$dimElpiY
#   .schnappsEnv$dimElpiCol <- input$dimElpiCol
#   .schnappsEnv$elpiSeed <- input$elpiSeed
#   .schnappsEnv$ElpiMethod <- input$ElpiMethod
#   .schnappsEnv$elpinReps <- input$elpinReps
#   .schnappsEnv$elpiNumNodes <- input$elpiNumNodes
#   .schnappsEnv$elpiProbPoint <- input$elpiProbPoint
#   .schnappsEnv$elpiEndNode <- input$elpiEndNode
#   .schnappsEnv$elpiStartNode <- input$elpiStartNode
#   .schnappsEnv$elpi_num_permutations <- input$elpi_num_permutations
#   .schnappsEnv$elpi_ntree <- input$elpi_ntree
#   .schnappsEnv$elpi_ntree_perm <- input$elpi_ntree_perm
#   .schnappsEnv$elpi_nGenes <- input$elpi_nGenes
# })

## observe traj_endpoints ----
observe(label = "ob18", {
  endpoints <- traj_endpoints()
  
  # browser()
  # Can use character(0) to remove all choices
  if (is.null(endpoints)) {
    return(NULL)
  }
  # # save endpoint in global variable to make sure that we don't update unnecessarily
  if (!is.null(.schnappsEnv$elpiEndpoints) &
      length(endpoints) == length(.schnappsEnv$elpiEndpoints) &
      all(sort(endpoints) == sort(.schnappsEnv$elpiEndpoints))
  ) {
    return(NULL)
  }
  if (DEBUG) cat(file = stderr(), "observeProj 18 started.\n")
  
  .schnappsEnv$elpiEndpoints <- endpoints
  # Can also set the label and select items
  updateSelectInput(session,
                    inputId = "elpiStartNode",
                    choices = endpoints,
                    selected = endpoints[1]
  )
  updateSelectInput(session,
                    inputId = "elpiEndNode",
                    choices = endpoints,
                    selected = endpoints[length(endpoints)]
  )
})

## observe projections ----
observe(label = "ob19", {
  projections <- Elpi_projections()
  projFactors = projFactors()
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  if (DEBUG) cat(file = stderr(), "observeProj 19 started.\n")
  
  # Can also set the label and select items
  updateSelectInput(session, "dimElpiX",
                    choices = colnames(projections),
                    selected = defaultValue("dimElpiX","tsne1")
  )
  
  updateSelectInput(session, "dimElpiY",
                    choices = colnames(projections),
                    selected = defaultValue("dimElpiY","tsne2")
  )
  updateSelectInput(session, "dimElpiCol",
                    choices = projFactors,
                    selected = defaultValue("dimElpiCol","dbCluster")
  )
  updateNumericInput(session, "elpiSeed",
                     value = defaultValue("elpiSeed",9)
  )
  updateSelectInput(session, "dimElpi",
                    selected = defaultValue("dimElpi","components")
  )
  updateSelectInput(session, "ElpiMethod",
                    selected = defaultValue("ElpiMethod","computeElasticPrincipalTree")
  )
  updateNumericInput(session, "elpiNumNodes",
                     value = defaultValue("elpiNumNodes",20)
  )
  updateNumericInput(session, "elpinReps",
                     value = defaultValue("elpinReps",1)
  )
  updateNumericInput(session, "elpiProbPoint",
                     value = defaultValue("elpiProbPoint",0.6)
  )
  updateNumericInput(session, "elpi_num_permutations",
                     value = defaultValue("elpi_num_permutations",3)
  )
  updateNumericInput(session, "elpi_ntree",
                     value = defaultValue("elpi_ntree",10000)
  )
  updateNumericInput(session, "elpi_ntree_perm",
                     value = defaultValue("elpi_ntree_perm",1000)
  )
  updateNumericInput(session, "elpi_nGenes",
                     value = defaultValue("elpi_nGenes",50)
  )
})
## observeProj ----
# update projections: add new projections
observe(label = "ob20", {
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "observeProj")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "observeProj")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("observeProj", id = "observeProj", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "observeProj 20 started.\n")
  
  startNode <- input$elpiStartNode
  endNode <- input$elpiEndNode
  elpimode <- input$ElpiMethod
  psTime <- traj_getPseudotime()
  scEx_log <- Elpi_scEx_log()
  scExNames <- isolate(colnames(scEx_log()))
  isolate({
    prjs <- sessionProjections$prjs
  })
  
  if (is.null(scEx_log) || is.null(psTime) || elpimode == "computeElasticPrincipalCircle") {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/observeProj.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/observeProj.RData")
  cn <- paste0("traj_", startNode, "_", endNode)
  # to avoid overwriting
  # TODO but what happens if we want to overwrite!!!
  # 
  if (cn %in% colnames(prjs)) {
    return(NULL)
  }
  # browser()
  if (ncol(prjs) > 0) {
    # make sure we are working with the correct cells. This might change when cells were removed.
    # prjs <- prjs[scExNames, , drop = FALSE]
    # didn't find a way to easily overwrite columns
    
    # if (cn %in% colnames(prjs)) {
    prjs = as.data.frame(prjs)
    prjs[, cn] = -1
    prjs[colnames(scEx_log), cn] <- psTime$Pt
    # } else {
    #   prjs <- base::cbind(prjs[colnames(scEx_log),], psTime$Pt, deparse.level = 0)
    #   colnames(prjs)[ncol(prjs)] <- cn
    # }
  } else {
    prjs <- data.frame(row.names = scExNames)
    prjs[, cn] = -1
    prjs[colnames(scEx_log),cn] <- psTime$Pt
    # colnames(prjs)[ncol(prjs)] <- cn
  }
  sessionProjections$prjs <- prjs
})

## set scorpiuseParameters on button pressed----
observeEvent(input$updatetScorpiusParameters,{
  # only react on this value
  cat(file = stderr(), paste("hit button updatetScorpiusParameters\n"))
  input$updatetScorpiusParameters
  # scorpiuseParameters$dimX = isolate(input$dimScorpiusX)
  # scorpiuseParameters$dimY = isolate(input$dimScorpiusY)
  # scorpiuseParameters$scInput = isolate(scorpiusInput())
  # scorpiuseParameters$scorpMaxGenes = isolate(input$scorpMaxGenes)
  # scorpiuseParameters$scorpRepeat = isolate(input$scorpRepeat)
  
  setRedGreenButton(
    vars = list(
      c("scorpDimX", isolate(input$dimScorpiusX)),
      c("scorpDimY", isolate(input$dimScorpiusY)),
      c("scorpMaxGenes", isolate(input$scorpMaxGenes)),
      c("scorpRepeat", isolate(input$scorpRepeat)),
      c("scorpInFile", isolate(input$trajInputFile))
    ),
    button = "updatetScorpiusParameters"
  )
  .schnappsEnv$defaultValues[["dimScorpiusX"]] <- isolate(input$dimScorpiusX)
  .schnappsEnv$defaultValues[["dimScorpiusY"]] <- isolate(input$dimScorpiusY)
  .schnappsEnv$defaultValues[["dimScorpiusCol"]] <- isolate(input$dimScorpiusCol)
  .schnappsEnv$defaultValues[["scorpMaxGenes"]] <- isolate(input$scorpMaxGenes)
  .schnappsEnv$defaultValues[["scorpRepeat"]] <- isolate(input$scorpRepeat)
  
})


## scropius_trajectory_plot ----
# The output type has to be in line with the tablist item. I.e. plotOutput in this case
output$scropius_trajectory_plot <- renderPlot({
  if (DEBUG) cat(file = stderr(), "scropius_trajectory_plot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scropius_trajectory_plot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scropius_trajectory_plot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scropius_trajectory_plot", id = "scropius_trajectory_plot", duration = NULL)
  }
  # sdi = Scorpius_dataInput() # variable not used but display depends on this => causes endless loop??
  traj <- scorpiusTrajectory()
  projections <- isolate(scorpius_projections())
  space <- isolate(scorpiusSpace())
  # upI <- updateScorpiusInput() # needed to update input
  dimX <- input$dimScorpiusX
  dimY <- input$dimScorpiusY
  dimCol <- input$dimScorpiusCol
  sampCol <- sampleCols$colPal
  ccols <- clusterCols$colPal
  
  # doCalc <- input$scorpiusCalc
  
  if (is.null(projections) ) {
    if (DEBUG) cat(file = stderr(), "scropius_trajectory_plot:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scropius_trajectory_plot.RData", list = c(".schnappsEnv", ls()))
  }
  # cp =load(file="~/SCHNAPPsDebug/scropius_trajectory_plot.RData")
  
  vChanged = valuesChanged(parameters = c(
    "scorpDimX", "scorpDimY",
    "scorpMaxGenes", "scorpRepeat" , "scorpInFile"
  ))
  
  prj = projections[,dimCol]
  mycolPal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(
    n = 12, name =
      "Paired"
  ))(length(levels(prj)))
  names(mycolPal) = levels(prj)
  if (dimCol == "sampleNames") {
    mycolPal <- sampCol
  }
  if (dimCol == "dbCluster") {
    mycolPal <- ccols
  }
  
  
  if (is.null(space) |is.null(traj) | vChanged) {
    if (vChanged) {
      cat(file = stderr(), "scropius Values changed\n")
    }
    p1 <- ggplot2::ggplot(projections, aes_string(dimX, dimY)) + 
      ggplot2::geom_point(colour = mycolPal[projections[,dimCol]]) + 
      theme_classic()
    return(p1)
  } else {
  
  # space <- projections[, c(dimX, dimY)]
  require(SCORPIUSbj)
  # traj <- SCORPIUSbj::infer_trajectory(space)
  # dimCol="CELLTYPES"
  colnames(traj) <- c("Comp1", "Comp2", "time")
  p1 = draw_trajectory_plot(space, progression_group = projections[rownames(space), dimCol], 
                       progression_group_palette = mycolPal[projections[,dimCol]],
                       path = as.matrix(traj[, 1:2]))
  }
  return(p1)
  
})

callModule(tableSelectionServer, "scorpiusTableMod", scorpiusModulesTable)
# selected clusters heatmap module

scorpiusHeatmapPlotReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusHeatmapPlotReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusHeatmapPlotReactive")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scorpiusHeatmapPlotReactive")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusHeatmapPlotReactive", id = "scorpiusHeatmapPlotReactive", duration = NULL)
  }
  
  # upI <- updateScorpiusInput() # needed to update input
  scEx = Scorpius_scEx_log()
  projections <- isolate(scorpius_projections())
  traj <- scorpiusTrajectory()
  expr_sel <- scorpiusExpSel()
  modules <- scorpiusModules()
  sampCol <- isolate(sampleCols$colPal)
  ccols <- isolate(clusterCols$colPal)
  
  dimCol <- input$dimScorpiusCol
  # doCalc <- input$scorpiusCalc
  pixelratio <- session$clientData$pixelratio
  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  
  # browser()
  if (is.null(projections) | is.null(modules) | is.null(expr_sel) | is.null(traj)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusHeatmapPlot:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusHeatmapPlot.RData", list = c(ls()))
  }
  # cp=load(file="~/SCHNAPPsDebug/scorpiusHeatmapPlot.RData")
  
  if (is.null(pixelratio)) pixelratio <- 1
  if (is.null(width)) {
    width <- 96 * 7
  } # 7x7 inch output
  if (is.null(height)) {
    height <- 96 * 7
  }
  annCols <- list(
    "sampleNames" = sampCol,
    "dbCluster" = ccols
  )
  
  outfile <- paste0(tempdir(), "/heatmapScorpius", base::sample(1:10000, 1), ".png")
  cat(file = stderr(), paste("saving to: ", outfile, "\n"))
  colnames(expr_sel$expr_sel) = rowData(scEx[colnames(expr_sel$expr_sel),])$symbol
  # modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = F)
  retVal <- drawTrajectoryHeatmap(x = expr_sel$expr_sel, time = traj$time, 
                                  progression_group = projections[rownames(expr_sel$expr_sel), dimCol], 
                                  modules = modules, 
                                  annotation_colors = annCols,
                                  show_labels_row = TRUE, show_labels_col = FALSE,scale_features = TRUE,
                                  filename = normalizePath(outfile, mustWork = FALSE)
  )
  
  exportTestValues(scorpiusHeatmapPlotReactive = {
    retVal
  })
  return(retVal)
})

callModule(
  pHeatMapModule,
  "scorpiusHeatmapPlotModule",
  scorpiusHeatmapPlotReactive
)

output$downLoadTraj <- downloadHandler(
  filename = paste0("scorpiusTraj.", Sys.Date(), ".csv"),
  content = function(file) {
    if (DEBUG) cat(file = stderr(), "downLoadTraj started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "downLoadTraj")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "downLoadTraj")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("downLoadTraj", id = "downLoadTraj", duration = NULL)
    }
    
    traj <- scorpiusTrajectory()
    if (is.null(traj)) {
      return(NULL)
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/downLoadTraj.RData", list = c(ls()))
    }
    # load(file="~/SCHNAPPsDebug/downLoadTraj.RData")
    write.csv(traj, file)
  }
)






# Elpi Graph -------------------------------------------------------------
## obeserver for elpi button -
## observe: cellNameTable_rows_selected =====
observe(label = "ob_elpiCalcParams", {
  if (DEBUG) cat(file = stderr(), "observe elpiCalc\n")
  input$elpiCalc
  setRedGreenButtonCurrent(
    vars = list(
      c("elpinReps", input$elpinReps),
      c("elpiNumNodes", input$elpiNumNodes),
      c("elpiProbPoint", input$elpiProbPoint),
      c("ElpiMethod", input$ElpiMethod),
      c("dimElpi", input$dimElpi),
      c("elpiSeed", input$elpiSeed),
      c("dimElpi", input$dimElpi),
      c("dimElpiX", input$dimElpiX),
      c("dimElpiY", input$dimElpiY)
    )
  )
  
  updateButtonColor(buttonName = "elpiCalc", parameters = c(
    "elpinReps","elpiNumNodes","elpiProbPoint","ElpiMethod","dimElpi",
    "elpiSeed","dimElpi","dimElpiX","dimElpiY"
  ))
})

## elpiHeatmapPlotReactive ----
elpiHeatmapPlotReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "elpiHeatmapPlotReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpiHeatmapPlotReactive")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "elpiHeatmapPlotReactive")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("elpiHeatmapPlotReactive", id = "elpiHeatmapPlotReactive", duration = NULL)
  }
  
  # upI <- updateScorpiusInput() # needed to update input
  scEx = Elpi_scEx()
  clicked <- input$elpiCalc
  projections <- Elpi_projections()
  psTime <- traj_getPseudotime()
  expr_sel <- traj_elpi_gimp()
  modules <- traj_elpi_modules()
  
  dimCol <- isolate(input$dimElpiCol)
  pixelratio <- session$clientData$pixelratio
  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  
  
  
  if (is.null(projections) | is.null(modules) | is.null(expr_sel) | is.null(psTime)) {
    if (.schnappsEnv$DEBUG) cat(file = stderr(), paste("elpiHeatmapPlot:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/elpiHeatmapPlotReactive.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/elpiHeatmapPlotReactive.RData")
  # browser()
  if (is.null(pixelratio)) pixelratio <- 1
  if (is.null(width)) {
    width <- 96 * 7
  } # 7x7 inch output
  if (is.null(height)) {
    height <- 96 * 7
  }
  colnames(expr_sel$expr_sel) = rowData(scEx[colnames(expr_sel$expr_sel),])$symbol
  
  outfile <- paste0(tempdir(), "/heatmapScorpius", base::sample(1:10000, 1), ".png")
  cat(file = stderr(), paste("saving to: ", outfile, "\n"))
  outfile = NULL
  pst <- psTime$Pt[which(!is.na(psTime$Pt))]
  # modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = F)
  retVal <- drawTrajectoryHeatmap(x =expr_sel$expr_sel,
                                  time = pst, progression_group = projections[rownames(expr_sel$expr_sel), dimCol], 
                                  modules = modules, show_labels_col = FALSE, show_labels_row = TRUE,
                                  filename = normalizePath(outfile, mustWork = FALSE)
  )
  
  exportTestValues(scorpiusHeatmapPlotReactive = {
    retVal
  })
  return(retVal)
})


## elpiHeatmapPlotModule -----
callModule(
  pHeatMapModule,
  "elpiHeatmapPlotModule",
  elpiHeatmapPlotReactive
)

#
# output$elpi_heatmap <- renderPlot({
#   start.time <- base::Sys.time()
#   on.exit({
#     printTimeEnd(start.time, "traj_getPseudotime")
#     if (!is.null(getDefaultReactiveDomain()))
#       removeNotification(id = "traj_getPseudotime")
#   })
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("traj_getPseudotime", id = "traj_getPseudotime", duration = NULL)
#   }
#   if (.schnappsEnv$DEBUG) cat(file = stderr(), "traj_getPseudotime started.\n")
#   scEx_log <- scEx_log()
#   projections <- projections()
#   TreeEPG <- elpiGraphCompute()
#   elpimode <- input$ElpiMethod
#   tree_data <- elpiTreeData()
#   tragetPath <- traj_tragetPath()
#   gene_sel <- traj_elpi_gimp()
#   modules <-  traj_elpi_modules()
#   psTime = traj_getPseudotime()
#
#   if (is.null(gene_sel) || is.null(scEx_log) || is.null(TreeEPG) || elpimode=="computeElasticPrincipalCircle" || is.null(psTime)) {
#     return(NULL)
#   }
#   if (.schnappsEnv$DEBUGSAVE) {
#     save(file = "~/SCHNAPPsDebug/traj_getPseudotime.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file="~/SCHNAPPsDebug/traj_getPseudotime.RData")
#
#   ## Select most important genes (set ntree to at least 10000!)
#   # gene_sel <- geneImport[1:50,]
#   gene_sel = gene_sel$gene_sel
#   expr_sel <- t(as.matrix(assays(scEx_log)[[1]][gene_sel$gene,which(!is.na(psTime$Pt))]))
#
#   pst = psTime$Pt[which(!is.na(psTime$Pt))]
#
#
#   p <- SCORPIUSbj::draw_trajectory_heatmap(x = expr_sel, time = pst, progression_group = projections$dbCluster[which(!is.na(psTime$Pt))] ,
#                                          modules=modules, show_labels_row = TRUE)
#
# })

## elpi_plot ----
output$elpi_plot <- renderPlot({
  if (DEBUG) cat(file = stderr(), "elpi_plot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpi_plot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "elpi_plot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("elpi_plot", id = "elpi_plot", duration = NULL)
  }
  
  doCalc <- input$elpiCalc
  
  tree_data <- elpiTreeData()
  cep <- elpiGraphCompute()
  PointLabel <- elpiPointLabel()
  projections <- Elpi_projections()
  dimX <- input$dimElpiX
  dimY <- input$dimElpiY
  dimCol <- input$dimElpiCol
  sampCol <- sampleCols$colPal
  ccols <- clusterCols$colPal
  
  vChanged = valuesChanged(parameters = c(
    "elpinReps","elpiNumNodes","elpiProbPoint","ElpiMethod",
    "elpiSeed","dimElpi","dimElpiX","dimElpiY"
  ))
  
  
  if (is.null(projections)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpi_plot.RData", list = c(base::ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/elpi_plot.RData")

  prj = projections[,dimCol]
  mycolPal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(
    n = 12, name =
      "Paired"
  ))(length(levels(prj)))
  
  if (dimCol == "sampleNames") {
    mycolPal <- sampCol
  }
  if (dimCol == "dbCluster") {
    mycolPal <- ccols
  }
  
  
    if (is.null(cep) | vChanged) {
    if (vChanged) {
      cat(file = stderr(), "elpi Values changed\n")
    }
    require(ggplot2)
    p1 <- ggplot(projections, aes_string(dimX, dimY, colour = dimCol)) + geom_point()
    p1 <- ggplot2::ggplot(projections, aes_string(dimX, dimY)) + 
      ggplot2::geom_point(colour = mycolPal[projections[,dimCol]]) + 
      theme_classic()
    return(p1)
  }
  
  if (is.null(tree_data) | is.null(cep)) {
    return(NULL)
  }
  
  
  require(ggrepel)
  require(igraph)
  
  NodeLabs <- 1:nrow(cep[[length(cep)]]$NodePositions)
  NodeLabs[degree(ConstructGraph(cep[[length(cep)]])) != 1] <- NA
  
  p <- PlotPG(
    X = tree_data, TargetPG = cep[[length(cep)]],
    NodeLabels = NodeLabs,
    LabMult = 5, PointSize = NA, p.alpha = .5,
    GroupsLab = projections[rownames(tree_data),dimCol]
  )
  
  q <- ggplot_build(p[[1]])
  q$data[[1]]$colour = mycolPal[projections[,dimCol]]
  q <- ggplot_gtable(q)
  p[[1]] = ggplotify::as.ggplot(q)
  

  
  # p <- PlotPG(X = tree_data, TargetPG = cep[[length(cep)]], GroupsLab = PointLabel, p.alpha = 0.9)
  # p[[1]] <- p[[1]] + geom_label_repel(
  #   data = plyr::ddply(p[[1]]$data, ~Group, summarise, meanA = mean(PCA), meanB = mean(PCB)),
  #   aes(x = meanA, y = meanB, label = Group),
  #   vjust = 1
  # )
  p
})

# output$elpi_moduleHeatmap <- renderPlot({
#   draw_trajectory_heatmap(expr_sel, pst, group_name, modules)
# })
## elpiTableMod ----
callModule(tableSelectionServer, "elpiTableMod", elpiModulesTable)

## elpi_histo ----
output$elpi_histo <- renderPlot({
  if (DEBUG) cat(file = stderr(), "elpi_histo started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpi_histo")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "elpi_histo")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("elpi_histo", id = "elpi_histo", duration = NULL)
  }
  
  PointLabel <- elpiPointLabel()
  if (is.null(PointLabel)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpi_histo.RData", list = c(base::ls()))
  }
  # load("~/SCHNAPPsDebug/elpi_histo.RData")
  
  barplot(table(PointLabel), las = 2, ylab = "Number of points")
})






# Tempora -----------------------------------------------------------------
## tempora_plot ====
output$tempora_plot <- renderPlot({
  if (DEBUG) cat(file = stderr(), "tempora_plot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "tempora_plot")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "tempora_plot")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tempora_plot", id = "tempora_plot", duration = NULL)
  }
  
  temporaObj <- temporaTrajectory()
  
  if(is.null(temporaObj)) {
    if (DEBUG) cat(file = stderr(), paste("tempora_plot:NULL\n"))
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/tempora_plot.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/tempora_plot.RData")
  # temporaObj = to
  
  af = PlotTrajectory
  # remove env because it is too big
  environment(af) = new.env(parent = emptyenv())
  
  .schnappsEnv[["TRAJ_tempora_plot"]] <- list(plotFunc = af,
                                              object = temporaObj
  )
  
  PlotTrajectory(temporaObj)
  
})

## tempora_screeplot ----
output$tempora_screeplot <- renderPlot({
  if (DEBUG) cat(file = stderr(), "tempora_screeplot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "tempora_screeplot")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "tempora_screeplot")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tempora_screeplot", id = "tempora_screeplot", duration = NULL)
  }
  
  temporaObj <- temporaPWProfiles()
  
  if(is.null(temporaObj)) {
    if (DEBUG) cat(file = stderr(), paste("tempora_screeplot:NULL\n"))
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/tempora_screeplot.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/tempora_screeplot.RData")
  
  
  af = screeplot
  # remove env because it is too big
  environment(af) = new.env(parent = emptyenv())
  
  .schnappsEnv[["TRAJ_tempora_screeplot"]] <- list(plotFunc = af,
                                                   x = temporaObj@cluster.pathways.dr,
                                                   npcs=25, 
                                                   type="lines",
                                                   main="PCA on pathway enrichment analysis result"
  )
  
  
  
  screeplot(temporaObj@cluster.pathways.dr, npcs=25, type="lines", main="PCA on pathway enrichment analysis result")
  
})

## temporaSelectedGOs ----
output$temporaSelectedGOs <- renderPlot({
  
  if (DEBUG) cat(file = stderr(), "temporaSelectedGOs started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaSelectedGOs")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "temporaSelectedGOs")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaSelectedGOs", id = "temporaSelectedGOs", duration = NULL)
  }
  
  temporaObj <- temporaIdentifyVaryingPWs()
  selectedGOs <- selectedGOs()
  if(is.null(temporaObj) | is.null(selectedGOs)) {
    if (DEBUG) cat(file = stderr(), paste("temporaSelectedGOs:NULL\n"))
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaSelectedGOs.RData", list = c(ls(), ".schnappsEnv"))
  }
  # cp = load(file="~/SCHNAPPsDebug/temporaSelectedGOs.RData")
  selectedGOs =  as.numeric(strsplit(selectedGOs,"," )[[1]])
  # max (temporaObj@varying.pws)
  if (length(selectedGOs)<1){
    return(NULL)
  }
  temporaObj@varying.pws = temporaObj@varying.pws[selectedGOs]
  # object=temporaObj
  # temporaObj@varying.pws = temporaObj@varying.pws[order(temporaObj@varying.pws)[1:10]]
  
  af = PlotVaryingPWs
  # remove env because it is too big
  environment(af) = new.env(parent = emptyenv())
  
  .schnappsEnv[["TRAJ_temporaSelectedGOs"]] <- list(plotFunc = af,
                                                    object = temporaObj
  )
  
  # browser()
  PlotVaryingPWs(temporaObj)
  
})

## install ggnetwork ----
if (!"ggnetwork" %in% installed.packages()){
  remotes::install_github("briatte/ggnetwork")
}

## tempora2dPlot ----
output$tempora2dPlot <- plotly::renderPlotly({
  if (DEBUG) cat(file = stderr(), "temporaScorpius2dPlot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaScorpius2dPlot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "temporaScorpius2dPlot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaScorpius2dPlot", id = "temporaScorpius2dPlot", duration = NULL)
  }
  
  temporaObj <- temporaTrajectory()
  projections <- isolate(projections())
  
  dimX <- input$dimTemporaX
  dimY <- input$dimTemporaY
  dimCol <- input$temporaCluster
  # doCalc <- input$scorpiusCalc
  
  if (is.null(projections) | is.null(temporaObj) ) {
    if (DEBUG) cat(file = stderr(), "temporaScorpius2dPlot:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaScorpius2dPlot.RData", list = c(ls()))
  }
  # cp =load(file="~/SCHNAPPsDebug/temporaScorpius2dPlot.RData")
  
  p2  = tempora2DPlotFunc(temporaObj, projections, dimX, dimY, dimCol)
  if(is.null(p2)){
    return(NULL)
  }
  
  af = tempora2DPlotFunc
  # remove env because it is too big
  environment(af) = new.env(parent = emptyenv())
  
  .schnappsEnv[["TRAJ_tempora2dPlot"]] <- list(plotFunc = af,
                                               temporaObj = temporaObj,
                                               projections = projections,
                                               dimX = dimX, 
                                               dimY = dimY,
                                               dimCol = dimCol
  )
  
  ggplotly(p2)
  # p2 
})


## temporaGOTableMod ----
selectedGOs <- callModule(tableSelectionServer, "temporaGOTableMod", temporaPvalModulesTable)


output$coE_temporaPWgenes = renderText({
  if (DEBUG) cat(file = stderr(), "coE_temporaPWgenes started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_temporaPWgenes")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "coE_temporaPWgenes")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_temporaPWgenes", id = "coE_temporaPWgenes", duration = NULL)
  }
  
  temporaObj <- temporaIdentifyVaryingPWs()
  sGOs = as.numeric(selectedGOs())
  tgi = temporaGeneIds()
  
  if (is.null(temporaObj) ) {
    if (DEBUG) cat(file = stderr(), paste("coE_temporaPWgenes:NULL\n"))
    return(NULL)
  }
  if (length(sGOs)<1) {
    if (DEBUG) cat(file = stderr(), paste("no pathway selected\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_temporaPWgenes.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/coE_temporaPWgenes.RData")
  goName = names(temporaObj@varying.pws)[sGOs]
  
  nn = grep(paste0("^",goName,"$"), names(tgi))
  if(is.null(nn)) return(NULL)
  if(length(nn) < 1) return(NULL)
  if(length(nn)>1){return(paste("more than one found:  ", names(tgi[nn]),"__",tgi[nn],  sep="", collapse=", "))}
  return(paste(tgi[[nn ]],  sep="", collapse=", "))
  
})
## observe dimX/Y ----

observe({
  .schnappsEnv$defaultValues[["dimTemporaX"]] = input$dimTemporaX
  .schnappsEnv$defaultValues[["dimTemporaY"]] = input$dimTemporaY
  
})

observe({
  projections = projections()
  if (is.null(projections)) {
    return(NULL)
  }
  
  updateSelectInput(session, "dimTemporaX",
                    choices = colnames(projections),
                    selected = defaultValue("dimTemporaX","tsne1")
  )
  updateSelectInput(session, "dimTemporaY",
                    choices = colnames(projections),
                    selected = defaultValue("dimTemporaY","tsne2")
  )
})

## observe temporaCluster/Factor ----
observe({
  projFactors = projFactors()
  if(is.null(projFactors)) {
    return(NULL)
  }
  
  updateSelectInput(session, "temporaCluster",
                    choices = projFactors,
                    selected = defaultValue("temporaCluster",projFactors[1])
  )
  updateSelectInput(session, "temporaFactor",
                    choices = projFactors,
                    selected = defaultValue("temporaFactor",projFactors[1])
  )
  
})

## observe temporaLevels ----
observe({
  projections <- projections()
  temporaFactor = input$temporaFactor
  selectedCells <- isolate(Tempora_dataInput()) #DE_Exp_dataInput
  if(is.null(selectedCells)) return(NULL)
  cellNs <- selectedCells$cellNames()
  
  
  if (DEBUG) cat(file = stderr(), paste("observe input$temporaFactor ", temporaFactor, ".\n"))
  if (is.null(projections)) {
    if (DEBUG) cat(file = stderr(), paste("temporaImport:NULL\n"))
    return(NULL)
  }
  projections = projections[cellNs,]
  cdat = projections
  if (! temporaFactor %in% colnames(cdat)) {
    return(NULL)
  }
  # save(file = "~/SCHNAPPsDebug/temporaobserveFactor.RData", list = c(ls()))
  # cp = load("~/SCHNAPPsDebug/temporaobserveFactor.RData")
  updateSelectInput(session, "temporaLevels",
                    choices = unique(cdat[,temporaFactor]),
                    selected = defaultValue("temporaFactor","tsne1"))
  
})


# saveToHistory ----
# save to history  observer ----
# save2Hist_tempora2dPlot
observe(label = "save2Hist_tempora2dPlot", {
  clicked  = input$save2Hist_tempora2dPlot
  if (DEBUG) cat(file = stderr(), "observe input$tempora2dPlot \n")
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "save2Hist_tempora2dPlot")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("save2Hist_tempora2dPlot", id = "save2Hist_tempora2dPlot", duration = NULL)
  }
  if (is.null(clicked)) return()
  if (clicked < 1) return()
  
  add2history(type = "save", input = isolate( reactiveValuesToList(input)), 
              comment = paste0("# tempora2dPlot\n",
                               "require(Tempora)\n",
                               "fun = plotData$plotData$plotFunc\n", 
                               "environment(fun) = environment()\n",
                               "print(do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)]))\n"
              ),
              plotData = .schnappsEnv[["TRAJ_tempora2dPlot"]])
  
})
# save2Hist_temporaSelectedGOs
observe(label = "save2Hist_temporaSelectedGOs", {
  clicked  = input$save2Hist_temporaSelectedGOs
  if (DEBUG) cat(file = stderr(), "observe input$temporaSelectedGOs \n")
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "save2Hist_temporaSelectedGOs")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("save2Hist_temporaSelectedGOs", id = "save2Hist_temporaSelectedGOs", duration = NULL)
  }
  if (is.null(clicked)) return()
  if (clicked < 1) return()
  
  add2history(type = "save", input = isolate( reactiveValuesToList(input)), 
              comment = paste0("# temporaSelectedGOs\n",
                               "require(Tempora)\n",
                               "fun = plotData$plotData$plotFunc\n", 
                               "environment(fun) = environment()\n",
                               "print(do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)]))\n"
              ),
              plotData = .schnappsEnv[["TRAJ_temporaSelectedGOs"]])
  
})
# save2Hist_tempora_screeplot
observe(label = "save2Hist_tempora_screeplot", {
  clicked  = input$save2Hist_tempora_screeplot
  if (DEBUG) cat(file = stderr(), "observe input$save2Hist_tempora_screeplot \n")
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "save2Hist_tempora_screeplot")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("save2Hist_tempora_screeplot", id = "save2Hist_tempora_screeplot", duration = NULL)
  }
  if (is.null(clicked)) return()
  if (clicked < 1) return()
  
  add2history(type = "save", input = isolate( reactiveValuesToList(input)), 
              comment = paste0("# tempora_screeplot\n",
                               "require(Tempora)\n",
                               "print(do.call(\"screeplot\",plotData$plotData[2:length(plotData$plotData)]))\n"
              ),
              plotData = .schnappsEnv[["TRAJ_tempora_screeplot"]])
  
})
# save2Hist_tempora_plot
observe(label = "save2Hist_tempora_plot", {
  clicked  = input$save2Hist_tempora_plot
  if (DEBUG) cat(file = stderr(), "observe input$save2Hist_tempora_plot \n")
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "save2Hist_tempora_plot")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("save2Hist_tempora_plot", id = "save2Hist_tempora_plot", duration = NULL)
  }
  if (is.null(clicked)) return()
  if (clicked < 1) return()
  
  add2history(type = "save", input = isolate( reactiveValuesToList(input)), 
              comment = paste0("# save2Hist_tempora_plot\n",
                               "require(Tempora)\n",
                               "require(igraph)\n",
                               "fun = plotData$plotData$plotFunc\n", 
                               "environment(fun) = environment()\n",
                               "tmp = do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)])\n"
              ),
              plotData = .schnappsEnv[["TRAJ_tempora_plot"]])
  
})
