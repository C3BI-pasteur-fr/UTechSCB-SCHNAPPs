require(ElPiGraph.R)
require(Tempora)
require(S4Vectors)
require(SingleCellExperiment)

# retVal <- drawTrajectoryHeatmap(x=expr_sel, time=traj$time, progression_group=projections[, dimCol], modules,
#                                 filename = normalizePath(outfile, mustWork = FALSE)
# )

# drawTrajectoryHeatmap ----
drawTrajectoryHeatmap <- function(x, time, progression_group = NULL, modules = NULL,
                                  show_labels_row = FALSE, show_labels_col = FALSE, scale_features = TRUE,
                                  ...) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    cat(sQuote("x"), " must be a numeric matrix or data frame",file = stderr())
    return(NULL)
  }
  if (!is.vector(time) || !is.numeric(time)) {
    cat(sQuote("time"), " must be a numeric vector",file = stderr())
    return(NULL)
  }
  if (nrow(x) != length(time)) {
    cat(
      sQuote("time"), " must have one value for each row in ",
      sQuote("x"),file = stderr()
    )
    return(NULL)
  }
  if ((!is.null(progression_group) && !is.vector(progression_group) &&
       !is.factor(progression_group)) || (!is.null(progression_group) &&
                                          length(progression_group) != nrow(x))) {
    cat(sQuote("progression_group"), " must be a vector or a factor of length nrow(x)",file = stderr())
    return(NULL)
  }
  if (is.null(rownames(x))) {
    rownames(x) <- paste("Row ", seq_len(nrow(x)))
  }
  col_ann <- data.frame(row.names = rownames(x), Time = time)
  x_part <- x[order(time), , drop = FALSE]
  if (scale_features) {
    x_part <- scale_quantile(x_part)
  }
  x_part <- t(x_part)
  gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
  }
  if (exists("annotation_colors")) {
    ann_col <- annotation_colors
  } else {
    ann_col <- list(Time = RColorBrewer::brewer.pal(5, "RdGy"))
  }
  
  if (!is.null(progression_group)) {
    if (!is.factor(progression_group)) {
      progression_group <- factor(progression_group)
    }
    col_ann$Progression <- progression_group
    num_progressions <- length(levels(progression_group))
    progression_cols <- if (num_progressions <= 9) {
      RColorBrewer::brewer.pal(num_progressions, "Set1")
    }
    else {
      gg_color_hue(num_progressions)
    }
    ann_col$Progression <- stats::setNames(
      progression_cols,
      levels(progression_group)
    )
  }
  labels_row <- if (!show_labels_row) {
    rep("", nrow(x_part))
  } else {
    NULL
  }
  labels_col <- if (!show_labels_col) {
    rep("", ncol(x_part))
  } else {
    NULL
  }
  if (!is.null(modules)) {
    x_part <- x_part[modules$symbol, ]
    gaps_row <- which(modules$module[-1] != modules$module[-length(modules$module)])
    cluster_rows <- F
  }
  else {
    gaps_row <- NULL
    cluster_rows <- T
  }
  # return(list(data = x_part, cluster_cols = F, cluster_rows = cluster_rows,
  #        annotation_col = col_ann, annotation_colors = ann_col,
  #        gaps_row = gaps_row, labels_row = labels_row, labels_col = labels_col,
  #        ...))
  # pheatmap::pheatmap(x_part, cluster_cols = F, cluster_rows = cluster_rows,
  #                    annotation_col = col_ann, annotation_colors = ann_col,
  #                    gaps_row = gaps_row, labels_row = labels_row, labels_col = labels_col,
  #                    ...)
  list(
    mat = x_part, cluster_cols = F, cluster_rows = cluster_rows,
    annotation_col = col_ann, annotation_colors = ann_col,
    gaps_row = gaps_row, labels_row = labels_row, labels_col = labels_col
  )
}

# scorpiuseParameters ----
# herewith we control that scorpius is only run if the button is pressed.
# scorpiuseParameters <- reactiveValues()
Scorpius_dataInput <- callModule(
  cellSelectionModule,
  "Scorpius_dataInput"
)


Scorpius_scEx_log <- reactive({
  scEx_log = scEx_log()
  if(is.null(scEx_log)) return(NULL)
  # browser()
  selectedCells <- isolate(Scorpius_dataInput()) #DE_Exp_dataInput
  if(is.null(selectedCells)) return(NULL)
  cellNs <- selectedCells$cellNames()
  if(length(cellNs)<1) return(NULL)
  scEx_log = scEx_log[,cellNs]
})


scorpius_projections <- reactive({
  projections <- projections()
  if (is.null(projections)) return(NULL)
  selectedCells <- isolate(Scorpius_dataInput()) #DE_Exp_dataInput
  if(is.null(selectedCells)) return(NULL)
  cellNs <- selectedCells$cellNames()
  if(length(cellNs)<1) return(NULL)
  projections[cellNs,]
})
# scorpiusInput ----
# read input space from file
scorpiusInput <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusInput started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusInput")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusInput")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusInput", id = "scorpiusInput", duration = NULL)
  }
  
  inFile <- input$trajInputFile
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusInput.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/scorpiusInput.RData")
  if (is.null(inFile)) {
    return(NULL)
  }
  if (!file.exists(inFile$datapath)) {
    return(NULL)
  }
  traj <- read.csv(file = inFile$datapath)
  if (colnames(traj) == c("path.Comp1", "path.Comp2", "time")) {
    return(traj)
  } else {
    warning("file not correct")
    return(NULL)
  }
})

# scorpiusSpace ----
# 2D space in which trajectory is calculated
# return projections or loaded space.
scorpiusSpace <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusSpace started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusSpace")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusSpace")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusSpace", id = "scorpiusSpace", duration = NULL)
  }
  
  
  doCalc <- input$updatetScorpiusParameters
  dimX <- isolate(input$dimScorpiusX)
  dimY <- isolate(input$dimScorpiusY)
  scInput <- isolate(scorpiusInput())
  projections <- scorpius_projections()
  
  if (!is.null(scInput)) {
    return(scInput[, c(1, 2)])
  }
  if ( is.null(projections)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusSpace:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusSpace.RData", list = c(ls()))
  }
  # cp=load(file="~/SCHNAPPsDebug/scorpiusSpace.RData")
  
  space <- projections[, c(dimX, dimY)]
  
  return(space)
})


# scorpiusTrajectory ----
# calculate trajectory
scorpiusTrajectory <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusTrajectory started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusTrajectory")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusTrajectory")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusTrajectory", id = "scorpiusTrajectory", duration = NULL)
  }
  
  # only execute if button is pressed
  clicked = input$updatetScorpiusParameters
  cat(file = stderr(), paste("scorpiusTrajectory clicked:", clicked,"\n"))
  # browser()
  # create dependancy on button and return if not pressed once
  if (input$updatetScorpiusParameters == 0) {
    # when coming from history
    if(!is.null(.schnappsEnv$react.scorpiusTrajectory))
      return(.schnappsEnv$react.scorpiusTrajectory$path)
    return(NULL)
  }
  
  space <- scorpiusSpace()
  # doCalc <- input$scorpiusCalc
  scInput <- scorpiusInput()
  
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusTrajectory.RData", list = c(ls()))
  }
  # cp=load(file="~/SCHNAPPsDebug/scorpiusTrajectory.RData")
  if (!is.null(scInput)) {
    return(scInput)
  }
  if (is.null(space)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusTrajectory:NULL\n"))
    return(NULL)
  }
  if (nrow(space)<10 | ncol(space)<2) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusTrajectory:NULL; need more samples/columns\n"))
    return(NULL)
  }
  require(SCORPIUSbj)
  dig = digest(space, algo = "sha256")
  if(!is.null(.schnappsEnv$react.scorpiusTrajectory) & length(.schnappsEnv$react.scorpiusTrajectory)==2){
    if (dig == .schnappsEnv$react.scorpiusTrajectory[[1]]){
      return(.schnappsEnv$react.scorpiusTrajectory$path)
    }
  }
  
  traj <- SCORPIUSbj::infer_trajectory(space,thresh = 0.00001)
  traj$path = data.frame(traj$path)
  traj$path$idx <- 1:nrow(traj$path)
  traj$path <- traj$path[order(traj$time),] 
  rownames(traj$path) = names(sort(traj$time))
  traj$path$time = sort(traj$time)
  traj$path = traj$path[order(traj$path$idx),]
  traj$path = traj$path[, -3]
  # rownames(traj$path) == names(traj$time)
  
  .schnappsEnv$react.scorpiusTrajectory = list(digest = dig, path = traj$path)
  return(traj$path)
})

# traj$path[1:10,]
# traj$path[order(traj$path$idx)[1:10],]

# scorpiusExpSel ----
scorpiusExpSel <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusExpSel started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusExpSel")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusExpSel")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusExpSel", id = "scorpiusExpSel", duration = NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification( id = "scorpiusExpSelWARNING")
  }
  
  
  # create dependancy on button and return if not pressed once
  if (input$updatetScorpiusParameters == 0) {
    # when coming from history
    if(!is.null(.schnappsEnv$react.scorpiusExpSel))
      return(.schnappsEnv$react.scorpiusExpSel[[2]])
    return(NULL)
  }
  scEx_log <- isolate(Scorpius_scEx_log())
  traj <- isolate(scorpiusTrajectory())
  # doCalc <- input$scorpiusCalc
  scorpMaxGenes <- isolate(input$scorpMaxGenes)
  scorpRepeat <- isolate(input$scorpRepeat)
  
  if (is.null(traj) | is.null(scEx_log) | is.null(traj)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusExpSel:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusExpSel.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/scorpiusExpSel.RData")
  # cellsNotFound <- colnames(assays(scEx_log)[[1]])[!colnames(assays(scEx_log)[[1]]) %in% rownames(traj)]
  dig = digest(list(assays(scEx_log)[[1]][,rownames(traj)], traj$time, scorpRepeat), algo = "sha256")
  if(!is.null(.schnappsEnv$react.scorpiusExpSel))
    if(dig == .schnappsEnv$react.scorpiusExpSel[[1]]) {
      return(.schnappsEnv$react.scorpiusExpSel[[2]])
    }
  expression <- t(as.matrix(assays(scEx_log)[[1]][,rownames(traj)]))
  gimp <- gene_importances(expression[rownames(traj),], traj$time, num_permutations = scorpRepeat, 
                           num_threads = detectCores())
  maxRow <- min(scorpMaxGenes, nrow(gimp))
  gene_sel <- gimp[1:maxRow, ]
  expr_sel <- expression[, gene_sel$gene]
  # dfTmp = data.frame(matrix(0,nrow = length(cellsNotFound), ncol = ncol(expr_sel)))
  # rownames(dfTmp) = cellsNotFound
  # colnames(dfTmp) = colnames(expr_sel)
  # retVal = rbind(expr_sel,dfTmp)
  .schnappsEnv$react.scorpiusExpSel = list(dig, list(expr_sel = expr_sel, gene_sel = gene_sel))
  return(list(expr_sel = expr_sel, gene_sel = gene_sel))
})

# scorpiusModules ----
scorpiusModules <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusModules started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusModules")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusModules")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusModules", id = "scorpiusModules", duration = NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification( id = "scorpiusModulesWARNING")
  }
  # browser()
  
  
  # browser()
  scEx_log <- isolate(Scorpius_scEx_log())
  # projections = projections()
  # space <- scorpiusSpace()
  traj <- scorpiusTrajectory()
  expr_sel <- scorpiusExpSel()
  
  if (is.null(traj) | is.null(scEx_log)| is.null(expr_sel)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusModules:NULL\n"))
    return(NULL)
  }
  
  dig = digest(list(expr_sel$expr_sel, traj$time), algo = "sha256") 
  
  # create dependancy on button and return if not pressed once
  if (input$updatetScorpiusParameters == 0) {
    # when coming from history
    if(!is.null(.schnappsEnv$react.scorpiusModules))
      return(.schnappsEnv$react.scorpiusModules[[2]])
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusModules.RData", list = c(ls()))
  }
  # cp =load(file="~/SCHNAPPsDebug/scorpiusModules.RData")
  
  # space = projections[,c(dimX, dimY)]
  # traj <- infer_trajectory(space)
  # expression = as.matrix(exprs(scEx_log))
  # gimp <- gene_importances(t(expression), traj$time, num_permutations = 0, num_threads = 8)
  # gene_sel <- gimp[1:50,]
  # expr_sel <- t(expression)[,gene_sel$gene]
  
  # modules <- extract_modules(scale_quantile(expr_sel$expr_sel), traj$time, verbose = T)
  if(!is.null(.schnappsEnv$react.scorpiusModules))
    if(dig == .schnappsEnv$react.scorpiusModules[[1]]) {
      return(.schnappsEnv$react.scorpiusModules[[2]])
    }
  modules = tryCatch ({extract_modules(scale_quantile(expr_sel$expr_sel), traj$time, verbose = T)},
                      error = function(e) {
                        cat(file = stderr(), paste("ERROR: extract_modules",e))
                        save(file = "~/SCHNAPPsDebug/scorpiusModulesError.RData", list = c(ls(),"expr_sel", "traj"))
                        if (!is.null(getDefaultReactiveDomain())) {
                          showNotification("GSEABase::getGmt not a valid file", 
                                           id = "scorpius.ERROR", type = "error", duration = NULL)
                        }
                        return(NULL)
                      })
  if(is.null(modules)) {
    return(NULL)
  }
  
  
  modules <- as.data.frame(modules)
  fd <- rowData(scEx_log)
  modules$symbol <- fd[modules$feature, "symbol"]
  rownames(modules) <- make.unique(as.character(modules$symbol, sep = "___"))
  
  .schnappsEnv$react.scorpiusModules = list(digest = dig, modules = modules)
  
  return(modules)
})

# scorpiusModulesTable ----
scorpiusModulesTable <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusModulesTable started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusModulesTable")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusModulesTable")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusModulesTable", id = "scorpiusModulesTable", duration = NULL)
  }
  
  
  # create dependancy on button and return if not pressed once
  if (input$updatetScorpiusParameters == 0) {
    # when coming from history
    if(!is.null(.schnappsEnv$react.scorpiusModulesTable))
      return(.schnappsEnv$react.scorpiusModulesTable[[2]])
    return(NULL)
  }
  
  scEx_log <- isolate(Scorpius_scEx_log())
  # projections = projections()
  # space <- scorpiusSpace()
  traj <- isolate(scorpiusTrajectory())
  expr_sel <- isolate(scorpiusExpSel())
  modules <- isolate(scorpiusModules())
  
  # scorpiusModules = scorpiusModules()
  # upI <- updateScorpiusInput() # needed to update input
  dimX <- isolate(input$dimScorpiusX)
  dimY <- isolate(input$dimScorpiusY)
  # dimCol = input$dimScorpiusCol
  # doCalc <- input$scorpiusCalc
  dig = digest(list(scEx_log, traj, expr_sel, modules, dimX, dimY), algo = "sha256")
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusModulesTable.RData", list = c(ls()))
  }
  if (is.null(traj) | is.null(scEx_log)| is.null(expr_sel)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusModules:NULL\n"))
    return(NULL)
  }
  # cp= load(file="~/SCHNAPPsDebug/scorpiusModulesTable.RData")
  # space = projections[,c(dimX, dimY)]
  # traj <- infer_trajectory(space)
  # expression = as.matrix(exprs(scEx_log))
  # gimp <- gene_importances(t(expression), traj$time, num_permutations = 0, num_threads = 8)
  # gene_sel <- gimp[1:50,]
  # expr_sel <- t(expression)[,gene_sel$gene]
  
  gene_selDF <- as.data.frame(expr_sel$gene_sel)
  rownames(gene_selDF) = gene_selDF[,1]
  gene_selDF = gene_selDF[,-1]
  if(is.null(modules))
    retVal = gene_selDF[modules$feature,]
  else
    retVal = cbind(modules,gene_selDF[modules$feature,])
  .schnappsEnv$react.scorpiusModulesTable = list(dig = dig, retVal)
  return(retVal)
})



# --------------------------
# Elpi Graph
# --------------------------
# TODO change to projections input

Elpi_dataInput <- callModule(
  cellSelectionModule,
  "Elpi_dataInput"
)
Elpi_scEx_log <- reactive({
  scEx_log = scEx_log()
  if(is.null(scEx_log)) return(NULL)
  
  selectedCells <- isolate(Elpi_dataInput()) #DE_Exp_dataInput
  if(is.null(selectedCells)) return(NULL)
  # browser()
  cellNs <- selectedCells$cellNames()
  if(is.null(cellNs)) return(NULL)
  if(length(cellNs)<1) return(NULL)
  scEx_log = scEx_log[,cellNs]
  return(scEx_log)
})
Elpi_scEx <- reactive({
  scEx = scEx()
  if(is.null(scEx)) return(NULL)
  
  selectedCells <- isolate(Elpi_dataInput()) #DE_Exp_dataInput
  if(is.null(selectedCells)) return(NULL)
  cellNs <- selectedCells$cellNames()
  if(length(cellNs)<1) return(NULL)
  scEx = scEx[,cellNs]
  return(scEx)
})

Elpi_projections <- reactive({
  projections <- projections()
  if (is.null(projections)) return(NULL)
  selectedCells <- isolate(Elpi_dataInput()) #DE_Exp_dataInput
  if(is.null(selectedCells)) return(NULL)
  cellNs <- selectedCells$cellNames()
  if(length(cellNs)<1) return(NULL)
  projections[cellNs,]
})

# elpiTreeData ----
elpiTreeData <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "elpiTreeData\n")
  }
  
  clicked <- input$elpiCalc
  scEx <- Elpi_scEx()
  projections <- isolate(Elpi_projections())
  dimElpi <- isolate(input$dimElpi)
  dim1 <- input$dimElpiX
  dim2 <- input$dimElpiY
  if (is.null(scEx)) {
    cat(file = stderr(), "--- elpiTreeData: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpiTreeData.RData", list = c(base::ls()))
  }
  # cp =load(file = "~/SCHNAPPsDebug/elpiTreeData.RData")
  if (dimElpi == "elpiPCA") {
    return(t(as.matrix(assays(scEx)[[1]])))
  }
  if (!all(c(dim1, dim2) %in% colnames(projections))){
    cat(file = stderr(), paste("--- dims not in projctions: NULL:", c(dim1, dim2), "\n"))
    return(NULL)
    
  }
  if (dimElpi == "components") {
    if (!dim1 %in% colnames(projections) | !dim2 %in% colnames(projections)) {
      return(NULL)
    }
    return(as.matrix(projections[,c(dim1,dim2)]))
  }
  cat(file = stderr(), "elpiTreeData should not happen\n")
  t(assays(scEx)[[1]])
})

# traj_endpoints ----
traj_endpoints <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "traj_endpoints")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "traj_endpoints")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("traj_endpoints", id = "traj_endpoints", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "traj_endpoints started.\n")
  
  clicked <- input$elpiCalc
  TreeEPG <- elpiGraphCompute()
  scEx_log <- Elpi_scEx_log()
  
  projections <- isolate(Elpi_projections())
  elpimode <- isolate(input$ElpiMethod)
  seed <- isolate(input$elpiSeed)
  
  if (is.null(projections) || is.null(scEx_log) || is.null(TreeEPG) || elpimode=="computeElasticPrincipalCircle") {
    cat(file = stderr(), "--- traj_endpoints: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_endpoints.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/traj_endpoints.RData")
  set.seed(seed = seed)
  Tree_Graph <- ElPiGraph.R::ConstructGraph(TreeEPG[[length(TreeEPG)]])
  Tree_e2e <- ElPiGraph.R::GetSubGraph(Net = Tree_Graph, Structure = 'end2end', Circular = T,Nodes = 30,KeepEnds = T)
  # get all end-points:
  endPoints = unique(c(sapply(Tree_e2e, function(x) x[1]), sapply(Tree_e2e, function(x) x[length(x)])))
  return(endPoints)
})


# updateElpiInput ----
# update the start/end positions with the possible end point values
# traj_getPseudotime ----
traj_getPseudotime <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "traj_getPseudotime")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "traj_getPseudotime")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("traj_getPseudotime", id = "traj_getPseudotime", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "traj_getPseudotime started.\n")
  clicked <- input$elpiCalc
  scEx_log <- Elpi_scEx_log()
  # scEx_log_sha <- scEx_log_sha()
  TreeEPG <- elpiGraphCompute()
  # TreeEPG_sha <- TreeEPG_sha()
  
  tree_data <- elpiTreeData()
  tragetPath <- traj_tragetPath()
  seed <- isolate(input$elpiSeed)
  elpimode <- isolate(input$ElpiMethod)
  
  if (is.null(scEx_log) || is.null(TreeEPG) || elpimode=="computeElasticPrincipalCircle" || length(tragetPath) == 0) {
    cat(file = stderr(), "--- traj_getPseudotime: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_getPseudotime.RData", list = c(ls()))
  }
  # cp =load(file="~/SCHNAPPsDebug/good/traj_getPseudotime.RData")
  
  set.seed(seed)
  # computing a Partition structure
  PartStruct <- ElPiGraph.R::PartitionData(X = tree_data, NodePositions = TreeEPG[[length(TreeEPG)]]$NodePositions)
  
  # projection structure
  ProjStruct <- ElPiGraph.R::project_point_onto_graph(X = tree_data,
                                                      NodePositions = TreeEPG[[length(TreeEPG)]]$NodePositions,
                                                      Edges = TreeEPG[[length(TreeEPG)]]$Edges$Edges,
                                                      Partition = PartStruct$Partition)
  psTime = ElPiGraph.R::getPseudotime(ProjStruct = ProjStruct, NodeSeq = names(tragetPath[[1]]))
  
  # writeShaCache(moduleName = "traj_getPseudotime", 
  #               moduleParameters = list(scEx_log_sha, TreeEPG,
  #                                       tree_data, targetPathSha,
  #                                       elpimode, seed),
  #               retVal = psTime,
  #               status = "finished",
  #               message = "")
  return(psTime)
  
})

# traj_targetPathSha <- reactive({
#   start.time <- base::Sys.time()
#   on.exit({
#     printTimeEnd(start.time, "traj_targetPathSha")
#     if (!is.null(getDefaultReactiveDomain()))
#       removeNotification(id = "traj_targetPathSha")
#   })
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("traj_targetPathSha", id = "traj_targetPathSha", duration = NULL)
#   }
#   if (DEBUG) cat(file = stderr(), "traj_getPseudotime started.\n")
#   # browser()
#   targetPath <- traj_tragetPath()
#   if(is.null(targetPath) | length(targetPath)==0)
#     return("")
#   sha1(as.character(targetPath[[1]]))
# })

## traj_elpi_modules -----
traj_elpi_modules <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "traj_elpi_modules")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "traj_elpi_modules")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("traj_elpi_modules", id = "traj_elpi_modules", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "traj_elpi_modules started.\n")
  
  clicked <- input$elpiCalc
  scEx_log <- Elpi_scEx_log()
  TreeEPG <- elpiGraphCompute()
  gene_sel <- traj_elpi_gimp()
  
  elpimode <- isolate(input$ElpiMethod)
  seed <- isolate(input$elpiSeed)
  
  # add a button to check result.
  # This button will invalidate this reactive and restart checking of cache.
  # input$elpi_modules_check
  
  if (clicked == 0) {
    # when coming from history
    if(!is.null(.schnappsEnv$react.traj_elpi_modules) & length(.schnappsEnv$react.traj_elpi_modules)==2)
      return(.schnappsEnv$react.traj_elpi_modules[[2]])
    return(NULL)
  }
  
  
  if (is.null(scEx_log) | is.null(gene_sel) | elpimode=="computeElasticPrincipalCircle" ) {
    cat(file = stderr(), "--- traj_elpi_modules: NULL\n")
    return(NULL)
  }
  if (nrow(gene_sel$gene_sel)<1) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_elpi_modules.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/traj_elpi_modules.RData")
  
  gene_sel = gene_sel$gene_sel
  set.seed(seed)
  expr_sel <- t(as.matrix(assays(scEx_log)[[1]][gene_sel$gene,]))
  
  ## Group the genes into modules and visualise the modules in a heatmap
  # group_name should be dbCluster or other selectable option
  dig = digest(list(expr_sel), algo = "sha256")
  
  if(!is.null(.schnappsEnv$react.traj_elpi_modules) & length(.schnappsEnv$react.traj_elpi_modules)==2){
    if (dig == .schnappsEnv$react.traj_elpi_modules[[1]]){
      return(.schnappsEnv$react.traj_elpi_modules[[2]])
    }
  }
  
  modules <- SCORPIUSbj::extract_modules(SCORPIUS::scale_quantile(expr_sel))
  modules <- as.data.frame(modules)
  fd <- rowData(scEx_log)
  modules$symbol <- fd[modules$feature, "symbol"]
  rownames(modules) <- make.unique(as.character(modules$symbol, sep = "___"))
  
  .schnappsEnv$react.traj_elpi_modules = list(digest = dig, retVal = modules)
  
  return(modules)
})


# traj_elpi_gimp -----
traj_elpi_gimp <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "traj_elpi_gimp")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "traj_elpi_gimp")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("traj_elpi_gimp", id = "traj_elpi_gimp", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "traj_elpi_gimp started.\n")
  scEx_log <- Elpi_scEx_log()
  projections <- Elpi_projections()
  TreeEPG <- elpiGraphCompute()
  psTime = traj_getPseudotime()
  clicked <- input$elpiCalc
  
  elpimode <- isolate(input$ElpiMethod)
  num_permutations <- input$elpi_num_permutations
  ntree <- input$elpi_ntree
  ntree_perm <- input$elpi_ntree_perm
  nGenes <- input$elpi_nGenes
  seed <- isolate(input$elpiSeed)
  
  # create dependancy on button and return if not pressed once
  if (clicked == 0) {
    # when coming from history
    if(!is.null(.schnappsEnv$react.traj_elpi_gimp) & length(.schnappsEnv$traj_elpi_gimp)==2)
      return(.schnappsEnv$react.traj_elpi_gimp[[2]])
    return(NULL)
  }
  
  
  if (is.null(scEx_log) | is.null(psTime) | elpimode=="computeElasticPrincipalCircle") {
    cat(file = stderr(), "--- traj_elpi_gimp: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_elpi_gimp.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/traj_elpi_gimp.RData")
  
  set.seed(seed)
  require(SCORPIUSbj)
  logCounts <- as.matrix(assays(scEx_log)[[1]][,which(!is.na(psTime$Pt))])
  pst = psTime$Pt[which(!is.na(psTime$Pt))]
  dig = digest(list(logCounts, pst, num_permutations, ntree, ntree_perm), algo = "sha256")
  if(!is.null(.schnappsEnv$react.traj_elpi_gimp) & length(.schnappsEnv$react.traj_elpi_gimp)==2){
    if (dig == .schnappsEnv$react.traj_elpi_gimp[[1]]){
      return(.schnappsEnv$react.traj_elpi_gimp[[2]])
    }
  }
  geneImport <- SCORPIUSbj::gene_importances(t(logCounts), pst, num_permutations = num_permutations, ntree = ntree,
                                             ntree_perm = ntree_perm, mtry = ncol(logCounts) * 0.01, num_threads = detectCores()-1)
  gene_sel <- geneImport[1:nGenes,]
  expr_sel <- t(logCounts)[, gene_sel$gene]
  
  .schnappsEnv$react.traj_elpi_gimp = list(digest = dig, retVal = list(expr_sel = expr_sel, gene_sel = gene_sel))
  
  return(list(expr_sel = expr_sel, gene_sel = gene_sel))
})
# 
# # observeProj ----
# # update projections
# observe({
#   start.time <- base::Sys.time()
#   on.exit({
#     printTimeEnd(start.time, "observeProj")
#     if (!is.null(getDefaultReactiveDomain()))
#       removeNotification(id = "observeProj")
#   })
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("observeProj", id = "observeProj", duration = NULL)
#   }
#   if (DEBUG) cat(file = stderr(), "observeProj started.\n")
#   
#   startNode <- input$elpiStartNode
#   endNode <- input$elpiEndNode
#   elpimode <- input$ElpiMethod
#   psTime = traj_getPseudotime()
#   scEx_log <- scEx_log()
#   isolate({
#     prjs <- sessionProjections$prjs
#   })
#   
#   if (is.null(scEx_log) || is.null(psTime) || elpimode=="computeElasticPrincipalCircle") {
#     return(NULL)
#   }
#   if (.schnappsEnv$DEBUGSAVE) {
#     save(file = "~/SCHNAPPsDebug/observeProj.RData", list = c(ls()))
#   }
#   # load(file="~/SCHNAPPsDebug/observeProj.RData")
#   cn = paste0("traj_", startNode, "_", endNode)
#   if (cn %in% colnames(prjs)) {
#     return(NULL)
#   }
#   # browser()
#   if (ncol(prjs) > 0) {
#     # make sure we are working with the correct cells. This might change when cells were removed.
#     prjs = prjs[colnames(scEx_log),,drop=FALSE]
#     # didn't find a way to easily overwrite columns
#     
#     if (cn %in% colnames(prjs)) {
#       prjs[, cn] <- psTime$Pt
#     } else {
#       prjs <- base::cbind(prjs, psTime$Pt, deparse.level = 0)
#       colnames(prjs)[ncol(prjs)] <- cn
#     }
#     sessionProjections$prjs <- prjs
#   } else {
# 
#         prjs <- data.frame(cn = psTime$Pt )
#     rownames(prjs) = colnames(scEx_log)
#     colnames(prjs)[ncol(prjs)] = cn
#     sessionProjections$prjs = prjs
#   }
# })

#traj_tragetPath ----
traj_tragetPath <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "traj_tragetPath")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "traj_tragetPath")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("traj_tragetPath", id = "traj_tragetPath", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "traj_tragetPath started.\n")
  scEx_log <- Elpi_scEx_log()
  projections <- Elpi_projections()
  TreeEPG <- elpiGraphCompute()
  clicked <- input$elpiCalc
  startNode <- input$elpiStartNode
  endNode <- input$elpiEndNode
  
  elpimode <- isolate(input$ElpiMethod)
  seed <- isolate(input$elpiSeed)
  
  if (is.null(scEx_log) || is.null(TreeEPG) || elpimode=="computeElasticPrincipalCircle") {
    cat(file = stderr(), "--- traj_tragetPath: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_tragetPath.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/traj_tragetPath.RData")
  set.seed(seed)
  Tree_Graph <- ElPiGraph.R::ConstructGraph(TreeEPG[[length(TreeEPG)]])
  Tree_e2e <- ElPiGraph.R::GetSubGraph(Net = Tree_Graph, Structure = 'end2end')
  
  retVal <- Tree_e2e[sapply(Tree_e2e, function(x){any(x[1] %in% c(startNode,endNode)) & any(x[length(x)] %in% c(startNode,endNode)) })]
  
})

# elpiGraphCompute ----
elpiGraphCompute <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpiGraphCompute")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "elpiGraphCompute")
  })
  if (DEBUG) {
    cat(file = stderr(), "elpiGraphCompute\n")
  }
  
  clicked <- input$elpiCalc
  tree_data <- isolate(elpiTreeData())
  
  nReps <- isolate(input$elpinReps) # 1-50
  NumNodes <- isolate(input$elpiNumNodes) # 10 - 100
  ProbPoint <- isolate(input$elpiProbPoint) # 0.1-1.0
  method <- isolate(input$ElpiMethod)
  dimUse <-  isolate(input$dimElpi)
  seed <- isolate(input$elpiSeed)
  
  require(parallel)
  
  if (is.null(tree_data)) {
    cat(file = stderr(), "--- elpiGraphCompute: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpiCalc.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load("~/SCHNAPPsDebug/elpiCalc.RData")
  
  set.seed(seed)
  if (dimUse == "elpiPCA") {
    elpiDoPCA = TRUE
    elipCenter = TRUE
  } else {
    elpiDoPCA = FALSE
    elipCenter = FALSE
    
  }
  cep <- do.call(method, list(
    X = tree_data,
    NumNodes = NumNodes,
    drawAccuracyComplexity = F,
    nReps = nReps, # bootstrapping
    ProbPoint = ProbPoint, # bootstrapping
    drawPCAView = F,
    drawEnergy = F,
    Do_PCA = elpiDoPCA,
    CenterData = elipCenter, 
    n.cores = detectCores() - 1
  ))
  
  setRedGreenButton(
    vars = list(
      c("elpinReps", isolate(input$elpinReps)),
      c("elpiNumNodes", isolate(input$elpiNumNodes)),
      c("elpiProbPoint", isolate(input$elpiProbPoint)),
      c("ElpiMethod", isolate(input$ElpiMethod)),
      c("dimElpi", isolate(input$dimElpi)),
      c("elpiSeed", isolate(input$elpiSeed)),
      c("dimElpi", isolate(input$dimElpi)),
      c("dimElpiX", isolate(input$dimElpiX)),
      c("dimElpiY", isolate(input$dimElpiY))
    ),
    button = "elpiCalc"
  )
  
  .schnappsEnv$defaultValues[["dimElpi"]] <- input$dimElpi
  .schnappsEnv$defaultValues[["dimElpiX"]] <- input$dimElpiX
  .schnappsEnv$defaultValues[["dimElpiY"]] <- input$dimElpiY
  .schnappsEnv$defaultValues[["dimElpiCol"]] <- input$dimElpiCol
  .schnappsEnv$defaultValues[["elpiSeed"]] <- input$elpiSeed
  .schnappsEnv$defaultValues[["ElpiMethod"]] <- input$ElpiMethod
  .schnappsEnv$defaultValues[["elpinReps"]] <- input$elpinReps
  .schnappsEnv$defaultValues[["elpiNumNodes"]] <- input$elpiNumNodes
  .schnappsEnv$defaultValues[["elpiProbPoint"]] <- input$elpiProbPoint
  .schnappsEnv$defaultValues[["elpiEndNode"]] <- input$elpiEndNode
  .schnappsEnv$defaultValues[["elpiStartNode"]] <- input$elpiStartNode
  .schnappsEnv$defaultValues[["elpi_num_permutations"]] <- input$elpi_num_permutations
  .schnappsEnv$defaultValues[["elpi_ntree"]] <- input$elpi_ntree
  .schnappsEnv$defaultValues[["elpi_ntree_perm"]] <- input$elpi_ntree_perm
  .schnappsEnv$defaultValues[["elpi_nGenes"]] <- input$elpi_nGenes
  
  return(cep)
})

# elpiGraphConstruct ----
elpiGraphConstruct <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "elpiGraphConstruct\n")
  }
  clicked <- input$elpiCalc
  cep <- elpiGraphCompute()
  tree_data <- elpiTreeData()
  
  seed <- isolate(input$elpiSeed)
  
  if (is.null(cep)) {
    cat(file = stderr(), "--- elpiGraphConstruct: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpiConstruct.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/elpiConstruct.RData")
  
  set.seed(seed = seed)
  Tree_Graph <- ElPiGraph.R::ConstructGraph(PrintGraph = cep[[length(cep)]])
  Tree_Brches <- ElPiGraph.R::GetSubGraph(Net = Tree_Graph, Structure = "branches")
  PartStruct <- ElPiGraph.R::PartitionData(X = tree_data, NodePositions = cep[[length(cep)]]$NodePositions)
  
  list(
    Tree_Graph = Tree_Graph,
    Tree_e2e = GetSubGraph(Net = Tree_Graph, Structure = "end2end"),
    Tree_Brches = Tree_Brches,
    Tree_BrBrPt = GetSubGraph(Net = Tree_Graph, Structure = "branches&bpoints"),
    PartStruct = PartStruct,
    PtInBr = lapply(Tree_Brches, function(x) {
      which(PartStruct$Partition %in% x)
    })
  )
})

# elpiPointLabel ----
elpiPointLabel <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "elpiPointLabel\n")
  }
  
  clicked <- input$elpiCalc
  elpiGraphConstruct <- elpiGraphConstruct()
  if (is.null(elpiGraphConstruct)) {
    cat(file = stderr(), "--- elpiPointLabel: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpiPointLabel.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/elpiPointLabel.RData")
  PartStruct <- elpiGraphConstruct$PartStruct
  Tree_BrBrPt <- elpiGraphConstruct$Tree_BrBrPt
  
  
  PointLabel <- rep("", length(PartStruct$Partition))
  
  for (i in 1:length(Tree_BrBrPt)) {
    PointLabel[PartStruct$Partition %in% Tree_BrBrPt[[i]]] <- names(Tree_BrBrPt)[i]
  }
  return(PointLabel)
})

# elpiModulesTable ----
elpiModulesTable <- reactive({
  if (DEBUG) cat(file = stderr(), "elpiModulesTable started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpiModulesTable")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "elpiModulesTable")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("elpiModulesTable", id = "elpiModulesTable", duration = NULL)
  }
  clicked <- input$elpiCalc
  
  scEx_log <- isolate(Elpi_scEx_log())
  traj <- traj_getPseudotime()
  expr_sel <- traj_elpi_gimp()
  modules <- traj_elpi_modules()
  
  dimX <- isolate(input$dimElpiX)
  dimY <- isolate(input$dimElpiY)
  
  if (is.null(scEx_log) | is.null(expr_sel)) {
    if (DEBUG) cat(file = stderr(), paste("elpiModulesTable:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/elpiModulesTable.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/elpiModulesTable.RData")
  
  gene_selDF <- as.data.frame(expr_sel$gene_sel)
  rownames(gene_selDF) = gene_selDF[,1]
  gene_selDF = gene_selDF[,-1]
  return(cbind(modules,gene_selDF[modules$feature,]))
})


# TEMPORA ----------
# temporaImport ----
Tempora_dataInput <- callModule(
  cellSelectionModule,
  "Tempora_dataInput"
)
Tempora_scEx_log <- reactive({
  scEx_log = scEx_log()
  if(is.null(scEx_log)) return(NULL)
  
  selectedCells <- isolate(Tempora_dataInput()) #DE_Exp_dataInput
  cellNs <- selectedCells$cellNames()
  if(length(cellNs)<1) return(NULL)
  scEx_log = scEx_log[,cellNs]
  return(scEx_log)
})
Tempora_projections <- reactive({
  projections <- projections()
  if (is.null(projections)) return(NULL)
  selectedCells <- isolate(Tempora_dataInput()) #DE_Exp_dataInput
  cellNs <- selectedCells$cellNames()
  if(length(cellNs)<1) return(NULL)
  projections[cellNs,]
})

temporaImport <- reactive({
  if (DEBUG) cat(file = stderr(), "temporaImport started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaImport")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "temporaImport")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaImport", id = "temporaImport", duration = NULL)
  }
  clicked <- input$updatetTemporaParameters
  
  scEx_log <- isolate(Tempora_scEx_log())
  projections <- isolate(Tempora_projections())
  tCluster <- isolate(input$temporaCluster) # clusters to use
  tFactor <- isolate(input$temporaFactor) # time variable
  tLevels <- isolate(input$temporaLevels) # time points
  # selectedCells <- isolate(Tempora_dataInput()) #DE_Exp_dataInput
  # cellNs <- selectedCells$cellNames()
  # if(length(cellNs)<1) return(NULL)
  
  if (is.null(scEx_log) | is.null(projections)) {
    if (DEBUG) cat(file = stderr(), paste("temporaImport:NULL\n"))
    return(NULL)
  }
  
  
  if (!all(c(tCluster,tFactor) %in% colnames(projections))){
    if (DEBUG) cat(file = stderr(), paste("tCluster,tFactor:NULL\n"))
    return(NULL)
  }
  # if (!all(tLevels %in% levels(projections[,tFactor]))){
  #   if (DEBUG) cat(file = stderr(), paste("tLevels:NULL\n"))
  #   return(NULL)
  # }
  if(length(tLevels)<1){
    if (DEBUG) cat(file = stderr(), paste("tLevels:empty\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaImport.RData", list = c(ls()),compress = F)
  }
  # cp = load(file="~/debug/temporaImport.RData")
  # scEx_log = scEx_log[,cellNs]
  colData(scEx_log) <- S4Vectors::DataFrame(projections[rownames(colData(scEx_log)),])
  # HACK
  # TODO
  # should be a variable, here we are forcing symbol but could be anything.
  rownames(scEx_log) = rowData(scEx_log)[rownames(scEx_log), "symbol"]
  # # tempora hack 
  # suppressMessages(
  #   setMethod("getMD","SingleCellExperiment",
  #             function(x) data.frame(SingleCellExperiment::colData(x)))
  # )
  tcl = table(colData(scEx_log)[, tCluster])
  for( tIdx in which(tcl <2) ) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(paste("temporaWarning: ", names(tcl)[[tIdx]]),
                       id = "temporaImport",
                       type = "warning", duration = NULL)
    }
    # warning that we are removing things
    scEx_log[,!colData(scEx_log)[,tCluster] == tIdx]
    # to remove empty levels
    colData(scEx_log)[,tCluster] = factor(colData(scEx_log)[,tCluster])
  }
  # We only work with the levels that are given. 
  scEx_log = scEx_log[,colData(scEx_log)[,tFactor] %in% tLevels]
  #in case there are levels where there are no data we have to update this
  colData(scEx_log)[,tFactor] = factor(colData(scEx_log)[,tFactor])
  colData(scEx_log)[,tCluster] = factor(colData(scEx_log)[,tCluster])
  temporaObj <- ImportSeuratObject(seuratobj = scEx_log, 
                                   assayType = "logcounts",
                                   clusters = tCluster,
                                   timepoints = tFactor,
                                   timepoint_order = tLevels,
                                   cluster_labels = levels(colData(scEx_log)[,tCluster])
  )
  
  return(temporaObj)
})


temporaGeneIds = reactive({
  gmt_path = input$temporaGMTFile
  gs1 = tryCatch ({GSEABase::getGmt(gmt_path$datapath)},
                  error = function(e) {
                    cat(file = stderr(), "GSEABase::getGmt")
                    if (!is.null(getDefaultReactiveDomain())) {
                      showNotification("GSEABase::getGmt not a valid file", 
                                       id = "temporaPWProfiles", type = "error", duration = NULL)
                    }
                    return(NULL)
                  })
  if(!is.null(gs1)) {
    return(GSEABase::geneIds(gs1))
  }
  return(NULL)
})

# temporaPWProfiles ----
temporaPWProfiles <- reactive({
  if (DEBUG) cat(file = stderr(), "temporaPWProfiles started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaPWProfiles")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "temporaPWProfiles")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaPWProfiles", id = "temporaPWProfiles", duration = NULL)
  }
  BPPARAM=bpparam()
  temporaObj <- temporaImport()
  gmt_path = isolate(input$temporaGMTFile)
  
  if(is.null(temporaObj) | is.null(gmt_path)) {
    if (DEBUG) cat(file = stderr(), paste("temporaPWProfiles:NULL\n"))
    return(NULL)
  }
  
  min.sz = isolate(input$temporaMinSz)
  max.sz = isolate(input$temporaMaxSz)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaPWProfiles.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/temporaPWProfiles.RData")
  
  if (!file.exists(gmt_path$datapath)) {
    if (DEBUG) cat(file = stderr(), paste("gmt_path:NULL\n"))
    return(NULL)
  }
  
  good = tryCatch ({GSEABase::getGmt(gmt_path$datapath)},
                   error = function(e) {
                     cat(file = stderr(), "GSEABase::getGmt")
                     if (!is.null(getDefaultReactiveDomain())) {
                       showNotification("GSEABase::getGmt not a valid file", 
                                        id = "temporaPWProfiles", type = "error", duration = NULL)
                     }
                     return(NULL)
                   })
  if(is.null(good)) {
    if (DEBUG) cat(file = stderr(), paste("GSEABase::getGmt:NULL\n"))
    return(NULL)
  }
  # DEBUG
  # gmt_path$datapath = "~/Rstudio/SCHNAPPsContributions/Mouse_GOBP_AllPathways_no_GO_iea_February_05_2021_symbol.gmt"
  temporaObj <- CalculatePWProfiles(temporaObj, 
                                    gmt_path = gmt_path$datapath,
                                    method="gsva", 
                                    min.sz = min.sz, 
                                    max.sz = max.sz, 
                                    
                                    parallel.sz = BiocParallel::bpnworkers(BPPARAM))
  # to = temporaObj
  
  return(temporaObj)
})

# temporaTrajectory ----
temporaTrajectory <- reactive({
  if (DEBUG) cat(file = stderr(), "temporaTrajectory started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaTrajectory")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "temporaTrajectory")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaTrajectory", id = "temporaTrajectory", duration = NULL)
  }
  
  temporaObj <- temporaPWProfiles()
  
  if(is.null(temporaObj)) {
    if (DEBUG) cat(file = stderr(), paste("temporaTrajectory:NULL\n"))
    return(NULL)
  }
  
  n_pcs = isolate(input$temporaNPCs)
  difference_threshold = isolate(input$temporaDiff_thresh)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaTrajectory.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/temporaTrajectory.RData")
  # temporaObj = to
  temporaObj = tryCatch ({
    BuildTrajectory(temporaObj, 
                    n_pcs = n_pcs, 
                    difference_threshold = difference_threshold)
  },
  error = function(e) {
    cat(file = stderr(), "tempora BuildTrajectory")
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("tempora BuildTrajectory", 
                       id = "temporaBuildTrajectory", type = "error", duration = NULL)
    }
    return(NULL)
  })
  # to = temporaObj
  
  .schnappsEnv$defaultValues[["temporaCluster"]] = isolate(input$temporaCluster)
  .schnappsEnv$defaultValues[["temporaFactor"]] = isolate(input$temporaFactor)
  .schnappsEnv$defaultValues[["temporaLevels"]] = isolate(input$temporaLevels)
  .schnappsEnv$defaultValues[["temporaMinSz"]] = isolate(input$temporaMinSz)
  .schnappsEnv$defaultValues[["temporaMaxSz"]] = isolate(input$temporaMaxSz)
  .schnappsEnv$defaultValues[["temporaNPCs"]] = isolate(input$temporaNPCs)
  .schnappsEnv$defaultValues[["temporaDiff_thresh"]] = isolate(input$temporaDiff_thresh)
  .schnappsEnv$defaultValues[["temporaPval_thresh"]] = isolate(input$temporaPval_thresh)
  
  
  return(temporaObj)
})



# temporaIdentifyVaryingPWs ----
temporaIdentifyVaryingPWs <- reactive({
  if (DEBUG) cat(file = stderr(), "temporaIdentifyVaryingPWs started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaIdentifyVaryingPWs")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "temporaIdentifyVaryingPWs")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaIdentifyVaryingPWs", id = "temporaIdentifyVaryingPWs", duration = NULL)
  }
  
  
  # create dependancy on button and return if not pressed once
  if (input$updatetTemporaParameters == 0) {
    return(NULL)
  }
  
  temporaObj <- temporaTrajectory()
  temporaPval_thresh <- input$temporaPval_thresh
  if (is.null(temporaObj) ) {
    if (DEBUG) cat(file = stderr(), paste("temporaIdentifyVaryingPWs:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaIdentifyVaryingPWs.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/temporaIdentifyVaryingPWs.RData")
  # object = to
  
  #Fit GAMs on pathway enrichment profile
  temporaObj <- IdentifyVaryingPWsParallel(object = temporaObj, pval_threshold = temporaPval_thresh)
  # temporaObj <- IdentifyVaryingPWs(object = temporaObj, pval_threshold = temporaPval_thresh)
  
  
  
  return(temporaObj)
})

# temporaPvalModulesTable ----

temporaPvalModulesTable <- reactive({
  if (DEBUG) cat(file = stderr(), "temporaPvalModulesTable started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaPvalModulesTable")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "temporaPvalModulesTable")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaPvalModulesTable", id = "temporaPvalModulesTable", duration = NULL)
  }
  
  
  # create dependancy on button and return if not pressed once
  if (input$updatetTemporaParameters == 0) {
    return(NULL)
  }
  
  temporaObj <- temporaIdentifyVaryingPWs()
  
  if (is.null(temporaObj) ) {
    if (DEBUG) cat(file = stderr(), paste("temporaPvalModulesTable:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaPvalModulesTable.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/temporaPvalModulesTable.RData")
  
  outTable = data.frame(goTerm = names(temporaObj@varying.pws))
  outTable$pValues = temporaObj@varying.pws
  
  return(outTable)
  
})


# tempora2DPlotFunc ----
tempora2DPlotFunc <- function(temporaObj, projections, dimX, dimY, dimCol) {
  require(ggnetwork)
  require(ggplot2)
  require(ggrepel)
  require(network)
  space <- projections[, c(dimX, dimY)]
  require(SCORPIUSbj)
  
  
  # temporaObj@cluster.metadata
  
  data = projections[,c(dimX, dimY, dimCol)]
  if(!all(unlist(lapply(1:2,  FUN = function(x) is.numeric(data[,x]))))) return (NULL)
  
  mean.points <- aggregate(data[, 1:2], list(data[,3]), mean)
  
  # construct a network
  nCl = length(levels(projections[,dimCol]))
  net = matrix(0, nCl, nCl) 
  rownames(net) = levels(projections[,dimCol])
  colnames(net) = levels(projections[,dimCol])
  temporaObj@trajectory$from = as.character( temporaObj@trajectory$from)
  temporaObj@trajectory$to = as.character( temporaObj@trajectory$to)
  apply(temporaObj@trajectory, 1, FUN = function(x){
    net[x[[1]],x[[2]]] <<- 1
    if (! x[[5]] == "unidirectional") {
      net[x[[2]],x[[1]]] <<- 1
    }
  } )
  # traj <- SCORPIUSbj::infer_trajectory(space)
  n = as.network(net)
  # 
  # gnn = ggnetwork(n, layout = as.matrix(mean.points[,2:3]))
  # # rownames(temporaObj@cluster.metadata) = temporaObj@cluster.metadata$Id
  # 
  # # gnn$vertex.names = temporaObj@cluster.metadata[gnn$vertex.names,"label"]
  # 
  # ggplot(gnn) + geom_edges(aes(x,y,xend = xend,yend=yend),color = "black") + 
  #   geom_nodetext(aes(x,y,label = vertex.names ),
  #                                     fontface = "bold")
  
  dat = ggnetwork(n, layout = as.matrix(mean.points[,2:3]))
  colnames(dat)
  colnames(data) <- c("x", "y", "vertex.names")
  data$xend = NA
  data$yend = NA
  scaleBJ = function(x,y) scale(x, center = min(y), scale = diff(range(y)))
  data$x = scaleBJ(data$x,mean.points[,2])
  data$y = scaleBJ(data$y,mean.points[,3])
  dat = rbind(dat,data)
  p2 = ggplot(dat) +
    geom_point(aes(x, y, col=vertex.names )) +
    geom_edges(data = subset(dat, !is.na(xend)), 
               aes(x, y, xend = xend, yend = yend), 
               color="black",
               alpha = 1, 
               arrow = arrow(length = unit(16, "pt"), type = "closed")) +
    geom_nodes(data = subset(dat, !is.na(xend)), 
               aes(x, y, col=vertex.names), 
               size = 5, color = "white") +
    geom_nodetext(data = subset(dat, !is.na(xend)), 
                  aes(x, y, label = vertex.names), 
                  fontface = "bold") +
    theme_blank() 
  
  p2
}
