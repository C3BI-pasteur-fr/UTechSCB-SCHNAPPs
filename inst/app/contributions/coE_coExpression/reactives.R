# coExpression reactives.R

# heatmapFunc ---------------------------------
# used by both selection and all to create input for heatmap module
coE_heatmapFunc <- function(featureData, scEx_matrix, projections, genesin, cells, sampCol, ccols) {
  if (DEBUG) cat(file = stderr(), "coE_heatmapFunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_heatmapFunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_heatmapFunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_heatmapFunc", id = "coE_heatmapFunc", duration = NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_heatmapFunc.RData", list = c(ls()))
  }
  # cp=load(file = "~/SCHNAPPsDebug/coE_heatmapFunc.RData")
  # browser()
  #  create parameters used for pheatmap module
  genesin <- geneName2Index(genesin, featureData)
  if (is.null(genesin) | is.null(cells)) {
    return(NULL)
  }
  if (length(genesin) < 2) {
    return(NULL)
  }
  if (is.list(cells)) {
    cells = unlist(cells)
  }
  expression <- scEx_matrix[genesin, cells]
  
  validate(need(
    is.na(sum(expression)) != TRUE,
    "Gene symbol incorrect or genes not expressed"
  ))
  
  projections <- projections[order(as.numeric(as.character(projections$dbCluster))), ]
  
  if (!("sampleNames" %in% colnames(projections))) {
    projections$sample <- 1
  }
  annotation <- data.frame(projections[cells, c("dbCluster", "sampleNames")])
  rownames(annotation) <- colnames(expression)
  colnames(annotation) <- c("Cluster", "sampleNames")
  
  # For high-res displays, this will be greater than 1
  # pixelratio <- session$clientData$pixelratio
  pixelratio <- 1
  if (is.null(pixelratio)) pixelratio <- 1
  # width <- session$clientData$output_plot_width
  # height <- session$clientData$output_plot_height
  width <- NULL
  height <- NULL
  if (is.null(width)) {
    width <- 96 * 7
  } # 7x7 inch output
  if (is.null(height)) {
    height <- 96 * 7
  }
  
  # nonZeroRows <- which(Matrix::rowSums(expression) > 0)
  
  annCols <- list(
    "sampleNames" = sampCol,
    "dbCluster" = ccols
  )
  
  retVal <- tryCatch({list(
    # mat = expression[nonZeroRows, order(annotation[, 1], annotation[, 2])],
    mat = expression[, order(annotation[, 1], annotation[, 2])],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    # scale = "row",
    fontsize_row = 14,
    # labels_col = colnames(expression),
    labels_row = featureData[rownames(expression), "symbol"],
    show_rownames = TRUE,
    annotation_col = annotation,
    show_colnames = FALSE,
    annotation_legend = TRUE,
    # breaks = seq(minBreak, maxBreak, by = stepBreak),
    # filename = 'test.png',
    # filename = normalizePath(outfile),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(
      n = 6, name =
        "RdBu"
    )))(6),
    annotation_colors = annCols
  )},
  error = function(x) {return(NULL)}
  )
  # TODO Error handling
  if (is.null(retVal) ) return(NULL)
  
  
  # print debugging information on the console
  printTimeEnd(start.time, "inputData")
  # for automated shiny testing using shinytest
  # exportTestValues(coE_heatmapFunc = {
  #   retVal
  # })
  
  # this is what is run in the module
  # do.call(TRONCO::pheatmap, retVal)
  return(retVal)
}


# coE_dotPlot_GeneSetsFunc ----
coE_dotPlot_GeneSets <- function(projections = projections,
                                 scEx_log = scEx_log,
                                 clusters = clusters,
                                 geneSets = geneSets,
                                 gmtData = gmtData,
                                 summarize = T,
                                 col = "RdBu",
                                 col.min = col.min,
                                 col.max = col.max,
                                 dot.min = dot.min,
                                 dot.scale = dot.scale,
                                 scale.by = scale.by
){
  # replace "_", with "."
  newRN = stringr::str_replace_all(rownames(scEx_log),"_",".")
  rownames(scEx_log) = newRN
  seurDat <- CreateSeuratObject(
    counts = assays(scEx_log)[[1]]
  )
  Idents(seurDat) <- as.factor(projections[,clusters])
  featureDat = rowData(scEx_log)
  if(length(geneSets) == 1){
    features = geneName2Index(paste(gmtData[[geneSets]]$genes,collapse = ", "), featureDat)
  } else {
    FUN = function(x){
      geneName2Index(paste(gmtData[[x]]$genes,collapse = ", "), featureDat)
    }
    features = lapply(geneSets, FUN = FUN)
    names(features) = geneSets
  }
  # handle duplicated gene names
  ulFeatures = unlist(features)
  FUN = function(x,g){
    # cat(file = stderr(), g)
    if(length(which(x==g)>0)) x = x[-which(x==g)]
    x
  }
  for(dupGene in ulFeatures[which(ulFeatures %>% duplicated())]){
    features =lapply(features,FUN = FUN,g=dupGene)
    features$common =c(features$common, dupGene)
  }
  p = DotPlot(seurDat, 
              assay="RNA", 
              features = features, 
              cols = col,
              col.min = col.min,
              col.max = col.max,
              dot.min = dot.min,
              dot.scale = dot.scale,
              idents = NULL,
              group.by = NULL,
              split.by = NULL,
              cluster.idents = FALSE,
              scale = TRUE,
              scale.by = scale.by,
              scale.min = NA,
              scale.max = NA
  )
  labFun <- function(breakval){featureDat[breakval,"symbol"]}
  ylabFun <- function(val){stringr::str_replace(val,"SeuratProject_","bla")}
  p= p + scale_x_discrete(labels = labFun) +
    scale_y_discrete(labels = ylabFun) +
    ylab(clusters) +
    theme(axis.text.x = element_text(angle = 25, vjust = 0.5),
          strip.text.x = element_text(angle=5))
  ggplotly(p)
  p
}

# coE_dotPlot_GeneSetsModuleScore ----
# same as coE_dotPlot_GeneSets, only that the X labels are already set inside the function and adds ModuleScore.
coE_dotPlot_GeneSetsModuleScore <- function(projections = projections,
                                            scEx_log = scEx_log,
                                            clusters = clusters,
                                            geneSets = geneSets,
                                            gmtData = gmtData,
                                            summarize = T,
                                            col = "RdBu",
                                            col.min = col.min,
                                            col.max = col.max,
                                            dot.min = dot.min,
                                            dot.scale = dot.scale,
                                            scale.by = scale.by
){
  require(Seurat)
  require(SingleCellExperiment)
  # save(file = "~/SCHNAPPsDebug/coE_dotPlot_GeneSetsModuleScoreF.RData", list = c(ls()))
  # cp = load("~/SCHNAPPsDebug/coE_dotPlot_GeneSetsModuleScoreF.RData")
  # replace "_", with "."
  newRN = stringr::str_replace_all(rownames(scEx_log),"_",".")
  rownames(scEx_log) = newRN
  seurDat <- CreateSeuratObject(
    counts = assays(scEx_log)[[1]]
  )
  Idents(seurDat) <- as.factor(projections[,clusters])
  featureDat = rowData(scEx_log)
  if(length(geneSets) == 1){
    features = geneName2Index(paste(gmtData[[geneSets[1]]]$genes,collapse = ", "), featureDat) %>% list()
    names(features) = gmtData[[geneSets[1]]]$name
  } else {
    FUN = function(x){
      geneName2Index(paste(gmtData[[x]]$genes,collapse = ", "), featureDat)
    }
    features = lapply(geneSets, FUN = FUN)
    names(features) = geneSets
  }
  # handle duplicated gene names
  ulFeatures = unlist(features)
  FUN = function(x,g){
    # cat(file = stderr(), g)
    if(length(which(x==g)>0)) x = x[-which(x==g)]
    x
  }
  for(dupGene in ulFeatures[which(ulFeatures %>% duplicated())]){
    features =lapply(features,FUN = FUN,g=dupGene)
    features$common =c(features$common, dupGene)
  }
  p = DotPlotwithModuleScore(seurDat, 
                             assay="RNA", 
                             features = features,
                             featureDat = featureDat,
                             cols = col,
                             col.min = col.min,
                             col.max = col.max,
                             dot.min = dot.min,
                             dot.scale = dot.scale,
                             idents = NULL,
                             clusters = clusters,
                             group.by = NULL,
                             split.by = NULL,
                             cluster.idents = FALSE,
                             scale = TRUE,
                             scale.by = scale.by,
                             scale.min = NA,
                             scale.max = NA
  )
  p
}


# coE_heatmapSelectedReactive ----
# reactive function for selected heatmap
coE_heatmapSelectedReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_heatmapSelectedReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_heatmapSelectedReactive")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_heatmapSelectedReactive")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_heatmapSelectedReactive", id = "coE_heatmapSelectedReactive", duration = NULL)
  }
  # deepDebug()
  scEx_log <- scEx_log()
  projections <- projections()
  clicked <- input$updateHeatMapSelectedParameters
  genesin <- isolate(input$coE_heatmapselected_geneids)
  sc <- isolate(coE_selctedCluster())
  # scCL <- sc$cluster # "1" "2" "3" "4"
  # scCL <- isolate(levels(projections$dbCluster))
  scCells <- isolate(sc$selectedCells()) # [1] "AAACCTGAGACACTAA-1" "AAACCTGAGACGACGT-1" "AAACCTGAGTCAAGCG-1" "AAACCTGCAAGAAGAG-1" "AAACCTGCAGACAGGT-1" "AAACCTGCATACAGCT-1" "AAACCTGGTGTGACCC-1"
  sampCol <- isolate(projectionColors$sampleNames)
  ccols <- isolate(projectionColors$dbCluster)
  # coE_heatmapSelectedModuleShow <- input$coE_heatmapSelectedModuleShow
  
  if (is.null(scEx_log) ||is.null(scCells) || length(scCells) == 0 ||
      is.null(projections)) {
    # output$coE_heatmapNull = renderUI(tags$h3(tags$span(style="color:red", "please select some cells")))
    return(NULL)
    # return(
    #   list(
    #     src = "empty.png",
    #     contentType = "image/png",
    #     width = 96,
    #     height = 96,
    #     alt = "heatmap should be here"
    #   )
    # )
  }
  # else {
  # output$coE_heatmapNull = NULL
  # }
  
  if (is.null(scCells) || length(scCells) == 0) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("No cells selected", id = "coE_heatmapSelectedReactiveProbl", type = "error", duration = 10)
    }
  }
  
  scEx_matrix <- assays(scEx_log)[[1]]
  featureData <- rowData(scEx_log)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/selectedHeatmap.RData", list = c(ls()))
  }
  # load(file = "~/SCHNAPPsDebug/selectedHeatmap.RData")
  
  retval <- coE_heatmapFunc(featureData, scEx_matrix, projections, genesin,
                            cells = scCells, sampCol = sampCol, ccols = ccols
  )
  
  setRedGreenButton(
    vars = list(
      c("coE_heatmapselected_geneids", isolate(input$coE_heatmapselected_geneids)),
      c("coE_heatmapselected_cells", isolate(coE_selctedCluster()$selectedCells())),
      c("coE_heatmapselected_sampcolPal", isolate(projectionColors$sampleNames)),
      c("coE_heatmapselected_cluscolPal", isolate(projectionColors$dbCluster))
    ),
    button = "updateHeatMapSelectedParameters"
  )
  
  exportTestValues(coE_heatmapSelectedReactive = {
    retVal
  })
  return(retval)
})

# coE_topExpGenesTable ----
#' coE_topExpGenesTable
#' in coexpressionSelected tab, showing the table of top expressed genes for a given
#' selection
#' coEtgPerc = genes shown have to be expressed in at least X % of cells
#' coEtgMinExpr = genes shown have at least to X UMIs expressed
coE_topExpGenesTable <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_topExpGenesTable started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_topExpGenesTable")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_topExpGenesTable")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_topExpGenesTable", id = "coE_topExpGenesTable", duration = NULL)
  }
  
  scEx_log <- scEx_log()
  projections <- projections()
  clicked <- input$updateMinExprSelectedParameters
  coEtgPerc <- isolate(input$coEtgPerc)
  coEtgminExpr <- isolate(input$coEtgMinExpr)
  sc <- isolate(coE_selctedCluster())
  # scCL <- sc$cluster
  # scCL <- isolate(levels(projections$dbCluster))
  scCells <- isolate(sc$selectedCells())
  # coEtgMinExprShow <- input$coEtgMinExprShow
  
  if (is.null(scEx_log) || is.null(scCells)) {
    if (DEBUG) if (is.null(scEx_log)) cat(file = stderr(), "coE_topExpGenesTable scEx_log null.\n")
    if (DEBUG) if (is.null(scCells)) cat(file = stderr(), "coE_topExpGenesTable scCells null.\n")
    return(NULL)
  }
  if (is.null(scCells) || length(scCells) == 0) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("No cells selected", id = "coE_topExpGenesTableProbl", type = "error", duration = 10)
    }
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/output_coE_topExpGenes.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/output_coE_topExpGenes.RData")
  
  featureData <- rowData(scEx_log)
  # we only work on cells that have been selected
  mat <- assays(scEx_log)[[1]][, scCells]
  # only genes that express at least coEtgminExpr UMIs
  mat[mat < coEtgminExpr] <- 0
  # only genes that are expressed in coEtgPerc or more cells
  allexpressed <- Matrix::rowSums(mat > 0) / length(scCells) * 100 >= coEtgPerc
  mat <- mat[allexpressed, ]
  
  cv <- function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  }
  matCV <- apply(mat, 1, cv)
  # top.genes <- as.data.frame(exprs(scEx_log))
  maxRows <- min(nrow(mat), 200)
  top.genesOrder <- order(matCV, decreasing = TRUE)[1:maxRows]
  retVal <- NULL
  if (dim(mat)[1] > 0) {
    mat <- mat[top.genesOrder, ]
    fd <- featureData[rownames(mat), c("symbol", "Description")]
    matCV <- matCV[rownames(mat)]
    fd <- cbind(fd, matCV)
    colnames(fd) <- c("gene", "description", "CV")
    # since we are returning a table to be plotted, we convert to regular table (non-sparse)
    outMat <- cbind2(fd, as.matrix(mat))
    rownames(outMat) <- make.unique(as.character(outMat$gene), sep = "_#_")
    retVal <- as.data.frame(outMat)
  }
  
  setRedGreenButton(
    vars = list(
      c("coEtgPerc", isolate(input$coEtgPerc)),
      c("coEtgMinExpr", isolate(input$coEtgMinExpr)),
      c("coE_heatmapselected_cells", isolate(coE_selctedCluster()$selectedCells()))
    ),
    button = "updateMinExprSelectedParameters"
  )
  
  exportTestValues(coE_topExpGenesTable = {
    retVal
  })
  return(retVal)
})


scranFindMarkerFullReactiveTable <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_scranFindMarkerTableReact started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_scranFindMarkerTableReact")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_scranFindMarkerTableReact")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_scranFindMarkerTableReact", id = "coE_scranFindMarkerTableReact", duration = NULL)
  }
  # deepDebug()
  scEx_log <- isolate(scEx_log())
  projections <- isolate(projections())
  direction <- isolate(input$coE_direction)
  lfc <- isolate(input$coE_lfc)
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "pca:NULL because scEx_log is null\n")
    }
    return(NULL)
  }
  clicked <- input$scranFindMarkerApply
  wmarkers <- tryCatch({
    scran::findMarkers(scEx_log, 
                       projections$dbCluster,
                       direction = direction,
                       lfc = lfc)
    
    
    
  }, error = function(e) {
    return(NULL)
  })
  updateSelectInput(session = session, inputId = "coE_scranFindMarkerCluster",
                    choices = levels(projections$dbCluster))
  
  
  return(wmarkers)
})

# coE_scranFindMarkerTableReact ----
coE_scranFindMarkerTableReact <- reactive({
  
  markerlist = scranFindMarkerFullReactiveTable()
  selectedCluster = input$coE_scranFindMarkerCluster
  projections = projections()
  req(markerlist)
  if(! selectedCluster %in% levels(projections$dbCluster)) return(NULL)
  
  return(markerlist[[selectedCluster]] %>% as.data.frame())
})


# coE_topExpCCTable ----
#' coE_topExpCCTable
#' in coE_topExpCCTable tab we show a table correlation coefficients and associated p-values
#' using the Hmisc::rcorr function
#' coEtgPerc = genes shown have to be expressed in at least X % of cells
#' coEtgMinExpr = genes shown have at least to X UMIs expressed
coE_topExpCCTable <- reactive({
  if (!"Hmisc" %in% rownames(installed.packages())) {
    showNotification("Please install Hmisc", id = "hmiscError", type = "error", duration = NULL)
    return(NULL)
  }
  require("Hmisc")
  if (DEBUG) cat(file = stderr(), "coE_topExpCCTable started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_topExpCCTable")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_topExpCCTable")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_topExpCCTable", id = "coE_topExpCCTable", duration = NULL)
  }
  # deepDebug()
  scEx_log <- scEx_log()
  projections <- projections()
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "coE_topExpCCTable: scEx_log:NULL\n")
    }
    return(NULL)
  }
  clicked <- input$updatetopCCGenesSelectedParameters
  
  genesin <- isolate(input$coE_heatmapselected_geneids)
  # coEtgPerc <- input$coEtgPerc
  # coEtgminExpr <- input$coEtgMinExpr
  sc <- isolate(coE_selctedCluster())
  # scCL <- sc$cluster
  # scCL <- levels(projections$dbCluster)
  scCells <- isolate(sc$selectedCells())
  # coE_topCCGenesShow <- input$coE_topCCGenesShow
  
  if (is.null(scCells) || length(scCells) == 0) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("No cells selected", id = "coE_topExpCCTableProbl", type = "error", duration = 10)
      return(NULL)
    }
  }
  
  featureData <- rowData(scEx_log)
  genesin <- geneName2Index(genesin, featureData)
  if (is.null(genesin)) {
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_topExpCCTable.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/coE_topExpCCTable.RData")
  
  # get numeric columns from projections.
  nums <- unlist(lapply(projections, is.numeric))
  numProje <- projections[, nums]
  # colnames(numProje)
  genesin <- unique(genesin)
  scCells <- scCells[scCells %in% colnames(assays(scEx_log)[[1]])]
  # we only work on cells that have been selected
  mat <- assays(scEx_log)[[1]][genesin, scCells, drop = FALSE]
  # only genes that express at least coEtgminExpr UMIs
  # mat[mat < coEtgminExpr] <- 0
  # only genes that are expressed in coEtgPerc or more cells
  # allexpressed <- Matrix::rowSums(mat > 0) / length(scCells) * 100 >= coEtgPerc
  # mat <- mat[allexpressed, ]
  
  if (length(mat) == 0) {
    return(NULL)
  }
  rownames(mat) <- featureData[rownames(mat), "symbol"]
  mat <- mat[!Matrix::rowSums(mat) == 0, , drop = FALSE]
  numProje <- t(numProje)[, colnames(mat), drop = FALSE]
  corrInput <- as.matrix(rbind(numProje, mat))
  # rownames(res2$r)
  res2 <- rcorr(t(corrInput))
  
  flatMat <- flattenCorrMatrix(res2$r, res2$P)
  
  retVal <- as.data.frame(flatMat[, ])
  # duplicated(retVal$row)
  rownames(retVal) <- make.unique(paste0(retVal$row, "_##_", retVal$column), sep = "_#_")
  
  setRedGreenButton(
    vars = list(
      c("coE_heatmapselected_geneids", isolate(input$coE_heatmapselected_geneids)),
      c("coE_heatmapselected_cells", isolate(coE_selctedCluster()$selectedCells()))
    ),
    button = "updatetopCCGenesSelectedParameters"
  )
  
  
  exportTestValues(coE_topExpCCTable = {
    retVal
  })
  return(retVal)
})

#' coE_geneGrp_vioFunc
#' generates a ggplot object with a violin plot
#' optionally creates all combinations
coE_geneGrp_vioFunc <- function(genesin, projections, scEx, featureData, minMaxExpr = c(-1,1),
                                dbCluster, coE_showPermutations = FALSE, 
                                projectionColors, showExpression = FALSE, coE_scale = "count") {
  
  
  if (DEBUG) cat(file = stderr(), "coE_geneGrp_vioFunc started.\n")
  start.time <- base::Sys.time()
  require(BiocParallel)
  on.exit({
    printTimeEnd(start.time, "coE_geneGrp_vioFunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_geneGrp_vioFunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_geneGrp_vioFunc", id = "coE_geneGrp_vioFunc", duration = NULL)
  }
  
  suppressMessages(require(gtools))
  suppressMessages(require(stringr))
  
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")[[1]]
  
  map <- rownames(featureData[which(toupper(featureData$symbol) %in% genesin), ])
  pc = projectionColors 
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_geneGrp_vioFunc.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/coE_geneGrp_vioFunc.RData")
  
  if (coE_showPermutations & showExpression) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "Please use only one of show expression and permuationas",
        id = "nogtools",
        type = "error",
        duration = NULL
      )
    }
    return(NULL)
  }
  if (length(map) == 0) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "no genes found",
        id = "heatmapWarning",
        type = "warning",
        duration = 20
      )
    }
    return(NULL)
  }
  
  
  if(showExpression) {
    expression <- Matrix::colSums(assays(scEx)[[1]][map, , drop = F])
    ylabText <- "sum of normalized read counts"
    
  } else {
    expression <- Matrix::colSums(assays(scEx)[[1]][map, , drop = F] >= minMaxExpr[1] & 
                                    assays(scEx)[[1]][map, , drop = F] <= minMaxExpr[2])
    ylabText <- "number genes from list"
  }
  projections <- cbind(projections, coExpVal = expression)
  permsNames <- as.character(1:max(expression))
  
  
  
  
  # only meaningful if not showing expression values
  if (coE_showPermutations) {
    if ("gtools" %in% rownames(installed.packages())) {
      perms <- rep("", length(expression))
      ylabText <- "Combinations"
      xPerm <- length(genesin)
      if (xPerm > 5) {
        xPerm <- 5
        warning("reducing number of combinations to 5")
      }
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification(
          "show permutations Reactive domain null",
          id = "heatmapWarning",
          type = "warning",
          duration = NULL
        )
      }
      x <- bplapply(1:xPerm, FUN = function(r) finner(xPerm, r, genesin, featureData, scEx, perms, minMaxExpr))
      
      for (idx in 1:length(x)) {
        perms <- combinePermutations(perms, x[[idx]])
      }
      perms <- factor(perms)
      permsNames <- levels(perms)
      permsNum <- unlist(lapply(strsplit(permsNames, "\\+"), length))
      perms <- factor(as.character(perms), levels = permsNames[order(permsNum)])
      permsNames <- str_wrap(levels(perms))
      perms <- as.integer(perms)
      projections[,'coExpVal'] <- perms
    } else {
      if (!"gtools" %in% rownames(installed.packages())) {
        showNotification(
          "please install gtools",
          id = "nogtools",
          type = "error",
          duration = NULL
        )
      }
      cat(file = stderr(), "Please install gtools: install.packages('gtools')")
      return(NULL)
    }
  }
  
  # browser()
  prj <- factor(projections[, dbCluster])
  if(dbCluster %in% names(pc)){
    mycolPal = pc[[dbCluster]]
  }else{
     mycolPal <- colorRampPalette(RColorBrewer::brewer.pal(
      n = 6, name =
        "RdYlBu"
    ))(length(levels(prj)))
  }
  # if (dbCluster == "sampleNames") {
  #   mycolPal <- sampCol
  # }
  # if (dbCluster == "dbCluster") {
  #   mycolPal <- ccols
  # }
  
  
  # p1 <- projections %>% plotly::plot_ly(
  #     x = prj,
  #     y = ~coExpVal,
  #      split = prj,
  #     type = 'violin',
  #     box = list(
  #       visible = T
  #     ),
  #     meanline = list(
  #       visible = T
  #     ), color=prj, colors = ccols
  #   ) %>%
  #     layout(
  #       xaxis = list(
  #         title = dbCluster
  #       ),
  #       yaxis = list(
  #         title = ylabText,
  #         zeroline = F,
  #         ticknames = permsNames
  #       ),
  #       annotations = list(y = permsNames, yref = "y")
  #     )
  #
  # p1
  #
  
  if(showExpression) {
    scY = NULL
  } else {
    scY = scale_y_continuous(breaks = 1:length(permsNames), labels = str_wrap(permsNames))
  }
  # browser()
  p1 <-
    ggplot(projections, aes(prj, coExpVal,
                            fill = .data[[dbCluster]]
    )) +
    geom_violin(scale = coE_scale) +
    scale_fill_manual(values = mycolPal, aesthetics = "fill") +
    stat_summary( # plot the centered dots
      fun = median,
      geom = "point",
      size = 5,
      color = "black"
    ) +
    stat_summary(fun.data = n_fun, geom = "text") +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 60,
        size = 12,
        vjust = 0.5
      ),
      axis.text.y = element_text(size = 10),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 12),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      legend.position = "right"
    ) +
    xlab(dbCluster) + ylab(ylabText) + 
    scY
  
  # p1 <- ggplotly(p1)
  # p1 + NULL
  
  
  
  return(p1)
}

######## grouped plotly version

#' coE_geneGrp_vioFunc
#' generates a ggplot object with a violin plot
#' optionally creates all combinations
coE_geneGrp_vioFunc2 <- function(genesin, projections, scEx, featureData, minMaxExpr = c(-1,1),
                                 dbCluster, projectionColors ) {
  if (DEBUG) cat(file = stderr(), "coE_geneGrp_vioFunc2 started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_geneGrp_vioFunc2")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_geneGrp_vioFunc2")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_geneGrp_vioFunc2", id = "coE_geneGrp_vioFunc2", duration = NULL)
  }
  
  suppressMessages(require(gtools))
  suppressMessages(require(stringr))
  
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")[[1]]
  pc = projectionColors
  map <- rownames(featureData[which(toupper(featureData$symbol) %in% genesin), ])
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_geneGrp_vioFunc2.RData", list = c(ls()))
  }
  # cp =load(file="~/SCHNAPPsDebug/coE_geneGrp_vioFunc2.RData")
  
  if (length(map) == 0) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "no genes found",
        id = "heatmapWarning",
        type = "warning",
        duration = 20
      )
    }
    return(NULL)
  }
  
  # expression <- Matrix::colSums(assays(scEx)[[1]][map, , drop = F] >= minExpr)
  expression <- Matrix::colSums(assays(scEx)[[1]][map, , drop = F] >= minMaxExpr[1] & 
                                  assays(scEx)[[1]][map, , drop = F] <= minMaxExpr[2])
  ylabText <- "number of genes from list"
  
  # coExpVal = number of cells with exprssion over minExpr
  projections <- cbind(projections, coExpVal = expression * 1.0)
  permsNames <- as.character(1:max(expression))
  
  #first coordinate
  prj <- factor(projections[, dbCluster[1]])
  if(dbCluster[1] %in% names(pc)){
    mycolPal = pc[[dbCluster[1]]]
  }else{
    mycolPal <- colorRampPalette(RColorBrewer::brewer.pal(
      n = 6, name =
        "RdYlBu"
    ))(length(levels(prj)))
  }
  
  # mycolPal <- colorRampPalette(RColorBrewer::brewer.pal(
  #   n = 6, name =
  #     "RdYlBu"
  # ))(length(levels(prj)))
  # 
  # if (dbCluster[1] == "sampleNames") {
  #   mycolPal <- sampCol
  #   names(mycolPal) = names(sampCol)
  # }
  # if (dbCluster[1] == "dbCluster") {
  #   mycolPal <- ccols
  #   names(mycolPal) = names(ccols)
  # }
  
  
  
  if (length(dbCluster) == 2 ) {
    prj2 = factor(projections[, dbCluster[2]])
    #second coordinate
    # prj2 <- factor(projections[, dbCluster[2]])
    if(dbCluster[2] %in% names(pc)){
      mycolPal2 = pc[[dbCluster[2]]]
    }else{
      mycolPal2 <- colorRampPalette(RColorBrewer::brewer.pal(
        n = 6, name =
          "RdYlBu"
      ))(length(levels(prj)))
    }
    # mycolPal2 <- colorRampPalette(RColorBrewer::brewer.pal(
    #   n = 6, name =
    #     "RdYlBu"
    # ))(length(levels(prj2)))
    # names(mycolPal2) = levels(prj2)
    # 
    # if (dbCluster[2] == "sampleNames") {
    #   mycolPal2 <- sampCol
    #   names(mycolPal2) = names(sampCol)
    # }
    # if (dbCluster[2] == "dbCluster") {
    #   mycolPal2 <- ccols
    #   names(mycolPal2) = names(ccols)
    # }
  } else {
    if (length(dbCluster) == 1 ) {
      prj2 = factor(rep(1,nrow(projections)))
      mycolPal2 = list('1'="#2D96FA")
    } else {
      # should not happen
    }
  }
  
  showLegend <- rep(F, length(levels(projections[, dbCluster[1]])))
  showLegend[1] = T
  
  if(lapply(levels(projections[, dbCluster[1]]), FUN = function(x)!is.na(suppressWarnings(as.numeric(x)))) %>% unlist() %>% any() ){
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("Some values can be interpreted as numericals and will be removed", id = "coE_geneGrp_vioFunc2Warn", type = "warning", duration = NULL)
    } 
  }
  

  p1 <- projections %>% plotly::plot_ly(type = 'violin', colors=mycolPal2)
  for (lv2 in levels(prj2)) {
    for (lv in levels(projections[, dbCluster[1]])) {
      p1 <- p1 %>% plotly::add_trace(
      x = {
        # cat(file = stderr(), projections[prj2 == lv2 & projections[dbCluster[1]] == lv, dbCluster[1]] %>% as.character() )
        projections[prj2 == lv2 & projections[dbCluster[1]] == lv, dbCluster[1]] %>% as.character() 
        },
      y = projections$coExpVal[prj2 == lv2 & projections[dbCluster[1]] == lv],
      legendgroup = lv,
      scalegroup = lv,
      showlegend = showLegend[which(lv == levels(projections[, dbCluster[1]]))],
      name = lv2,
      visible=T,
      box = list(
        visible = T
      ),
      meanline = list(
        visible = T
      )
      ,line = list(color=mycolPal2[[lv2]])
      , fillcolor = mycolPal2[[lv2]]
    )
    # print(p1)
    }
  }
  
  p1 <- p1  %>%
    plotly::layout(
      xaxis = list(
        title = paste(dbCluster, collapse = " + ")
      ),
      yaxis = list(
        title = ylabText,
        zeroline = F
        # ,
        # ticknames = permsNames
      )
      ,
      # violingap = 0,violingroupgap = 1,
      violinmode = 'group'
      # ,
      # annotations = list(y = permsNames, yref = "y")
    )
  
  # p1 <- p1 %>% 
  
  # p1
  
  
  
  # p1 <-
  #   ggplot(projections, aes_string(prj, "coExpVal",
  #                                  fill = factor(projections[, dbCluster])
  #   )) +
  #   geom_violin(scale = "count") +
  #   scale_fill_manual(values = mycolPal, aesthetics = "fill") +
  #   stat_summary( # plot the centered dots
  #     fun = median,
  #     geom = "point",
  #     size = 5,
  #     color = "black"
  #   ) +
  #   stat_summary(fun.data = n_fun, geom = "text") +
  #   theme_bw() +
  #   theme(
  #     axis.text.x = element_text(
  #       angle = 60,
  #       size = 12,
  #       vjust = 0.5
  #     ),
  #     axis.text.y = element_text(size = 10),
  #     strip.text.x = element_text(size = 16),
  #     strip.text.y = element_text(size = 12),
  #     axis.title.x = element_text(face = "bold", size = 16),
  #     axis.title.y = element_text(face = "bold", size = 16),
  #     legend.position = "right"
  #   ) +
  #   xlab(dbCluster) + ylab(ylabText) +
  #   scale_y_continuous(breaks = 1:length(permsNames), labels = str_wrap(permsNames))
  # 
  # # p1 <- ggplotly(p1)
  return(p1)
}





# save to history violoin observer 2----
observe(label = "save2histVio", {
  clicked  = input$save2HistVio
  if (DEBUG) cat(file = stderr(), "observe input$save2HistVio \n")
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
              comment = paste("# violin plot\n",
                              "fun = plotData$plotData$plotFunc\n", 
                              "environment(fun) = environment()\n",
                              "plotData$plotData$outfile=NULL\n",
                              "print(do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)]))\n"
              ),
              plotData = .schnappsEnv[["coE_geneGrp_vio_plot"]])
  
})

# save to history violoin observer ----
observe(label = "save2histVio2", {
  clicked  = input$save2HistVio2
  if (DEBUG) cat(file = stderr(), "observe input$save2HistVio2 \n")
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
              comment = paste("# violin plot\n",
                              "fun = plotData$plotData$plotFunc\n", 
                              "environment(fun) = environment()\n",
                              "plotData$plotData$outfile=NULL\n",
                              "print(do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)]))\n"
              ),
              plotData = .schnappsEnv[["coE_geneGrp_vio_plot2"]])
  
})


# coE_updateInputXviolinPlot----
#' coE_updateInputXviolinPlot
#' Update x/y axis selection possibilities for violin plot
#' could probably be an observer, but it works like this as well...
# .schnappsEnv$coE_vioGrp <- "sampleNames"


coE_updateInputXviolinPlot <- observe({
  if (DEBUG) cat(file = stderr(), "coE_updateInputXviolinPlot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_updateInputXviolinPlot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_updateInputXviolinPlot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_updateInputXviolinPlot", id = "coE_updateInputXviolinPlot", duration = NULL)
  }
  
  tsneData <- projections()
  projFactors <- projFactors()
  # Can use character(0) to remove all choices
  if (is.null(tsneData)) {
    return(NULL)
  }
  updateSelectInput(
    session,
    "coE_dimension_xVioiGrp",
    choices = projFactors,
    selected = .schnappsEnv$coE_dimension_xVioiGrp
  )
  updateSelectInput(
    session,
    "coE_dimension_xVioiGrp2",
    choices = projFactors,
    selected = .schnappsEnv$coE_dimension_xVioiGrp2
  )
})


observe(label = "obs_coE_heatmap_geneids", x= {
  .schnappsEnv$defaultValues[["coE_heatmap_geneids"]] = input$coE_heatmap_geneids
})

# coE_heatmapReactive -------
# reactive for module pHeatMapModule
# for all clusters menu item
coE_heatmapReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_heatmapReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_heatmapReactive")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_heatmapReactive")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_heatmapReactive", id = "coE_heatmapReactive", duration = NULL)
  }
  
  scEx_log <- scEx_log()
  projections <- projections()
  genesin <- input$coE_heatmap_geneids
  sampCol <- projectionColors$sampleNames
  ccols <- projectionColors$dbCluster
  projections <- projections()
  
  if (is.null(scEx_log) | is.null(projections)) {
    return(list(
      src = "empty.png",
      contentType = "image/png",
      width = 96,
      height = 96,
      alt = "heatmap should be here"
    ))
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/heatmap.RData", list = c(ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/heatmap.RData")
  
  featureData <- rowData(scEx_log)
  if (genesin == "") return(NULL)
  
  scEx_matrix <- as.matrix(assays(scEx_log)[["logcounts"]])
  retVal <- coE_heatmapFunc(
    featureData = featureData, scEx_matrix = scEx_matrix,
    projections = projections, genesin = genesin, cells = colnames(scEx_matrix),
    sampCol = sampCol, ccols = ccols
  )
  
  exportTestValues(coE_heatmapReactive = {
    retVal
  })
  return(retVal)
})


## alluvialPlotFunc ----

alluvialPlotFunc <- function(dat, alluiv1, alluiv2) {
  gg = ggplot(as.data.frame(dat),
              aes(  axis1 = .data[[alluiv1]], axis2 = .data[[alluiv2]])) +
    geom_alluvium(aes(fill = .data[[alluiv1]]), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", infer.label = TRUE) +
    scale_x_discrete(limits = c(alluiv1, alluiv2), expand = c(.05, .05)) +
    # scale_fill_brewer(type = "qual", palette = "Set1") +
    ggtitle(paste("Alluvial plot of ", alluiv1, "and", alluiv2))
}

## save2HistAlluvial observer ---- 
observe(label = "save2HistAlluvial", {
  clicked  = input$save2HistAlluvial
  if (DEBUG) cat(file = stderr(), "observe save2HistAlluvial \n")
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
              comment = paste0("# Alluvial plot\n",
                               "fun = plotData$plotData$plotFunc\n", 
                               "environment(fun) = environment()\n",
                               "print(do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)]))\n"
              ),
              plotData = .schnappsEnv[["coE_alluvialPlot"]])
  
})


