# heatmapFunc ---------------------------------
# used by both selection and all to create input for heatmap module
coE_heatmapFunc <- function(featureData, scEx_matrix, projections, genesin, cells, sampCol, ccols) {
  if (DEBUG) cat(file = stderr(), "coE_heatmapFunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_heatmapFunc")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "coE_heatmapFunc")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_heatmapFunc", id = "coE_heatmapFunc", duration = NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_heatmapFunc.RData", list = c(ls(envir = globalenv()), ls()))
  }
  # load(file = "~/SCHNAPPsDebug/coE_heatmapFunc.RData")
  
  #  create parameters used for pheatmap module
  genesin <- geneName2Index(genesin, featureData)
  if (is.null(genesin) | is.null(cells)) {
    return(NULL)
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
  
  nonZeroRows <- which(Matrix::rowSums(expression) > 0)
  
  annCols <- list("sampleNames" = sampCol,
                  "Cluster" = ccols)
  
  retVal <- list(
    mat = expression[nonZeroRows, order(annotation[, 1], annotation[, 2])],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    # scale = "row",
    fontsize_row = 14,
    labels_col = colnames(expression),
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
  )
  
  # print debugging information on the console
  printTimeEnd(start.time, "inputData")
  # for automated shiny testing using shinytest
  exportTestValues(coE_heatmapFunc = {
    retVal
  })
  
  # this is what is run in the module
  # do.call(TRONCO::pheatmap, retVal)
  return(retVal)
}


# coE_heatmapSelectedReactive ----
# reactive function for selected heatmap
coE_heatmapSelectedReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_heatmapSelectedReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_heatmapSelectedReactive")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "coE_heatmapSelectedReactive")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_heatmapSelectedReactive", id = "coE_heatmapSelectedReactive", duration = NULL)
  }
  
  scEx_log <- scEx_log()
  projections <- projections()
  genesin <- input$coE_heatmapselected_geneids
  sc <- coE_selctedCluster()
  # scCL <- sc$cluster # "1" "2" "3" "4"
  scCL <- levels(projections$dbCluster)
  scCells <- sc$selectedCells() # [1] "AAACCTGAGACACTAA-1" "AAACCTGAGACGACGT-1" "AAACCTGAGTCAAGCG-1" "AAACCTGCAAGAAGAG-1" "AAACCTGCAGACAGGT-1" "AAACCTGCATACAGCT-1" "AAACCTGGTGTGACCC-1"
  sampCol <- sampleCols$colPal
  ccols <- clusterCols$colPal
  
  if (is.null(scEx_log) ||
      is.null(projections) || is.null(scCells) || length(scCells) == 0) {
    # output$coE_heatmapNull = renderUI(tags$h3(tags$span(style="color:red", "please select some cells")))
    return(
      list(
        src = "empty.png",
        contentType = "image/png",
        width = 96,
        height = 96,
        alt = "heatmap should be here"
      )
    )
  } 
  # else {
    # output$coE_heatmapNull = NULL
  # }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/selectedHeatmap.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/selectedHeatmap.RData")
  
  
  scEx_matrix <- assays(scEx_log)[[1]]
  featureData <- rowData(scEx_log)
  
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/selectedHeatmap.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/selectedHeatmap.RData")
  
  retval <- coE_heatmapFunc(featureData, scEx_matrix, projections, genesin,
                            cells = scCells, sampCol = sampCol, ccols = ccols
  )

  exportTestValues(coE_heatmapSelectedReactive = {retVal})  
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
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "coE_topExpGenesTable")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_topExpGenesTable", id = "coE_topExpGenesTable", duration = NULL)
  }
  
  scEx_log <- scEx_log()
  projections <- projections()
  coEtgPerc <- input$coEtgPerc
  coEtgminExpr <- input$coEtgMinExpr
  sc <- coE_selctedCluster()
  # scCL <- sc$cluster
  scCL <- levels(projections$dbCluster)
  scCells <- sc$selectedCells()
  
  if (is.null(scEx_log) || is.null(scCells)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/output_coE_topExpGenes.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/output_coE_topExpGenes.RData")
  
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
  retVal = NULL
  if (dim(mat)[1] > 0) {
    mat <- mat[top.genesOrder, ]
    fd <- featureData[rownames(mat), c("symbol", "Description")]
    matCV <- matCV[rownames(mat)]
    fd <- cbind(fd, matCV)
    colnames(fd) <- c("gene", "description", "CV")
    # since we are returning a table to be plotted, we convert to regular table (non-sparse)
    outMat <- cbind2(fd, as.matrix(mat))
    rownames(outMat) <- make.unique(as.character(outMat$gene), sep = "___")
    retVal = as.data.frame(outMat)
  } 
  
  exportTestValues(coE_topExpGenesTable = {retVal})  
  return(retVal)
})



#' coE_geneGrp_vioFunc
#' generates a ggplot object with a violin plot
#' optionally creates all permutations.
coE_geneGrp_vioFunc <- function(genesin, projections, scEx, featureData, minExpr = 1,
                            dbCluster, coE_showPermutations = FALSE, sampCol, ccols) {
  if (DEBUG) cat(file = stderr(), "coE_geneGrp_vioFunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_geneGrp_vioFunc")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "coE_geneGrp_vioFunc")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_geneGrp_vioFunc", id = "coE_geneGrp_vioFunc", duration = NULL)
  }
  
  require(gtools)
  require(stringr)
  
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")[[1]]
  
  map <- rownames(featureData[which(toupper(featureData$symbol) %in% genesin), ])
  
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_geneGrp_vioFunc.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/coE_geneGrp_vioFunc.RData")
  
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
  
  expression <- Matrix::colSums(assays(scEx)[[1]][map, ] >= minExpr)
  ylabText <- "number genes from list"
  
  if (coE_showPermutations) {
    perms <- rep("", length(expression))
    ylabText <- "Permutations"
    xPerm <- length(genesin)
    if (xPerm > 10) {
      xPerm <- 10
      warning("reducing number of permutations to 10")
    }
    # cat(file = stderr(), paste("===violin-", vIdx,"-", difftime(Sys.time(), start.time, units = "min"), " min\n")); vIdx = vIdx+1;start.time <- Sys.time()
    for (r in 1:xPerm) {
      comb <- combinations(xPerm, r, genesin)
      for (cIdx in 1:nrow(comb)) {
        map <-
          rownames(featureData[which(toupper(featureData$symbol) %in% comb[cIdx, ]), ])
        # permIdx <- Matrix::colSums(exprs(gbm[map, ]) >= minExpr) == length(comb[cIdx, ])
        
        permIdx <- Matrix::colSums(assays(scEx)[[1]][map, , drop = FALSE] >= minExpr) == length(comb[cIdx, ])
        perms[permIdx] <- paste0(comb[cIdx, ], collapse = "+")
      }
    }
    # cat(file = stderr(), paste("===violin-", vIdx,"-", difftime(Sys.time(), start.time, units = "min"), " min\n")); vIdx = vIdx+1;start.time <- Sys.time()
    perms <- factor(perms)
    permsNames <- levels(perms)
    permsNum <- unlist(lapply(strsplit(permsNames, "\\+"), length))
    perms <- factor(as.character(perms), levels = permsNames[order(permsNum)])
    permsNames <- str_wrap(levels(perms))
    perms <- as.integer(perms)
    projections <- cbind(projections, coExpVal = perms)
  } else {
    projections <- cbind(projections, coExpVal = expression)
    permsNames <- as.character(1:max(expression))
  }
  
  prj <- factor(projections[, dbCluster])
  mycolPal <- colorRampPalette(RColorBrewer::brewer.pal(
    n = 6, name =
      "RdYlBu"
  ))(length(levels(prj)))
  
  if (dbCluster == "sampleNames") {
    mycolPal <- sampCol
  }
  if (dbCluster == "dbCluster") {
    mycolPal <- ccols
  }
  
  
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
    
    
    p1 <-
    ggplot(projections, aes_string(prj, "coExpVal",
                                   fill = factor(projections[, dbCluster])
    )) +
    geom_violin(scale = "count") +
    scale_fill_manual(values = mycolPal, aesthetics = "fill") +
    stat_summary( # plot the centered dots
      fun.y = median,
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
    xlab(dbCluster) +

    scale_y_continuous(breaks = 1:length(permsNames), labels = str_wrap(permsNames)) +
    ylab(ylabText)
   # p1 <- ggplotly(p1)
  return(p1)
}

#' coE_somFunction
#' iData = expression matrix, rows = genes
#' cluster genes in SOM 
#' returns genes associated together within a som-node
coE_somFunction <- function(iData, nSom, geneName) {
  if (DEBUG) cat(file = stderr(), "coE_somFunction started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_somFunction")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "coE_somFunction")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_somFunction", id = "coE_somFunction", duration = NULL)
  }

    require(kohonen)
  require(Rsomoclu)
  if (sum(geneName %in% rownames(iData)) == 0) return(NULL)
  res2 <- Rsomoclu.train(
    input_data = iData,
    nSomX = nSom, nSomY = nSom,
    nEpoch = 10,
    radius0 = 0,
    radiusN = 0,
    radiusCooling = "linear",
    mapType = "planar",
    gridType = "rectangular",
    scale0 = 1,
    scaleN = 0.01,
    scaleCooling = "linear"
  )
  
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_somFunction.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/coE_somFunction.RData")
  rownames(res2$globalBmus) <- make.unique(as.character(rownames(iData)), sep = "___")
  simGenes <- rownames(res2$globalBmus)[which(res2$globalBmus[, 1] == res2$globalBmus[geneName, 1] &
                                                res2$globalBmus[, 2] == res2$globalBmus[geneName, 2])]
  return(simGenes)
}

# coE_heatmapSOMReactive ----
#' coE_heatmapSOMReactive
#' calculates a self organizing map (SOM) on the expression data and identifies genes
#' that cluster together with a gene of interest
# TODO: check that we are using only raw counts and not normalized
coE_heatmapSOMReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_heatmapSOMReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_heatmapSOMReactive")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "coE_heatmapSOMReactive")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_heatmapSOMReactive", id = "coE_heatmapSOMReactive", duration = NULL)
  }
  
  if (!is.null(getDefaultReactiveDomain()))
    removeNotification(id = "heatmapWarning")

  scEx_log <- scEx_log()
  projections <- projections()
  genesin <- input$coE_geneSOM
  nSOM <- input$coE_dimSOM
  sampCol <- sampleCols$colPal
  ccols <- clusterCols$colPal
  
  if (is.null(scEx_log)) {
    return(
      list(
        src = "empty.png",
        contentType = "image/png",
        width = 96,
        height = 96,
        alt = "heatmap should be here"
      )
    )
  }

  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_heatmapSOMReactive.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/coE_heatmapSOMReactive.RData")

  scEx_matrix <- as.matrix(assays(scEx_log)[[1]])
  featureData <- rowData(scEx_log)
  # go from readable gene name to ENSG number
  genesin <- geneName2Index(genesin, featureData)

  geneNames <- coE_somFunction(iData = scEx_matrix, nSom = nSOM, geneName = genesin)

  # plot the genes found
  output$coE_somGenes <- renderText({
    paste(featureData[geneNames, "symbol"], collapse = ", ", sep = ",")
  })
  if (length(geneNames) < 2) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "no additional gene found. reduce size of SOM",
        id = "heatmapWarning",
        type = "warning",
        duration = 20
      )
    }
    return(NULL)
  }

  # create variables for heatmap module
  annotation <- data.frame(projections[, c("dbCluster", "sampleNames")])
  rownames(annotation) <- rownames(projections)
  colnames(annotation) <- c("Cluster", "sampleNames")
  annCols <- list("sampleNames" = sampCol,
                  "Cluster" = ccols)
  
  retVal <- list(
    mat = scEx_matrix[geneNames, ,drop = FALSE],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    # scale = "row",
    fontsize_row = 14,
    labels_row = featureData[geneNames, "symbol"],
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
  )
  # system.time(do.call(TRONCO::pheatmap, retVal))
  
  exportTestValues(coE_heatmapSOMReactive = {retVal })  
  return(retVal)
})

# coE_updateInputXviolinPlot----
#' coE_updateInputXviolinPlot
#' Update x/y axis selection possibilities for violin plot
#' could probably be an observer, but it works like this as well...
coE_vioGrp <<- "sampleNames"
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: coE_dimension_xVioiGrp\n"))
  coE_vioGrp <<- input$coE_dimension_xVioiGrp
})

coE_updateInputXviolinPlot <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_updateInputXviolinPlot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_updateInputXviolinPlot")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "coE_updateInputXviolinPlot")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_updateInputXviolinPlot", id = "coE_updateInputXviolinPlot", duration = NULL)
  }

    tsneData <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(tsneData)) {
    return(NULL)
  }
  
  # Can also set the label and select items
  # updateSelectInput(
  #   session,
  #   "dimension_x3",
  #   choices = colnames(tsneData),
  #   selected = colnames(tsneData)[1]
  # )
  # updateSelectInput(
  #   session,
  #   "dimension_y3",
  #   choices = colnames(tsneData),
  #   selected = colnames(tsneData)[2]
  # )
  
  coln <- colnames(tsneData)
  choices <- c()
  for (cn in coln) {
    if (length(levels(as.factor(tsneData[, cn]))) < 20) {
      choices <- c(choices, cn)
    }
  }
  if (length(choices) == 0) {
    choices <- c("no valid columns")
  }
  updateSelectInput(
    session,
    "coE_dimension_xVioiGrp",
    choices = choices,
    selected = coE_vioGrp
  )
})


# coE_heatmapReactive -------
# reactive for module pHeatMapModule
# for all clusters menu item
coE_heatmapReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_heatmapReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_heatmapReactive")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "coE_heatmapReactive")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_heatmapReactive", id = "coE_heatmapReactive", duration = NULL)
  }
  
  scEx_log <- scEx_log()
  projections <- projections()
  genesin <- input$coE_heatmap_geneids
  sampCol = sampleCols$colPal
  ccols <- clusterCols$colPal
  
  if (is.null(scEx_log) | is.null(projections)) {
    return(list(
      src = "empty.png",
      contentType = "image/png",
      width = 96,
      height = 96,
      alt = "heatmap should be here"
    ))
  }
  
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/heatmap.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/heatmap.RData")

  featureData <- rowData(scEx_log)
  scEx_matrix <- as.matrix(assays(scEx_log)[["logcounts"]])
  retVal <- coE_heatmapFunc(
    featureData = featureData, scEx_matrix = scEx_matrix,
    projections = projections, genesin = genesin, cells = colnames(scEx_matrix),
    sampCol = sampCol, ccols = ccols
  )
  
  exportTestValues(coE_heatmapReactive = {retVal})  
  return(retVal)
})
