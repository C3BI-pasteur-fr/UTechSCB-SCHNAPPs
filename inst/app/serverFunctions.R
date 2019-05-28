library(magrittr)

printTimeEnd <- function(start.time, messtr) {
  end.time <- base::Sys.time()
  if (DEBUG){
    cat(file = stderr(), paste("---", messtr,":done", difftime(end.time, start.time, units = "min"), " min\n"))  
  }
}



geneName2Index <- function(g_id, featureData) {
  if (DEBUG) cat(file = stderr(), "geneName2Index started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "geneName2Index")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "geneName2Index")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("geneName2Index", id = "geneName2Index", duration = NULL)
  }
  
  if(is.null(g_id)){
    return(NULL)
  }
  
  g_id <- toupper(g_id)
  g_id <- gsub(" ", "", g_id, fixed = TRUE)
  g_id <- strsplit(g_id, ",")
  g_id <- g_id[[1]]
  
  notFound <- g_id[!g_id %in% toupper(featureData$symbol)]
  if (length(featureData$symbol) == length(notFound)) {
    # in case there is only one gene that is not available.
    notFound <- g_id
  }
  if (length(notFound) > 0) {
    if (DEBUG) {
      cat(file = stderr(), paste("gene names not found: ", notFound, "\n"))
    }
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(paste("following genes were not found", notFound, collapse = " "),
                       id = "moduleNotFound", type = "warning",
                       duration = 20
      )
    }
  }
  
  geneid <- rownames(featureData[which(toupper(featureData$symbol) %in% toupper(g_id)), ])

  return(geneid)
}


updateProjectionsWithUmiCount <- function(dimX, dimY, geneNames, geneNames2 = NULL, scEx, projections) {
  featureData = rowData(scEx)
  # if ((dimY == "UmiCountPerGenes") | (dimX == "UmiCountPerGenes")) {
  geneNames <- geneName2Index(geneNames, featureData)
  # if (length(geneNames) > 0) {
  if ((length(geneNames) > 0) && (length(geneNames[[1]]) > 0)) {
    if (length(geneNames) == 1) {
      projections$UmiCountPerGenes <- assays(scEx)[[1]][geneNames, ]
    } else {
      projections$UmiCountPerGenes <- Matrix::colSums(assays(scEx)[[1]][geneNames, ])
    }
  } else {
    projections$UmiCountPerGenes <- 0
  }
  # }
  # }
  # if ((dimY == "UmiCountPerGenes2") | (dimX == "UmiCountPerGenes2")) {
  geneNames <- geneName2Index(geneNames2, featureData)
  # if (length(geneNames) > 0) {
  if ((length(geneNames) > 0) && (length(geneNames[[1]]) > 0)) {
    if (length(geneNames) == 1) {
      projections$UmiCountPerGenes2 <- assays(scEx)[[1]][geneNames, ]
    } else {
      projections$UmiCountPerGenes2 <- Matrix::colSums(assays(scEx)[[1]][geneNames, ])
    }
  } else {
    projections$UmiCountPerGenes2 <- 0
  }
  # }
  # }
  return(projections)
}


# append to heavyCalculations
append2list <- function(myHeavyCalculations, heavyCalculations) {
  for (hc in myHeavyCalculations) {
    if (length(hc) == 2 & is.character(hc[1]) & is.character(hc[2])) {
      heavyCalculations[[length(heavyCalculations) + 1]] <- hc
    } else {
      stop(paste("myHeavyCalculations is malformatted\nhc:", hc, "\nmyHeavyCalculations:", myHeavyCalculations, "\n"))
    }
  }
  return(heavyCalculations)
}

#### plot2Dprojection ----------------
# used in moduleServer and reports
plot2Dprojection <- function(scEx_log, projections, g_id, featureData,
                             geneNames, geneNames2, dimX, dimY, clId, grpN, legend.position, grpNs,
                             logx = FALSE, logy = FALSE, divXBy = "None", divYBy = "None", dimCol = "Gene.count",
                             colors = NULL) {
  geneid <- geneName2Index(g_id, featureData)
  
  if (length(geneid) == 0) {
    return(NULL)
  }
  # if (length(geneid) == 1) {
  #   expression <- exprs(scEx_log)[geneid, ,drop=FALSE]
  # } else {
  expression <- Matrix::colSums(assays(scEx_log)[[1]][geneid, , drop = FALSE])
  # }
  validate(need(is.na(sum(expression)) != TRUE, ""))
  # if (length(geneid) == 1) {
  #   expression <- exprs(scEx_log)[geneid, ]
  # } else {
  #   expression <- Matrix::colSums(exprs(scEx_log)[geneid, ])
  # }
  # validate(need(is.na(sum(expression)) != TRUE, ""))
  
  # geneid <- geneName2Index(geneNames, featureData)
  projections <- updateProjectionsWithUmiCount(
    dimX = dimX, dimY = dimY,
    geneNames = geneNames,
    geneNames2 = geneNames2,
    scEx = scEx_log, projections = projections
  )
  
  
  projections <- cbind(projections, expression)
  names(projections)[ncol(projections)] <- "exprs"
  
  if (DEBUG) {
    cat(file = stderr(), paste("output$sCA_dge_plot1:---", clId[1], "---\n"))
  }
  subsetData <- subset(projections, dbCluster %in% clId)
  # subsetData$dbCluster = factor(subsetData$dbCluster)
  # if there are more than 18 samples ggplot cannot handle different shapes and we ignore the
  # sample information
  if (length(as.numeric(as.factor(subsetData$sample))) > 18) {
    subsetData$shape <- as.factor(1)
  } else {
    subsetData$shape <- as.numeric(as.factor(subsetData$sample))
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/clusterPlot.RData", list = c(ls(), "legend.position", ls(envir = globalenv())))
    cat(file = stderr(), paste("plot2Dprojection saving done.\n"))
  }
  # load(file="~/SCHNAPPsDebug/clusterPlot.RData")
  if (nrow(subsetData) == 0) return(NULL)
  # subsetData$shape = as.factor(1)
  gtitle <- paste(g_id, clId, sep = "-Cluster", collapse = " ")
  if (nchar(gtitle) > 50) {
    gtitle <- paste(substr(gtitle, 1, 50), "...")
  }
  
  require(plotly)
  f <- list(
    family = "Courier New, monospace",
    size = 18,
    color = "#7f7f7f"
  )
  if (divXBy != "None") {
    subsetData[, dimX] <- subsetData[, dimX] / subsetData[, divXBy]
  }
  if (divYBy != "None") {
    subsetData[, dimY] <- subsetData[, dimY] / subsetData[, divYBy]
  }
  
  typeX <- typeY <- "linear"
  if (logx) {
    typeX <- "log"
  }
  if (logy) {
    typeY <- "log"
  }
  if (is.factor(subsetData[, dimX])) {
    typeX <- NULL
  }
  if (is.factor(subsetData[, dimY])) {
    typeY <- NULL
  }
  xAxis <- list(
    title = dimX,
    titlefont = f,
    type = typeX
  )
  yAxis <- list(
    title = dimY,
    titlefont = f,
    type = typeY
  )
  
  
  if (is.factor(subsetData[,dimX])) {
    subsetData[,dimX] = as.character(subsetData[,dimX])
  }
  if (is.factor(subsetData[,dimY])) {
    subsetData[,dimY] = as.character(subsetData[,dimY])
  }
  # dimCol = "Gene.count"
  # dimCol = "sampleNames"
  # subsetData$"__key__" = rownames(subsetData)
  p1 <- plotly::plot_ly(data = subsetData, source = "subset",
                        key = rownames(subsetData)) %>%
    add_trace(x = ~ get(dimX)
              ,y = ~ get(dimY),  
              type = "scatter" ,mode = "markers"
              ,text = ~ paste(1:nrow(subsetData), " ", rownames(subsetData), "<br />", subsetData$exprs)
              ,color = ~ get(dimCol)
              ,colors = colors
              ,showlegend =  TRUE
    ) %>% 
    layout(
      xaxis = xAxis,
      yaxis = yAxis,
      title = gtitle,
      dragmode = "select"
    )     
  if( is.factor(subsetData[,dimCol]) ) {
    
  } else {
    p1 = colorbar(p1, title = dimCol)
  }
  
  selectedCells <- NULL
  if (length(grpN) > 0) {
    if (length(grpNs[rownames(subsetData), grpN]) > 0 & sum(grpNs[rownames(subsetData), grpN], na.rm = TRUE) > 0) {
      grpNSub <- grpNs[rownames(subsetData), ]
      selectedCells <- rownames(grpNSub[grpNSub[, grpN], ])
    }
  }
  if (!is.null(selectedCells)) {
    shape <- rep("a", nrow(subsetData))
    selRows <- which(rownames(subsetData) %in% selectedCells)
    shape[selRows] <- "b"
    x1 <- subsetData[selectedCells, dimX, drop = FALSE]
    y1 <- subsetData[selectedCells, dimY, drop = FALSE]
    p1 <- p1 %>%
      add_trace(data = subsetData[selectedCells, ],
        x = x1[, 1], y = y1[, 1],
        marker = list(
          color = rep("red", nrow(x1)),
          size = 5
        ),
        text = ~ paste(
          rownames(subsetData[selectedCells, ]),
          "<br />", subsetData$exprs[selectedCells]
        ),
        type = "scatter",
        mode = "markers",
        name = "selected",
        key = rownames(subsetData[selectedCells, ])
      )
  }
  p1
}


# functions should go in external file

n_fun <- function(x) {
  return(data.frame(y = -0.5, label = paste0(length(x), "\ncells")))
}

diffLRT <- function(x, y, xmin = 1) {
  lrtX <- bimodLikData(x)
  lrtY <- bimodLikData(y)
  lrtZ <- bimodLikData(c(x, y))
  lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  return(pchisq(lrt_diff, 3, lower.tail = F))
}

bimodLikData <- function(x, xmin = 0) {
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  xal <- minmax(length(x2) / length(x), min = 1e-05, max = (1 - 1e-05))
  likA <- length(x1) * log(1 - xal)
  mysd <- sd(x2)
  if (length(x2) < 2) {
    mysd <- 1
  }
  likB <- length(x2) * log(xal) + sum(dnorm(x2, mean(x2), mysd, log = TRUE))
  return(likA + likB)
}

ainb <- function(a, b) {
  a2 <- a[a %in% b]
  return(a2)
}

minmax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}
set.ifnull <- function(x, y) {
  if (is.null(x)) {
    return(y)
  }
  return(x)
}

expMean <- function(x) {
  return(log(mean(exp(x) - 1) + 1))
}


DiffExpTest <- function(expression, cells.1, cells.2, genes.use = NULL, print.bar = TRUE) {
  genes.use <- set.ifnull(genes.use, rownames(expression))
  p_val <- unlist(lapply(genes.use, function(x) diffLRT(as.numeric(expression[x, cells.1]), as.numeric(expression[x, cells.2]))))
  to.return <- data.frame(p_val, row.names = genes.use)
  return(to.return)
}


heatmapPlotFromModule <- function(heatmapData, moduleName, input, projections) {
  
  addColNames <- input[[paste0(moduleName, "-ColNames")]]
  orderColNames <- input[[paste0(moduleName, "-orderNames")]]
  moreOptions <- input[[paste0(moduleName, "-moreOptions")]]
  colTree <- input[[paste0(moduleName, "-showColTree")]]
  
  if (is.null(heatmapData) | is.null(projections) | is.null(heatmapData$mat)) {
    return(NULL)
  }
  
  heatmapData$filename <- NULL
  
  if (length(addColNames) > 0 & moreOptions) {
    heatmapData$annotation_col <- projections[rownames(heatmapData$annotation_col), addColNames, drop = FALSE]
  }
  if (sum(orderColNames %in% colnames(projections)) > 0 & moreOptions) {
    heatmapData$cluster_cols <- FALSE
    colN <- rownames(psych::dfOrder(projections, orderColNames))
    colN <- colN[colN %in% colnames(heatmapData$mat)]
    heatmapData$mat <- heatmapData$mat[, colN, drop = FALSE]
    # return()
  }
  if (moreOptions) {
    heatmapData$cluster_cols = colTree
  }
  heatmapData$fontsize <- 14
  system.time(do.call(TRONCO::pheatmap, heatmapData))
  
}

# twoDplotFromModule ----
#' function to be used in markdown docs to ease the plotting of the clusterServer module
twoDplotFromModule <- function(twoDData, moduleName, input, projections, g_id, legend.position = "none") {

    grpNs <- groupNames$namesDF
  grpN <- make.names(input$groupName)
  
  dimY <- input[[paste0(moduleName,"-dimension_y")]]
  dimX <- input[[paste0(moduleName,"-dimension_x")]]
  dimCol <- input[[paste0(moduleName,"-dimension_col")]]
  clId <- input[[paste0(moduleName,"-clusters")]]

  geneNames <- input[[paste0(moduleName,"-geneIds")]]
  geneNames2 <- input[[paste0(moduleName,"-geneIds2")]]
  logx <- input[[paste0(moduleName,"-logX")]]
  logy <- input[[paste0(moduleName,"-logY")]]
  divXBy <- input[[paste0(moduleName,"-devideXBy")]]
  divYBy <- input[[paste0(moduleName,"-devideYBy")]]
  scols <- sampleCols$colPal
  ccols <- clusterCols$colPal
  
  
  if (is.null(scEx_log) | is.null(scEx_log) | is.null(projections)) {
    if (DEBUG) cat(file = stderr(), paste("output$clusterPlot:NULL\n"))
    return(NULL)
  }
  
  featureData <- rowData(scEx_log)
  if (is.null(g_id) || nchar(g_id) == 0) {
    g_id <- featureData$symbol
  }
  if (is.null(logx)) logx <- FALSE
  if (is.null(logy)) logy <- FALSE
  if (is.null(divXBy)) divXBy <- "None"
  if (is.null(divYBy)) divYBy <- "None"
  
 
  if (dimCol == "sampleNames") {
    myColors <- scols
  } else {
    myColors <- NULL
  }
  if (dimCol == "dbCluster") {
    myColors <- ccols
  }
  
  p1 <- plot2Dprojection(scEx_log, projections, g_id, featureData, geneNames,
                         geneNames2, dimX, dimY, clId, grpN, legend.position,
                         grpNs = grpNs, logx, logy, divXBy, divYBy, dimCol, colors = myColors
  )
  return(p1)
  
}
