suppressMessages(library(magrittr))
require(digest)
require(psychTools)
# printTimeEnd ----
printTimeEnd <- function(start.time, messtr) {
  end.time <- base::Sys.time()
  if (DEBUG) {
    cat(file = stderr(), paste("---", messtr, ":done", difftime(end.time, start.time, units = "min"), " min\n"))
  }
}


# some comments removed because they cause too much traffic ----
geneName2Index <- function(g_id, featureData) {
  # if (DEBUG) cat(file = stderr(), "geneName2Index started.\n")
  # start.time <- base::Sys.time()
  on.exit({
    # printTimeEnd(start.time, "geneName2Index")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "geneName2Index")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("geneName2Index", id = "geneName2Index", duration = NULL)
  }
  
  if (is.null(g_id)) {
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
  
  geneid <- unique(rownames(featureData[which(toupper(featureData$symbol) %in% toupper(g_id)), ]))
  
  return(geneid)
}

# updateProjectionsWithUmiCount ----
updateProjectionsWithUmiCount <- function(dimX, dimY, geneNames, geneNames2 = NULL, scEx, projections) {
  featureData <- rowData(scEx)
  # if ((dimY == "UmiCountPerGenes") | (dimX == "UmiCountPerGenes")) {
  geneNames <- geneName2Index(g_id = geneNames, featureData = featureData)
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


# append to heavyCalculations ----
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
  # in case it is called not from within schnapps
  if(!exists(".schnappsEnv")) {
    .schnappsEnv <- new.env(parent=emptyenv())
    .schnappsEnv$DEBUGSAVE = FALSE
  } else {
    if (is.null(.schnappsEnv$DEBUGSAVE)) {
      .schnappsEnv$DEBUGSAVE = FALSE
    }
  }
  
  
  # geneid <- geneName2Index(g_id, featureData)
  
  # if (length(geneid) == 0) {
  #   return(NULL)
  # }
  # if (length(geneid) == 1) {
  #   expression <- exprs(scEx_log)[geneid, ,drop=FALSE]
  # } else {
  # expression <- Matrix::colSums(assays(scEx_log)[[1]][geneid, rownames(projections), drop = FALSE])
  # }
  # validate(need(is.na(sum(expression)) != TRUE, ""))
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
  
  # histogram as y and cellDensity as color is not allowed
  
  if (dimY == "histogram") {
    if (!all(c(dimX, dimCol) %in% colnames(projections))) {
      return(NULL)
    }
  } else {
    # need to do proper checking of possibilities
    # removing dimCol for now
    if (!all(c(dimX, dimY) %in% colnames(projections))) {
      return(NULL)
    }
  }
  
  # projections <- cbind(projections, expression)
  # names(projections)[ncol(projections)] <- "exprs"
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/clusterPlot.RData", 
         list = c(ls(), "legend.position")
    )
    cat(file = stderr(), paste("plot2Dprojection saving done.\n"))
  }
  # load(file="~/SCHNAPPsDebug/clusterPlot.RData")
  if (DEBUG) {
    cat(file = stderr(), paste("output$sCA_dge_plot1:---", clId[1], "---\n"))
  }
  if ("dbCluster" %in% clId) {
    subsetData <- subset(projections, dbCluster %in% clId)
  } else {
    subsetData <- projections
  }
  #ensure that the highest values are plotted last.
  if (dimCol %in% colnames(subsetData)){
    if (is.numeric(subsetData[,dimCol])){
      subsetData <- subsetData[order(subsetData[,dimCol]),]
    }
  }
  # subsetData$dbCluster = factor(subsetData$dbCluster)
  # if there are more than 18 samples ggplot cannot handle different shapes and we ignore the
  # sample information
  if (length(as.numeric(as.factor(subsetData$sampleNames))) > 18) {
    subsetData$shape <- as.factor(1)
  } else {
    subsetData$shape <- as.numeric(as.factor(subsetData$sampleNames))
  }
  if (nrow(subsetData) == 0) {
    return(NULL)
  }
  # subsetData$shape = as.factor(1)
  gtitle <- paste(g_id, collapse = " ")
  if (nchar(gtitle) > 50) {
    gtitle <- paste(substr(gtitle, 1, 50), "...")
  }
  
  suppressMessages(require(plotly))
  f <- list(
    family = "Courier New, monospace",
    size = 18,
    color = "#7f7f7f"
  )
  if (divXBy != "None") {
    subsetData[, dimX] <- subsetData[, dimX] / subsetData[, divXBy]
  }
  if (divYBy != "None" & dimY != "histogram") {
    subsetData[, dimY] <- subsetData[, dimY] / subsetData[, divYBy]
  }
  
  typeX <- typeY <- "linear"
  if (logx) {
    typeX <- "log"
  }
  if (logy) {
    typeY <- "log"
  }
  if (is.factor(subsetData[, dimX]) | is.logical(subsetData[, dimX])) {
    typeX <- NULL
  }
  if (dimY != "histogram"){
    if (is.factor(subsetData[, dimY]) | is.logical(subsetData[, dimY])) {
      typeY <- NULL
    }
  } else {
    typeX = NULL
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
  if (dimX == "barcode") {
    subsetData$"__dimXorder" <- rank(subsetData[, dimY])
    dimX <- "__dimXorder"
    if (dimY == "histogram"){
      # Error message
      return(NULL)
    }
  }
  # save(file = "~/SCHNAPPsDebug/2dplot.RData", list = ls())
  # load("~/SCHNAPPsDebug/2dplot.RData")
  
  if (dimY != "histogram"){
    if (is.factor(subsetData[, dimX]) | is.logical(subsetData[, dimX])) {
      subsetData[, dimX] <- as.character(subsetData[, dimX])
    }
    if (is.factor(subsetData[, dimY]) | is.logical(subsetData[, dimY])) {
      subsetData[, dimY] <- as.character(subsetData[, dimY])
    }
  } else {
    # if (is.factor(subsetData[, dimX]) | is.logical(subsetData[, dimX])) {
    #   # barchart
    #   # subsetData[, dimX] <- as.character(subsetData[, dimX])
    # } else {
    # histogram
    
    if (is(subsetData[,dimCol], "factor")) {
      one_plot <- function(d) {
        add_histogram(x = ~d[,dimX])  
        # %>%
        #   add_annotations(
        #     ~unique(clarity), x = 0.5, y = 1, 
        #     xref = "paper", yref = "paper", showarrow = FALSE
        #   )
      }
      p1 <- subsetData %>% split(.[dimCol])
      
      # %>% lapply(one_plot)
      p = plot_ly(alpha = 0.6)
      for (le in 1:length(p1)) {
        if(nrow(p1[[le]])>0){
          
          p = add_histogram(p, x = p1[[le]][,dimX], name = levels(subsetData[,dimCol])[le])
          print(le)
        }
      }
      p = p %>% layout(barmode = "stack")
      p
      # plot_ly(alpha = 0.6) %>% lapply(p1, one_plot) %>%
      #   layout(barmode = "stack")
      # plot_ly(alpha = 0.6) %>% add_histogram(x=~p1[[3]][,dimX]) %>% add_histogram(x=~p1[[2]][,dimX]) %>%
      #   layout(barmode = "stack")
      
      # %>% p1 %>% layout(barmode = "overlay")
    } else {
      
      
      
      p <- plot_ly( x=~subsetData[, dimX], type = "histogram") %>%
        layout(
          xaxis = xAxis,
          yaxis = yAxis,
          title = gtitle,
          dragmode = "select"
        )
    }
    return (p)
    # %>%
    #   layout(yaxis=list(type='linear'))
    # }
    
  }
  
  if (dimCol == "cellDensity") {
    subsetData$cellDensity <- get_density(subsetData[,dimX], subsetData[,dimY], n = 100)
  }
  
  # dimCol = "Gene.count"
  # dimCol = "sampleNames"
  # subsetData$"__key__" = rownames(subsetData)
  
  p1 <- plotly::plot_ly(
    data = subsetData, source = "subset",
    key = rownames(subsetData)
  ) %>%
    add_trace(
      x = ~ get(dimX),
      y = ~ get(dimY),
      type = "scatter", mode = "markers",
      # text = ~ paste(1:nrow(subsetData), " ", rownames(subsetData), "<br />", subsetData$exprs),
      text = ~ paste(1:nrow(subsetData), " ", rownames(subsetData)),
      color = ~ get(dimCol),
      colors = colors,
      showlegend = TRUE
    ) %>%
    layout(
      xaxis = xAxis,
      yaxis = yAxis,
      title = gtitle,
      dragmode = "select"
    )
  
  
  if (is.factor(subsetData[, dimCol])) {
    
  } else {
    p1 <- colorbar(p1, title = dimCol)
  }
  
  selectedCells <- NULL
  if (length(grpN) > 0) {
    if (length(grpNs[rownames(subsetData), grpN] == "TRUE") > 0 & sum(grpNs[rownames(subsetData), grpN] == "TRUE", na.rm = TRUE) > 0) {
      grpNSub <- grpNs[rownames(subsetData), ]
      selectedCells <- rownames(grpNSub[grpNSub[, grpN] == "TRUE", ])
    }
  }
  if (!is.null(selectedCells)) {
    shape <- rep("a", nrow(subsetData))
    selRows <- which(rownames(subsetData) %in% selectedCells)
    shape[selRows] <- "b"
    x1 <- subsetData[selectedCells, dimX, drop = FALSE]
    y1 <- subsetData[selectedCells, dimY, drop = FALSE]
    p1 <- p1 %>%
      add_trace(
        data = subsetData[selectedCells, ],
        x = x1[, 1], y = y1[, 1],
        marker = list(
          color = rep("red", nrow(x1)),
          size = 5
        ),
        # text = ~ paste(rownames(subsetData[selectedCells, ]), "<br />", subsetData$exprs[selectedCells]),
        text = ~ paste(rownames(subsetData[selectedCells, ])),
        type = "scatter",
        mode = "markers",
        name = "selected",
        key = rownames(subsetData[selectedCells, ])
      )
  }
  # add density on y-axis if sorted by barcode
  if (dimX == "__dimXorder") {
    density <- stats::density(subsetData[, dimY], na.rm = T)
    library(MASS)
    fit <- fitdistr(subsetData[!is.na(subsetData[, dimY]), dimY], "normal", na.rm = T)
    hline <- function(y, dash) {
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y, line = list(dash = dash))
    }
    p1 <- p1 %>% layout(shapes = list(
      hline(fit$estimate[1], "dash"),
      # hline(fit$estimate[1] + fit$estimate[2]),
      # hline(fit$estimate[1] - fit$estimate[2]),
      hline(fit$estimate[1] + 3 * fit$estimate[2], "dot"),
      hline(fit$estimate[1] - 3 * fit$estimate[2], "dot")
    ))
  }
  p1
}


# functions should go in external file
# n_fun ----
n_fun <- function(x) {
  return(data.frame(y = -0.5, label = paste0(length(x), "\ncells")))
}
#' diffLRT ----
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

expMean <- function(x, normFactor = 1) {
  if (is.null(normFactor)){
    normFactor = 1
  }
  return(log(mean(exp(x/normFactor) - 1) + 1)*normFactor)
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
  # moreOptions <- input[[paste0(moduleName, "-moreOptions")]]
  colTree <- input[[paste0(moduleName, "-showColTree")]]
  
  if (is.null(heatmapData) | is.null(projections) | is.null(heatmapData$mat)) {
    return(NULL)
  }
  
  heatmapData$filename <- NULL
  
  # if (length(addColNames) > 0 & moreOptions) {
  if (length(addColNames) > 0) {
    heatmapData$annotation_col <- projections[rownames(heatmapData$annotation_col), addColNames, drop = FALSE]
  }
  # if (sum(orderColNames %in% colnames(projections)) > 0 & moreOptions) {
  if (sum(orderColNames %in% colnames(projections)) > 0) {
    heatmapData$cluster_cols <- FALSE
    colN <- rownames(dfOrder(projections, orderColNames))
    colN <- colN[colN %in% colnames(heatmapData$mat)]
    heatmapData$mat <- heatmapData$mat[, colN, drop = FALSE]
    # return()
  }
  # if (moreOptions) {
  heatmapData$cluster_cols <- colTree
  # }
  heatmapData$fontsize <- 14
  system.time(do.call(TRONCO::pheatmap, heatmapData))
}

# twoDplotFromModule ----
#' function to be used in markdown docs to ease the plotting of the clusterServer module
# TODO relies on reactive groupNames, should be a variable! Same goes for input$groupName!
twoDplotFromModule <- function(twoDData, moduleName, input, projections, g_id, legend.position = "none") {
  grpNs <- groupNames$namesDF
  grpN <- make.names(input$groupName, unique = TRUE)
  
  dimY <- input[[paste0(moduleName, "-dimension_y")]]
  dimX <- input[[paste0(moduleName, "-dimension_x")]]
  dimCol <- input[[paste0(moduleName, "-dimension_col")]]
  clId <- input[[paste0(moduleName, "-clusters")]]
  
  geneNames <- input[[paste0(moduleName, "-geneIds")]]
  geneNames2 <- input[[paste0(moduleName, "-geneIds2")]]
  logx <- input[[paste0(moduleName, "-logX")]]
  logy <- input[[paste0(moduleName, "-logY")]]
  divXBy <- input[[paste0(moduleName, "-divideXBy")]]
  divYBy <- input[[paste0(moduleName, "-divideYBy")]]
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

###### caching ----

# return values:    status = ["finished", "running", "error", "new" ]
#                   message = running message (or error)
#                   retVal = return value
checkShaCache <- function(moduleName = "traj_elpi_modules",
                          moduleParameters = list(
                            scEx_log, projections, TreeEPG,
                            elpimode, gene_sel, seed
                          )) {
  require(digest)
  retVal <- NULL
  message <- ""
  status <- "new"
  
  shaStr <- ""
  idStr <- paste0(moduleName, getshaStr(moduleParameters), collapse = "_")
  infile <- paste0("schnappsCache/", moduleName, "_", idStr, ".RData")
  if (!file.exists("schnappsCache")) {
    dir.create("schnappsCache")
  }
  if (!dir.exists("schnappsCache")) {
    error("cache directory is not a directory")
  }
  # nothing has been done so far
  if (!file.exists(infile)) {
    return(list(status = status, message = message, retVal = retVal))
  }
  # now we know that there is a result, check if it is the right one.
  message(paste("\nloading cache", idStr, "\n"))
  cp <- load(infile)
  if (!all(c("message", "retVal", "status") %in% cp)) {
    message <- "not all required values returned by cache function"
    return(list(status = status, message = message, retVal = retVal))
  }
  # check if number of objects is the same
  # moduleParameters
  # + message
  # + retVal
  # + status
  # + list(moduleParameters)
  if (length(cp) != 4) {
    # now we have problem: same sha but different input
  }
  # all values are set in writeShaCache
  return(list(status = status, message = message, retVal = retVal))
}

# writeShaCache ---
# write return values to cache directory
writeShaCache <- function(moduleName = "",
                          moduleParameters = list(),
                          retVal = NULL,
                          status = "error",
                          message = "") {
  require(digest)
  idStr <- paste0(moduleName, getshaStr(moduleParameters), collapse = "_")
  outfile <- paste0("schnappsCache/", moduleName, "_", idStr, ".RData")
  if (!file.exists("schnappsCache")) {
    dir.create("schnappsCache")
  }
  if (!dir.exists("schnappsCache")) {
    error("cache directory is not a directory")
  }
  save(file = outfile, list = c(
    "message", "retVal", "status", "moduleParameters"
  ))
}


getshaStr <- function(moduleParameters) {
  shaStr <- ""
  idx <- 1
  for (md in moduleParameters) {
    # print(class(md))
    # print(idx); idx=idx+1
    shaStr <- paste(
      shaStr,
      tryCatch(sha1(md, digits = 14),
               warning = function(x) {
                 # print("warning")
                 # print(x)
                 # print(idx)
                 return(
                   sha1(capture.output(str(md, vec.len = 40, digits.d = 14, nchar.max = 1400000, list.len = 100)))
                 )
               },
               error = function(x) {
                 if (class(md) == "SingleCellExperiment") {
                   return(sha1(as.matrix(assays(md)[[1]])))
                 } else {
                   print(idx)
                   return(
                     sha1(capture.output(str(md, vec.len = 40, digits.d = 14, nchar.max = 1400000, list.len = 100)))
                   )
                 }
               }
      )
    )
  }
  return(sha1(shaStr))
}

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}



# # recHistory ----
# # record history in env
# # needs pdftk https://www.pdflabs.com/tools/pdftk-server/
# # only save to history file if variable historyFile in schnappsEnv is set
# if (!all(c("pdftools", "gridExtra", "png") %in% rownames(installed.packages()))) {
#   recHistory <- function(...) {
#     return(NULL)
#   }
# } else {
#   require(pdftools)
#   recHistory <- function(name, plot1, envir = .schnappsEnv) {
#     if (!exists("historyFile", envir = envir)) {
#       return(NULL)
#     }
#     if (!exists("history", envir = .schnappsEnv)) {
#       .schnappsEnv$history <- list()
#     }
#     name <- paste(name, date())
#     tmpF <- tempfile(fileext = ".pdf")
#     cat(file = stderr(), paste0("history tmp File: ", tmpF, "\n"))
#     # save(file = "~/SCHNAPPsDebug/save2History2.RData", list = c(ls()))
#     # cp =load(file="~/SCHNAPPsDebug/save2History2.RData")
#     clP <- class(plot1)
#     cat(file = stderr(), paste0("class: ", clP[1], "\n"))
#     # here we create a PDF file for a given plot that is then combined later
#     created <- FALSE
#     switch(clP[1],
#       "plotly" = {
#         cat(file = stderr(), paste0("plotly\n"))
#         plot1 <- plot1 %>% layout(title = name)
#         if ("plotly" %in% class(plot1)) {
#           # requires orca bing installed (https://github.com/plotly/orca#installation)
#           withr::with_dir(dirname(tmpF), plotly::orca(p = plot1, file = basename(tmpF)))
#         }
#         created <- TRUE
#       },
#       "character" = {
#         # in case this is a link to a file:
#         cat(file = stderr(), paste0("character\n"))
#         if (file.exists(plot1)) {
#           if (tools::file_ext(plot1) == "png") {
#             pdf(tmpF)
#             img <- png::readPNG(plot1)
#             plot(1:2, type = "n")
#             rasterImage(img, 1.2, 1.27, 1.8, 1.73, interpolate = FALSE)
#             dev.off()
#           }
#           created <- TRUE
#         }
#       },
#       "datatables" = {
#         # # // this takes too long
#         # cat(file = stderr(), paste0("datatables\n"))
#         # save(file = "~/SCHNAPPsDebug/save2History2.RData", list = c(ls()))
#         # # cp =load(file="~/SCHNAPPsDebug/save2History2.RData")
#         #
#         # pdf(tmpF)
#         # if (nrow(img) > 20) {
#         #   maxrow = 20
#         # } else {
#         #   maxrow = nrow(plot1)
#         # }
#         # gridExtra::grid.table(img[maxrow],)
#         # dev.off()
#         # created = TRUE
#       }
#     )
# 
#     if (!created) {
#       return(FALSE)
#     }
# 
#     if (file.exists(.schnappsEnv$historyFile)) {
#       tmpF2 <- tempfile(fileext = ".pdf")
#       file.copy(.schnappsEnv$historyFile, tmpF2)
#       tryCatch(
#         pdf_combine(c(tmpF2, tmpF), output = .schnappsEnv$historyFile),
#         error = function(x) {
#           cat(file = stderr(), paste0("problem while combining PDF files:", x, "\n"))
#         }
#       )
#     } else {
#       file.copy(tmpF, .schnappsEnv$historyFile)
#     }
#     return(TRUE)
#     # pdf(file = tmpF,onefile = TRUE)
#     # ggsave(filename = tmpF, plot = plot1, device = pdf())
#     # dev.off()
#   }
# }

#' addColData
#' adds empty column to SingleCellExperiment object
#' the columns from scEx that are not present allScEx_log will be added to the latter one
#' return singleCellObject with added columns
addColData <- function(allScEx_log, scEx) {
  cd1 <- colnames(colData(scEx))
  cd2 <- colnames(colData(allScEx_log))
  
  for (cc in setdiff(cd1, cd2)) {
    lv <- NULL
    nCol <- NULL
    switch(class(colData(scEx)[, cc]),
           "factor" = {
             levels(colData(scEx)[, cc]) <- c(levels(colData(scEx)[, cc]), "NA")
             lv <- levels(colData(scEx)[, cc])
             nCol <- data.frame(factor(rep("NA", nrow(colData(allScEx_log))), levels = lv))
             colnames(nCol) <- cc
           },
           "character" = {
             nCol <- data.frame(cc = rep("NA", nrow(colData(allScEx_log))), stringsAsFactors = F)
             colnames(nCol) <- cc
           },
           "numeric" = {
             nCol <- data.frame(cc = rep(0, nrow(colData(allScEx_log))))
             colnames(nCol) <- cc
           }
    )
    colData(allScEx_log) <- cbind(colData(allScEx_log), nCol)
  }
  return(allScEx_log)
}

# ( = "updatetsneParameters",  = c("calculated_gQC_tsneDim", "calculated_gQC_tsnePerplexity",
#                                                                       "calculated_gQC_tsneTheta", "calculated_gQC_tsneSeed"))

# valuesChanged ----
valuesChanged <- function(parameters) {
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/valuesChanged.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/valuesChanged.RData')
  modified <- FALSE
  for (var in parameters) {
    oldVar <- paste0("calculated_", var)
    currVar <- var
    if (!exists(oldVar, envir = .schnappsEnv) | !exists(currVar, envir = .schnappsEnv)) {
      cat(file = stderr(), green("modified1\n"))
      cat(file = stderr(), paste("var:", var, "old: ", exists(oldVar, envir = .schnappsEnv), " new:", exists(currVar, envir = .schnappsEnv)))
      modified <- TRUE
      next()
    }
    if (is.null(get(oldVar, envir = .schnappsEnv)) | is.null(get(currVar, envir = .schnappsEnv))) {
      if (!(is.null(get(oldVar, envir = .schnappsEnv)) & is.null(get(currVar, envir = .schnappsEnv)))) {
        cat(file = stderr(), "modified2\n")
        modified <- TRUE
      }
    } else {
      if (is.na(get(oldVar, envir = .schnappsEnv)) | is.na(get(currVar, envir = .schnappsEnv))){
        if (!(is.na(get(oldVar, envir = .schnappsEnv)) & is.na(get(currVar, envir = .schnappsEnv)))){
          # browser()
          cat(file = stderr(), "modified3a\n")
          modified <- TRUE
        }
      } else
        if (!get(oldVar, envir = .schnappsEnv) == get(currVar, envir = .schnappsEnv)) {
          # browser()
          cat(file = stderr(), "modified3\n")
          cat(file = stderr(), oldVar)
          cat(file = stderr(), get(oldVar, envir = .schnappsEnv))
          cat(file = stderr(), get(currVar, envir = .schnappsEnv))
          modified <- TRUE
        }
    }
  }
  if (!modified) {
    cat(file = stderr(), "\n\nnot modified\n\n\n")
  } else {
    cat(file = stderr(), "\n\nmodified4\n\n\n")
  }
  return(modified)
}


# updateButtonColor ----
updateButtonColor <- function(buttonName, parameters) {
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/updateButtonColor.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/updateButtonColor.RData')
  modified <- valuesChanged(parameters)
  if (!modified) {
    removeClass(buttonName, "red")
    addClass(buttonName, "green")
  } else {
    removeClass(buttonName, "green")
    addClass(buttonName, "red")
  }
}

# setRedGreenButton ----
setRedGreenButton <- function(vars = list(), button = "") {
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/setRedGreenButton.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/setRedGreenButton.RData')
  for (v1 in vars) {
    # if (is.null(v1[2])) {
    #   assign(paste0("calculated_", v1[1]), "NULL", envir = .schnappsEnv)
    # }else {
    assign(paste0("calculated_", v1[1]), paste(v1[-1], collapse = "; "), envir = .schnappsEnv)
    # }
  }
  addClass(button, "green")
}

# setRedGreenButtonCurrent ----
setRedGreenButtonCurrent <- function(vars = list()) {
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/setRedGreenButtonCurrent.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/setRedGreenButtonCurrent.RData')
  for (v1 in vars) {
    # if (is.null(v1[2])) {
    #   assign(paste0("", v1[1]), "NULL", envir = .schnappsEnv)
    # }else {
    assign(paste0("", v1[1]), paste(v1[-1], collapse = "; "), envir = .schnappsEnv)
    # }
  }
}





reportFunction <- function(tmpPrjFile) {
  return(NULL)
}

# getReactEnv ----

getReactEnv <- function(DEBUG) {
  report.env <- new.env()
  # translate reactiveValues to lists
  # this way they can be saved
  rectVals <- c()
  isolate({
    for (var in c(names(globalenv()), names(parent.env(environment())))) {
      if (DEBUG) cat(file = stderr(), paste("var: ", var, "---", class(get(var))[1], "\n"))
      if (var == "reacativeReport") {
        next()
      }
      if (class(get(var))[1] == "reactivevalues") {
        if (DEBUG) cat(file = stderr(), paste("is reactiveValue: ", var, "\n"))
        rectVals <- c(rectVals, var)
        assign(var, reactiveValuesToList(get(var)), envir = report.env)
      } else if (class(get(var))[1] == "reactiveExpr") {
        if (DEBUG) {
          cat(
            file = stderr(),
            paste("is reactiveExpr: ", var, "--", class(get(var))[1], "\n")
          )
        }
        # if ( var == "coE_selctedCluster")
        # browser()
        rectVals <- c(rectVals, var)
        tempVar <- tryCatch(eval(parse(text = paste0(
          "\`", var, "\`()"
        ))),
        error = function(e) {
          browser()
          cat(file = stderr(), paste("error var", var, ":(", e, ")\n"))
          e
        },
        warning = function(e) {
          cat(file = stderr(), paste("warning with var", var, ":(", e, ")\n"))
          e
        }
        )
        assign(var, tempVar, envir = report.env)
        # for modules we have to take care of return values
        # this has to be done manually (for the moment)
        # and is only required for clusterServer
        if (class(report.env[[var]])[1] == "reactivevalues") {
          if (all(c("selectedCells") %in% names(report.env[[var]]))) {
            # cat(
            #   file = stderr(),
            #   paste(
            #     "is reactivevalues2: ",
            #     paste0(var, "-cluster"),
            #     "\n"
            #   )
            # )
            # if( paste0(var,"-cluster") == "coE_selctedCluster-cluster")
            #   browser()
            # assign(paste0(var, "-cluster"),
            # eval(report.env[[var]][["cluster"]]),
            # envir = report.env
            # )
            tempVar <- report.env[[var]][["selectedCells"]]
            assign(paste0(var, "-selectedCells"),
                   eval(parse(text = "tempVar()")),
                   envir = report.env
            )
          }
        }
      }
    }
    
    assign("input", reactiveValuesToList(get("input")), envir = report.env)
  })
  return(report.env)
}

# updateButtonUI = function(input, name, variables){
#   # browser()
#   if (is.null(input[[name]])){
#     cat(file = stderr(), paste("\ncreating button: ", name,"\n"))
#     return(renderUI({actionButton(inputId = name, label = "apply changes", width = '80%')}))
#   }
#
#   renderUI({
#     # inp <- isolate(reactiveValuesToList(get("input")))
#     # save(file = "~/SCHNAPPsDebug/render.RData", list = c(ls(), "name", "variables", ".schnappsEnv", "inp", "session"))
#     # cp =load(file = "~/SCHNAPPsDebug/render.RData")
#     # browser()
#     if (DEBUG) cat(file = stderr(), "updateButtonUI\n")
#     modified = FALSE
#     # input$updatetsneParameters
#     # input$gQC_tsneDim
#     # input$gQC_tsnePerplexity
#     # input$gQC_tsneTheta
#     # input$gQC_tsneSeed
#     # variables = c("gQC_tsneDim", "gQC_tsnePerplexity", "gQC_tsneTheta", "gQC_tsneSeed"  )
#
#     # create button if not exists
#
#     for (var in variables) {
#       oldVar = paste0("calculated_", var)
#       currVar = var
#       if (!exists(oldVar, envir = .schnappsEnv)  | !exists(currVar, envir = .schnappsEnv)) {
#         # cat(file = stderr(), "modified1\n")
#         modified = TRUE
#         next()
#       }
#       save(file = paste0("~/SCHNAPPsDebug/updateButtonUI", var,".RData"), list = c(
#         "var", "variables", ".schnappsEnv", "name"
#       ))
#       cat(file = stderr(), paste("\noldVar: ", oldVar,"\n"))
#       cat(file = stderr(), paste("noldVar: ", get(oldVar, envir = .schnappsEnv),"\n"))
#       cat(file = stderr(), paste("currVar: ", currVar ,"\n"))
#       cat(file = stderr(), paste("currVar: ", get(currVar, envir = .schnappsEnv),"\n"))
#       if(is.null(get(oldVar, envir = .schnappsEnv)) | is.null(get(currVar, envir = .schnappsEnv))) {
#         if (!(is.null(get(oldVar, envir = .schnappsEnv)) & is.null(get(currVar, envir = .schnappsEnv)))) {
#           modified = TRUE
#         }
#       } else {
#         if (!get(oldVar, envir = .schnappsEnv) == get(currVar, envir = .schnappsEnv)) {
#           # cat(file = stderr(), "modified12\n")
#           modified = TRUE
#         }
#       }
#       # observe("input$gQC_tsnePerplexity", quoted = T)
#     }
#     # ob = quote("input$gQC_tsnePerplexity")
#     # observe(ob , quoted = TRUE)
#     # observe(quote(paste0("input$",currVar)), quoted = F)
#     if (!modified) {
#       # cat(file = stderr(), "\n\ntsne not modified\n\n\n")
#       addClass(name, "green")
#       # updateActionButton(inputId = name, label = "apply changes", width = '80%',
#       #              style = "color: #fff; background-color: #00b300; border-color: #2e6da4")
#     } else {
#       # cat(file = stderr(), "\n\ntsne modified\n\n\n")
#       removeClass(name, "green")
#       # updateActionButton(inputId = name, label = "apply changes", width = '80%',
#       #              style = "color: #fff; background-color: #cc0000; border-color: #2e6da4")
#     }
#   })
# }


# add2history ----

add2history <- function(type, comment = "", input = input, ...) {
  if (!exists("historyPath", envir = .schnappsEnv)) {
    # if this variable is not set we are not saving
    return(NULL)
  }
  
  varnames <- lapply(substitute(list(...))[-1], deparse)
  arg <- list(...)
  if(is.null(arg[[1]])) return(NULL)
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/add2history.RData", list = c(ls()))
  }
  # cp =load(file='~/SCHNAPPsDebug/add2history.RData')
  if (type == "text") {
    cat(file = stderr(), paste0("history text: \n"))
    assign(names(varnames[1]), arg[1])
    line <- paste0(
      "\n", get(names(varnames[1])), "\n"
    )
    write(line, file = .schnappsEnv$historyFile, append = TRUE)
    
  }
  
  if (type == "save") {
    # browser()
    tfile <- tempfile(pattern = paste0(names(varnames[1]), "."), tmpdir = .schnappsEnv$historyPath, fileext = ".RData")
    assign(names(varnames[1]), arg[1])
    save(file = tfile, list = c(names(varnames[1]), "input"))
    # the load is commented out because it is not used at the moment and only takes time to load
    line <- paste0(
      "```{R}\n#load ", names(varnames[1]), "\n#load(file = \"", basename(tfile),
      "\")\n```\n"
    )
    write(line, file = .schnappsEnv$historyFile, append = TRUE)
  }
  
  if (type == "renderPlotly") {
    
    tfile <- tempfile(pattern = paste0(names(varnames[1]), "."), tmpdir = ".", fileext = ".png")
    assign(names(varnames[1]), arg[1])
    # save(file = tfile, list = c(names(varnames[1])))
    tryCatch({
      # orca(plotData$plotData, file = tfile, format = "png")
      withr::with_dir(normalizePath(.schnappsEnv$historyPath), orca(plotData$plotData, file = tfile, format = "png"))
      line <- paste0(
        # "```{R}\n#load ", names(varnames[1]), "\nload(file = \"", basename(tfile),
        # "\")\nhtmltools::tagList(", names(varnames[1]), ")\n```\n",
        "\n![](",basename(tfile),")\n\n"
      )
      write(line, file = .schnappsEnv$historyFile, append = TRUE)
    })
  }
  
  if (type == "tronco") {
    # browser()
    tfile <- tempfile(pattern = paste0(names(varnames[1]), "."), tmpdir = .schnappsEnv$historyPath, fileext = ".RData")
    assign(names(varnames[1]), arg[[1]])
    # report.env <- getReactEnv(DEBUG = .schnappsEnv$DEBUG)
    save(file = tfile, list = c(names(varnames[1]), "input"))
    
    line <- paste0(
      "```{R}\n#load ", names(varnames[1]), "\nload(file = \"", basename(tfile),"\")\n",
      "\n", names(varnames[1]) ,"$filename <- NULL \n",
      "\ndo.call(TRONCO::pheatmap, ", names(varnames[1]), ")\n```\n"
    )
    write(line, file = .schnappsEnv$historyFile, append = TRUE)
  }
  
  if (type == "renderPlot") {
    tfile <- tempfile(pattern = paste0(names(varnames[1]), "."), tmpdir = .schnappsEnv$historyPath, fileext = ".RData")
    assign(names(varnames[1]), arg[[1]])
    # report.env <- getReactEnv(DEBUG = .schnappsEnv$DEBUG)
    save(file = tfile, list = c(names(varnames[1]), "input"))
    
    line <- paste0(
      "```{R}\n#load ", names(varnames[1]), "\nload(file = \"", basename(tfile),"\")\n",
      "\n", names(varnames[1]), "\n```\n"
    )
    write(line, file = .schnappsEnv$historyFile, append = TRUE)
    
  }
  
  if (type == "renderDT") {
    tfile <- tempfile(pattern = paste0(names(varnames[1]), "."), tmpdir = .schnappsEnv$historyPath, fileext = ".RData")
    assign(names(varnames[1]), arg[[1]])
    # report.env <- getReactEnv(DEBUG = .schnappsEnv$DEBUG)
    save(file = tfile, list = c(names(varnames[1]), "input"))
    
    line <- paste0(
      "```{R}\n#load ", names(varnames[1]), "\nload(file = \"", basename(tfile),"\")\n",
      "\n", names(varnames[1]), "\n```\n"
    )
    write(line, file = .schnappsEnv$historyFile, append = TRUE)
    
  }
}




# Get density of points in 2 dimensions. ----
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# ns <- session$ns
# heatmapData <- pheatmapList()
# addColNames <- input$ColNames
# orderColNames <- input$orderNames
# colTree <- input$showColTree
# scale <- input$normRow
# myns <- ns("pHeatMap")
# save2History <- input$save2History
# pWidth <- input$heatmapWidth
# pHeight <- input$heatmapHeight
# colPal <- input$colPal
# minMaxVal <- input$heatmapMinMaxValue
# # maxVal <- input$heatmapMaxValue
# 
# proje <- projections()



heatmapModuleFunction <- function(
  heatmapData = NULL,
  addColNames = "sampleNames",
  orderColNames = c(), 
  colTree = FALSE,
  scale = "none",
  pWidth = 300, 
  pHeight = 900,
  colPal= "none",
  minMaxVal = c(min(heatmapData$mat), max(heatmapData$mat)),
  proje = NULL,
  outfile = NULL
) {
  
  
  
  if (is.null(heatmapData) | is.null(proje) | is.null(heatmapData$mat)) {
    
    return(list(
      src = "empty.png",
      contentType = "image/png",
      width = 96,
      height = 96,
      alt = "pHeatMapPlot should be here (null)"
    ))
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/heatmapModuleFunction.RData", list = c(ls()))
  }
  # load(file = "~/SCHNAPPsDebug/heatmapModuleFunction.RData")
  
  if (is.null(pWidth)) {
    pWidth <- 800
    pHeight <- 300
  }
  if (!is(heatmapData$mat, "matrix") & !is(heatmapData$mat, "Matrix")) {
    cat(file = stderr(), "!!!!! output$pHeatMapModule:mat is not a matrix\n")
    heatmapData$mat <- as.matrix(t(heatmapData$mat))
  }
  if (is(heatmapData$mat, "sparseMatrix")) {
    heatmapData$mat = as.matrix(heatmapData$mat)
  }
  heatmapData$mat[is.nan(heatmapData$mat)] <- 0.0
  if (is.null(scale)) {
    heatmapData$scale <- "none"
  } else {
    heatmapData$scale <- scale
    if (scale == "row") {
      rmRows = -which(apply(heatmapData$mat, 1, sd) == 0)
      if(length(rmRows) > 1)
        heatmapData$mat <- heatmapData$mat[rmRows, ]
    }
    if (scale == "column") {
      rmCols = -which(apply(heatmapData$mat, 2, sd) == 0)
      if (length(rmCols) > 1)
        heatmapData$mat <- heatmapData$mat[, rmCols]
    }
  }
  
  heatmapData$filename <- outfile
  # heatmapData$filename = NULL
  # if (length(addColNames) > 0 & moreOptions) {
  if (length(addColNames) > 0) {
    heatmapData$annotation_col <- proje[rownames(heatmapData$annotation_col), addColNames, drop = FALSE]
  }
  # if (length(addColNames) > 0 & moreOptions) {
  if (length(addColNames) > 0) {
    heatmapData$annotation_col <- proje[rownames(heatmapData$annotation_col), addColNames, drop = FALSE]
  }
  # if (sum(orderColNames %in% colnames(proje)) > 0 & moreOptions) {
  if (sum(orderColNames %in% colnames(proje)) > 0) {
    heatmapData$cluster_cols <- FALSE
    colN <- rownames(psychTools::dfOrder(proje, orderColNames))
    matN <- colnames(heatmapData$mat)
    if (is.null(matN)) {
      matN <- names(heatmapData$mat)
    }
    colN <- colN[colN %in% matN]
    heatmapData$mat <- heatmapData$mat[, colN, drop = FALSE]
    # return()
  }
  # if (moreOptions) {
  heatmapData$cluster_cols <- colTree
  # }
  # orgMat = heatmapData$mat
  
  # heatmapData$mat = orgMat
  # system.time(do.call(pheatmap, heatmapData))
  # heatmapData$mat = as(orgMat, "dgTMatrix")
  heatmapData$fontsize <- 14
  # heatmapData$fontsize_row = 18
  # heatmapData$filename=NULL
  if (nrow(heatmapData$mat) > 1000) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "more than 1000 row in heatmap. This can be very slow to display. Only showing first 1000 rows",
        id = "pHeatMapPlotWARNING",
        type = "warning",
        duration = 20
      )
    }
    heatmapData$mat <- heatmapData$mat[1:1000, ]
    heatmapData$gaps_row <- heatmapData$gaps_row[heatmapData$gaps_row < 1000]
  }
  if (nrow(heatmapData$mat) < 2) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "Less than two rows to display",
        id = "pHeatMapPlotWARNING",
        type = "warning",
        duration = 20
      )
    }
    return(list(
      src = "empty.png",
      contentType = "image/png",
      width = 96,
      height = 96,
      alt = "pHeatMapPlot should be here (no rows)"
    ))
  }
  heatmapData$width <- pWidth / 72
  heatmapData$height <- pHeight / 72
  
  if (colPal == "none") {
    # use the supplied colors
  } else {
    heatmapData$color <- colorRampPalette(rev(brewer.pal(
      n = 7, name =
        colPal
    )))(100)
  }
  
  heatmapData$mat[heatmapData$mat <= minMaxVal[1]] = minMaxVal[1]
  heatmapData$mat[heatmapData$mat >= minMaxVal[2]] = minMaxVal[2]
  
  
  do.call(TRONCO::pheatmap, heatmapData)
  
  
  # pixelratio <- session$clientData$pixelratio
  # if (is.null(pixelratio)) pixelratio <- 1
  # width <- session$clientData$output_plot_width
  # height <- session$clientData$output_plot_height
  # if (is.null(width)) {
  #   width <- 96 * 7
  # } # 7x7 inch output
  # if (is.null(height)) {
  #   height <- 96 * 7
  # }
  outfilePH <- outfile
  return(list(
    src = outfilePH,
    contentType = "image/png",
    width = paste0(pWidth, "px"),
    height = paste0(pHeight, "px"),
    alt = "heatmap should be here"
  ))
}

