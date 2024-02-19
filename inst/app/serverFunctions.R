if (!exists(".schnappsEnv")) {
  .schnappsEnv <- new.env(parent = emptyenv())
  .schnappsEnv$DEBUG = TRUE
}

if (.schnappsEnv$DEBUG) {
  cat(file = stderr(), "started severFunctions.R\n")
}

suppressMessages(library(magrittr))
suppressMessages(require(digest))
suppressMessages(require(psychTools))
suppressMessages(require(tidyr))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(InteractiveComplexHeatmap))
library(dendsort)
library(MASS)


### Try catch from extended examples ----
#' tryCatch with Warning Extraction
#' 
#' A wrapper around the tryCatch function that captures both the value and warnings, allowing for the extraction of warning messages.
#' 
#' @param expr The expression to be evaluated within the tryCatch block.
#' 
#' @return A list with two elements: "value" containing the result of the expression evaluation (or the error if an error occurs),
#' and "warning" containing the warning message (if any) or NULL if no warning occurred.
#' 
#' @export tryCatch.W.E
#' 

tryCatch.W.E <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}

# Print Time End ----
#' Print Time End ----
#' 
#' Prints the elapsed time and a completion message based on the provided start time.
#' 
#' @param start.time The starting time, typically obtained using `Sys.time()` before the execution of a task.
#' @param messtr A character string containing the completion message.
#' 
#' @return The function prints the elapsed time and a completion message to the standard error stream.
#' 
#' @export printTimeEnd
#' 

printTimeEnd <- function(start.time, messtr) {
  require(hms)
  end.time <- base::Sys.time()
  if (.schnappsEnv$DEBUG & !is.null(start.time)) {
    duration = difftime(end.time, start.time, units = "min")
    if(.schnappsEnv$DEBUGSAVE){
      save(file = "~/SCHNAPPsDebug/printTimeEnd.RData",
           list = c(ls())
      )
    }
    # cp = load("~/SCHNAPPsDebug/printTimeEnd.RData")
    cat(file = stderr(), paste("---", hms::round_hms(hms::as_hms(duration),0.25), "--- done:", messtr, "\n"))
  }
}
#' Plot High Expression Genes
#' 
#' Generates a plot of the highest expression genes using the `plotHighestExprs` function from the 'scater' package.
#' 
#' @param scaterReads A 'SingleCellExperiment' object or any object with a 'counts' assay containing expression values.
#' @param n The number of highest expression genes to plot.
#' @param scols A vector of colors for the cells.
#' 
#' @return A ggplot2 plot displaying the highest expression genes.
#' 
#' @export pltHighExp
#' 
pltHighExp <- function( scaterReads, n, scols) {
  # since we are saving without the environment, we need require scater for storing
  require(scater)
  p1 <- scater::plotHighestExprs(scaterReads, colour_cells_by = "sampleNames", n = n) + scale_color_manual(values= scols)
  p1
}

#' Convert gene names to indices in featureData
#'
#' This function takes in a vector of gene names and a featureData object and 
#' returns the indices of the gene names in the featureData object. It performs 
#' case-insensitive comparisons to match the gene names. If any of the gene names 
#' are not found in the featureData, a warning notification is displayed.
#'
#' @param g_id A comma separated string of gene names to be converted to indices.
#' @param featureData A data frame or matrix containing gene information, with gene 
#'   names stored in the 'symbol' column.
#'
#' @return A vector of indices corresponding to the input gene names in the 
#'   featureData object.
#'
#' @examples
#' # Create a featureData object 'features' with gene names in the 'symbol' column
#' features <- data.frame(symbol = c("GeneA", "GeneB", "GeneC"))
#' 
#' # Convert gene names to indices
#' gene_indices <- geneName2Index(c("GeneA", "GeneC"), features)
#' gene_indices
#'
#' @export
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
    }else{
      warning(paste("following genes were not found", notFound, collapse = " "))
    }
  }
  
  
  # we need to compare case insensitive and then use case-sensitive rownames in the order of g_id
  featureIdx = seq(nrow(featureData))
  names(featureIdx) = toupper(featureData$symbol)
  geneid <- rownames(featureData)[featureIdx[intersect(toupper(g_id) , toupper(featureData$symbol) )]]
  # unique(rownames(featureData[which(toupper(featureData$symbol) %in% toupper(g_id)), ]))
  
  return(geneid)
}

# updateProjectionsWithUmiCount ----
#' Update Projections with UMI Count
#'
#' This function updates the `projections` object with UMI count information 
#' for specified genes from a `SingleCellExperiment` object. It handles both 
#' single genes and multiple genes, calculating the sum of UMI counts across 
#' multiple genes if necessary.
#'
#' @param geneNames A vector of gene names or indices for which UMI counts 
#'                  are to be retrieved.
#' @param geneNames2 An optional second vector of gene names or indices for 
#'                   an alternative set of UMI counts.
#' @param scEx A `SingleCellExperiment` object containing the UMI count data.
#' @param projections A list or data structure where the UMI count information 
#'                    will be stored.
#'
#' @return Returns the updated `projections` object with added UMI count data.
#'
#' @examples
#' # Assuming scEx is a valid SingleCellExperiment object 
#' # and projections is a valid projections object
#' gene_list = c("Gene1", "Gene2")
#' updateProjectionsWithUmiCount(gene_list, NULL, scEx, projections)
#'
#' @export
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Matrix colSums
updateProjectionsWithUmiCount <- function(
    geneNames, geneNames2 = NULL, scEx, projections) {
  featureData <- rowData(scEx)
  geneNames <- geneName2Index(g_id = geneNames, featureData = featureData)
  if ((length(geneNames) > 0) && (length(geneNames[[1]]) > 0)) {
    #TODO add drop=F to remove if statement
    if (length(geneNames) == 1) {
      projections$UmiCountPerGenes <- assays(scEx)[[1]][geneNames, ]
    } else {
      projections$UmiCountPerGenes <- Matrix::colSums(assays(scEx)[[1]][geneNames, ])
    }
  } else {
    projections$UmiCountPerGenes <- 0
  }
  geneNames <- geneName2Index(geneNames2, featureData)
  if ((length(geneNames) > 0) && (length(geneNames[[1]]) > 0)) {
    if (length(geneNames) == 1) {
      projections$UmiCountPerGenes2 <- assays(scEx)[[1]][geneNames, ]
    } else {
      projections$UmiCountPerGenes2 <- Matrix::colSums(assays(scEx)[[1]][geneNames, ])
    }
  } else {
    projections$UmiCountPerGenes2 <- 0
  }
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

# plot2Dprojection ----------------
# used in moduleServer and reports
plot2Dprojection <- function(scEx_log, projections, g_id, featureData,
                             geneNames, geneNames2, dimX, dimY, clId, grpN, legend.position, grpNs,
                             logx = FALSE, logy = FALSE, divXBy = "None", divYBy = "None", dimCol = "sampleNames",
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
  # 
  # In the module the input table is coming from the projections (tData) that is a parameter
  # to the module. Since this can be a subset we need to subset the scEx as well.
  if (is.null(projections)) return(NULL)
  projections <- updateProjectionsWithUmiCount(
    geneNames = geneNames,
    geneNames2 = geneNames2,
    scEx = scEx_log[, rownames(projections)], 
    projections = projections
  )
  
  # histogram as y and cellDensity as color is not allowed
  if (dimCol == "") return(NULL)
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
  # cp = load(file="~/SCHNAPPsDebug/clusterPlot.RData")
  if (DEBUG) {
    cat(file = stderr(), paste("output$sCA_dge_plot1:---", clId[1], "---\n"))
  }
  if ("dbCluster" %in% clId) {
    subsetData <- subset(projections, dbCluster %in% clId)
  } else {
    subsetData <- projections
  }
  
  # Clean data from NA
  for (colNa in c(dimX, dimY, dimCol)) {
    if(! colNa %in% colnames(subsetData)) next()
    naRows = which(is.na(subsetData[,colNa]))
    if (length(naRows) > 0) {
      if (is.factor(subsetData[,colNa])) {
        if (! "NA" %in% levels(subsetData[,colNa])){
          warning(paste("we found NA in factor", colNa,"\n"))
          levels(subsetData[,colNa]) = c(levels(subsetData[,colNa]), "NA")
          subsetData[naRows,colNa] = "NA"
        }
      } else {
        if (is.character(subsetData[,colNa])) {
          warning("we setting NA to character\n")
          subsetData[naRows,colNa] = "NA"
        } else {
          if (is.numeric(subsetData[,colNa])) { # remove lines
            warning("we are removing rows due to NA\n")
            subsetData = subsetData[-naRows,]
          }
        }
      }
    }
  }
  #ensure that the highest values are plotted last.
  if (dimCol %in% colnames(subsetData)){
    if (is.numeric(subsetData[,dimCol])){
      subsetData <- subsetData[order(subsetData[,dimCol]),]
    }
  }
  # # subsetData$dbCluster = factor(subsetData$dbCluster)
  # # if there are more than 18 samples ggplot cannot handle different shapes and we ignore the
  # # sample information
  # if (length(as.numeric(as.factor(subsetData$sampleNames))) > 18) {
  #   subsetData$shape <- as.factor(1)
  # } else {
  #   subsetData$shape <- as.numeric(as.factor(subsetData$sampleNames))
  # }
  if (nrow(subsetData) == 0) {
    return(NULL)
  }
  # subsetData$shape = as.factor(1)
  gtitle <- paste(g_id, collapse = " ")
  if (nchar(gtitle) > 50) {
    gtitle <- paste(substr(gtitle, 1, 50), "...")
    gtitle <- ""
  }
  
  suppressMessages(require(plotly))
  f <- list(
    family = "Courier New, monospace",
    size = 18,
    color = "#7f7f7f"
  )
  if (divXBy != "None") {
    subsetData[, dimX] <- subsetData[, dimX] / subsetData[, divXBy]
    subsetData = subsetData[!is.infinite(subsetData[, dimX]), ]
  }
  if (divYBy != "None" & divYBy != "normByCol" & dimY != "histogram") {
    subsetData[, dimY] <- subsetData[, dimY] / subsetData[, divYBy]
    subsetData = subsetData[!is.infinite(subsetData[, dimY]), ]
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
  dimxFact=FALSE
  dimyFact=FALSE
  dimyNum=FALSE
  if (dimY != "histogram"){
    if (is.numeric(subsetData[, dimY]) ) {
      dimyNum = TRUE
    }
    if (is.factor(subsetData[, dimX]) | is.logical(subsetData[, dimX])) {
      subsetData[, dimX] <- as.character(subsetData[, dimX])
      dimxFact = TRUE
    }
    if (is.factor(subsetData[, dimY]) | is.logical(subsetData[, dimY])) {
      subsetData[, dimY] <- as.character(subsetData[, dimY])
      dimyFact = TRUE
    }
  } else {
    # if (is.factor(subsetData[, dimX]) | is.logical(subsetData[, dimX])) {
    #   # barchart
    #   # subsetData[, dimX] <- as.character(subsetData[, dimX])
    # } else {
    
    #
    # histogram
    #
    
    if (is(subsetData[,dimCol], "factor")) {
      # Special case: histogram with color as factor and divYBy == "normByCol"
      histnorm = NULL
      barmode = "stack"
      if (divYBy == "normByCol" & dimY == "histogram"){
        histnorm = "probability"
        barmode = "group"
      }
      
      p1 <- subsetData %>% split(.[dimCol])
      
      # %>% lapply(one_plot)
      p = plot_ly(alpha = 0.6)
      for (le in 1:length(p1)) {
        if(nrow(p1[[le]])>0){
          mcol = NULL
          if (!is.null(colors)) {
            mcol = list(color = colors[names(p1)[le]])
          }
          p = add_histogram(p, x = p1[[le]][,dimX], 
                            name = levels(subsetData[,dimCol])[le],
                            marker = mcol,
                            histnorm = histnorm)
          print(le)
        }
      }
      p = p %>% layout(barmode = barmode)
      # p
      return (p)
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
  
  # 
  # violoin plot (X = factor, Y = numeric)
  # if color factor with two levels split
  # 
  
  if (dimxFact & dimyNum) {
    p1 <- plotly::plot_ly(
      data = subsetData, 
      source = "subset",
      key = rownames(subsetData)
    )
    if (is.factor(subsetData[,dimCol])) {
      p1 = subsetData %>% 
        plotly::plot_ly( 
          source = "subset",
          type = 'violin') 
      for ( lv in levels(subsetData[,dimCol])) {
        if (lv %in% names(colors)){
          col = I(colors[lv])
        } else {
          col = NULL
        }
        p1 = p1 %>%
          add_trace(
            x = formula(paste0("~get(dimX)[subsetData[,dimCol] == \"", lv, "\"]")),
            y = formula(paste0("~get(dimY)[subsetData[,dimCol] == \"", lv, "\"]")),
            # key = rownames(formula(paste0("~get(dimCol)[subsetData[,dimCol] == \"", lv, "\"]"))),
            name = lv,
            points = 'all',
            pointpos = 0,
            scalegroup = lv,
            legendgroup = lv,
            # alignmentgroup = lv,
            marker = list(size = 1),
            unselected = list ( 
              marker = list(size = 1)),
            selected = list ( 
              marker = list(size = 6)),
            box = list(
              visible = T
            ),
            meanline = list(
              visible = T
            ),
            color = col
          )
        # p1
        
      }
      p1 = p1 %>% layout(
        xaxis = xAxis,
        yaxis = yAxis,
        title = gtitle,
        dragmode = "select",
        violinmode = 'group'
      ) 
      # p1
      return (p1)
    }
    # ignore color
    p1 = p1 %>%
      add_trace(
        x = ~ get(dimX),
        y = ~ get(dimY),
        type = "violin",
        # text = ~ paste(1:nrow(subsetData), " ", rownames(subsetData), "<br />", subsetData$exprs),
        text = ~ paste(1:nrow(subsetData), " ", rownames(subsetData)),
        box = list(
          visible = T
        ),
        meanline = list(
          visible = T
        ),
        showlegend = TRUE
      ) %>%
      layout(
        xaxis = xAxis,
        yaxis = yAxis,
        title = gtitle,
        dragmode = "select"
      )
    return(p1)
  }
  
  
  if (dimCol == "cellDensity") {
    subsetData$cellDensity <- get_density(subsetData[,dimX], subsetData[,dimY], n = 100)
  }
  
  # dimCol = "Gene.count"
  # dimCol = "sampleNames"
  # subsetData$"__key__" = rownames(subsetData)
  markerInfo = NULL
  
  maxPointSize = 10
  if (dimxFact & dimyFact) {
    tb = as.data.frame(table(subsetData[,c(dimX, dimY)]))
    subsetData$rownames = rownames(subsetData)
    subsetData = merge(y=tb, x= subsetData, by=c(dimX, dimY))
    rownames(subsetData) =  subsetData$rownames
    # subsetData$dotSize = subsetData[,dimX] == tb[,dimX] & subsetData[,dimY] == tb[,dimY]
    # marker = list(size = ~Gap, opacity = 0.5))
    markerInfo = list(size=~Freq, sizes = c(100, 500))
    maxPointSize = 50
    textDisp = apply(subsetData[,c(dimX, dimY, "Freq")], 1, FUN = function(x) glue::glue("({x[1]}, {x[2]}): {x[3]}"))
  } else {
    subsetData$Freq = 1
    textDisp = paste(1:nrow(subsetData), " ", rownames(subsetData))
  }
  
  
  p1 <-
    plotly::plot_ly(
      data = subsetData, source = "subset",
      key = rownames(subsetData)
    ) %>%
    add_trace(
      x = ~ get(dimX),
      y = ~ get(dimY),
      type = "scatter", mode = "markers",
      # text = ~ paste(1:nrow(subsetData), " ", rownames(subsetData), "<br />", subsetData$exprs),
      text = ~ textDisp,
      color = ~ get(dimCol),
      colors = colors,
      size = ~ Freq,
      sizes = c(1,maxPointSize),
      marker = list(sizemod = 'diameter'),
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
    if (length(grpNs[rownames(subsetData), grpN] == "TRUE") > 0 & 
        sum(grpNs[rownames(subsetData), grpN] == "TRUE", na.rm = TRUE) > 0) {
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
    if(is.logical(x1)) {
      x1 = as.factor(x1)
    }
    if(is.logical(y1)) {
      y1 = as.factor(y1)
    }
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
    if (!is.numeric(subsetData[, dimY])) {
      cat(file=stderr(), paste("\nWARNING: ", dimY, "is not numeric\n\n"))
      return(NULL)
    }
    density <- stats::density(subsetData[, dimY], na.rm = T)
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
  # parallel done
  p_val <- bplapply(genes.use, function(x) diffLRT(as.numeric(expression[x, cells.1]), as.numeric(expression[x, cells.2]))) %>% unlist()
  to.return <- data.frame(p_val, row.names = genes.use)
  return(to.return)
}


heatmapPlotFromModule <- function(heatmapData, moduleName, input, projections) {
  addColNames <- input[[paste0(moduleName, "-ColNames")]]
  orderColNames <- input[[paste0(moduleName, "-orderNames")]]
  # moreOptions <- input[[paste0(moduleName, "-moreOptions")]]
  colTree <- FALSE
  if (input[[paste0(moduleName, "-sortingCols")]] == "dendrogram")
    colTree <- TRUE
  
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
  # scols <- projectionColors$sampleNames
  # ccols <- projectionColors$dbCluster
  scols <- projectionColors[["sampleNames"]]
  ccols <- projectionColors[["dbCluster"]]
  
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
#           # https://github.com/plotly/orca/releases
#           # https://github.com/plotly/orca
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
           },
           "integer" = {
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
      # cat(file = stderr(), green("modified1\n"))
      # cat(file = stderr(), paste("var:", var, "old: ", exists(oldVar, envir = .schnappsEnv), " new:", exists(currVar, envir = .schnappsEnv)))
      modified <- TRUE
      next()
    }
    if (is.null(get(oldVar, envir = .schnappsEnv)) | is.null(get(currVar, envir = .schnappsEnv))) {
      if (!(is.null(get(oldVar, envir = .schnappsEnv)) & is.null(get(currVar, envir = .schnappsEnv)))) {
        # cat(file = stderr(), "modified2\n")
        modified <- TRUE
      }
    } else {
      if (is.na(get(oldVar, envir = .schnappsEnv)) | is.na(get(currVar, envir = .schnappsEnv))){
        if (!(is.na(get(oldVar, envir = .schnappsEnv)) & is.na(get(currVar, envir = .schnappsEnv)))){
          # deepDebug()
          # cat(file = stderr(), "modified3a\n")
          modified <- TRUE
        }
      } else
        if (!get(oldVar, envir = .schnappsEnv) == get(currVar, envir = .schnappsEnv)) {
          # deepDebug()
          # cat(file = stderr(), "modified3\n")
          # cat(file = stderr(), oldVar)
          # cat(file = stderr(), get(oldVar, envir = .schnappsEnv))
          # cat(file = stderr(), get(currVar, envir = .schnappsEnv))
          modified <- TRUE
        }
    }
  }
  # if (!modified) {
  #   cat(file = stderr(), "\n\nnot modified\n\n\n")
  # } else {
  #   cat(file = stderr(), "\n\nmodified4\n\n\n")
  # }
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
    shinyjs::removeClass(buttonName, "red")
    shinyjs::addClass(buttonName, "green")
  } else {
    shinyjs::removeClass(buttonName, "green")
    shinyjs::addClass(buttonName, "red")
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
  shinyjs::addClass(button, "green")
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
        # deepDebug()
        rectVals <- c(rectVals, var)
        tempVar <- tryCatch(eval(parse(text = paste0(
          "\`", var, "\`()"
        ))),
        error = function(e) {
          # deepDebug()
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
            #   deepDebug()
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
#   # deepDebug()
#   if (is.null(input[[name]])){
#     cat(file = stderr(), paste("\ncreating button: ", name,"\n"))
#     return(renderUI({actionButton(inputId = name, label = "apply changes", width = '80%')}))
#   }
#
#   renderUI({
#     # inp <- isolate(reactiveValuesToList(get("input")))
#     # save(file = "~/SCHNAPPsDebug/render.RData", list = c(ls(), "name", "variables", ".schnappsEnv", "inp", "session"))
#     # cp =load(file = "~/SCHNAPPsDebug/render.RData")
#     # deepDebug()
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
  if (DEBUG) cat(file = stderr(), paste0("add2history: ", type, "\n"))
  if (!base::exists("historyPath", envir = .schnappsEnv)) {
    # if this variable is not set we are not saving
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "history path not set, not saving",
        id = "noHistoryWARNING",
        type = "warning",
        duration = 20
      )
    }
    return(NULL)
  }
  if (is.null(.schnappsEnv$historyPath)) {
    # if this variable is not set we are not saving
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "history path NULL, not saving",
        id = "noHistoryWARNING",
        type = "warning",
        duration = 20
      )
    }
    return(NULL)
  }
  # browser()
  
  # new feature to prevent saving during start-up
  # modules don't have access to the global input...
  if("startingUp" %in% names(input)){
    if(input$startingUp=="1") {
      # browser()
      return(NULL)
    } 
  }
  # deepDebug()
  sessionI = sessionInfo()
  defaultValues = .schnappsEnv$defaultValues
  dvFile = paste0(.schnappsEnv$historyPath, "/defaultValues.RData")
  if (DEBUG) cat(file = stderr(), paste0("add2history: saving default Values to", dvFile, "\n"))
  save(file=dvFile, list = c("defaultValues"))
  # browser()
  varnames <- lapply(substitute(list(...))[-1], deparse)
  arg <- list(...)
  inp = input
  schnappsEnv = .schnappsEnv
  if(is.null(arg[[1]])) {
    if (DEBUG) cat(file = stderr(), paste0("add2history: arg[[1]] is null\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/add2history.RData", list = c(ls(), ".schnappsEnv"), compress = F)
  }
  # cp =load(file='~/SCHNAPPsDebug/add2history.RData')
  if (type == "text") {
    # cat(file = stderr(), paste0("history text: \n"))
    assign(names(varnames[1]), arg[1])
    line <- paste0(
      "\n", get(names(varnames[1])), "\n"
    )
    write(line, file = .schnappsEnv$historyFile, append = TRUE)
    if (DEBUG) cat(file = stderr(), paste0("add2history: saving text to", .schnappsEnv$historyFile, "\n"))
    
  }
  
  # deepDebug()
  if (type == "save") {
    tfile <- tempfile(pattern = paste0(names(varnames[1]), "."), tmpdir = .schnappsEnv$historyPath, fileext = ".RData")
    if (file.exists(tfile)){cat(file = stderr(), paste("file exists, should not happen", tfile, "\n\n"))}
    assign(names(varnames[1]), arg[1])
    save(file = tfile, list = c(names(varnames[1]), "inp","schnappsEnv", "sessionI"), compress = F)
    if (DEBUG) cat(file = stderr(), paste0("add2history: inp, schnappsEnv to", tfile, "\n"))
    # the load is commented out because it is not used at the moment and only takes time to load
    if(comment %in% c("scEx", "scEx_log")) {
      commentOutLoad = "#"
    } else {
      commentOutLoad = ""
    }
    line <- paste0(
      "```{R, fig.cap = \"\", fig.width=20,fig.height=7}\n#",date(),"\n#load ", names(varnames[1]), "\n", commentOutLoad,"load(file = \"", basename(tfile),
      "\")\n#", comment, "\n```\n"
    )
    write(line, file = .schnappsEnv$historyFile, append = TRUE)
    if (DEBUG) cat(file = stderr(), paste0("add2history: saving stuff to ", .schnappsEnv$historyFile, "\n"))
  }
  
  if (type == "renderPlotly") {
    
    tfile <- tempfile(pattern = paste0(names(varnames[1]), "."), tmpdir = ".", fileext = ".png")
    assign(names(varnames[1]), arg[1])
    # save(file = tfile, list = c(names(varnames[1])))
    tryCatch({
      # orca(plotData$plotData, file = tfile, format = "png")
      withr::with_dir(normalizePath(.schnappsEnv$historyPath), orca(plotData$plotData, file = tfile, format = "png"))
      line <- paste0(
        # "```{R, fig.cap = \"\"}\n#", comment, "\n#load ", names(varnames[1]), "\nload(file = \"", basename(tfile),
        # "\")\nhtmltools::tagList(", names(varnames[1]), ")\n```\n",
        "\n", date(), "\n![](",basename(tfile),")\n\n"
      )
      write(line, file = .schnappsEnv$historyFile, append = TRUE)
      if (DEBUG) cat(file = stderr(), paste0("add2history: saving renderPlotly to", .schnappsEnv$historyFile, "\n"))
    }, error = function(w){
      cat(file = stderr(),paste("problem with orca:",w,"\n"))
    })
  }
  
  if (type == "tronco") {
    # deepDebug()
    tfile <- tempfile(pattern = paste0(names(varnames[1]), "."), tmpdir = .schnappsEnv$historyPath, fileext = ".RData")
    assign(names(varnames[1]), arg[[1]])
    # report.env <- getReactEnv(DEBUG = .schnappsEnv$DEBUG)
    save(file = tfile, list = c(names(varnames[1]), "inp", "schnappsEnv", "sessionI"), compress = F)
    if (DEBUG) cat(file = stderr(), paste0("add2history: tronco data to", tfile, "\n"))
    
    line <- paste0(
      "```{R, fig.cap = \"\", fig.width=20,fig.height=7}\n#",date(),"\n#", comment, "\n#load ", names(varnames[1]), "\nload(file = \"", basename(tfile),"\")\n",
      "\n", names(varnames[1]) ,"$filename <- NULL \n",
      "\ndo.call(TRONCO::pheatmap, ", names(varnames[1]), ")\n```\n"
    )
    write(line, file = .schnappsEnv$historyFile, append = TRUE)
    if (DEBUG) cat(file = stderr(), paste0("add2history: saving tronco to", .schnappsEnv$historyFile, "\n"))
  }
  
  if (type == "renderPlot") {
    tfile <- tempfile(pattern = paste0(names(varnames[1]), "."), tmpdir = .schnappsEnv$historyPath, fileext = ".RData")
    assign(names(varnames[1]), arg[[1]])
    # report.env <- getReactEnv(DEBUG = .schnappsEnv$DEBUG)
    save(file = tfile, list = c(names(varnames[1]), "inp", "schnappsEnv", "sessionI"), compress = F)
    if (DEBUG) cat(file = stderr(), paste0("add2history: saving renderPlot to", tfile, "\n"))
    
    line <- paste0(
      "```{R, fig.cap = \"\", fig.width=20,fig.height=7}\n#",date(),"\n#", comment, "\n#load ", names(varnames[1]), "\nload(file = \"", basename(tfile),"\")\n",
      "\n", names(varnames[1]), "\n```\n"
    )
    write(line, file = .schnappsEnv$historyFile, append = TRUE)
    if (DEBUG) cat(file = stderr(), paste0("add2history: saving renderPlot to", .schnappsEnv$historyFile, "\n"))
    
  }
  
  if (type == "renderDT") {
    tfile <- tempfile(pattern = paste0(names(varnames[1]), "."), tmpdir = .schnappsEnv$historyPath, fileext = ".RData")
    assign(names(varnames[1]), arg[[1]])
    # report.env <- getReactEnv(DEBUG = .schnappsEnv$DEBUG)
    save(file = tfile, list = c(names(varnames[1]), "inp", "schnappsEnv", "sessionI"), compress = F)
    if (DEBUG) cat(file = stderr(), paste0("add2history: saving renderDT to", tfile, "\n"))
    
    line <- paste0(
      "```{R}\n#",date(),"\n#", comment, "\n#load ", names(varnames[1]), "\nload(file = \"", basename(tfile),"\")\n",
      "\n", names(varnames[1]), "\n```\n"
    )
    write(line, file = .schnappsEnv$historyFile, append = TRUE)
    if (DEBUG) cat(file = stderr(), paste0("add2history: saving renderDT to", .schnappsEnv$historyFile, "\n"))
  }
  if (DEBUG) cat(file = stderr(), paste0("add2history done\n"))
  
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

# heatmapModuleFunction =====

heatmapModuleFunction <- function(
    heatmapData = NULL,
    addColNames = "sampleNames",
    orderColNames = c(), 
    sortingCols = "dendrogram",
    sortingRows = "dendrogram",
    scale = "none",
    colPal= "none",
    minMaxVal = NULL,
    proje = NULL,
    outfile = NULL,
    heatmapCellGrp = 5
) {
  
  colTree = FALSE
  if (is.null(sortingCols)) return(NULL)
  if (sortingCols == "dendrogram") colTree = TRUE
  
  if (is.null(heatmapData) | is.null(proje) | is.null(heatmapData$mat)) {
    return(NULL)
    # return(list(
    #   src = "empty.png",
    #   contentType = "image/png",
    #   width = 96,
    #   height = 96,
    #   alt = "pHeatMapPlot should be here (null)"
    # ))
  }
  if(is.null(minMaxVal)) minMaxVal = c(min(heatmapData$mat), max(heatmapData$mat))
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/heatmapModuleFunction.RData", list = c(ls()))
  }
  # cp =load(file = "~/SCHNAPPsDebug/heatmapModuleFunction.RData")
  # require(heatmaply)
  # if (is.null(pWidth)) {
  #   pWidth <- 800
  #   pHeight <- 300
  # }
  if (!is(heatmapData$mat, "matrix") & !is(heatmapData$mat, "Matrix")) {
    cat(file = stderr(), "!!!!! output$pHeatMapModule:mat is not a matrix\n")
    heatmapData$mat <- as.matrix(t(heatmapData$mat))
  }
  if (is(heatmapData$mat, "sparseMatrix")) {
    heatmapData$mat = as.matrix(heatmapData$mat)
  }
  
  if (sortingRows == "dendrogram") {
    heatmapData$cluster_rows = T
  }
  if (sortingRows == "list") {
    heatmapData$cluster_rows = F
  }
  # should be handled by plotly
  # heatmapData$mat[is.nan(heatmapData$mat)] <- 0.0
  # 
  # to check if this is handled better in plotly
  if (is.null(scale)) heatmapData$scale <- "none"
  heatmapData$name = "Expr."
  switch(scale,
         none = {heatmapData$scale <- "none"},
         row = {
           heatmapData$scale <- scale
           heatmapData$name = "cent.+scaled Expr."
           rmRows = -which(apply(heatmapData$mat, 1, sd) == 0)
           if(length(rmRows) > 1){
             if (.schnappsEnv$DEBUG) {
               cat(file = stderr(), "!!!!! output$pHeatMapModule:mat removing row due to sd = 0 for some genes\n")
             }
             heatmapData$mat <- heatmapData$mat[rmRows, ]
           }
         },
         column = {
           heatmapData$scale <- scale
           heatmapData$name = "cent.+scaled Expr."
           rmCols = -which(apply(heatmapData$mat, 2, sd) == 0)
           if (length(rmCols) > 1){
             if (.schnappsEnv$DEBUG) {
               cat(file = stderr(), "!!!!! output$pHeatMapModule:mat removing cols due to sd = 0 for some genes\n")
             }
             heatmapData$mat <- heatmapData$mat[, rmCols]
           }
         },
         row_order = {
           heatmapData$scale <- "none"
           heatmapData$name = "ranked Expr."
           
           # the "[]" are needed to preserve the row/column names
           heatmapData$mat[] <- t(apply(heatmapData$mat,1,FUN = function(x)rank(x, ties.method = "min")))
         },
         col_order = {
           heatmapData$scale <- "none"
           heatmapData$name = "ranked Expr."
           heatmapData$mat[] <- apply(heatmapData$mat,2,FUN = function(x)rank(x, ties.method = "min"))
         }
  )
  
  
  heatmapData$annotation_col = data.frame(row.names = colnames(heatmapData$mat))
  # heatmapData$filename <- outfile
  # heatmapData$filename = NULL
  # if (length(addColNames) > 0 & moreOptions) {
  addColNames = addColNames[addColNames %in% colnames(proje)]
  if (length(addColNames) > 0) {
    heatmapData$annotation_col <- proje[rownames(heatmapData$annotation_col), addColNames, drop = FALSE]
  }
  # if (sum(orderColNames %in% colnames(proje)) > 0 & moreOptions) {
  orderColNames = orderColNames[orderColNames %in% colnames(proje)]
  if (sum(orderColNames %in% colnames(proje)) > 0) {
    heatmapData$cluster_cols <- FALSE
    proje2 = proje[,orderColNames,drop = F]
    to.replace <- sapply(proje2, is.logical)
    proje2[,to.replace] <- lapply(proje2[,to.replace,drop = F], as.numeric)
    proje[,names(to.replace)[to.replace]] <- proje2[,names(to.replace)[to.replace],drop = F]
    if(ncol(proje)<2) return(NULL) # this cannot happen, but is a requirement fo dfOrder to work properly
    colN <- rownames(psychTools::dfOrder(proje[,,drop=F], orderColNames))
    matN <- colnames(heatmapData$mat)
    if (is.null(matN)) {
      matN <- names(heatmapData$mat)
    }
    colN <- colN[colN %in% matN]
    heatmapData$mat <- heatmapData$mat[, colN, drop = FALSE]
    heatmapData$annotation_col = heatmapData$annotation_col[colN,,drop = FALSE]
    
    # return()
  }
  # if (moreOptions) {
  heatmapData$cluster_cols <- colTree
  # }
  # orgMat = heatmapData$mat
  
  # heatmapData$mat = orgMat
  # system.time(do.call(pheatmap, heatmapData))
  # heatmapData$mat = as(orgMat, "TsparseMatrix")
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
  # heatmapData$width <- pWidth / 72
  # heatmapData$height <- pHeight / 72
  
  if (colPal == "none") {
    # use the supplied colors
  } else {
    heatmapData$color <- colorRampPalette(rev(brewer.pal(
      n = 7, name = colPal
    )))(100)
  }
  
  heatmapData$mat[heatmapData$mat <= minMaxVal[1]] = minMaxVal[1]
  heatmapData$mat[heatmapData$mat >= minMaxVal[2]] = minMaxVal[2]
  
  # working with large trees is not working
  # cluster_within_group and other dendextend function are recursive and crash
  # thus we only provide the cut 
  
  if (colTree) {
    # retVal = do.call(pheatmap::pheatmap, heatmapData)
    # retVal$tree_col
    callback = function(hc, ...){dendsort(hc)}
    heatmapData$clustering_callback = callback
    require(ComplexHeatmap)
    my_hclust <- hclust(dist(t(heatmapData$mat)), method = "complete")
    
    my_gene_col <- cutree(tree =  my_hclust, k = heatmapCellGrp)
    # while(max(table(my_gene_col)) > 1000){
    #   heatmapCellGrp = heatmapCellGrp + 10
    #   my_gene_col <- cutree(tree =  my_hclust, k = heatmapCellGrp)
    # }
    heatmapData$annotation_col[names(my_gene_col), "DendrogramCluster"] = factor(my_gene_col)
    heatmapData$mat = heatmapData$mat[,rownames(heatmapData$annotation_col)]
    heatmapData$cutree_cols = heatmapCellGrp
    if (ncol(heatmapData$mat)>2000){
      sortedCells = names(sort(my_gene_col))
      heatmapData$mat = heatmapData$mat[,sortedCells, drop = FALSE]
      heatmapData$annotation_col = heatmapData$annotation_col[sortedCells,,drop = FALSE]
      heatmapData$cluster_cols = F
      heatmapData$col_split = NULL
      heatmapData$border = NULL
      heatmapData$border_color = NA
      heatmapData$clustering_callback = NA
      my_gene_col = my_gene_col[sortedCells]
      heatmapData$gaps_col = which(my_gene_col[-1] != my_gene_col[-length(my_gene_col)])
    }
    
  }
  
  set.seed(1) # to make clustering reproducible
  # run_draw is needed to interprete mouse click events
  # drawback is that it creates a local plot (outside of Shiny)
  heatmapData$run_draw = T 
  if (sortingCols == "gene (click)") heatmapData$run_draw = T
  
  # make sure only values that are plotted are in the annotation
  heatmapData$annotation_col = droplevels(droplevels(heatmapData$annotation_col))
  
  
  # names(heatmapData)
  # orgHeatmapData = heatmapData
  # 
  # heatmapData = orgHeatmapData
  # heatmapData$scale = NULL
  # heatmapData$legend = T
  
  # heatmapData$legend_breaks = c(1,2,3,5)
  
  # heatmapData$run_draw = T
  # heatmapData$heatmap_legend_param = list(title="foo")
  # cutM =cut(heatmapData$mat,6)
  # heatmapData$mat = as.numeric(cutM)
  # heatmapData$legend_breaks = levels(cutM)
  # heatmapData
  # do.call(ComplexHeatmap::pheatmap, heatmapData)
  # 
  
  # if modules are available don't do cluster rows
  if(!is.null(heatmapData$gaps_row)) heatmapData$cluster_rows = FALSE
  ht_opt$message = F
  if(length(heatmapData$annotation_colors)>0){
    heatmapData$annotation_colors = heatmapData$annotation_colors[which(!unlist(lapply(heatmapData$annotation_colors,is.null)))]
  }
  # col_fun =  c("blue", "white", "red","green")
  # heatmapData$annotation_colors$PC2=col_fun
  # heatmapData$annotation_colors$PC2
  # 
  # hande missing values in annotation_col
  # heatmapData$annotation_col <- 
  heatmapData$annotation_col = heatmapData$annotation_col %>% mutate_if(is.numeric, function(x) ifelse( is.na(x), min(x, na.rm = T),x))
  
  
  retVal = tryCatch(
    do.call(ComplexHeatmap::pheatmap, heatmapData),
    # do.call(TRONCO::pheatmap, heatmapData),
    error = function(e){
      cat (file = stderr(), paste(e))
      save(file = "~/SCHNAPPsDebug/heatmapError.rdata", list = ls(parent.frame(n=4)))
      # cp = load("~/SCHNAPPsDebug/heatmapError.rdata")
      # browser()
      return(NULL)
    }
  )
  
  return(retVal)
}


# this has to be called after the cache setup... moved to server.R
# heatmapModuleFunction_m = memoise::memoise(heatmapModuleFunction,cache=do.call(cachem::cache_disk,.schnappsEnv$cacheDir))
heatmapModuleFunction_m = heatmapModuleFunction

# consolidateScEx ----

consolidateScEx <-
  function(scEx, projections, scEx_log, pca, tsne) {
    # save(file = "~/SCHNAPPsDebug/consolidate.RData", list = c(ls(), "myProjections"))
    # load(file = "~/SCHNAPPsDebug/consolidate.RData")
    commCells <- colnames(scEx)
    if(!is.null(scEx_log)){
      commCells <- base::intersect(colnames(scEx), colnames(scEx_log))
      commGenes <- base::intersect(rownames(scEx), rownames(scEx_log))
      scEx <- scEx[commGenes, commCells]
      assays(scEx)[["logcounts"]] <- assays(scEx_log)[[1]][commGenes, commCells]
    }
    # what about UMAP??? others? => they are considered as projections not as reducedDims
    if(!is.null(pca)){
      reducedDims(scEx) <- SimpleList(PCA = pca$x[commCells, ], TSNE = tsne[commCells, ])
    }
    if(!is.null(projections)){
      for (name in colnames(projections)) {
        colData(scEx)[[name]] <- projections[commCells, name]
      }
    }
    # colData(scEx)[["before.Filter"]] <- projections[commCells, "before.filter"]
    # colData(scEx)[["dbCluster"]] <- projections[commCells, "dbCluster"]
    # colData(scEx)[["UmiCountPerGenes"]] <- projections[commCells, "UmiCountPerGenes"]
    # colData(scEx)[["UmiCountPerGenes2"]] <- projections[commCells, "UmiCountPerGenes2"]
    
    return(scEx)
  }


# loadLiteData ----

loadLiteData <- function(fileName = NULL) {
  if (is.null(fileName)) return(NULL)
  # fileName = "~/Rstudio/UTechSCB-SCHNAPPs/data/scExLite.RData"
  cp = load(fileName)
  
  # The data has to be stored in scEx as it come from save from the main SCHNAPPs app
  if (!all(c("scEx", "pcaReact") %in% cp)) {
    # this is to cope with a change of reactive name
    if (all(c("scEx", "pca") %in% cp)) {
      pcaReact = pca
    } else {
      return(NULL)
    }
  }
  
  projections = colData(scEx)
  dbCluster = projections$dbCluster
  counts = scEx
  assays(counts)[["logcounts"]] = NULL
  logcounts = scEx
  assays(logcounts)[["counts"]] = NULL
  if (!exists("ccol", inherits = F)){
    # cluster colors
    inCols <- list()
    lev <- levels(dbCluster)
    # inCols <- allowedColors[1:length(lev)]
    inCols <-  allowedColors[rep(1:length(allowedColors),ceiling(length(lev) / length(allowedColors)))[1:length(lev)]]
    names(inCols) <- lev
    ccol <- unlist(inCols)
  }
  projections$sampleNames = factor(projections$sampleNames)
  if (!exists("scol", inherits = F)) {
    sampNames = levels(projections$sampleNames)
    scol <-  rev(allowedColors[rep(1:length(allowedColors),ceiling(length(lev) / length(allowedColors)))[1:length(lev)]])
    # scol <- rev(allowedColors)[seq_along(sampNames)]
    names(scol) <- sampNames
  }
  # pcaReact = reducedDims(scEx)[["PCA"]] # now stored separately
  returnList = list(scEx = counts, scEx_log = logcounts, pcaReact = pcaReact, projections = projections, dbCluster = dbCluster, clusterCol = ccol, sampleCol = scol)
  for (va in cp[!cp %in% c("scEx", "pcaReact", "ccol", "scol")]) {
    # .schnappsEnv = global envirnment for schnapps
    #  - projectionFunctions = list of 2 entries each: 1st: display name, 2nd reactive name (defined in a reactive.R)
    
    returnList[[va]] = get(va)
  }
  
  return(returnList)
}

# defaultValue ---
# get a value that has been supplied as a parameter or return the default value val.
defaultValue <- function(param = "coEtgMinExpr", val ) {
  # if (DEBUG) cat(file = stderr(), paste( "defaultValue : ",param, " val: ", str(val), "\n"))
  # deepDebug()
  # if(param == "Scorpius_dataInput-Mod_clusterPP")
  #   browser()
  if (!exists(".schnappsEnv")) return (val)
  if (exists(envir = .schnappsEnv, x = "defaultValues")) {
    if ( param %in% names(.schnappsEnv$defaultValues)) {
      # if (.schnappsEnv$DEBUG) cat(file = stderr(), paste( "value: ", str(.schnappsEnv$defaultValues[[param]]), "\n"))
      return(.schnappsEnv$defaultValues[[param]])
    }
  } 
  return(val)
}

# dheader ----
# it is used in ui.R and ui-lite.R and is still being developped
dheader <- function() {
  shinydashboardPlus::dashboardHeader(
    controlbarIcon = "Wkfl",
    title = paste("SCHNAPPs", packageVersion("SCHNAPPs")),
    # title = tags$a(tags$img(src='www/images/logo.32.32.png'),  paste("SCHNAPPs", packageVersion("SCHNAPPs"))),
    shinydashboard::dropdownMenu(type = "task", icon = icon("fas fa-question"),badgeStatus = NULL,
                                 headerText = "Help",
                                 # notificationItem(text =  actionButton("menuTour", label = "short Tour", icon = icon("fas fa-directions")),
                                 #                  icon = icon("", verify_fa = FALSE)
                                 # ),
                                 notificationItem(text =  "online documentation",
                                                  href="https://c3bi-pasteur-fr.github.io/UTechSCB-SCHNAPPs/", 
                                                  icon("fas fa-book-medical")
                                 )
    ), 
    shinydashboard::dropdownMenu(type = "task", icon = icon("fas fa-info"),badgeStatus = NULL,
                                 headerText = "About",
                                 notificationItem(text =  actionButton("AboutApp", label = "about SCHNAPPs"),
                                                  icon = icon("", verify_fa = FALSE)
                                 )
    )
  )
}

# boxWhelp ----
# box with help
# modified box function that places a question mark with a botton at the far right
# 
# boxWhelp <- function (..., title = NULL, footer = NULL, status = NULL, solidHeader = FALSE, 
#                       background = NULL, width = 6, height = NULL, collapsible = FALSE, 
#                       collapsed = FALSE, helpID = NULL) 
# {
#   boxClass <- "box"
#   if (solidHeader || !is.null(background)) {
#     boxClass <- paste(boxClass, "box-solid")
#   }
#   if (!is.null(status)) {
#     # validateStatus(status)
#     boxClass <- paste0(boxClass, " box-", status)
#   }
#   if (collapsible && collapsed) {
#     boxClass <- paste(boxClass, "collapsed-box")
#   }
#   if (!is.null(background)) {
#     # validateColor(background)
#     boxClass <- paste0(boxClass, " bg-", background)
#   }
#   style <- NULL
#   if (!is.null(height)) {
#     style <- paste0("height: ", validateCssUnit(height))
#   }
#   titleTag <- NULL
#   if (!is.null(title)) {
#     titleTag <- h3(class = "box-title", title)
#   }
#   helpTag <- NULL
#   if (!is.null(helpID)) {
#     helpTag <- actionButton(inputId = helpID, label = "", icon = icon("fas fa-question")
#     )
#   }
#   collapseTag <- NULL
#   if (collapsible) {
#     buttonStatus <- status 
#     if (is.null(buttonStatus)) buttonStatus = "default"
#     collapseIcon <- if (collapsed) 
#       "plus"
#     else "minus"
#     collapseTag <- div(class = "box-tools pull-right", helpTag, tags$button(class = paste0("btn btn-box-tool"), 
#                                                                             `data-widget` = "collapse", shiny::icon(collapseIcon)))
#   } else {
#     if (!is.null(helpTag))
#       collapseTag <- div(class = "box-tools pull-right", 
#                          helpTag
#       )
#   }
#   
#   headerTag <- NULL
#   if (!is.null(titleTag) || !is.null(collapseTag)) {
#     headerTag <- div(class = "box-header", titleTag, collapseTag)
#   }
#   div(class = if (!is.null(width)) 
#     paste0("col-sm-", width), div(class = boxClass, style = if (!is.null(style)) 
#       style, headerTag, div(class = "box-body", ...), if (!is.null(footer)) 
#         div(class = "box-footer", footer)))
# }

# setId ----
# gives an area an ID to be referenced by introjs
# to be able to use %>%
setId = function(inp, id) {return(tags$div(id=id, inp))}

### used in coE - violin plot - combinations
# combinePermutations ----
combinePermutations <- function(perm1, perm2) {
  perms <- rep("", length(perm1))
  for (cIdx in 1:length(perm1)) {
    if (perm2[cIdx] == "") {
      perms[cIdx] <- perm1[cIdx]
    } else {
      perms[cIdx] <- perm2[cIdx]
    }
  }
  perms
}

# finner
finner <- function(xPerm, r, genesin, featureData, scEx_log, perms, minMaxExpr) {
  comb <- gtools::combinations(xPerm, r, genesin)
  for (cIdx in 1:nrow(comb)) {
    map <-
      rownames(featureData[which(toupper(featureData$symbol) %in% comb[cIdx, ]), ])
    # permIdx <- Matrix::colSums(exprs(gbm[map, ]) >= minExpr) == length(comb[cIdx, ])
    
    permIdx <- Matrix::colSums(
      SummarizedExperiment::assays(scEx_log)[[1]][map, , drop = FALSE] >=  minMaxExpr[1] &
        SummarizedExperiment::assays(scEx_log)[[1]][map, , drop = FALSE] <=  minMaxExpr[2]) == length(comb[cIdx, 1:ncol(comb)])
    perms[permIdx] <- paste0(comb[cIdx, 1:ncol(comb)], collapse = "+")
  }
  perms
}



# only return if there is no variable in env.
checkAllowed <- function(x, env=.schnappsEnv) {
  if(!is(x, "shiny.tag")) {
    return(x)
  }
  id = paste0("allowFunctionality--", x$attribs$id)
  if (exists(id, envir = env)) {
    if (!get(id, envir = env)) {
      return("")
    }
  }
  return(x)
}

### History
# createHistory ---

createHistory <- function(.schnappsEnv) {
  
  .schnappsEnv$historyPath = paste0(.schnappsEnv$historyPath, "/hist_",format(Sys.time(), "%Y-%b-%d.%H.%M"))
  if (!dir.exists(.schnappsEnv$historyPath)){
    dir.create(.schnappsEnv$historyPath, recursive = T)
  }  
  if (!exists("historyFile", envir = .schnappsEnv)) {
    .schnappsEnv$historyFile = paste0("history.",format(Sys.time(), "%Y-%b-%d.%H.%M"),".Rmd")
  }
  if (is.null(.schnappsEnv$historyFile)) {
    .schnappsEnv$historyFile = "history2.Rmd"
  }
  .schnappsEnv$historyFile <- paste0(.schnappsEnv$historyPath, .Platform$file.sep, basename(.schnappsEnv$historyFile))
  line=paste0("---
title: \"history\"
output:
  bookdown::html_document2:
    toc: true
    toc_depth: 3
    toc_float: true
    number_sections: true
    code_folding: hide
---
  
```{r setup, include=FALSE}
  knitr::opts_chunk$set(echo = TRUE)
  suppressMessages(require(shiny))
  suppressMessages(require(shinyTree))
  suppressMessages(require(tibble))
  suppressMessages(require(plotly))
  suppressMessages(require(shinythemes))
  suppressMessages(require(ggplot2))
  suppressMessages(require(ggalluvial))
  suppressMessages(require(DT))
  suppressMessages(require(pheatmap))
  # suppressMessages(require(threejs))
  suppressMessages(require(RColorBrewer))
  # suppressMessages(require(mclust))
  suppressMessages(require(reshape2))
  suppressMessages(require(knitr))
  suppressMessages(require(shinyWidgets))
  suppressMessages(require(scater))
  suppressMessages(require(kohonen))
  # suppressMessages(require(Rsomoclu))
  suppressMessages(require(SingleCellExperiment))
  suppressMessages(require(Matrix))
  suppressMessages(require(colourpicker))
  # suppressMessages(require(shinytest))
  suppressMessages(require(ggalluvial))
  suppressMessages(require(scran))
  suppressMessages(require(BiocSingular))
  suppressMessages(require(dplyr))

  if (\"debugme\" %in% rownames(installed.packages())) {
    suppressMessages(require(debugme))
  }
  if (\"gtools\" %in% rownames(installed.packages())) {
    suppressMessages(require(gtools))
  }
  if (\"kableExtra\" %in% rownames(installed.packages())) {
    suppressMessages(require(kableExtra))
  }
  if (\"reactlog\" %in% rownames(installed.packages())) {
    suppressMessages(require(reactlog))
  }
  suppressMessages(require(Seurat))
  .schnappsEnv = list()
  .schnappsEnv$DEBUGSAVE = FALSE
  source(system.file(\"app\", \"serverFunctions.R\", package = \"SCHNAPPs\"))
  DEBUG=FALSE
```\n\n" )
  write(line,file=.schnappsEnv$historyFile,append=FALSE)
}

sc_textInput <- function(...){
  arg <- list(...)
  .schnappsEnv$textInputlist = unique(c(.schnappsEnv$textInputlist, arg[[1]][1]))
  return(textInput(...))
}
sc_numericInput <- function(...){
  arg <- list(...)
  .schnappsEnv$numericInputList = unique(c(.schnappsEnv$numericInputList, arg[[1]][1]))
  return(numericInput(...))
}
sc_checkboxInput <- function(...){
  arg <- list(...)
  # browser()
  # cat(file = stderr(), ..., "\n")
  .schnappsEnv$checkboxInputList = unique(c(.schnappsEnv$checkboxInputList, arg[[1]][1]))
  return(checkboxInput(...))
}
sc_selectInput <- function(...){
  arg <- list(...)
  .schnappsEnv$selectInputList = unique(c(.schnappsEnv$selectInputList, arg[[1]][1]))
  return(selectInput(...))
}
sc_radioButtons <- function(...){
  arg <- list(...)
  .schnappsEnv$radioButtonsList = unique(c(.schnappsEnv$radioButtonsList, arg[[1]][1]))
  return(radioButtons(...))
}
sc_selectizeInput <- function(...){
  arg <- list(...)
  .schnappsEnv$selectizeInputList = unique(c(.schnappsEnv$selectizeInputList, arg[[1]][1]))
  return(selectizeInput(...))
}



# textAreaInput => used for projections rename not useful
# orderInput => used for rearrange levels not useful

deepDebug <-function(){
  # callingFun = as.list(sys.call(-1))[[1]]
  # cat(file = stderr(), paste("called from :", callingFun, "\n"))
  # deepDebug()
  # lobstr::cst()
  # traceback()
  # # trace()
  # recover()
  # # lobstr::sxp()
}

## loadInput ---
loadInput = function(inp, session, input){
  # deepDebug()
  deepDebug()
  if (is.null(getDefaultReactiveDomain())) {
    cat(file = stderr(), "I need to be in a reactive context.")
    stop()
  }
  
  for (id in .schnappsEnv$textInputlist){
    if(!id %in% names(inp)) {
      cat(file = stderr(), paste("not found:", id,"\n"))
      next()
    }
    updateTextInput(session = session, inputId = id, value = inp[[id]])
  }
  for (id in .schnappsEnv$numericInputList){
    if(!id %in% names(inp)) {
      cat(file = stderr(), paste("not found numericInputList:", id,"\n"))
      next()
    }
    updateNumericInput(session = session, inputId = id, value = inp[[id]])
  }
  for (id in .schnappsEnv$checkboxInputList){
    if(!id %in% names(inp)) {
      cat(file = stderr(), paste("not found checkboxInputList:", id,"\n"))
      next()
    }
    updateCheckboxInput(session = session, inputId = id, value = inp[[id]])
  }
  for (id in .schnappsEnv$selectInputList){
    if(!id %in% names(inp)) {
      cat(file = stderr(), paste("not found checkboxInputList:", id,"\n"))
      next()
    }
    updateSelectInput(session = session, inputId = id, selected = inp[[id]])
  }
  for (id in .schnappsEnv$radioButtonsList){
    deepDebug()
    inp[[id]]
    if(!id %in% names(inp)) {
      cat(file = stderr(), paste("not found checkboxInputList:", id,"\n"))
      next()
    }
    updateRadioButtons(session = session, inputId = id, selected = inp[[id]])
  }
  for (id in .schnappsEnv$selectizeInputList){
    inp[[id]]
    if(!id %in% names(inp)) {
      cat(file = stderr(), paste("not found checkboxInputList:", id,"\n"))
      next()
    }
    updateSelectizeInput(session = session, inputId = id, selected = inp[[id]])
  }
  
  # updateRadioButtons(session = session, inputId = "whichscLog", selected = "calcLog")
  
  for(inName in names(inp)){
  }
}

debugControl <- function( name = "cellSelectionModule", list = c(ls())){
  if(.schnappsEnv[["DEBUG"]]){
    cat(file = stderr(), paste("debugControl: ", name, "\n"))
  }
  envir = parent.frame(n=1)
  listObj = lapply(list, FUN=function(x) get(x, envir = envir))
  names(listObj) = list
  saveVal = FALSE
  if(!is.null(.schnappsEnv[["DEBUGSAVE"]])) {
    saveVal = .schnappsEnv[["DEBUGSAVE"]]
  }
  saveVal = .schnappsEnv$DEBUGSAVE
  if(!is.null(.schnappsEnv[[paste0("DEBUGSAVE_",name)]])) {
    saveVal = .schnappsEnv[[paste0("DEBUGSAVE_",name)]]
  }
  if (saveVal ) {
    browser()
    save(file = paste0("~/SCHNAPPsDebug/", name, ".RData"), list = list)
  }
  if(.schnappsEnv[["DEBUG"]]){
    cat(file = stderr(), paste("debugControl: ", name, "done\n"))
  }
}


# add2Projections ----
#' add a column to the newPrjs oject
#' 
#' This functions makes sure that all rows from the original input data are used, the order of rows is correct
#' that there are no name collisions in the columns
#' 
#' @param newPrjs is projectionsTable$newProjections, a reactive that holds user defined projections
#' this is added to the predefined projections (projections())
#' 
#' @param proj2Add a data frame that should be added to newPrjs
#' 
#' @param acn output from allCellNames(), a reactive, (needed if newPrjs is empty)
#' 
#' @param projections the complete projections data frame from projections() reactive
#' 
#' @return newPrjs that can replace projectionsTable$newProjections
#' 
#' 
add2Projections <- function(proj2Add, acn, newPrjs, projections){
  if (is.null(proj2Add)) {
    return(NULL)
  }
  if(!is.data.frame(newPrjs))
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/add2Projections.RData",
           list = c( ls())
      )
    }
  # cp=  load(file="~/SCHNAPPsDebug/gQC_renameLevButton.RData")
  which(colnames(proj2Add) %in% colnames(projections)) 
  if(is.null(
    tryCatch({
      if (ncol(newPrjs) == 0) {
        newPrjs = data.frame(row.names = acn)
      }
      # 2DO: what happes if newPrjs has more rows than original? How could this happen?
      newPrjs <- dplyr::full_join(
        tibble::rownames_to_column(newPrjs), 
        tibble::rownames_to_column(proj2Add), 
        by='rowname')
      rownames(newPrjs) = newPrjs[,1]
      newPrjs = newPrjs[,-1]
    }, error=function(w){
      cat(file = stderr(), paste("something went wrong during add2Projections", w,"\n"))
      if (!is.null(getDefaultReactiveDomain()))
        showNotification("problem with names", id = "renameProbl", duration = NULL, type = "error")
      return(NULL)
    }))) return(NULL)
  # what happens if NA is introduced by the join
}



if (.schnappsEnv$DEBUG) {
  cat(file = stderr(), "end severFunctions.R\n")
}
