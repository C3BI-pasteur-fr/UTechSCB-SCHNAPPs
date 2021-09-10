suppressMessages(require(ggplot2))
# deProjTable ----
# for 2D plot in ExpressionPanel
deProjTable <- reactive({
  projections <- projections()
  selectedCells <- DE_Exp_dataInput() #DE_Exp_dataInput
  cellNs <- selectedCells$cellNames()
  
  if (DEBUG) cat(file = stderr(), "observeEvent: input$DE_clusterPP\n")
  # Can use character(0) to remove all choices
  if (is.null(projections) | is.null(cellNs)) {
    return(NULL)
  }
  # browser()
  cellNs <- cellNs[cellNs %in% rownames(projections)]
  projections[cellNs, ]
})

.schnappsEnv$coE_PPGrp <- "sampleNames"
observe(label = "DE_PPGrp", {
  if (DEBUG) cat(file = stderr(), paste0("observe: DE_PPGrp\n"))
  .schnappsEnv$DE_PPGrp <- input$DE_PPGrp
})
.schnappsEnv$coE_PPSelection <- "1"
observe(label = "DE_clusterPP", {
  if (DEBUG) cat(file = stderr(), paste0("observe: DE_clusterPP\n"))
  .schnappsEnv$DE_clusterPP <- input$DE_clusterPP
})

# DE_updateInputPPt ====
DE_updateInputPPt <- reactive({
  if (DEBUG) cat(file = stderr(), "DE_updateInputPPt started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_updateInputPPt")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_updateInputPPt")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_updateInputPPt", id = "DE_updateInputPPt", duration = NULL)
  }
  tsneData <- projections()
  projFactors <- projFactors()
  
  # Can use character(0) to remove all choices
  if (is.null(tsneData)) {
    return(NULL)
  }
  # save(file = "~/SCHNAPPsDebug/DE_updateInputPPt.Rdata", list = c(ls()))
  # load(file = "~/SCHNAPPsDebug/DE_updateInputPPt.Rdata")
  
  # coln <- colnames(tsneData)
  # choices <- c()
  # for (cn in coln) {
  #   if (length(levels(as.factor(tsneData[, cn]))) < 50) {
  #     choices <- c(choices, cn)
  #   }
  # }
  # if (length(choices) == 0) {
  #   choices <- c("no valid columns")
  # }
  updateSelectInput(
    session,
    "DE_clusterPP",
    choices = projFactors,
    selected = .schnappsEnv$DE_PPGrp
  )
})

observeEvent(label = "DE_clusterPP", input$DE_clusterPP,{
  projections <- projections()
  if (DEBUG) cat(file = stderr(), "observeEvent: input$DE_clusterPP\n")
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  if(!input$DE_clusterPP %in% colnames(projections)) return(NULL)
  choicesVal = levels(projections[, input$DE_clusterPP])
  updateSelectInput(
    session,
    "DE_PPGrp",
    choices = choicesVal,
    selected = .schnappsEnv$DE_clusterPP
  )
  
})




observe(label = "save2HistScater", {
  clicked  = input$save2HistScater
  if (DEBUG) cat(file = stderr(), "observe input$save2HistScater \n")
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
              comment = paste0("scater plot\n", 
                               "fun = plotData$plotData$plotFunc\n", 
                               "environment(fun) = environment()\n",
                               "print(do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)]))\n"), 
              plotData = .schnappsEnv[["DE_scaterPNG"]])
  
})

observe(label = "save2HistPanel", {
  clicked  = input$save2HistPanel
  if (DEBUG) cat(file = stderr(), "observe input$save2HistPanel \n")
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
              comment = paste0("# Panel plot\n",
                               "fun = plotData$plotData$plotFunc\n", 
                               "environment(fun) = environment()\n",
                               "print(do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)]))\n"
              ),
              plotData = .schnappsEnv[["DE_panelPlot"]])
})



# DE_scaterPNG ----
#' DE_scaterPNG
#' reactive to plot highest expressed genes
#' take quite some time to compute, but since we normally don't need it
#' it is not in the heavyCalculations list.
#' TODO
#' maybe in a future version there can be a button to enable caclulations
DE_scaterPNG <- reactive({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_scaterPNG")
    }
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_scaterPNG", id = "DE_scaterPNG", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "DE_scaterPNG\n")
  
  clicked <- input$runScater
  cat(file = stderr(), paste("DE_scaterPNG", clicked, "\n"))
  # takes too long, commenting out for course
  # if (is.null(.schnappsEnv$scaterRan)){
  #   .schnappsEnv$scaterRan = 0
  #   return(list(
  #     src = "",
  #     contentType = "image/png",
  #     width = 10,
  #     height = 10,
  #     alt = "Scater plot will be here when 'apply changes' is checked"
  #   ))
  # }
  if (clicked < 1) {
    return(list(
      src = "",
      contentType = "image/png",
      width = 10,
      height = 10,
      alt = "Scater plot will be here when 'apply changes' is checked"
    ))
  }
  scaterReads <- isolate(scaterReads())
  if (is.null(scaterReads)) {
    return(list(
      src = "",
      contentType = "image/png",
      width = 10,
      height = 10,
      alt = "Scater plot will be here when 'apply changes' is clicked"
    ))
  }
  scols <- isolate(sampleCols$colPal)
  
  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scater.RData", list = c(ls()))
  }
  # cp=load(file='~/SCHNAPPsDebug/scater.RData')
  
  # calculations
  if (is.null(width)) {
    width <- 96 * 7
  }
  if (is.null(height)) {
    height <- 96 * 7
  }
  
  myPNGwidth <- width / 96
  myPNGheight <- height / 96
  
  outfile <- paste0(tempdir(), "/scaterPlot.png")
  # outfile <- paste0("~/SCHNAPPsDebug",'/scaterPlot.png')
  if (DEBUG) cat(file = stderr(), paste("output file: ", outfile, "\n"))
  if (DEBUG) cat(file = stderr(), paste("output file normalized: ", normalizePath(outfile, mustWork = FALSE), "\n"))
  n <- min(nrow(scaterReads), 50)
  
  rownames(scaterReads) <- rowData(scaterReads)$symbol
  p1 = pltHighExp( scaterReads, n, scols) 
  
  tryCatch(
    ggsave(file = normalizePath(outfile, mustWork = FALSE), plot = p1, width = myPNGwidth, height = myPNGheight, units = "in"),
    error = function(e) {
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("Problem saving ggplot", type = "warning", duration = NULL)
      }
      return(NULL)
    }
  )
  retVal <- list(
    src = normalizePath(outfile, mustWork = FALSE),
    contentType = "image/png",
    width = width,
    height = height,
    alt = "Scater plot should be here"
  )
  # end calculation
  af = pltHighExp
  # remove env because it is too big
  environment(af) = new.env(parent = emptyenv())
  
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
  
  printTimeEnd(start.time, "DE_scaterPNG")
  exportTestValues(DE_scaterPNG = {
    retVal
  })
  return(retVal)
})

# DE_dataExpltSNEPlot ---
#' DE_dataExpltSNEPlot
#' plot 3D tSNE in DataExploration - Expression
#' Here, only the expression of a gene or gene list is shown, compared to the other tSNE plot
#' in General QC - tSNE
#' DEPRICATED
DE_dataExpltSNEPlot <- function(scEx_log, g_id, projections,x,y,z) {
  if (is.null(scEx_log)) {
    return(NULL)
  }
  featureData <- rowData(scEx_log)
  geneid <- geneName2Index(g_id, featureData)
  if (length(geneid) == 0) {
    return(NULL)
  }
  if (length(geneid) == 1) {
    expression <- assays(scEx_log)[[1]][geneid, ]
  } else {
    expression <- Matrix::colSums(assays(scEx_log)[[1]][geneid, ])
  }
  
  validate(need(
    is.na(sum(expression)) != TRUE,
    "Gene symbol incorrect or gene not expressed"
  ))
  
  projections <- cbind(projections, expression)
  names(projections)[ncol(projections)] <- "values"
  if (!all(c("tsne1", "tsne2", "tsne3") %in% colnames(projections))) {
    showNotification("some tsne projections are not available.",
                     id = "DE_tsne_pltERROR",
                     duration = NULL, type = "error"
    )
  }
  
  p <-
    plotly::plot_ly(
      projections,
      x = ~get(x),
      y = ~tsne2,
      z = ~tsne3,
      type = "scatter3d",
      hoverinfo = "text",
      text = ~ paste(
        1:nrow(projections), " ", rownames(projections), "<br />",
        "Cluster:", as.numeric(as.character(projections$dbCluster))
      ),
      # text = paste("Cluster:", as.numeric(as.character(projections$dbCluster))),
      mode = "markers",
      source = "C",
      marker = list(
        size = 2,
        line = list(width = 0),
        color = ~values
        # ,colors = "Greens"
      )
    )
  layout(p, title = paste(featureData[geneid, "symbol"], collapse = ", "))
}

DE_geneViolinFunc <- function(scEx_log, g_id, projections, ccols) {
  if (is.null(scEx_log)) {
    return(NULL)
  }
  featureData <- rowData(scEx_log)
  geneid <- geneName2Index(g_id, featureData)
  
  if (length(geneid) == 0) {
    return(NULL)
  }
  
  # expression <- exprs(scEx_log)[geneid, ]
  if (length(geneid) == 1) {
    expression <- assays(scEx_log)[[1]][geneid, ]
  } else {
    expression <- Matrix::colSums(assays(scEx_log)[[1]][geneid, ])
  }
  
  projections <- cbind(projections, expression)
  names(projections)[length(projections)] <- "values"
  
  p1 <-
    ggplot(projections, aes(factor(dbCluster), values, fill = factor(dbCluster))) +
    geom_violin(scale = "width") +
    scale_color_manual(values = ccols) +
    scale_fill_manual(values = ccols, aesthetics = "fill") +
    
    stat_summary(
      fun = median,
      geom = "point",
      size = 5,
      color = "black"
    ) +
    stat_summary(fun.data = n_fun, geom = "text") +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        size = 12,
        vjust = 0.5
      ),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 14),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      legend.position = "none"
    ) +
    xlab("Cluster") +
    ylab("Expression") +
    ggtitle(paste(featureData[geneid, "symbol"], collapse = ", "))
  
  return(p1)
}

# DE_selectedCells <- reactiveValues(
#   selectedCells <- ""
# )


# panelplotFunc ----
# scEx_log singlecell Experiment object
# projections as used in schnapps
# genesin gene names to be plotted
# dimx4, dimy4 dimensions to be plotted on
# sameScale True/False
# nCol number of columns for final plot
# sampdes header for plot
# cellNs cell names to be used


panelPlotFunc <- function(scEx_log, projections, genesin, dimx4, dimy4, sameScale, nCol, sampdesc, cellNs,
                          lowCol = "blue", highCol = "red", midCol = "white",
                          midFunc = function(x){(max(x)-min(x))/2}) {
  
  featureData <- rowData(scEx_log)
  # featureData$symbol = toupper(featureData$symbol)
  genesin <- genesin[which(genesin %in% toupper(featureData$symbol))]
  if (length(genesin) < 1) {
    return(NULL)
  }
  par(mfrow = c(ceiling(length(genesin) / 4), 4), mai = c(0., .3, .3, .3))
  rbPal <- colorRampPalette(c("#f0f0f0", "red"))
  ylim <- c(min(projections[, dimy4]), max(projections[, dimy4]))
  if (is(projections[, dimx4], "factor") & dimy4 == "UMI.count") {
    ymax <- 0
    for (i in 1:length(genesin)) {
      geneIdx <- which(toupper(featureData$symbol) == genesin[i])
      ymax <- max(ymax, max(Matrix::colSums(assays(scEx_log)[[1]][geneIdx, , drop = FALSE])))
    }
    ylim <- c(0, ymax)
    if (!sameScale) {
      ylim <- NULL
    }
  }
  plotList <- list()
  plotIdx <- 0
  # if (cl4 == "All") {
  #   for (i in 1:length(genesin)) {
  #     geneIdx <- which(toupper(featureData$symbol) == genesin[i])
  #     Col <- rbPal(10)[
  #       as.numeric(
  #         cut(
  #           as.numeric(
  #             assays(scEx_log)[[1]][
  #               rownames(featureData[geneIdx, ]),
  #               ]
  #           ),
  #           breaks = 10
  #         )
  #       )
  #       ]
  #     plotIdx = plotIdx +1
  #
  #     plotList[[plotIdx]] = ggplot(projections, aes_string(x=dimx4, y=dimy4))
  #     if (is(projections[, dimx4], "factor") & dimy4 == "UMI.count") {
  #       projections[, dimy4] <- Matrix::colSums(assays(scEx_log)[["logcounts"]][geneIdx, , drop = FALSE])
  #       plotList[[plotIdx]] = plotList[[plotIdx]] + geom_boxplot(show.legend = FALSE) + ggtitle(genesin[i])
  #     } else{
  #       plotList[[plotIdx]] = plotList[[plotIdx]] + geom_point(color = Col, show.legend = FALSE) + ggtitle(genesin[i])
  #
  #     }
  #     if (!is.null(ylim)) {
  #       plotList[[plotIdx]] = plotList[[plotIdx]]  + ylim(ylim)
  #     }
  #     # plot(projections[, dimx4], projections[, dimy4],
  #     #      col = Col, pch = 16, frame.plot = TRUE, ann = FALSE, ylim = ylim
  #     # )
  #     # title(genesin[i], line = -1.2, adj = 0.05, cex.main = 2)
  #     if (DEBUG) cat(file = stderr(), genesin[i])
  #   }
  # } else {
  printData = data.frame(gene = character(), dimx = numeric(), dimy = numeric(), col = numeric)
  # save(file = "~/SCHNAPPsDebug/DE_panelPlot.RData", list = c(ls()), compress = F)
  # cp = load(file="~/SCHNAPPsDebug/DE_panelPlot.RData")
  
  for (i in 1:length(genesin)) {
    geneIdx <- which(toupper(featureData$symbol) == genesin[i])
    subsetTSNE <- projections[cellNs, ] # only work the subset of cells
    
    # color for each cell based on the expression
    # this will always show the max value per cell
    # this should apply when sameScale = FALSE
    Col <- rbPal(10)[
      as.numeric(
        cut(
          as.numeric(
            assays(scEx_log)[[1]][
              rownames(featureData[geneIdx, ]),
            ]
          ),
          breaks = 10
        )
      )
    ]
    
    names(Col) <- rownames(projections)
    plotCol <- Col[rownames(subsetTSNE)]
    if (is(projections[, dimx4], "factor") & dimy4 == "UMI.count") {
      projections[, dimy4] <- Matrix::colSums(assays(scEx_log)[[1]][geneIdx, , drop = FALSE])
      subsetTSNE <- projections[cellNs, ]
    }
    colData = assays(scEx_log)[[1]][
      rownames(featureData[geneIdx, ]),cellNs
    ]
    if (!sameScale){
      colData = (colData -min(colData) )/ (max(colData) - min(colData)) 
    }
    printData = rbind(printData, data.frame(gene = genesin[i], 
                                            dimx = subsetTSNE[cellNs, dimx4] , 
                                            dimy = subsetTSNE[cellNs, dimy4],
                                            col = colData
    ))
    # plotIdx <- plotIdx + 1
    # plotList[[plotIdx]] <- ggplot(subsetTSNE, aes_string(x = dimx4, y = dimy4))
    # 
    # # TODO show.legend should be based on sameScale
    # if (is(subsetTSNE[, dimx4], "factor") & dimy4 == "UMI.count") {
    #   subsetTSNE[, dimy4] <- Matrix::colSums(assays(scEx_log)[[1]][geneIdx, rownames(subsetTSNE), drop = FALSE])
    #   plotList[[plotIdx]] <- plotList[[plotIdx]] + geom_boxplot(show.legend = FALSE) + ggtitle(genesin[i])
    # } else {
    #   plotList[[plotIdx]] <- plotList[[plotIdx]] + geom_point(color = plotCol, show.legend = FALSE) + ggtitle(genesin[i])
    # }
    # if (!is.null(ylim)) {
    #   plotList[[plotIdx]] <- plotList[[plotIdx]] + ylim(ylim)
    # }
    # plotList[[plotIdx]] <- plotList[[plotIdx]] 
    # plot(subsetTSNE[, dimx4], subsetTSNE[, dimy4],
    #      col = plotCol, pch = 16, frame.plot = TRUE,
    #      ann = FALSE, ylim = ylim
    # )
    # title(genesin[i], line = -1.2, adj = 0.05, cex.main = 2)
    # if (DEBUG) cat(file = stderr(), ppgrp)
  }
  # }
  # 
  require(cowplot)
  printData = printData[order(printData$col, decreasing = F),]
  printData$gene = factor(printData$gene, levels = unique(genesin))
  if (is(subsetTSNE[, dimx4], "factor") & dimy4 == "UMI.count") {
    retVal <- ggplot(printData, aes(dimx, dimy)) + geom_boxplot(show.legend = FALSE)
  } else {
    retVal <- ggplot(printData, aes(dimx, dimy)) + geom_point(aes(colour = col), alpha = 0.4) 
  }
  if (sameScale){
    retVal = retVal + 
      facet_wrap(vars(gene),ncol = nCol) 
  } else {
    retVal = retVal + facet_wrap(vars(gene),ncol = nCol,scales = "free")
  }
  retVal = retVal + xlab(dimx4) + 
    ylab(dimy4) 
  # +
  #   scale_colour_gradient2()
  
  retVal =  retVal + scale_color_gradient2(low = lowCol, high = highCol, mid = midCol, midpoint = midFunc(colData))
  
  
  require(ggpubr)
  # retVal <-
  #   ggarrange(
  #     plotlist = plotList, ncol = nCol, nrow = ceiling(length(plotList) / nCol),
  #     label.x = "test", legend = "right",
  #     common.legend = T
  #   )
  retVal <-
    annotate_figure(retVal,
                    top = text_grob(sampdesc)
    )
  return(retVal)
}



