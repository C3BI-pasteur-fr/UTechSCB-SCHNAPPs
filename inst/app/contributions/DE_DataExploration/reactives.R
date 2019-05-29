require(ggplot2)

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
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_scaterPNG")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_scaterPNG", id = "DE_scaterPNG", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "DE_scaterPNG\n")

  runScater <- input$runScater
  if (!runScater) {
    return(list(
      src = "",
      contentType = "image/png",
      width = 10,
      height = 10,
      alt = "Scater plot will be here when 'run scater' is checked"
    ))
  }
  scaterReads <- scaterReads()
  if (is.null(scaterReads)) {
    return(NULL)
  }


  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height

  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scater.Rmd", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/SCHNAPPsDebug/scater.Rmd')

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

  rownames(scaterReads) = rowData(scaterReads)$symbol
  p1 <- scater::plotHighestExprs(scaterReads, colour_cells_by = "log10_total_counts", n=n)
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
  
  printTimeEnd(start.time, "DE_scaterPNG")
  exportTestValues(DE_scaterPNG = {retVal})  
  return(retVal)
})

# DE_dataExpltSNEPlot ---
#' DE_dataExpltSNEPlot
#' plot 3D tSNE in DataExploration - Expression
#' Here, only the expression of a gene or gene list is shown, compared to the other tSNE plot 
#' in General QC - tSNE
DE_dataExpltSNEPlot <- function(scEx_log, g_id, projections) {
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
      x = ~tsne1,
      y = ~tsne2,
      z = ~tsne3,
      type = "scatter3d",
      hoverinfo = "text",
      text = ~ paste(1:nrow(projections), " ", rownames(projections), "<br />",
            "Cluster:", as.numeric(as.character(projections$dbCluster))),
      # text = paste("Cluster:", as.numeric(as.character(projections$dbCluster))),
      mode = "markers",
      marker = list(
        size = 2,
        line = list(width = 0),
        color = ~values,
        colors = "Greens"
      )
    )
  layout(p, title = paste(featureData[geneid, "symbol"], collapse = ", "))
}

DE_geneViolinFunc <- function(scEx_log, g_id, projections, ccols) {
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
      fun.y = median,
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
