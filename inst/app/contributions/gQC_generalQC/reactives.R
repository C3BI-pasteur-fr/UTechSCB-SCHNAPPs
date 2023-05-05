suppressMessages(require(uwot))
# suppressMessages(require(tidyverse))
suppressMessages(require(SingleCellExperiment))
require(plotly)

# here we define reactive values/variables

# save to history violoin observer ----
observe(label = "save2histumi", {
  clicked  = input$save2Histumi
  deepDebug()
  if (DEBUG) cat(file = stderr(), "observe input$save2histumi \n")
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
  
  add2history(type = "renderPlot", input = isolate( reactiveValuesToList(input)), comment = "UMI histogram",  
              plotData = .schnappsEnv[["gQC_plotUmiHist"]])
  
})


# save to history save2HistSample observer ----
observe({
  deepDebug()
  clicked  = input$save2HistSample
  if (DEBUG) cat(file = stderr(), "observe input$save2HistSample \n")
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
  add2history(type = "renderPlot", input = isolate( reactiveValuesToList(input)), comment = "Sample histogram",  
              plotData = .schnappsEnv[["gQC_plotSampleHist"]])
  
})

# save to history save2HistSample observer ----
observe(label = "save2histvar", {
  deepDebug()
  clicked  = input$save2Histvar  
  if (DEBUG) cat(file = stderr(), "observe input$save2histvar \n")
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
  # add2history(type = "renderPlot", input = isolate( reactiveValuesToList(input)), comment = "PC variance",  
  #             plotData = .schnappsEnv[["gQC_variancePCA"]])
  add2history(type = "save", input = isolate( reactiveValuesToList(input)), 
              comment = paste0("# PC variance\n",
                               "fun = plotData$plotData$plotFunc\n", 
                               "environment(fun) = environment()\n",
                               "print(do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)]))\n"),
              plotData = .schnappsEnv[["gQC_variancePCA"]])
  
})

plotHistVarPC <- function(df, pc, var) {
  ggplot(data = df,aes(x=pc, y=var)) + geom_bar(stat = "identity")
}  


# gQC_scaterReadsFunc ----
#' gQC_scaterReadsFunc
#' calculate the QC metrix and return updated singleCellExperiment object
#' TODO make the 200 a parameter
gQC_scaterReadsFunc <- function(scEx) {
  if (DEBUG) cat(file = stderr(), "gQC_scaterReadsFunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_scaterReadsFunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_scaterReadsFunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_scaterReadsFunc", id = "gQC_scaterReadsFunc", duration = NULL)
  }
  
  if (is(assays(scEx)[[1]], "TsparseMatrix")) {
    assays(scEx)[["counts"]] <- as(assays(scEx)[["counts"]], "CsparseMatrix")
  }
  
  # ercc <- rownames(scEx)[grepl("ERCC-", rownames(scEx), ignore.case = TRUE)]
  #
  # mt <- rownames(scEx)[grepl("^MT", rownames(scEx), ignore.case = TRUE)]
  
  cellMet <- scater::perCellQCMetrics(
    scEx
  )
  featureMet <- scater::perFeatureQCMetrics(
    scEx
  )
  filter_by_expr_features <- (cellMet$detected > 200)
  scEx$use <- (
    # sufficient features (genes)
    filter_by_expr_features
    # sufficient molecules counted
    # filter_by_total_counts &
    # sufficient endogenous RNA
    # filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    # filter_by_MT
  )
  
  return(scEx)
}

# scaterReads ----
#' scaterReads
#' singleCellExperiment object/reactive with QC metrix
scaterReads <- reactive({
  if (DEBUG) cat(file = stderr(), "scaterReads started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scaterReads")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scaterReads")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scaterReads", id = "scaterReads", duration = NULL)
  }
  
  scEx <- scEx()
  # scEx_log = scEx_log()
  if (is.null(scEx) ) {
    return(NULL)
  }
  featureData <- rowData(scEx)
  rownames(scEx) = featureData$symbol
  retVal <- gQC_scaterReadsFunc(scEx)
  
  exportTestValues(scaterReads = {
    str(retVal)
  })
  return(retVal)
}) 


# gQC_sampleHistFunc ----
#' gQC_sampleHistFunc
#' create a histogram from samples
gQC_sampleHistFunc <- function(sampleInf, scols) {
  if (DEBUG) cat(file = stderr(), "gQC_sampleHistFunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_sampleHistFunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_sampleHistFunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_sampleHistFunc", id = "gQC_sampleHistFunc", duration = NULL)
  }
  # deepDebug()
  
  # if (!is.factor(sampleInf)) {
  sampleInf = factor(sampleInf)
  # }
  
  plotly::plot_ly(x=sampleInf, 
                  type="histogram", 
                  marker = list(color = scols[levels(sampleInf)])) %>%
    layout(
      title = "Number of cells per sample",
      xaxis = list(title = "Samples"),
      yaxis = list(title = "Count"))
  # counts <- table(sampleInf)
  # df <- as.data.frame(counts)
  # ggplot(data = df,aes(x=sampleInf, y=Freq, fill=sampleInf)) + geom_bar(stat = "identity")  + 
  #   scale_color_manual(values=scols) 
  # barplot(counts,
  #         main = "histogram of number of cell per sample",
  #         xlab = "Samples",
  #         col=scols
  # )
}


# projectionTable -----
#' projectionTable
#' input for tableSelectionServer
#' presents all projections
projectionTable <- reactive({
  if (DEBUG) cat(file = stderr(), "projectionTable started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "projectionTable")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "projectionTable")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("projectionTable", id = "projectionTable", duration = NULL)
  }
  
  projections <- projections()
  
  if (is.null(projections)) {
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/projectionTable.RData", list = c(ls()))
  }
  # load(file = "~/SCHNAPPsDebug/projectionTable.RData")
  
  printTimeEnd(start.time, "projectionTable")
  exportTestValues(projectionTable = {
    projections
  })
  return(projections)
})

# tsne ----
#' tsne
#' reactive calculating the tSNE projections

# should be a module :
# https://shiny.rstudio.com/articles/modules.html
# https://stackoverflow.com/questions/43976128/create-a-reactive-function-outside-the-shiny-app
#
# I guess that modules are self-contained and one cannot retrieve information from the parent session
# this means that I cannot create rectivity from within the module to inputs that are outside.

# output$updatetsneParametersButton <- updateButtonUI(name = "updatetsneParameters",
#                                                     variables = c("gQC_tsneDim", "gQC_tsnePerplexity", "gQC_tsneTheta", "gQC_tsneSeed"  ) )
# updateButton(name = "updatetsneParameters",
#                                                   )


runTSNEclicked <- reactive({
  runme = input$updatetsneParameters
  if(runme>0) return(1)
  return(0)
})

tsne <- reactive({
  if (DEBUG) cat(file = stderr(), "tsne started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "tsne")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "tsne")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tsne", id = "tsne", duration = NULL)
  }
  
  pca <- pcaReact()
  # only recalculate when button is pressed.
  input$updatetsneParameters
  gQC_tsneDim <- isolate(input$gQC_tsneDim)
  gQC_tsnePerplexity <- isolate(input$gQC_tsnePerplexity)
  gQC_tsneTheta <- isolate(input$gQC_tsneTheta)
  gQC_tsneSeed <- isolate(input$gQC_tsneSeed)
  
  if (is.null(pca)) {
    if (DEBUG) cat(file = stderr(), "tsne: NULL\n")
    return(NULL)
  }
  
  retVal <- tsneFunc(pca, gQC_tsneDim, gQC_tsnePerplexity, gQC_tsneTheta, gQC_tsneSeed)
  if (is.null(tsne)) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(paste("Problem with tsne:", e),
                       id = "gQC_tsneWarning",
                       type = "error",
                       duration = NULL
      )
    }
    return(NULL)
  }
  
  setRedGreenButton(
    vars = list(
      c("gQC_tsneDim", gQC_tsneDim),
      c("gQC_tsnePerplexity", gQC_tsnePerplexity),
      c("gQC_tsneTheta", gQC_tsneTheta),
      c("gQC_tsneSeed", gQC_tsneSeed)
    ),
    button = "updatetsneParameters"
  )
  
  exportTestValues(tsne = {
    retVal
  })
  return(retVal)
}) %>% bindCache(  pcaReact(),
                   runTSNEclicked(),
                   isolate(input$gQC_tsneDim),
                   isolate(input$gQC_tsnePerplexity),
                   isolate(input$gQC_tsneTheta),
                   isolate(input$gQC_tsneSeed)
)

# tsneFunc ----
tsneFunc <- function(pca, gQC_tsneDim, gQC_tsnePerplexity, gQC_tsneTheta, gQC_tsneSeed) {
  if (DEBUG) cat(file = stderr(), "tsneFunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "tsneFunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "tsneFunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tsneFunc", id = "tsneFunc", duration = NULL)
  }
  
  set.seed(seed = gQC_tsneSeed)
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/tsne.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/tsne.RData')
  suppressMessages(require(parallel))
  suppressMessages(require(Rtsne))
  np <- dim(pca$x)[2]
  tsne <- tryCatch(
    {
      Rtsne::Rtsne(
        pca$x[, 1:np],
        pca = FALSE, dims = gQC_tsneDim,
        perplexity = gQC_tsnePerplexity,
        theta = gQC_tsneTheta,
        check_duplicates = FALSE, num_threads = detectCores()
      )
    },
    error = function(e) {
      return(NULL)
    }
  )
  if (is.null(tsne)) {
    return(NULL)
  }
  retVal <- data.frame(tsne$Y)
  colnames(retVal) <- paste0("tsne", c(1:ncol(retVal)))
  return(retVal)
}

activatedUMAP <- reactive({
  inp = input$activateUMAP
  if(inp>0) return(1)
  return(0)
})
# umapReact ----
#' umapReact
#' reactive for calculating UMAP projection
umapReact <- reactive({
  if (DEBUG) cat(file = stderr(), "umapReact started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "umapReact")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "umapReact")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("umapReact", id = "umapReact", duration = NULL)
  }
  
  # xaxis <- input$um_xaxis
  # yaxis <- input$um_yaxis
  # cellT <- input$um_ct
  # inputCT <- input$um_inputCT
  # sampleRatio <- as.numeric(input$um_sampleRatio)
  # sampleIds <- input$um_sampleIds
  # InlevelOrd <- input$um_levelOrd
  # UMAP1 <- input$um_umap1
  # UMAP2 <- input$um_umap2
  runUMAP <- input$activateUMAP
  scEx_log <- scEx_log()
  pca <- pcaReact()
  
  myseed <- isolate(input$gQC_um_randSeed)
  n_neighbors <- isolate(as.numeric(input$gQC_um_n_neighbors))
  n_components <- isolate(as.numeric(input$gQC_um_n_components))
  n_epochs <- isolate(as.numeric(input$gQC_um_n_epochs))
  # alpha <- isolate(as.numeric(input$um_alpha))
  init <- isolate(input$gQC_um_init)
  min_dist <- isolate(as.numeric(input$gQC_um_min_dist))
  set_op_mix_ratio <- isolate(as.numeric(input$gQC_um_set_op_mix_ratio))
  local_connectivity <- isolate(as.numeric(input$gQC_um_local_connectivity))
  bandwidth <- isolate(as.numeric(input$gQC_um_bandwidth))
  gamma <- isolate(as.numeric(input$um_gamma))
  negative_sample_rate <- isolate(as.numeric(input$gQC_um_negative_sample_rate))
  metric <- isolate(input$gQC_um_metric)
  spread <- isolate(as.numeric(input$gQC_um_spread))
  
  if (is.null(scEx_log)) {
    if (DEBUG) cat(file = stderr(), "output$umap_react:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/umap_react.RData", list = c(ls()))
  }
  # load("~/SCHNAPPsDebug/umap_react.RData")
  if (!runUMAP) {
    if (DEBUG) cat(file = stderr(), "output$umap_react:NULL\n")
    return(NULL)
  }
  umapData <- as.matrix(assays(scEx_log)[[1]])
  compCases <- complete.cases(umapData)
  
  # TODO it might be possible to reuse nearest neighbor information to speeed up recomputations
  # with eg. new seed
  
  set.seed(myseed)
  # embedding <- uwot::umap(t(as.matrix(assays(scEx_log)[[1]])),
  embedding <- uwot::umap(pca$x,
                          n_neighbors = n_neighbors,
                          n_components = n_components,
                          n_epochs = n_epochs,
                          # alpha = alpha,
                          init = init,
                          spread = spread,
                          min_dist = min_dist,
                          set_op_mix_ratio = set_op_mix_ratio,
                          local_connectivity = local_connectivity,
                          bandwidth = bandwidth,
                          # gamma = gamma,
                          negative_sample_rate = negative_sample_rate,
                          metric = metric,
                          n_threads = detectCores()
  )
  embedding <- as.data.frame(embedding)
  colnames(embedding) <- paste0("UMAP", 1:n_components)
  rownames(embedding) <- colnames(scEx_log)
  
  setRedGreenButton(
    vars = list(
      c("gQC_um_randSeed", myseed),
      c("gQC_um_n_neighbors", n_neighbors),
      c("gQC_um_n_components", n_components),
      c("gQC_um_n_epochs", n_epochs),
      # c("um_alpha", alpha),
      c("gQC_um_init", init),
      c("gQC_um_min_dist", min_dist),
      c("gQC_um_set_op_mix_ratio", set_op_mix_ratio),
      c("gQC_um_local_connectivity", local_connectivity),
      c("gQC_um_bandwidth", bandwidth),
      c("um_gamma", gamma),
      c("gQC_um_negative_sample_rate", negative_sample_rate),
      c("gQC_um_metric", metric),
      c("gQC_um_spread", spread)
    ),
    button = "activateUMAP"
  )
  updateSelectizeInput(session = session, inputId = "gQC_umap_main-dimension_x", selected = "UMAP1")
  updateSelectizeInput(session = session, inputId = "gQC_umap_main-dimension_y", selected = "UMAP2")
  updateSelectizeInput(session = session, inputId = "gQC_umap_main-dimension_col", selected = "dbCluster")
  .schnappsEnv$defaultValues[["gQC_umap_main-dimension_x"]] <- "UMAP1"
  .schnappsEnv$defaultValues[["gQC_umap_main-dimension_y"]] <- "UMAP2"
  .schnappsEnv$defaultValues[["gQC_umap_main-dimension_col"]] <- "dbCluster"
  .schnappsEnv[["gQC_umap_main-dimension_x"]] <- "UMAP1"
  .schnappsEnv[["gQC_umap_main-dimension_y"]] <- "UMAP2"
  .schnappsEnv[["gQC_umap_main-dimension_col"]] <- "dbCluster"
  
  exportTestValues(umapReact = {
    embedding
  })
  return(embedding)
})%>% bindCache(activatedUMAP(),
                isolate(input$gQC_um_n_neighbors),
                isolate(input$gQC_um_n_components),
                isolate(input$gQC_um_n_epochs),
                isolate(input$gQC_um_min_dist),
                isolate(input$gQC_um_set_op_mix_ratio),
                isolate(input$gQC_um_local_connectivity),
                isolate(input$gQC_um_bandwidth),
                isolate(input$um_gamma),
                isolate(input$gQC_um_negative_sample_rate),
                isolate(input$gQC_um_spread),
                isolate(input$gQC_um_metric),
                isolate(input$gQC_um_randSeed),
                isolate(input$gQC_um_init),
                scEx_log(),
                pcaReact()
)


# myProjections ----
myProjections <- list(
  c("tsne", "tsne"),
  c("dbCluster", "dbCluster"),
  c("umap", "umapReact")
)

# myHeavyCalculations 
# declare function as heavy
# myHeavyCalculations <- list(
#   c("scaterReads", "scaterReads"),
#   c("tsne", "tsne")
# )


#' tsnePlot ----
#' function that plots in 3D the tsne projection
tsnePlot <- function(projections, dimX, dimY, dimZ, dimCol, projColors) {
  if (DEBUG) cat(file = stderr(), "tsnePlot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "tsnePlot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "tsnePlot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tsnePlot", id = "tsnePlot", duration = NULL)
    removeNotification(id = "tsnePlotERROR")
  }
  
  projections <- as.data.frame(projections)
  if (!all(c(dimX, dimY, dimZ) %in% colnames(projections))) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("Selected projections not available. Did you run normalization?", id = "tsnePlotERROR", type = "error", duration = NULL)
    }
    return(NULL)
  }
  
  projections$dbCluster <- as.factor(projections$dbCluster)
  
  myColors = projColors[[dimCol]]
  # if (dimCol == "sampleNames") {
  #   myColors <- scols
  # } else {
  #   myColors <- NULL
  # }
  # if (dimCol == "dbCluster") {
  #   myColors <- ccols
  # }
  
  p <-
    plotly::plot_ly(
      projections,
      x = formula(paste("~ ", dimX)),
      y = formula(paste("~ ", dimY)),
      z = formula(paste("~ ", dimZ)),
      type = "scatter3d",
      color = formula(paste("~ ", dimCol)),
      colors = myColors,
      hoverinfo = "text",
      text = paste("Cluster:", as.numeric(as.character(projections$dbCluster))),
      mode = "markers",
      source = "B",
      marker =
        list(
          line = list(width = 0),
          size = rep(10, nrow(projections)),
          sizeref = 3
        )
    )
  return(p)
}

# cp = load("/Volumes/LaCie2022/downloads/project.2023-04-03.RData")
# scEx_log = scEx

### DoubletFinder related -----

if("DoubletFinder" %in% installed.packages()){
  require(DoubletFinder)
  require(Seurat)
  require(SingleCellExperiment)
  
  # estimation of multiplet rate from a Poisson distribution
  # - mu = probability of loading a cell
  # - n = number of droplets containing at least one cell
  # - N = total number of droplets
  # - P(c>=1) = 1 - exp(-mu)
  # - P(c>=2) = 1 - exp(-mu) - mu*exp(-mu)
  # - P(c>=1) = n/N -> mu = -log(1-n/N)
  # - multiplet rate = P(c>=2)/P(c>=1)
  # - note that N was estimated by fitting the 10X data 
  multiplet_rate = function(n_recovered, total_droplets=6.8e4) {
    # total droplet estimation from 10X data
    #n_recovered = seq(1000, 10000, by=1000)
    #f_doublets = c(0.8, 1.5, 2.3, 3.0, 3.8, 4.6, 5.3, 6.1, 6.8, 8.1)/100
    # total_droplets ~ 6.8e4
    mu = -log(1-n_recovered/total_droplets)
    return(1-mu*exp(-mu)/(1-exp(-mu)))
  }
  
  # pN -> number of artificial doublets generated (expressed as fraction of merged real-artifical data)
  # pK -> neighborhood size to estimate likelihood of being doublet (expressed as fraction of merged real-artifical data)
  # n_recovered -> number of cells recovered by CellRanger (before QC) -> estimate from #cells in srt object by default?
  find_doublets = function(scEx, dims=1:20, n_recovered=10000, pK=20/ncol(srt), pN=0.25) {
    ## pK Identification (no ground-truth)
    # a) author suggested
    #my_params = paramSweep_v3(srt, PCs = dims, sct = FALSE)
    #my_params = summarizeSweep(my_params, GT = FALSE)
    #my_pk = find.pK(my_params) #0.005
    # b) my suggestion -> 20 nearest neighors (default)
    
    # create Seurat object
    require(Seurat)
    require(DoubletFinder)
    srt = CreateSeuratObject(counts = assays(scEx)[["counts"]])
    
    # cellMeta <- colData(scEx_log)
    # meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
    # seurDat <- CreateSeuratObject(
    #   RNA = assays(scEx_log)[[1]],
    #   meta.data = meta.data
    # )
    # Mm_H_E18_filt@meta.data
    # 
    ## Homotypic Doublet Proportion Estimate
    homotypic_prop = modelHomotypic(Idents(srt))
    nExp_poi = round(multiplet_rate(n_recovered)*ncol(srt))
    nExp_poi_adj = round(nExp_poi*(1-homotypic_prop))
    
    
    ## Run DoubletFinder
    pann_name = paste0("pANN_",pN, "_", pK, "_", nExp_poi_adj)
    df_name = paste0("DF.classifications_",pN,"_", pK, "_", nExp_poi_adj)
    srt = doubletFinder_v3(srt, PCs = dims, pN = pN, pK = pK, nExp = nExp_poi_adj,sct = T)
    result = srt[[c(pann_name, df_name)]]
    # colnames(result) = c("DF.score", "DF.class")
    return(result)
  }
  
}