if (DEBUG) cat(file = stderr(), "outputs.R started.\n")
suppressMessages(require(shinyTree))
suppressMessages(require(stringr))
# require(rintrojs)
# SUMMARY STATS ----------------------------------------------------------------
base::source(paste0(packagePath, "/moduleServer.R"), local = TRUE)

DEBUGSAVE <- get(".SCHNAPPs_DEBUGSAVE", envir = .schnappsEnv)
# normalizationRadioButtonValue --------------------------------
# Parameters / normalization
output$normalizationRadioButtonValue <- renderPrint({
  input$normalizationRadioButton
})

library(profvis)
callModule(profvis_server, "profiler")

normaliztionParameters <- list(raw = "no Parameters needed")
# localContributionDir <- .SCHNAPPs_locContributionDir
parFiles <-
  dir(
    path = c(paste0(packagePath, "/contributions"), localContributionDir),
    pattern = "parameters.R",
    full.names = TRUE,
    recursive = TRUE
  )
for (fp in parFiles) {
  myNormalizationParameters <- list()
  source(fp, local = TRUE)
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/normalizationsParameters.RData",
         list = c("normaliztionParameters", ls())
    )
  }
  # load(file = '~/SCHNAPPsDebug/normalizationsParameters.RData')
  if (length(myNormalizationParameters) > 0) {
    for (li in 1:length(myNormalizationParameters)) {
      lVal <- myNormalizationParameters[[li]]
      if (length(lVal) > 0) {
        # if (DEBUG) {
        #   cat(
        #     file = stderr(),
        #     paste(
        #       "normalization Choice: ",
        #       names(myNormalizationParameters)[li],
        #       " ",
        #       lVal,
        #       "\n"
        #     )
        #   )
        #   cat(file = stderr(), paste(
        #     "class: ",
        #     class(myNormalizationParameters[[li]]),
        #     " ",
        #     lVal,
        #     "\n"
        #   ))
        # }
        oldNames <- names(normaliztionParameters)
        normaliztionParameters[[length(normaliztionParameters) + 1]] <-
          lVal
        names(normaliztionParameters) <-
          c(oldNames, names(myNormalizationParameters)[li])
      }
    }
  }
}


observe(label ="obs_pcaRank", x = {
  .schnappsEnv$defaultValues[["pcaRank"]] = input$pcaRank
})
observe(label ="obs_normalizationRadioButton", x = {
  .schnappsEnv$defaultValues[["normalizationRadioButton"]] = input$normalizationRadioButton
})
observe(label ="obs_cellSelectionComment", x = {
  .schnappsEnv$defaultValues[["cellSelectionComment"]] = input$cellSelectionComment
})
observe(label ="obs_cellsFiltersOut", x = {
  .schnappsEnv$defaultValues[["cellsFiltersOut"]] = input$cellsFiltersOut
})
observe(label ="obs_cellKeepOnly", x = {
  .schnappsEnv$defaultValues[["cellKeepOnly"]] = input$cellKeepOnly
})
observe(label ="obs_cellKeep", x = {
  .schnappsEnv$defaultValues[["cellKeep"]] = input$cellKeep
})
observe(label ="obs_cellPatternRM", x = {
  .schnappsEnv$defaultValues[["cellPatternRM"]] = input$cellPatternRM
})
observe(label ="obs_maxGenes", x = {
  .schnappsEnv$defaultValues[["maxGenes"]] = input$maxGenes
})
observe(label ="obs_minGenes", x = {
  .schnappsEnv$defaultValues[["minGenes"]] = input$minGenes
})
observe(label ="obs_minNonExpGenes", x = {
  .schnappsEnv$defaultValues[["minNonExpGenes"]] = input$minNonExpGenes
})
observe(label ="obs_minExpGenes", x = {
  .schnappsEnv$defaultValues[["minExpGenes"]] = input$minExpGenes
})
observe(label ="obs_genesKeep", x = {
  .schnappsEnv$defaultValues[["genesKeep"]] = input$genesKeep
})
observe(label ="obs_minGenesGS", x = {
  .schnappsEnv$defaultValues[["minGenesGS"]] = input$minGenesGS
})
observe(label ="obs_selectIds", x = {
  .schnappsEnv$defaultValues[["selectIds"]] = input$selectIds
})
observe(label ="obs_whichscLog", x = {
  .schnappsEnv$defaultValues[["whichscLog"]] = input$whichscLog
})
observe(label ="obs_subsampleNum", x = {
  .schnappsEnv$defaultValues[["subsampleNum"]] = input$subsampleNum
})
observe(label ="obs_sampleInput", x = {
  .schnappsEnv$defaultValues[["sampleInput"]] = input$sampleInput
})

observe(label ="obs_simlr_maxClust", x = {
  .schnappsEnv$defaultValues[["simlr_maxClust"]] = input$simlr_maxClust
})
observe(label ="obs_simlr_nClust", x = {
  .schnappsEnv$defaultValues[["simlr_nClust"]] = input$simlr_nClust
})
observe(label ="obs_snnType", x = {
  .schnappsEnv$defaultValues[["snnType"]] = input$snnType
})
observe(label ="obs_snnClusterSource", x = {
  .schnappsEnv$defaultValues[["snnClusterSource"]] = input$snnClusterSource
})
observe(label ="obs_geneSelectionClustering", x = {
  .schnappsEnv$defaultValues[["geneSelectionClustering"]] = input$geneSelectionClustering
})
observe(label ="obs_useRanks", x = {
  .schnappsEnv$defaultValues[["useRanks"]] = input$useRanks
})
observe(label ="obs_minClusterSize", x = {
  .schnappsEnv$defaultValues[["minClusterSize"]] = input$minClusterSize
})
observe(label ="obs_clusterMethod", x = {
  .schnappsEnv$defaultValues[["clusterMethod"]] = input$clusterMethod
})
observe(label ="obs_clusterSource", x = {
  .schnappsEnv$defaultValues[["clusterSource"]] = input$clusterSource
})
observe(label ="obs_seurClustresolution", x = {
  .schnappsEnv$defaultValues[["seurClustresolution"]] = input$seurClustresolution
})
observe(label ="obs_seurClustk.param", x = {
  .schnappsEnv$defaultValues[["seurClustk.param"]] = input$seurClustk.param
})
observe(label ="obs_seurClustDims", x = {
  .schnappsEnv$defaultValues[["seurClustDims"]] = input$seurClustDims
})
observe(label ="obs_tabsetCluster", x = {
  .schnappsEnv$defaultValues[["tabsetCluster"]] = input$tabsetCluster
})
observe(label ="obs_genesRMPCA", x = {
  .schnappsEnv$defaultValues[["genesRMPCA"]] = input$genesRMPCA
})
observe(label ="obs_genes4PCA", x = {
  .schnappsEnv$defaultValues[["genes4PCA"]] = input$genes4PCA
})
observe(label ="obs_useSeuratPCA", x = {
  .schnappsEnv$defaultValues[["useSeuratPCA"]] = input$useSeuratPCA
})
observe(label ="obs_hvgSelection", x = {
  .schnappsEnv$defaultValues[["hvgSelection"]] = input$hvgSelection
})
observe(label ="obs_pcaScale", x = {
  .schnappsEnv$defaultValues[["pcaScale"]] = input$pcaScale
})
observe(label ="obs_pcaN", x = {
  .schnappsEnv$defaultValues[["pcaN"]] = input$pcaN
})


# check gene names ----
observe({
  scEx = scEx()
  req(scEx)
  if(any(stringr::str_detect( rownames(scEx), "_"))){
    showNotification(
      "gene names contain '_', which will be replaced by Seurat by '.', which can cause artefacts",
      type = "error",
      duration = NULL
    )
  }
})


output$noLogWarning <- renderText({
  logCalc <- input$whichscLog
  if(logCalc != "calcLog") return("Warning normalization not being calculated due to input page selection of Compute normalizations?")
  return("")
})

# dimPlotPCA ----
# <- reactive({
output$dimPlotPCA <- renderPlot({
  if (DEBUG) {
    cat(file = stderr(), "dimPlotPCA started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "dimPlotPCA")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "dimPlotPCA")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("dimPlotPCA", id = "dimPlotPCA", duration = NULL)
  }
  
  input$updateDimPlot
  scEx_log <- isolate(scEx_log())
  scEx <- isolate(scEx())
  pca <- isolate(pcaReact())
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "dimPlotPCA:NULL\n")
    }
    return(0)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/dimPlotPCA.RData", list = c(ls()))
  }
  # cp =load(file='~/SCHNAPPsDebug/dimPlotPCA.RData')
  
  # return NuLL because it is not working correctly
  # return(NULL)
  
  scEx = scEx[rownames(pca$rotation),]
  scEx_log = scEx_log[rownames(pca$rotation),]
  
  cellMeta = colData(scEx_log)
  rData = rowData(scEx)
  meta.data = cellMeta[,"sampleNames", drop = FALSE]
  dat = assays(scEx)[[1]][rownames(scEx_log),]
  rownames(dat) = rData[rownames(scEx_log),"symbol"]
  rownames(pca$rotation) = rData[rownames(pca$rotation),"symbol"]
  seurDat <- CreateSeuratObject(
    counts = dat,
    meta.data = meta.data
  )
  
  # TODO use scEx_log
  logDat = assays(scEx_log)[[1]]
  rData = rowData(scEx_log)
  rownames(logDat) = rData$symbol
  seurDat@assays$RNA@data = as(logDat,"CsparseMatrix")
  # seurDat <- NormalizeData(seurDat, normalization.method = "LogNormalize", scale.factor = 10000)
  # seurDat <- FindVariableFeatures(seurDat, selection.method = "vst", nfeatures = 2000)
  
  # recalculating because createDimReducObject is not working
  all.genes <- rownames(seurDat)
  # seurDat <- ScaleData(seurDat, features = all.genes)
  # seurDat <- RunPCA(seurDat, features = VariableFeatures(object = seurDat))
  
  colnames(pca$x) = str_replace(colnames(pca$x), "PC", "PC_")
  ndim = min(15,ncol(pca$x))
  # pca.res = irlba(A=t(x=seurDat@assays$RNA@data), nv=50)
  # not working
  seurDat[["pca"]] = CreateDimReducObject(embeddings = pca$x[colnames(seurDat),], 
                                          loadings = pca$rotation, 
                                          stdev = pca$var_pcs, 
                                          key = "PC_", 
                                          assay = "RNA")
  # seurDat <- ProjectDim(object = seurDat, reduction = "pca", assay = "RNA")
  
  # DimPlot(seurDat, reduction = "pca")
  # seurDat <- ProjectDim(seurDat, reduction = "pca", assay = 'RNA')
  
  d = DimHeatmap(seurDat, dims = 1:ndim, slot = 'data',
                 balanced = TRUE, fast = TRUE, projected = FALSE, 
                 reduction = "pca")
  d
})


# normalizationsParametersDynamic -------------------------
output$normalizationsParametersDynamic <- renderUI({
  if (is.null(input$normalizationRadioButton)) {
    return(NULL)
  }
  selectedChoice <- input$normalizationRadioButton
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/normalizationsParametersDynamic.RData",
         list = c("normaliztionParameters", ls())
    )
  }
  # load(file = '~/SCHNAPPsDebug/normalizationsParametersDynamic.RData')
  do.call("switch",
          args = c(
            selectedChoice,
            normaliztionParameters,
            h3("no parameters provided")
          )
  )
})


cellSelectionValues <- reactiveVal(
  list(
    minExpGenes = defaultValue("minExpGenes", defaultValueRegExGene),
    minGenes = defaultValue("minGenes", 20),
    maxGenes = defaultValue("maxGenes", 1000000),
    cellPatternRM = defaultValue("cellPatternRM", ""),
    cellKeep = defaultValue("cellKeep", ""),
    cellKeepOnly = defaultValue("cellKeepOnly", ""),
    cellsFiltersOut = defaultValue("cellsFiltersOut", ""),
    minNonExpGenes = defaultValue("minNonExpGenes", "")
  )
)
geneSelectionValues <- reactiveVal(
  {
    list(
      selectIds = defaultValue("selectIds", "^MT-|^RP|^MRP"),
      geneListSelection = defaultValue("geneListSelection", NULL),
      minGenesGS = defaultValue("minGenesGS", 2),
      genesKeep = defaultValue("genesKeep", "")
    )
  }
)

observeEvent(
  label = "ob20",
  eventExpr = input$updateCellSelectionParameters,
  handlerExpr = {
    deepDebug()
    if (DEBUG) cat(file = stderr(), "observe updateCellSelectionParameters\n")
    cellSelectionValues(list(
      minExpGenes = input$minExpGenes,
      minGenes = input$minGenes,
      maxGenes = input$maxGenes,
      cellPatternRM = input$cellPatternRM,
      cellKeep = input$cellKeep,
      cellKeepOnly = input$cellKeepOnly,
      cellsFiltersOut = input$cellsFiltersOut,
      minNonExpGenes = input$minNonExpGenes
    ))
    if (DEBUG) cat(file = stderr(), "\nCellSelectionValues\n")
    updateButtonColor(buttonName = "updateCellSelectionParameters", parameters = c(
      "minExpGenes", "minGenes", "minNonExpGenes", 
      "maxGenes", "cellPatternRM", "cellKeep", "cellKeepOnly", "cellsFiltersOut"
    ))
    
  }
)

observe(label = "ob_cellSelection",
        {
          deepDebug()
          if (DEBUG) cat(file = stderr(), "observe ob_cellSelection\n")
          setRedGreenButtonCurrent(
            vars = list(
              c("minExpGenes", input$minExpGenes),
              c("minGenes", input$minGenes),
              c("maxGenes", input$maxGenes),
              c("cellPatternRM", input$cellPatternRM),
              c("cellKeep", input$cellKeep),
              c("cellKeepOnly", input$cellKeepOnly),
              c("cellsFiltersOut", input$cellsFiltersOut),
              c("minNonExpGenes", input$minNonExpGenes)
            )
          )
          
          updateButtonColor(buttonName = "updateCellSelectionParameters", parameters = c(
            "minExpGenes", "minGenes", "minNonExpGenes", 
            "maxGenes", "cellPatternRM", "cellKeep", "cellKeepOnly", "cellsFiltersOut"
          ))
        })

# observe: clustering Button ----
ob_clusteringParams <- observe(label = "ob_clusteringParams", {
  deepDebug()
  if (DEBUG) cat(file = stderr(), "observe ob_clusteringParams\n")
  
  # this happens when the lite version is used
  if (is.null(input$tabsetCluster)){
    ob_clusteringParams$destroy()
    return(NULL)
  }
  
  input$updateClusteringParameters
  whichClustering = isolate(input$tabsetCluster)
  req(whichClustering)
  
  if ( whichClustering == "scran_Cluster"){
    setRedGreenButtonCurrent(
      vars = list(
        c("seed", input$seed),
        c("useRanks", input$useRanks),
        c("clusterSource", clusterMethodReact$clusterSource),
        c("geneSelectionClustering", input$geneSelectionClustering),
        c("minClusterSize", input$minClusterSize),
        c("clusterMethod", clusterMethodReact$clusterMethod)
      )
    )
    
    updateButtonColor(buttonName = "updateClusteringParameters", parameters = c(
      "seed", "useRanks", "minClusterSize", "clusterMethod",
      "clusterSource", "geneSelectionClustering"
    ))
  }
})



observeEvent(
  label = "ob21",
  eventExpr = input$updateGeneSelectionParameters,
  handlerExpr = {
    deepDebug()
    geneSelectionValues(list(
      selectIds = input$selectIds,
      geneListSelection = input$geneListSelection,
      minGenesGS = input$minGenesGS,
      genesKeep = input$genesKeep
    ))
    if (DEBUG) cat(file = stderr(), "\ngeneSelectionValues\n")
    updateButtonColor(buttonName = "updateGeneSelectionParameters", parameters = c(
      "selectIds", "geneListSelection",
      "minGenesGS", "genesKeep"
    ))
    
  }
)

observe(label = "ob_geneSelection", 
        {
          deepDebug()
          if (DEBUG) cat(file = stderr(), "observe ob_geneSelection\n")
          setRedGreenButtonCurrent(
            vars = list(
              c("selectIds", input$selectIds),
              c("geneListSelection", input$geneListSelection),
              c("minGenesGS", input$minGenesGS),
              c("genesKeep", input$genesKeep)
            )
          )
          
          updateButtonColor(buttonName = "updateGeneSelectionParameters", parameters = c(
            "selectIds", "geneListSelection",
            "minGenesGS", "genesKeep"
          ))
          
        })

# summaryStatsSideBar -----------------------------
output$summaryStatsSideBar <- renderUI({
  if (DEBUG) {
    cat(file = stderr(), "output$summaryStatsSideBar\n")
  }
  scEx <- scEx()
  scEx_log <- scEx_log()
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "output$summaryStatsSideBar:NULL\n")
    }
    return(NULL)
  }
  # if (input$noStats) {
  #   if (DEBUG) {
  #     cat(file = stderr(), "output$summaryStatsSideBar:off\n")
  #   }
  #   return(NULL)
  # }
  annFile <- inputFile$annFile
  medianUMI <- medianUMI()
  medianENSG <- medianENSG()
  memoryUsed <- getMemoryUsed()
  infile <- inputFile$inFile
  normalizationRadioButton <- input$normalizationRadioButton
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/summaryStatsSideBar.RData",
         list = c("normaliztionParameters", ls())
    )
  }
  # load("~/SCHNAPPsDebug/summaryStatsSideBar.RData")
  line0 <- paste(infile, " _ ", annFile)
  line0a <- paste("Number of samples: ", length(levels(scEx$sampleNames)), sep = "\t")
  line1 <- paste("No. of cells: ", dim(scEx)[2], sep = "\t")
  line2 <- paste("No. of genes: ", dim(scEx)[1], sep = "\t")
  line1a <- paste("No. of cells (log): ", dim(scEx_log)[2], sep = "\t")
  line2a <- paste("No. of genes (log): ", dim(scEx_log)[1], sep = "\t")
  line3 <- paste("Median UMIs per cell: ", medianUMI, sep = "\t")
  line4 <-
    paste("Median Genes with min 1 UMI: ", medianENSG, sep = "\t")
  line5 <-
    paste("Total number of reads: ", sum(assays(scEx)[["counts"]]))
  line6 <- paste("Memory used:", memoryUsed)
  line7 <-
    paste("Normalization selected:", normalizationRadioButton)
  htmlOut <- paste0(
    "<br/>",
    "<br/>",
    "<br/>",
    "Summary statistics of this dataset:",
    "<br/>",
    "<br/>",
    line0,
    "<br/>",
    line0a,
    "<br/>",
    line1,
    "<br/>",
    line2,
    "<br/>",
    line1a,
    "<br/>",
    line2a,
    "<br/>",
    line3,
    "<br/>",
    line4,
    "<br/>",
    line5,
    "<br/>",
    # line6,
    # "<br/>",
    line7,
    "<br/>",
    "<br/>"
    
  )
  exportTestValues(summaryStatsSideBar = {
    htmlOut
  })
  
  HTML(htmlOut)
})

if ("shinyBS" %in% rownames(installed.packages())) {
  addPopover(
    session = session, id = "summaryStatsSideBar", title = "Data summary",
    content = "<ul><li>medium UMI: shows how many genes are  expressed in log2 space of normalized data</li> </ul> ",
    trigger = "click", options = list(container = "body")
  )
}
# Select Genes ----
# this is part of the basic functionality from this
# tools and thus, can stay in this file.
output$geneListSelection <- shinyTree::renderTree({
  geneLists
})

# selectedGenesTable ----
# ONOFF TAB RENDER TABLE ALL CELLS
# TODO module for DT this is part
# of the basic functionality from this tools and thus, can stay in this file.
# output$selectedGenesTable <- DT::renderDataTable({
#   if (DEBUG) {
#     cat(file = stderr(), "output$selectedGenesTable\n")
#   }
#   dataTables <- inputData()
#   useGenes <- useGenes()
#   useCells <- useCells()
#   minGenes <- input$minGenesGS
#
#   if (is.null(dataTables) | is.null(useGenes) | is.null(useCells)) {
#     return(NULL)
#   }
#   if (.schnappsEnv$DEBUGSAVE) {
#     save(file = "~/SCHNAPPsDebug/selectedGenesTable.RData",
#       list = c("normaliztionParameters", ls())
#     )
#   }
#   # load("~/SCHNAPPsDebug/selectedGenesTable.RData")
#
#   scEx <- assays(dataTables$scEx)[[1]]
#   fd <- rowData(dataTables$scEx)
#   dt = fd[useGenes,]
#   dt$rowSums <- Matrix::rowSums(scEx[useGenes, useCells])
#   dt$rowSamples <- Matrix::rowSums(scEx[useGenes, useCells] > 0)
#   # get the order of the frist two columns correct
#   firstCol = which(colnames(dt) == "symbol")
#   firstCol = c(firstCol, which(colnames(dt) == "Description"))
#   # those we created so we know they are there
#   firstCol = firstCol = c(firstCol,which (colnames(dt) %in% c("rowSums", "rowSamples")))
#   colOrder = c(firstCol, (1:ncol(dt))[-firstCol])
#   dt <- dt[, colOrder]
#   dt <- dt[dt$rowSums >= minGenes, ]
#   exportTestValues(selectedGenesTable = {
#     as.data.frame(dt)
#   })
#   DT::datatable(as.data.frame(dt),
#                 options = list(scrollX = TRUE))
# })

# removedGenesTable --------------------------
# TODO module for DT TODO move to were it belongs
# output$removedGenesTable <- DT::renderDataTable({
#   if (DEBUG) {
#     cat(file = stderr(), "output$removedGenesTable\n")
#   }
#   dataTables <- inputData()
#   useGenes <- useGenes()
#   useCells <- useCells()
#   minGenes <- input$minGenesGS
#
#     if (is.null(dataTables) | is.null(useGenes) | is.null(useCells)) {
#     return(NULL)
#   }
#   useGenes <- !useGenes
#
#   if (.schnappsEnv$DEBUGSAVE) {
#     save(file = "~/SCHNAPPsDebug/removedGenesTable.RData",
#       list = c("normaliztionParameters", ls())
#     )
#   }
#   # load("~/SCHNAPPsDebug/removedGenesTable.RData")
#   scEx <- assays(dataTables$scEx)[[1]]
#   fd <- rowData(dataTables$scEx)
#   dt <- fd[useGenes, c("symbol", "Description")]
#   dt$rowSums <- Matrix::rowSums(scEx[useGenes, useCells])
#   dt$rowSamples <- Matrix::rowSums(scEx[useGenes, useCells] > 0)
#
#   dt <- dt[dt$rowSums < minGenes, ]
#   exportTestValues(removedGenesTable = {
#     as.data.frame(dt)
#   })
#   DT::datatable(as.data.frame(dt))
# })

# gsSelectedGenes ---------------------------
# TODO module of DT with selected names above Print names of selected genes for gene
# selection above table

# gsSelectedGenesMod ----
callModule(
  tableSelectionServer,
  "gsSelectedGenesMod",
  gsSelectedGenesTable, caption = "Tables with used genes"
)

callModule(
  tableSelectionServer,
  "gsRMGenesMod",
  gsRMGenesTable, caption = "Tables with removed genes"
)

callModule(
  tableSelectionServer,
  "PCAloadingsMod",
  PCAloadingsTable, caption = "Tables with PCA loadings"
)

callModule(
  tableSelectionServer,
  "HVAinfoMod",
  HVAinfoTable, caption = "Tables with variable feature info"
)


# DEBUGSAVEstring ----
output$DEBUGSAVEstring <- renderText({
  if (DEBUG) {
    .schnappsEnv$DEBUGSAVE <- input$DEBUGSAVE
    DEBUGSAVE <- input$DEBUGSAVE
  } else {
    NULL
  }
})

# output$currentTabInfo <- renderText({
#   # deepDebug()
#   str(input$sideBarID)
# })

# output$save2Historystring <- renderText({
#   if (DEBUG) {
#     .schnappsEnv$saveHistorycheckbox <- input$save2History
#     saveHistorycheckbox <- input$save2History
#   } else {
#     NULL
#   }
# })

# cellSelectionMod ----
callModule(tableSelectionServer, 
           "cellSelectionMod", 
           inputSample, caption = "Table with input cells")

# normalizationResult ----
callModule(
  tableSelectionServer,
  "normalizationResult",
  scExLogMatrixDisplay, caption = "Tables with normalization results"
)

# descriptionOfWork ----
output$descriptOfWorkOutput <- renderPrint({
  input$descriptionOfWork
})

# # sampleColorSelection ----
# output$sampleColorSelection <- renderUI({
#   scEx <- scEx()
#   sampCol <- sampleCols$colPal
#   prFct = projFactors()
#   projections = projections()
#   
#   if (is.null(scEx)) {
#     return(NULL)
#   }
#   if (.schnappsEnv$DEBUGSAVE) {
#     save(file = "~/SCHNAPPsDebug/sampleColorSelection.RData",
#       list = c("normaliztionParameters", ls())
#     )
#   }
#   # cp = load("~/SCHNAPPsDebug/sampleColorSelection.RData")
#   
#   lev <- levels(colData(scEx)$sampleNames)
#   # cols <- gg_fill_hue(length(lev))
#   
#   # New IDs "colX1" so that it partly coincide with input$select...
#   lapply(seq_along(lev), function(i) {
#     colourpicker::colourInput(
#       inputId = paste0("sampleNamecol", lev[i]),
#       label = paste0("Choose colour for sample ", "\"", lev[i], "\""),
#       # value = "#762A83"
#       # ,
#       value = sampCol[i],
#       allowedCols = allowedColors,
#       palette = "limited"
#     )
#   })
# })
# sampleColorSelection ----
output$ColorSelection <- renderUI({
  scEx <- scEx()
  sampCol <- sampleCols$colPal
  prFct = projFactors()
  projections = projections()
  clusterCol <- clusterCols$colPal
  
  if (is.null(scEx)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sampleColorSelection.RData",
         list = c("normaliztionParameters", ls())
    )
  }
  # cp = load("~/SCHNAPPsDebug/sampleColorSelection.RData")
  
  lev <- levels(colData(scEx)$sampleNames)
  # cols <- gg_fill_hue(length(lev))
  lev1 <- levels(projections$dbCluster)
  lev2 <- levels(colData(scEx)$sampleNames)
  # deepDebug()
  tmpFun <- function(name = "Sample", value = "SampleColorPanel", lev = lev2, idStr = "sampleNamecol", sampCol, allowedColors){
    tabPanel(
      name, value = value,
      fluidRow(
        column(
          width = 6,
          lapply(seq_along(lev), function(i) {
            colourpicker::colourInput(
              inputId = paste0(idStr, lev[i]),
              label = paste0("Choose color for ",name,  "\"", lev[i], "\""),
              # value = "#762A83"
              # ,
              value = sampCol[i],
              allowedCols = allowedColors,
              palette = "limited"
            )
          })
        )))
  }
  
  tabs = list( 
    tmpFun(name = "Sample", value = "SampleColorPanel", lev = lev2, idStr = "sampleNamecol", sampCol, allowedColors),
    tmpFun(name = "Cluster", value = "ClusterColorPanel", lev = lev1, idStr = "clusterNamecol", clusterCol, allowedColors)
  )
  do.call(tabsetPanel, tabs)
})

# # clusterColorSelection ----
# output$clusterColorSelection <- renderUI({
#   scEx <- scEx()
#   projections <- projections()
#   clusterCol <- clusterCols$colPal
#   
#   if (is.null(scEx) || is.null(projections)) {
#     return(NULL)
#   }
#   if (.schnappsEnv$DEBUGSAVE) {
#     save(file = "~/SCHNAPPsDebug/clusterColorSelection.RData",
#       list = c("normaliztionParameters", ls())
#     )
#   }
#   # load("~/SCHNAPPsDebug/clusterColorSelection.RData")
#   
#   lev <- levels(projections$dbCluster)
#   # cols <- gg_fill_hue(length(lev))
#   
#   # New IDs "colX1" so that it partly coincide with input$select...
#   lapply(seq_along(lev), function(i) {
#     colourpicker::colourInput(
#       inputId = paste0("clusterNamecol", lev[i]),
#       label = paste0("Choose colour for cluster ", "\"", lev[i], "\""),
#       # value = "#762A83"
#       # ,
#       value = clusterCol[i],
#       allowedCols = allowedColors,
#       palette = "limited"
#     )
#   })
# })

# history store to file ----
#' 

# askComment <- function(failed = FALSE) {
#   modalDialog(
#     sc_textInput("HistComment", "add a comment", value = paste("created at ",date())),
#     footer = tagList(
#       modalButton("Cancel"),
#       actionButton("HistCommentok", "OK")
#     )
#   )
# }
# observeEvent(input$HistCommentok, {
#   if (DEBUG) {
#     cat(file = stderr(), "writing history.\n")
#   }
#   start.time <- base::Sys.time()
#   on.exit({
#     printTimeEnd(start.time, "HistCommentok")
#     if (!is.null(getDefaultReactiveDomain())) {
#       removeNotification(id = "HistCommentok")
#     }
#   })
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("writing history", id = "HistCommentok", duration = NULL)
#   }
#   
#   panelLinkHistory = list("coexpressionSelected" = "coE")
#   id <- input$sideBarID
#   cat(file = stderr(), paste0("observeEvent input$save2History\n"))
#   save(file = "~/SCHNAPPsDebug/save2History.RData", list = c(ls()))
#   # cp =load(file="~/SCHNAPPsDebug/save2History.RData")
#   lsS = ls(envir = .schnappsEnv)
#   for (pl in lsS[grep(paste0("^historyPlot-",panelLinkHistory[[id]]), lsS)]) {
#     cat(file = stderr(), paste0("writing to history: ",pl ,"\n"))
#     sp <- strsplit(  pl, "-" )[[1]]
#     recHistory(sp[[length(sp)]], .schnappsEnv[[pl]], envir = .schnappsEnv)
#     
#   }
#   
#   removeModal()
#   
# })
# 
# observeEvent(input$save2History, {
#   showModal(askComment())
# })

# observe: input$updateColors ----
observeEvent(
  label = "ob22",
  eventExpr = input$updateColors,
  handlerExpr = {
    deepDebug()
    cat(file = stderr(), paste0("observeEvent input$updateColors\n"))
    scExx <- scEx()
    projections <- projections()
    
    if (is.null(scExx) || is.null(projections)) {
      return(NULL)
    }
    # sample colors
    scols <- sampleCols$colPal
    
    inCols <- list()
    lev <- levels(colData(scExx)$sampleNames)
    
    inCols <- lapply(seq_along(lev), function(i) {
      input[[paste0("sampleNamecol", lev[i])]]
    })
    names(inCols) <- lev
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/updateColors.RData", list = c(ls()))
      cat(file = stderr(), paste0("observeEvent save done\n"))
    }
    # load(file="~/SCHNAPPsDebug/updateColors.RData")
    
    # isolate({
    sampleCols$colPal <- unlist(inCols)
    add2history(type = "save", input=isolate( reactiveValuesToList(input)), comment = "scol", scol = sampleCols$colPal)
    # })
    
    # cluster colors
    ccols <- clusterCols$colPal
    
    inCols <- list()
    lev <- levels(projections$dbCluster)
    
    inCols <- lapply(seq_along(lev), function(i) {
      input[[paste0("clusterNamecol", lev[i])]]
    })
    names(inCols) <- lev
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/updateColors2.RData", list = c(ls()))
      cat(file = stderr(), paste0("observeEvent 2 save done\n"))
    }
    # load(file="~/SCHNAPPsDebug/updateColors2.RData")
    
    # isolate({
    clusterCols$colPal <- unlist(inCols)
    add2history(type = "save", input=isolate( reactiveValuesToList(input)), comment = "ccol", ccol = clusterCols$colPal)
    
    # })
    setRedGreenButton(
      vars = list(
        c("sampleNamecol", sampleCols$colPal),
        c("clusterCols", clusterCols$colPal)
      ),
      button = "updateColors"
    )
  }
)

# observe: color selection----
observeEvent(eventExpr = input$updateColors, label = "ob_colorParams", {
  deepDebug()
  if (DEBUG) cat(file = stderr(), "observe color Vars\n")
  
  scEx <- scEx()
  projections <- projections()
  if (is.null(scEx) || is.null(projections)) {
    return(NULL)
  }
  
  lev <- levels(projections$dbCluster)
  ccols <- lapply(seq_along(lev), function(i) {
    input[[paste0("clusterNamecol", lev[i])]]
  })
  lev <- levels(colData(scEx)$sampleNames)
  scols <- lapply(seq_along(lev), function(i) {
    input[[paste0("sampleNamecol", lev[i])]]
  })
  setRedGreenButtonCurrent(
    vars = list(
      c("sampleNamecol", unlist(scols)),
      c("clusterCols", unlist(ccols))
    )
  )
  
  updateButtonColor(buttonName = "updateColors", parameters = c(
    "sampleNamecol", "clusterCols"
  ))
})

# Nclusters ----
output$Nclusters <- renderText({
  dbCluster <- dbCluster()
  if (is.null(dbCluster)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/Nclusters.RData", list = c(ls()))
    cat(file = stderr(), paste0("observeEvent save done\n"))
  }
  # load(file="~/SCHNAPPsDebug/Nclusters.RData")
  retVal <- paste(levels(dbCluster))
  exportTestValues(Nclusters = {
    retVal
  })
  return(retVal)
})

# download handler countscsv ----
output$countscsv <- downloadHandler(
  filename = paste0("counts.", Sys.Date(), ".csv"),
  content = function(file) {
    if (DEBUG) {
      cat(file = stderr(), "RDSsave started.\n")
    }
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "RDSsave")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "RDSsave")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("RDSsave", id = "RDSsave", duration = NULL)
    }
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "RDSsave")
    }
    
    scEx_log <- scEx_log()
    if (is.null(scEx_log)) {
      return(NULL)
    }
    write.csv(as.matrix(assays(scEx_log)[[1]]), file)
  }
)

# download RDS ----
output$RDSsave <- downloadHandler(
  filename = paste0("project.", Sys.Date(), ".RData"),
  content = function(file) {
    if (DEBUG) {
      cat(file = stderr(), "RDSsave started.\n")
    }
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "RDSsave")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "RDSsave")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("RDSsave", id = "RDSsave", duration = NULL)
    }
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "RDSsave")
    }
    
    
    # umaps???
    scEx <- scEx()
    projections <- projections()
    scEx_log <- scEx_log()
    pca <- pcaReact()
    # TODO should be taken from projections
    tsne <- tsne()
    ccol = clusterCols$colPal
    scol = sampleCols$colPal
    namesDF = groupNames$namesDF
    
    if (is.null(scEx)) {
      return(NULL)
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/RDSsave.RData", list = c(ls()))
    }
    # load(file='~/SCHNAPPsDebug/RDSsave.RData')
    deepDebug()
    scEx <- consolidateScEx(scEx, projections, scEx_log, pca, tsne)
    
    # we save the pca separately because I don't know how to store the rotation  otherwise.
    # mostly done to make the lite version work.
    
    saveList =  c("scEx" ,  "scol" , "ccol" , "namesDF")
    if(!is.null(pca)){
      saveList = c(saveList, "pca")
    }
    # deepDebug()
    # save projections that shouldn't be recalculated in lite version
    if (length(.schnappsEnv$projectionFunctions) > 0){
      for (idx in 1:length(.schnappsEnv$projectionFunctions) ){
        assign(.schnappsEnv$projectionFunctions[[idx]][2],
               eval(parse(text = paste0(.schnappsEnv$projectionFunctions[[idx]][2],"()"))))
        saveList = c(saveList, .schnappsEnv$projectionFunctions[[idx]][2])
      }
    }
    
    save(file = file, list = saveList)
    
    # write.csv(as.matrix(exprs(scEx)), file)
  }
)

# download Rmd ----
output$RmdSave <- downloadHandler(
  filename = "report.zip",
  content = function(outZipFile) {
    if (DEBUG) {
      cat(file = stderr(), "RmdSave started.\n")
    }
    start.time <- base::Sys.time()
    on.exit({
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "RmdSave")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("RmdSave", id = "RmdSave", duration = NULL)
    }
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "RmdSave")
    }
    
    if (is.null(.schnappsEnv$historyPath) ) {
      return(NULL)
    }
    # if (is.null(scEx)) {
    #   return(NULL)
    # }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/RmdSave.RData", list = c(ls(), ".schnappsEnv"))
    }
    # cp = load(file='~/SCHNAPPsDebug/RmdSave.RData')
    
    tempReport = .schnappsEnv$historyFile
    outReport = paste0(.schnappsEnv$historyPath, .Platform$file.sep,"report.html")
    # tempReport = "~/SCHNAPPsDebug/tmpReport.Rmd"
    # file.copy("contributions/gQC_generalQC//report.Rmd",
    #           '/var/folders/tf/jwlc7r3d48z7pkq0w38_v7t40000gp/T//RtmpTx4l4G/file1a6e471a698.Rmd', overwrite = TRUE)
    tryCatch(
      callr::r(
        function(inputFP, output_file, params, envir) {
          rmarkdown::render(
            input = inputFP,
            output_file = output_file,
            params = params,
            envir = envir
          )
        },
        args = list(
          inputFP = tempReport,
          output_file = outReport,
          # params = params,
          envir = new.env()
        ),
        stderr = stderr(),
        stdout = stderr()
      ),
      error = function(e) {
        cat(file = stderr(), paste("==== An error occurred during the creation of the report\n", e, "\n"))
      }
    )
    # file.copy(from = "contributions/sCA_subClusterAnalysis/report.Rmd",
    #           to = "/var/folders/_h/vtcnd09n2jdby90zkb6wyd740000gp/T//Rtmph1SRTE/file69aa37a47206.Rmd", overwrite = TRUE)
    # rmarkdown::render(input = tempReport, output_file = "report.html",
    #                   params = params, envir = new.env())
    
    # outZipFile <- paste0(tempdir(), .Platform$file.sep, "report.zip")
    
    # tDir <- paste0(.schnappsEnv$reportTempDir, .Platform$file.sep)
    # zippedReportFiles <- c(paste0(tDir, zippedReportFiles))
    zip(outZipFile, paste0(path.expand(.schnappsEnv$historyPath), .Platform$file.sep), flags = "-9jr")
    if (DEBUG) {
      end.time <- Sys.time()
      cat(
        file = stderr(),
        "===Report:done",
        difftime(end.time, start.time, units = "min"),
        "\n"
      )
    }
    return(outZipFile)
  }
)


# Report creation ------------------------------------------------------------------
output$report <- downloadHandler(
  filename = "report.zip",
  
  content = function(outZipFile) {
    outrepFile <- reacativeReport()
    file.copy(from = outrepFile, to = outZipFile)
  }
)

# dummy function to return NULL
returnNull <- function() {
  return(NULL)
}

# commented out because it is corrently not used
# # forceCalc -----# handling expensive calcualtions
# forceCalc <- shiny::observe({
#   if (DEBUG) cat(file = stderr(), paste0("observe: goCalc\n"))
#   go <- input$goCalc
#   start.time <- base::Sys.time()
#   if (go) {
#     isolate({
#       if (DEBUG) {
#         base::cat(file = stderr(), "forceCalc\n")
#       }
#       # list of output variable and function name
#
#       withProgress(message = "Performing heavy calculations", value = 0, {
#         n <- length(heavyCalculations)
#         for (calc in heavyCalculations) {
#           shiny::incProgress(1 / n, detail = base::paste("Creating ", calc[1]))
#           if (DEBUG) {
#             cat(file = stderr(), base::paste("forceCalc ", calc[1], "\n"))
#           }
#           assign(calc[1], eval(parse(text = base::paste0(
#             calc[2], "()"
#           ))))
#         }
#       })
#     })
#
#     printTimeEnd(start.time, "forceCalc")
#   }
# })

scranWarning <- function() {
  cat(file = stderr(), paste0("scranWarning\n"))
  modalDialog(
    span(
      "The parameters clusterSource=normData and/or clusterMethod=hclust can ",
      "can result in very long wait times (>6hrs). Do you really want to do this?"
    ),
    footer = tagList(
      actionButton("scranWarning_cancel", "Cancel"),
      actionButton("scranWarning_ok", "OK")
    )
  )
}





# handle long executions ----
observeEvent(
  label = "ob23",
  eventExpr = input$clusterMethod,
  handlerExpr = {
    deepDebug()
    if (DEBUG) cat(file = stderr(), paste0("observe: input$clusterMethod\n"))
    if (input$clusterMethod == "hclust") {
      showModal(scranWarning())
    } else {
      clusterMethodReact$clusterMethod <- "igraph"
    }
  }
)

observeEvent(
  label = "ob24",
  eventExpr = input$clusterSource,
  handlerExpr = {
    deepDebug()
    if (DEBUG) cat(file = stderr(), paste0("observe: input$clusterSource\n"))
    # if (input$clusterSource == "logcounts") {
    #   showModal(scranWarning())
    # } else {
    #   clusterMethodReact$clusterSource <- "counts"
    # }
    clusterMethodReact$clusterSource <- input$clusterSource
  }
)

observeEvent(
  label = "ob25",
  eventExpr = input$scranWarning_cancel,
  handlerExpr = {
    deepDebug()
    updateSelectInput(session, "clusterMethod",
                      selected = "igraph"
    )
    # updateSelectInput(session, "clusterSource",
    #                   selected = "counts"
    # )
    removeModal()
  }
)
observeEvent(
  label = "ob26",
  eventExpr = input$scranWarning_ok,
  handlerExpr = {
    deepDebug()
    if (input$clusterMethod == "hclust") {
      clusterMethodReact$clusterMethod <- "hclust"
    }
    # if (input$clusterSource == "normData") {
    #   clusterMethodReact$clusterSource <- "normData"
    # }
    removeModal()
  }
)



observe(label = "ob_pca",
        {
          deepDebug()
          if (DEBUG) cat(file = stderr(), "observe ob_pca\n")
          # out <- pcaReact()
          # if (is.null(out)) {
          #   .schnappsEnv$calculated_gQC_tsneDim <- "NA"
          # }
          input$updatePCAParameters
          
          setRedGreenButtonCurrent(
            vars = list(
              c("pcaRank", input$pcaRank),
              c("pcaN", input$pcaN),
              c("pcaCenter", input$pcaCenter),
              c("pcaScale", input$pcaScale),
              c("hvgSelection", input$hvgSelection),
              c("genes4PCA", input$genes4PCA)
            )
          )
          
          updateButtonColor(
            buttonName = "updatePCAParameters", 
            parameters = c(
              "pcaRank", "pcaN",
              "pcaCenter", "pcaScale", "genes4PCA"
            )
          )
          
        }
)

ob_clusterParams <- observe(label = "ob_clusterParams", {
  if (DEBUG) cat(file = stderr(), "observe ob_clusterParams\n")
  
  input$updateClusteringParameters
  tabsetCluster = input$tabsetCluster
  
  # this happens when the lite version is used
  if (is.null(tabsetCluster)){
    ob_clusterParams$destroy()
    return(NULL)
  }
  
  
  if (tabsetCluster == "seurat_Clustering") {
    setRedGreenButtonCurrent(
      vars = list(
        c("tabsetCluster", input$tabsetCluster),
        c("seurClustDims", input$seurClustDims),
        c("seurClustk.param", input$seurClustk.param),
        c("seurClustresolution", input$seurClustresolution)
      )
    )
    updateButtonColor(buttonName = "updateClusteringParameters", parameters = c(
      "seurClustDims", "seurClustk.param",
      "seurClustresolution", "tabsetCluster"
    ))
  }
  if (tabsetCluster == "scran_Cluster") {
    setRedGreenButtonCurrent(
      vars = list(
        c("useRanks", input$useRanks),
        c("clusterSource", clusterMethodReact$clusterSource),
        c("geneSelectionClustering", input$geneSelectionClustering),
        c("minClusterSize", input$minClusterSize),
        c("clusterMethod", input$clusterMethod),
        c("tabsetCluster", input$tabsetCluster)
      )
    )
    updateButtonColor(buttonName = "updateClusteringParameters", parameters = c(
      "useRanks", "clusterSource","geneSelectionClustering",
      "minClusterSize", "clusterMethod", "tabsetCluster"
    ))
  }
})

# about modal ----
observeEvent(input$AboutApp,{
  deepDebug()
  showModal(modalDialog(
    title = "About SCHNAPPs",
    tags$a(tags$b("Here comes the about text")),
    easyClose = TRUE,
    footer = NULL
  ))
})


inputHelpIJS <- tryCatch(read.delim(system.file("extdata", "inputHelpIJS.txt",package = "SCHNAPPs"), sep=";", stringsAsFactors = FALSE),
                         error = function(e) {
                           cat(file = stderr(), "There is an installation problem: inputHelpIJS.txt not in extdata of package SCHNAPPs.\n")
                           stop(e)
                         })
# inputHelpIJS<- read.delim("inst/extdata/inputHelpIJS.txt", sep=";", stringsAsFactors = FALSE)

observeEvent(input$inputHelp, {
  deepDebug()
  cat(file = stderr(), paste("inputHelp started\n"))
  cat(file = stderr(), apply(inputHelpIJS, 1, FUN = function(x) if(length(x)>0)cat(file = stderr(), paste(x, "\n"))))
  introjs(session,
          options = list(steps = inputHelpIJS)
  )
})

twoDselectedAddOptHelpIJS <- read.delim(system.file("extdata", "twoDselectedAddOptHelpIJS.txt",package = "SCHNAPPs"), sep=";", stringsAsFactors = FALSE)
# twoDselectedAddOptHelpIJS <- read.delim("inst/extdata/twoDselectedAddOptHelpIJS.txt", sep=";", stringsAsFactors = FALSE)
observeEvent(input$twoDselectedAddOpt, {
  cat(file = stderr(), paste("twoDselectedAddOpt started\n"))
  # cat(file = stderr(), apply(twoDselectedAddOptHelpIJS, 1, FUN = function(x) if(length(x)>0)cat(file = stderr(), paste(x, "\n"))))
  introjs(session,
          options = list(steps = twoDselectedAddOptHelpIJS,
                         "showBullets" = "false",
                         "showProgress" = "true",
                         "showStepNumbers" = "false",
                         "nextLabel" = "Next",
                         "prevLabel" = "Prev",
                         "skipLabel" = "Skip",
                         "highlightClass" = 'berndTest')
  )
})

# Heatmap for scran clustering ----
# All clusters heatmap ------
callModule(
  pHeatMapModule,
  "clusterBootstrap",
  clusterBootstrapReactive
)

source(paste0(packagePath, "/shortCuts.R"), local = TRUE)

if (DEBUG) {
  cat(file = stderr(), paste("end: outputs.R\n"))
}
