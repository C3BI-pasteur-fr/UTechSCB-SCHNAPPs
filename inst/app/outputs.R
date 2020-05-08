suppressMessages(require(shinyTree))

# SUMMARY STATS ----------------------------------------------------------------
source(paste0(packagePath, "/moduleServer.R"), local = TRUE)

DEBUGSAVE <- get(".SCHNAPPs_DEBUGSAVE", envir = .schnappsEnv)
# normalizationRadioButtonValue --------------------------------
# Parameters / normalization
output$normalizationRadioButtonValue <- renderPrint({
  input$normalizationRadioButton
})

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
    save(
      file = "~/SCHNAPPsDebug/normalizationsParameters.RData",
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

# normalizationsParametersDynamic -------------------------
output$normalizationsParametersDynamic <- renderUI({
  if (is.null(input$normalizationRadioButton)) {
    return(NULL)
  }
  selectedChoice <- input$normalizationRadioButton
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(
      file = "~/SCHNAPPsDebug/normalizationsParametersDynamic.RData",
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
    minExpGenes = defaultValueRegExGene,
    minGenes = 20,
    maxGenes = 1000000,
    cellPatternRM = "",
    cellKeep = "",
    cellKeepOnly = "",
    cellsFiltersOut = "",
    minNonExpGenes = ""
  )
)
geneSelectionValues <- reactiveVal(
  list(
    selectIds = "^MT-|^RP|^MRP",
    geneListSelection = NULL,
    minGenesGS = 2,
    genesKeep = ""
  )
)

observeEvent(
  label = "ob20",
  eventExpr = input$updateCellSelectionParameters,
  handlerExpr = {
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
  if (DEBUG) cat(file = stderr(), "observe ob_clusteringParams\n")
  
  # this happens when the lite version is used
  if (is.null(input$tabsetCluster)){
    ob_clusteringParams$destroy()
    return(NULL)
  }
  
  input$updateClusteringParameters
  whichClustering = isolate(input$tabsetCluster)
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
    save(
      file = "~/SCHNAPPsDebug/summaryStatsSideBar.RData",
      list = c("normaliztionParameters", ls())
    )
  }
  # load("~/SCHNAPPsDebug/summaryStatsSideBar.RData")
  line0 <- paste(infile, " _ ", annFile)
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
    "Summary statistics of this dataset:",
    "<br/>",
    "<br/>",
    line0,
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
    line6,
    "<br/>",
    line7
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
#     save(
#       file = "~/SCHNAPPsDebug/selectedGenesTable.RData",
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
#     save(
#       file = "~/SCHNAPPsDebug/removedGenesTable.RData",
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
  gsSelectedGenesTable
)

callModule(
  tableSelectionServer,
  "gsRMGenesMod",
  gsRMGenesTable
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
#   # browser()
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
callModule(tableSelectionServer, "cellSelectionMod", inputSample)

# normalizationResult ----
callModule(
  tableSelectionServer,
  "normalizationResult",
  scExLogMatrixDisplay
)

# descriptionOfWork ----
output$descriptOfWorkOutput <- renderPrint({
  input$descriptionOfWork
})

# sampleColorSelection ----
output$sampleColorSelection <- renderUI({
  scEx <- scEx()
  sampCol <- sampleCols$colPal
  
  if (is.null(scEx)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(
      file = "~/SCHNAPPsDebug/sampleColorSelection.RData",
      list = c("normaliztionParameters", ls())
    )
  }
  # load("~/SCHNAPPsDebug/sampleColorSelection.RData")
  
  lev <- levels(colData(scEx)$sampleNames)
  # cols <- gg_fill_hue(length(lev))
  
  # New IDs "colX1" so that it partly coincide with input$select...
  lapply(seq_along(lev), function(i) {
    colourpicker::colourInput(
      inputId = paste0("sampleNamecol", lev[i]),
      label = paste0("Choose colour for sample ", "\"", lev[i], "\""),
      # value = "#762A83"
      # ,
      value = sampCol[i],
      allowedCols = allowedColors,
      palette = "limited"
    )
  })
})

# clusterColorSelection ----
output$clusterColorSelection <- renderUI({
  scEx <- scEx()
  projections <- projections()
  clusterCol <- clusterCols$colPal
  
  if (is.null(scEx) || is.null(projections)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(
      file = "~/SCHNAPPsDebug/clusterColorSelection.RData",
      list = c("normaliztionParameters", ls())
    )
  }
  # load("~/SCHNAPPsDebug/clusterColorSelection.RData")
  
  lev <- levels(projections$dbCluster)
  # cols <- gg_fill_hue(length(lev))
  
  # New IDs "colX1" so that it partly coincide with input$select...
  lapply(seq_along(lev), function(i) {
    colourpicker::colourInput(
      inputId = paste0("clusterNamecol", lev[i]),
      label = paste0("Choose colour for cluster ", "\"", lev[i], "\""),
      # value = "#762A83"
      # ,
      value = clusterCol[i],
      allowedCols = allowedColors,
      palette = "limited"
    )
  })
})

# history store to file ----
#' 

# askComment <- function(failed = FALSE) {
#   modalDialog(
#     textInput("HistComment", "add a comment", value = paste("created at ",date())),
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
  if (DEBUG) cat(file = stderr(), "observe color Vars\n")
  
  scExx <- scEx()
  projections <- projections()
  if (is.null(scExx) || is.null(projections)) {
    return(NULL)
  }
  
  lev <- levels(projections$dbCluster)
  ccols <- lapply(seq_along(lev), function(i) {
    input[[paste0("clusterNamecol", lev[i])]]
  })
  lev <- levels(colData(scExx)$sampleNames)
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
    
    # TODO Warning if scEx_log is not set
    
    # umaps???
    scEx <- scEx()
    projections <- projections()
    scEx_log <- scEx_log()
    pca <- pca()
    # TODO should be taken from projections
    tsne <- tsne()
    ccol = clusterCols$colPal
    scol = sampleCols$colPal
    
    
    if (is.null(scEx)) {
      return(NULL)
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/RDSsave.RData", list = c(ls()))
    }
    # load(file='~/SCHNAPPsDebug/RDSsave.RData')
    
    scEx <- consolidateScEx(scEx, projections, scEx_log, pca, tsne)
    
    # we save the pca separately because I don't know how to store the rotation  otherwise.
    # mostly done to make the lite version work.
    
    saveList =  c("scEx" , "pca", "scol" , "ccol" )
# browser()
    # save projections that shouldn't be recalculated in lite version
    for (idx in 1:length(.schnappsEnv$projectionFunctions) ){
      assign(.schnappsEnv$projectionFunctions[[idx]][2], eval(parse(text = paste0(.schnappsEnv$projectionFunctions[[idx]][2],"()"))))
      saveList = c(saveList, .schnappsEnv$projectionFunctions[[idx]][2])
    }
    
    
    save(file = file, list = saveList)
    
    # write.csv(as.matrix(exprs(scEx)), file)
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
          if (DEBUG) cat(file = stderr(), "observe ob_pca\n")
          # out <- pca()
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


if (DEBUG) {
  cat(file = stderr(), paste("end: outputs.R\n"))
}
