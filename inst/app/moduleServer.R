if (DEBUG) {
  cat(file = stderr(), "\nloading Module server.\n")
}

# source(paste0(packagePath, "/reactives.R"), local = TRUE)
suppressMessages(library(psychTools))
suppressMessages(library(magrittr))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(InteractiveComplexHeatmap))
require(psychTools)

#' clusterServer
#'
#' server side shiny module function for printing a 2D represenation of cells
#' it uses the global projections object for plotting

#' @param gene_id name of the gene to be plotted (comma separated list, will be set to upper case)
#' @param tData projections, dataframe with cluster numbers
#' @param DEBUG whether or not to plot debugging messages on stderr
#' @param selectedCells cells that should be marked as a triangle
#' @param legend.position "none", ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#'
#' uses global reactives featueDataRact
#'                       log2cpm # for defining the color of the cells



#
# TODO parameter gene_id should be able to handle multiple gene_ids
# TODO coloring based on number of genes selected (if NULL=color by cluster NR)
#      handle coloring for binarized plot
# TODO potentially integrate gene_id selection into module?
clusterServer <- function(input, output, session,
                          tData,
                          gene_id = returnNull, # reactive
                          # selectedCells  = NULL,
                          legend.position = "none"
                          # ,
                          # defaultValues = c("tsne1", "tsne2")
) {
  ns <- session$ns
  subsetData <- NULL
  selectedGroupName <- ""
  groupName <- ""
  .schnappsEnv[[ns("addToGroupValue")]] <- FALSE
  
  
  # dim1 <- defaultValues[1]
  # dim2 <- defaultValues[2]
  # dim1 <- "PC1"
  # dim2 <- "PC2"
  # dimCol <- "Gene.count"
  divXBy <- "None"
  divYBy <- "None"
  # mod_cl1 <- ""
  # observe({
  #   if (DEBUG) cat(file = stderr(), paste0("observe: clusters\n"))
  #   mod_cl1 <<- input$clusters
  # })
  
  observe(label = "dimension_x", {
    # browser()
    if (DEBUG) cat(file = stderr(), paste0("observe: dimension_x\n"))
    .schnappsEnv[[ns('dimension_x')]] <- input$dimension_x
    .schnappsEnv$defaultValues[[ns('dimension_x')]] <- input$dimension_x
    
    .schnappsEnv[[ns('geneIds')]] <- input$geneIds
    .schnappsEnv$defaultValues[[ns('geneIds')]] <- input$geneIds
    .schnappsEnv[[ns('geneIds2')]] <- input$geneIds2
    .schnappsEnv$defaultValues[[ns('geneIds2')]] <- input$geneIds2
  })
  observe(label = "ob32", {
    if (DEBUG) cat(file = stderr(), paste0("observe: dimension_y\n"))
    .schnappsEnv[[ns('dimension_y')]]<- input$dimension_y
    .schnappsEnv$defaultValues[[ns('dimension_y')]]<- input$dimension_y
  })
  observe(label = "ob33", {
    if (DEBUG) cat(file = stderr(), paste0("observe: dimension_col\n"))
    .schnappsEnv[[ns('dimension_col')]] <- input$dimension_col
    .schnappsEnv$defaultValues[[ns('dimension_col')]] <- input$dimension_col
  })
  observe(label = "ob34", {
    if (DEBUG) cat(file = stderr(), paste0("observe: divideXBy\n"))
    .schnappsEnv[[ns("divideXBy")]] <- input$divideXBy
    .schnappsEnv$defaultValues[[ns("divideXBy")]] <- input$divideXBy
  })
  observe(label = "ob35", {
    if (DEBUG) cat(file = stderr(), paste0("observe: divideYBy\n"))
    .schnappsEnv[[ns("divideYBy")]] <- input$divideYBy
    .schnappsEnv$defaultValues[[ns("divideYBy")]] <- input$divideYBy
  })
  observe(label = "ob36", {
    if (DEBUG) cat(file = stderr(), paste0("observe: logX\n"))
    .schnappsEnv[[ns("logX")]] <- input$logX
    .schnappsEnv$defaultValues[[ns("logX")]] <- input$logX
  })
  observe(label = "ob37", {
    if (DEBUG) cat(file = stderr(), paste0("observe: logY\n"))
    .schnappsEnv[[ns("logY")]] <- input$logY
    .schnappsEnv$defaultValues[[ns("logY")]] <- input$logY
  })
  observe(label = "ob38", {
    if (DEBUG) cat(file = stderr(), paste0("observe: showCells\n"))
    .schnappsEnv[[ns("showCells")]] <- input$showCells
    .schnappsEnv$defaultValues[[ns("showCells")]] <- input$showCells
  })
  # clusterServer - observe input$addToGroup ----
  observe(label = "ob39", {
    # deepDebug()
    if (DEBUG) cat(file = stderr(), "observe input$addToGroup \n")
    if (!is.null(input$addToGroup)) {
      .schnappsEnv[[ns("addToGroup")]] <- input$addToGroup
      .schnappsEnv$defaultValues[[ns("addToGroup")]] <- input$addToGroup
    }
  })
  
  # clusterServer - observe input$groupNames ----
  # if we manually select a group from the list we update the groupName field
  # as this is the value that is actually being used for other calculations (except "plot" is selected)
  observe(label = "ob40", {
    if (DEBUG) cat(file = stderr(), "observe input$groupNames \n")
    if (!is.null(input$groupNames)) {
      if (input$groupNames == "plot") {
        
      } else {
        isolate({
          updateTextInput(
            session = session, 
            inputId = "groupName",
            value = input$groupNames
          )
        })
      }
    }
  })
  
  # observe save  2 history ----
  observe(label = "save2histMod2d", {
    # deepDebug()
    clicked <- input$save2Hist
    if (DEBUG) cat(file = stderr(), "observe input$save2Hist \n")
    myns <- session$ns("-")
    req(.schnappsEnv[[ns("historyPlot")]])
    start.time <- base::Sys.time()
    if (DEBUG) cat(file = stderr(), "cluster: save2Hist\n")
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "save2Hist")
      }
    )
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("save2Hist", id = "save2Hist", duration = NULL)
    }
    
    if (is.null(clicked)) {
      return()
    }
    if (clicked < 1) {
      return()
    }
    add2history(
      type = "save", input = c(isolate( reactiveValuesToList(input)), 
                               isolate( reactiveValuesToList(
                                 get("input", envir = parent.env(parent.env(environment())))))),
      comment = paste("# ",myns, "\n",
                      "fun = plotData$plotData$plotFunc\n", 
                      "environment(fun) = environment()\n",
                      "do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)])\n"
      ),
      plotData = .schnappsEnv[[ns("historyPlot")]]
    )
  })
  # clusterServer - updateInput ----
  # updateInput <-
  observe(label = "ob41", {
    if (DEBUG) cat(file = stderr(), paste0("updateInput\n"))
    projections <- projections()
    
    # Can use character(0) to remove all choices
    if (is.null(projections)) {
      return(NULL)
    }
    
    # Can also set the label and select items
    # if (is.null(mod_cl1) || mod_cl1 == "") mod_cl1 = levels(projections$dbCluster)
    # updateSelectInput(session, "clusters",
    #                   choices = levels(projections$dbCluster)
    #                   ,
    #                   selected = mod_cl1
    # )
    updateSelectInput(session, "dimension_x",
                      choices = c(colnames(projections), "UmiCountPerGenes", "UmiCountPerGenes2"),
                      selected = .schnappsEnv[[ns('dimension_x')]]
    )
    updateSelectInput(session, "dimension_y",
                      choices = c(colnames(projections), "histogram", "UmiCountPerGenes", "UmiCountPerGenes2"),
                      selected = .schnappsEnv[[ns('dimension_y')]]
    )
    updateSelectInput(session, "dimension_col",
                      choices = c(colnames(projections), "cellDensity", "UmiCountPerGenes", "UmiCountPerGenes2"),
                      selected = .schnappsEnv[[ns('dimension_col')]]
    )
    
    updateSelectInput(session, "divideXBy",
                      choices = c("None", colnames(projections), "UmiCountPerGenes", "UmiCountPerGenes2"),
                      selected = .schnappsEnv[[ns("divXBy")]]
    )
    updateSelectInput(session, "divideYBy",
                      choices = c("None", colnames(projections), "UmiCountPerGenes", "UmiCountPerGenes2", "normByCol"),
                      selected = .schnappsEnv[[ns("divYBy")]]
    )
    updateCheckboxInput(session, "logX",
                        value = .schnappsEnv[[ns("logX")]]
    )
    updateCheckboxInput(session, "logY",
                        value = .schnappsEnv[[ns("logY")]]
    )
    updateCheckboxInput(session, "showCells",
                        value = .schnappsEnv[[ns("showCells")]]
    )
  })
  
  grpNameDebounced <- reactive({
    make.names(input$groupName, unique = TRUE)
  }) %>% debounce(1000)
  
  # clusterServer - selectedCellNames ----
  selectedCellNames <- reactive({
    start.time <- base::Sys.time()
    if (DEBUG)
      cat(file = stderr(), "cluster: selectedCellNames\n")
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "selectedCellNames")
      }
    )
    
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("selectedCellNames", id = "selectedCellNames", duration = NULL)
    }
    projections <- projections()
    req(projections)
    # deepDebug()
    brushedPs <- plotly::event_data("plotly_selected", source = "subset")
    dimY <- input$dimension_y
    dimX <- input$dimension_x
    geneNames <- input$geneIds
    geneNames2 <- input$geneIds2
    scEx_log <- scEx_log()
    scEx <- scEx()
    namedGroup <- input$groupNames
    grpN <- grpNameDebounced()
    grpNs <- groupNames$namesDF
    
    if (is.null(projections) | is.null(brushedPs)) {
      if (DEBUG) cat(file = stderr(), "cluster: selectedCellNames: brush null\n")
      return(NULL)
    }
    # inpClusters <- input$clusters
    inpClusters <- levels(projections$dbCluster)
    
    if (.schnappsEnv$DEBUGSAVE) {
      if (DEBUG) cat(file = stderr(), "cluster: selectedCellNames: saving\n")
      save(file = "~/SCHNAPPsDebug/selectedCellNames.RData", list = c(ls(), "legend.position"))
    }
    # cp = load(file="~/SCHNAPPsDebug/selectedCellNames.RData")
    # deepDebug()
    
    # in case no normalization is done
    if (is.null(scEx_log)) {
      scEx_log <- scEx
    }
    featureData <- rowData(scEx)
    geneid <- geneName2Index(geneNames, featureData)
    projections <- updateProjectionsWithUmiCount(
      dimX = dimX, dimY = dimY,
      geneNames = geneNames,
      geneNames2 = geneNames2,
      scEx = scEx_log[, rownames(projections)], projections = projections
    )
    
    cells.names <- unlist(brushedPs$key)
    if (!is.null(namedGroup)) {
      if (namedGroup == "none") {
        if (DEBUG) cat(file = stderr(), "cluster: selectedCellNames: namedGroup none\n")
        return(NULL)
      }
      if (!namedGroup == "plot") {
        if (namedGroup %in% colnames(grpNs)) {
          grpNs = grpNs[unique(cells.names),]
          return(rownames(grpNs[grpNs[, namedGroup] == "TRUE", ,drop=FALSE]))
        } else {
          # TODO message about setting plot
          showNotification("Make sure 'the'group names' is set correctly ", id = "selectedCellNamesProbl", duration = 10)
          if (DEBUG) cat(file = stderr(), "cluster: selectedCellNames:some return null\n")
          return(NULL)
        }
      }
    }
    # deepDebug()
    # if(!is.null(inpClusters)){
    #   subsetData <- subset(projections, dbCluster %in% inpClusters)
    # }else{
    #   subsetData = projections
    # }
    # cells.names <- rownames(projections)[subset(brushedPs, curveNumber == 0)$pointNumber + 1]
    # cells.names <- rownames(projections)[subset(brushedPs)$pointNumber + 1]
    if (is.null(cells.names) )
      if (length(brushedPs) >0)
        if (nrow(brushedPs) > 0) { 
          cells.names = tryCatch(
            {
              rownames(projections[which(projections[,dimX] %in% brushedPs$x & 
                                           projections[,dimY] %in% brushedPs$y),])
            },
            error = function(e) {
              cat(file = stderr(), "\n\ncaught exception with dimx:", dimX, " dimy: ", dimY,  "\n\n")
              return(NULL)
            })
        }
    if(is.null(cells.names)) return(NULL)
    cells.names <- cells.names[cells.names %in% colnames(scEx_log)]
    cells.names <- unique(cells.names[!is.na(cells.names)])
    # if (DEBUG) {
    #   cat(file = stderr(), paste("curveNumbers:", unique(brushedPs$curveNumber), "\n"))
    # }
    printTimeEnd(start.time, "selectedCellNames")
    exportTestValues(selectedCellNames = {
      cells.names
    })
    if (DEBUG) cat(file = stderr(), "cluster: selectedCellNames: ", length(cells.names), "\n")
    
    return(cells.names)
  })
  
  # clusterServer - returnValues ----
  returnValues <- reactiveValues(
    # cluster = reactive(input$clusters),
    selectedCells = reactive({
      if (DEBUG) cat(file = stderr(), "reactiveValues: selectedCells.\n")
      start.time <- Sys.time()
      
      # moreOptions <- input$moreOptions
      retVal <- selectedCellNames()
      if (length(retVal) == 0) {
        if (DEBUG) cat(file = stderr(), paste("selectedCellNames is null\n"))
        retVal <- NULL
      }
      grpN <-grpNameDebounced()
      grpSelected <- make.names(input$groupNames, unique = TRUE)
      grpNs <- groupNames$namesDF
      if (length(grpN) == 0 | length(grpNs) == 0) {
        if (DEBUG) cat(file = stderr(), "reactiveValues: grpN empty\n")
        return(retVal)
      }
      # inpClusters <- input$clusters
      projections <- projections()
      dimY <- input$dimension_y
      dimX <- input$dimension_x
      geneNames <- input$geneIds
      geneNames2 <- input$geneIds2
      scEx_log <- scEx_log()
      scEx <- scEx()
      
      if (.schnappsEnv$DEBUGSAVE) {
        cat(file = stderr(), paste("selectedCell: saving\n"))
        base::save(file = "~/SCHNAPPsDebug/clusterServerreturnValues.RData", list = c(ls()))
      }
      # load(file="~/SCHNAPPsDebug/clusterServerreturnValues.RData")
      # in case no normalization is done:
      if (is.null(scEx_log)) {
        scEx_log <- scEx
      }
      featureData <- rowData(scEx)
      inpClusters <- levels(projections$dbCluster)
      # if (!is.null(projections) & moreOptions & !grpSelected == "plot") {
      if (!is.null(projections) & !grpSelected == "plot") {
        projections <- updateProjectionsWithUmiCount(
          dimX = dimX, dimY = dimY,
          geneNames = geneNames,
          geneNames2 = geneNames2,
          scEx = scEx_log[, rownames(projections)], projections = projections
        )
        if (!is.null(inpClusters)) {
          subsetData <- subset(projections, dbCluster %in% inpClusters)
        }else{
          subsetData = projections
        }
        grpSubset <- grpNs[rownames(subsetData), ]
        if (!grpN %in% colnames(grpSubset)) {
          if (!is.null(getDefaultReactiveDomain())) {
            # showing too often 
            # TODO check why and where
            # showNotification("group name is not available", id = "nogrpN", duration = NULL, type = "error")
          }
          return(NULL)
        }
        grpVal <- rownames(grpSubset[grpSubset[, grpN] == "TRUE", ])
        if (length(grpVal) > 0) {
          return(grpVal)
        }
      }
      
      retVal <- retVal[retVal %in% colnames(scEx)]
      printTimeEnd(start.time, "clusterServerReturnVal")
      exportTestValues(clusterServerReturnVal = {
        retVal
      })
      if (length(retVal) == 0) {
        return(NULL)
      }
      if (DEBUG) cat(file = stderr(), paste("reactiveValues: selectedCells:", length(retVal), " done\n"))
      
      return(retVal)
    })
  )
  
  # if (DEBUG) {
  #   cat(file = stderr(), paste("clusterServers", session$ns("clusters"), "\n"))
  # }
  
  # clusterServer - output$clusters ----
  # output$clusters <- renderUI({
  #   # if (DEBUG) cat(file = stderr(), paste("2observe: ns(input$clusters)", session$ns(input$clusters), "\n"))
  #   retVal <- NULL
  #   projections <- tData()
  #   upI <- updateInput() # needed to update input of this module
  #   ns <- session$ns
  #    if (DEBUG) cat(file = stderr(), paste("2observe: ns(mod_cl1)", ns(mod_cl1), "\n"))
  #   if (is.null(projections)) {
  #     HTML("Please load data first")
  #   } else {
  #     noOfClusters <- levels(as.factor(projections$dbCluster))
  #     # noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
  #     retVal <- selectizeInput(
  #       ns("clusters"),
  #       label = "Cluster",
  #       choices = noOfClusters,
  #       # selected = input$clusters, # not working because of stack, too slow and possible to create infinite loop
  #       selected = mod_cl1,
  #       multiple = TRUE
  #     )
  #   }
  #   retVal
  # })
  
  # clusterServer - output$clusterPlot ----
  output$clusterPlot <- plotly::renderPlotly({
    if (DEBUG) cat(file = stderr(), paste("Module: output$clusterPlot\n"))
    start.time <- base::Sys.time()
    
    # remove any notification on exit that we don't want
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "clusterPlot")
      }
    )
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("clusterPlot", id = "clusterPlot", duration = NULL)
    }
    
    scEx_log <- scEx_log()
    scEx <- scEx()
    tdata <- tData()
    projections <- projections()
    grpNs <- groupNames$namesDF
    grpN <- grpNameDebounced()
    pc = projectionColors %>% reactiveValuesToList()
    
    # returnValues$cluster <- input$clusters
    dimY <- input$dimension_y
    dimX <- input$dimension_x
    dimCol <- input$dimension_col
    g_id <- gene_id()
    geneNames <- input$geneIds
    geneNames2 <- input$geneIds2
    logx <- input$logX
    logy <- input$logY
    divXBy <- input$divideXBy
    divYBy <- input$divideYBy
    # scols <- projectionColors$sampleNames
    # ccols <- projectionColors$dbCluster
    # browser()
    if(!"sampleNames" %in% names(projectionColors)) return(NULL) # should not happen
    # scols <- projectionColors[["sampleNames"]]
    # ccols <- projectionColors[["dbCluster"]]
    # moreOptions <- input$moreOptions
    myns <- session$ns("-")
    save2History <- .schnappsEnv$saveHistorycheckbox
    if (is.null(save2History)) {
      save2History <- FALSE
    }
    if (is.null(scEx) | is.null(tdata)) {
      if (DEBUG) cat(file = stderr(), paste("output$clusterPlot:NULL\n"))
      .schnappsEnv[[ns("historyPlot")]] <- NULL
      return(NULL)
    }
    # in case the normalization is not done
    if (is.null(scEx_log)) {
      scEx_log <- scEx
    }
    # clId <- input$clusters
    clId <- levels(projections$dbCluster)
    
    featureData <- rowData(scEx_log)
    if (.schnappsEnv$DEBUGSAVE) {
      cat(file = stderr(), paste("cluster plot saving\n"))
      save(file = paste0("~/SCHNAPPsDebug/clusterPlot-", ns("-"), ".RData", collapse = "."),
           list = c(ls(), "legend.position")
      )
      cat(file = stderr(), paste("cluster plot saving done\n"))
    }
    
    # cp = load("~/SCHNAPPsDebug/clusterPlot-DE_expclusters--.RData");.schnappsEnv$DEBUGSAVE=FALSE
    if (is.null(g_id) || nchar(g_id) == 0) {
      g_id <- featureData$symbol
    }
    if (is.null(logx)) logx <- FALSE
    if (is.null(logy)) logy <- FALSE
    if (is.null(divXBy)) divXBy <- "None"
    if (is.null(divYBy)) divYBy <- "None"
    # TODO returns a table but is not used
    # updateProjectionsWithUmiCount(
    #   dimX = dimX, dimY = dimY,
    #   geneNames = geneNames,
    #   geneNames2 = geneNames2,
    #   scEx = scEx_log[, rownames(projections)], projections = tdata
    # )
    if(dimCol %in% names(pc)){
      myColors = pc[[dimCol]]
    }else {
      myColors <- NULL
    }
    
    # if (dimCol == "sampleNames") {
    #   myColors <- scols
    # } 
    # if (dimCol == "dbCluster") {
    #   myColors <- ccols
    # }
    # if (!moreOptions) grpN <- ""
    p1 <- plot2Dprojection(scEx_log,
                           projections = tdata, g_id, featureData, geneNames,
                           geneNames2, dimX, dimY, clId, grpN, legend.position,
                           grpNs = grpNs, logx, logy, divXBy, divYBy, dimCol, colors = myColors
    )
    
    # save p1 to .schnappsEnv for saving to history
    af = plot2Dprojection
    # remove env because it is too big
    environment(af) = new.env(parent = emptyenv())
    .schnappsEnv[[ns("historyPlot")]] <- list(plotFunc = af,
                                              scEx_log = scEx_log,
                                              projections = tdata,
                                              g_id = g_id, 
                                              featureData = featureData, 
                                              geneNames = geneNames, 
                                              geneNames2 = geneNames2, 
                                              dimX = dimX, 
                                              dimY = dimY, 
                                              clId = clId, 
                                              grpN = grpN, 
                                              legend.positio = legend.position,
                                              grpNs = grpNs, 
                                              logx = logx, 
                                              logy = logy, 
                                              divXBy = divXBy, 
                                              divYBy = divYBy, 
                                              dimCol = dimCol, 
                                              colors = myColors
    )
    
    # .schnappsEnv[[paste0("historyPlot-", myns)]] <- p1
    
    # add2history(type = "renderPlotly", plotData = p1, comment = paste(myns))
    # if (save2History) recHistory(myns, p1)
    # event_register(p1, 'plotly_selected')
    printTimeEnd(start.time, "clusterPlot")
    # deepDebug()
    # .schnappsEnv[[paste0()]] <- p
    exportTestValues(clusterPlot = {
      p1
    })
    suppressMessages(p1)
  })
  
  # observe({
  #   updateTextInput(session = session, inputId = "groupName",
  #                   value = input$groupNames)
  # })
  
  
  # clusterServer - visibleCellNames ----
  visibleCellNames <- reactive({
    if (DEBUG) cat(file = stderr(), "cluster: selectedCellNames\n")
    start.time <- base::Sys.time()
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "visibleCellNames")
      }
    )
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("visibleCellNames", id = "visibleCellNames", duration = NULL)
    }
    
    projections <- projections()
    if (is.null(projections)) {
      return(NULL)
    }
    # inpClusters <- input$clusters
    inpClusters <- levels(projections$dbCluster)
    if(!is.null(inpClusters)){
      subsetData <- subset(projections, dbCluster %in% inpClusters)
    }else{
      subsetData = projections
    }
    printTimeEnd(start.time, "visibleCellNames")
    exportTestValues(visibleCellNames = {
      subsetData
    })
    return(subsetData)
  })
  
  
  # clusterServer - observe input$changeGroups ----
  observe(label = "ob42", {
    # deepDebug()
    if (DEBUG) cat(file = stderr(), "observe input$changeGroups \n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "observe input$changeGroups")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "input-changeGroups")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("observe: changeGroups", id = "input-changeGroups", duration = NULL)
    }
    
    ns <- session$ns
    clicked <- input$changeGroups # action button
    if (clicked < 1) {
      if (DEBUG) cat(file = stderr(), "     input$changeGroups not clicked\n")
      return(NULL)
    }
    addToSelection <- isolate(input$addToGroup)
    inpGroupName = isolate(input$groupNames)
    # we isolate here because we only want to change if the button is clicked.
    # TODO what happens if new file is loaded??? => problem!
    isolate({
      # we should react on a changed filename, but that would imply calculating pca's etc directly after loading
      scEx <- scEx()
      acn = allCellNames()
      prjs <- sessionProjections$prjs
      # brushedPs <- plotly::event_data("plotly_selected", source = "subset")
      # scEx <- scEx()
      # inpClusters <- input$clusters
      
      # this used to be make.names(input$groupName) which created a column called "X"
      grpN <- input$groupName
      grpNs <- groupNames$namesDF
      cells.names <- selectedCellNames()
      visibleCells <- visibleCellNames()
      # 
      if (nrow(grpNs) == 0) {
        initializeGroupNames()
      }
      if (ncol(grpNs) == 0) {
        grpNs <- groupNames$namesDF
      }
      if (is.null(scEx)) {
        return(NULL)
      }
    })
    if (length(grpN) == 0 || nchar(grpN) == 0) {
      return(NULL)
    }
    # browser()
    if (.schnappsEnv$DEBUGSAVE) {
      cat(file = stderr(), "save: changeGroups\n")
      save(file = "~/SCHNAPPsDebug/changeGroups.RData", list = c(ls()))
      cat(file = stderr(), "done save: changeGroups\n")
      # deepDebug()
    }
    grpN =  make.names(grpN, unique = TRUE)
    # cp = load(file="~/SCHNAPPsDebug/changeGroups.RData")
    # in case the cell selection has changed
    grpNs <- grpNs[colnames(scEx), ]
    if (!grpN %in% colnames(grpNs)) {
      grpNs[, grpN] <- "FALSE"
    }
    grpNs[, grpN] = as.logical(grpNs[, grpN])
    if (!addToSelection) {
      grpNs[rownames(visibleCells), grpN] <- FALSE
    }
    if (class(cells.names) == "list") {
      cells.names = unlist(cells.names)
    }
    grpNs[cells.names, grpN] <- TRUE
    grpNs[, grpN] = as.factor(as.character(grpNs[, grpN]))
    # grpNs[, grpN] <- as.factor(grpNs[, grpN])
    # Set  reactive value
    # cat(file = stderr(), paste("DEBUG: ",colnames(grpNs)," \n"))
    # deepDebug()
    groupNames$namesDF <- grpNs
    updateSelectInput(session, "groupNames",
                      choices = c("plot", colnames(grpNs)),
                      # selected = grpN
                      selected = inpGroupName
    )
    updateTextInput(
      session = session, inputId = "groupName",
      value = grpN
    )
    tryCatch({
      if (ncol(prjs) > 0) {
        # We overwrite the common columns
        comColPrjs = which(colnames(prjs) %in% colnames(grpNs) )
        # didn't find a way to easily overwrite columns
        # we need to add 1 to comColPrjs because tibble added rownames column 
        prjs[rownames(grpNs), comColPrjs] = grpNs[, colnames(prjs)[comColPrjs]]
        
        newColsGrpns = which(!colnames(grpNs) %in% colnames(prjs))
        if(length(newColsGrpns) > 0) {
          # prj has all the rows from the very first input data set.
          # prjs <- tibble::rownames_to_column(prjs)
          prjs <- dplyr::left_join(
            tibble::rownames_to_column(prjs), 
            tibble::rownames_to_column(grpNs[,newColsGrpns, drop=FALSE]), 
            by='rowname')
          rownames(prjs) = prjs[,1]
          prjs = prjs[,-1]
        }
        
        # for (cn in colnames(grpNs)) {
        #   if (cn %in% colnames(prjs)) {
        #     prjs[, cn] <- grpNs[, cn]
        #   } else {
        #     prjs <- base::cbind(prjs, grpNs[, cn], deparse.level = 0)
        #     colnames(prjs)[ncol(prjs)] <- cn
        #   }
        # }
        
        # sessionProjections$prjs <- prjs
      } else {
        prjs = data.frame(row.names = acn)
        prjs <- dplyr::left_join(
          tibble::rownames_to_column(prjs), 
          tibble::rownames_to_column(grpNs), 
          by='rowname')
        rownames(prjs) = prjs[,1]
        prjs = prjs[,-1]
      }
      prjs[is.na(prjs)] <- "FALSE"
      prjs$all = "TRUE"
      if ('rowname' %in% colnames(prjs)) prjs = prjs [,-which(colnames(prjs)=='rowname')]
    },
    error = function(e){
      # deepDebug()
      cat(file = stderr(), paste("error in grp", e))
    })
    sessionProjections$prjs = prjs
    selectedGroupName <<- grpN
  })
  
  # clusterServer - output$nCellsVisibleSelected ----
  # display the number of cells that belong to the group, but only from the visible ones
  output$nCellsVisibleSelected <- renderText({
    if (DEBUG) cat(file = stderr(), "nCellsVisibleSelected.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "nCellsVisibleSelected")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "nCellsVisibleSelected")
      }
    })
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("nCellsVisibleSelected", id = "nCellsVisibleSelected", duration = NULL)
    }
    
    grpN <- grpNameDebounced()
    grpNs <- groupNames$namesDF
    cells.names <- selectedCellNames()
    visibleCells <- visibleCellNames()
    
    # inpClusters <- input$clusters
    projections <- projections()
    if (is.null(projections)) {
      return(NULL)
    }
    # deepDebug()
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/nCellsVisibleSelected.RData", list = c(ls()))
    }
    # cp = load(file="~/SCHNAPPsDebug/nCellsVisibleSelected.RData")
    inpClusters <- levels(projections$dbCluster)
    if(!is.null(inpClusters)){
      subsetData <- subset(projections, dbCluster %in% inpClusters)
    }else{
      subsetData = projections
    }
    # deepDebug()
    retVal <- paste("Number of visible cells in section ", 
                    sum(grpNs[rownames(subsetData), grpN] == "TRUE"), 
                    "\n",
                    "selected cells ", length(cells.names), "\n",
                    "selected & named: ",sum(grpNs[cells.names, grpN] == "TRUE"), "\n",
                    "visible cells ", nrow(visibleCells))
    if (DEBUG) cat(file = stderr(), retVal)
    
    exportTestValues(DummyReactive = {
      retVal
    })
    return(retVal)
  })
  
  # clusterServer - output$nCellsAllSelected ----
  # display the number of cells that belong to the group, including the cells from non visible clusters
  # output$nCellsAllSelected <- renderText({
  #   if (DEBUG) cat(file = stderr(), "nCellsAllSelected started.\n")
  #   start.time <- base::Sys.time()
  #   on.exit({
  #     printTimeEnd(start.time, "nCellsAllSelected")
  #     if (!is.null(getDefaultReactiveDomain())) {
  #       removeNotification(id = "nCellsAllSelected")
  #     }
  #   })
  #   # show in the app that this is running
  #   if (!is.null(getDefaultReactiveDomain())) {
  #     showNotification("nCellsAllSelected", id = "nCellsAllSelected", duration = NULL)
  #   }
  #
  #   grpNs <- groupNames$namesDF
  #   grpN <- make.names(input$groupName)
  #   val <- 0
  #   if (grpN %in% colnames(grpNs)) {
  #     val <- sum(grpNs[, grpN])
  #   }
  #   retVal <- paste("number of cells in group over all cells", val)
  #
  #   exportTestValues(DummyReactive = {
  #     retVal
  #   })
  #   return(retVal)
  # })
  
  
  # clusterServer - output$additionalOptions ----
  # TODO: not sure if we still need this after changing to box layout. Maybe it should be an observer?
  output$additionalOptions <- renderUI({
    if (DEBUG) cat(file = stderr(), "additionalOptions started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "additionalOptions")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "additionalOptions")
      }
    })
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("additionalOptions", id = "additionalOptions", duration = NULL)
    }
    
    ns <- session$ns
    # moreOptions <- input$moreOptions
    # groupNs <- groupNames$namesDF
    grpNs <- groupNames$namesDF
    projections <- projections()
    # if (!moreOptions) {
    #   if (DEBUG) cat(file = stderr(), "additionalOptions NULL\n")
    #   groupName <<- ""
    #   # this doesn't seem to work... cannot set to NULL after it has been initialized
    #   updateSelectInput(session, "groupNames",
    #     choices = c("plot", colnames(grpNs)),
    #     selected = "plot"
    #   )
    #   updateTextInput(
    #     session = session, inputId = "groupName",
    #     value = "none"
    #   )
    #   # updateTextInput(
    #   #   session = session, inputId = "groupName",
    #   #   value = "none"
    #   # )
    #
    #   return("")
    # }
    
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/additionalOptions.RData", list = c(ls()))
    }
    # load(file="~/SCHNAPPsDebug/additionalOptions.RData")
    
    
    return("")
  })
  
  # clusterServer - output$cellSelection ----
  output$cellSelection <- renderText({
    if (DEBUG)
      cat(file = stderr(), "cellSelection started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "cellSelection")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "cellSelection")
      }
    })
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("cellSelection", id = "cellSelection", duration = NULL)
    }
    
    ns <- session$ns
    projections <- projections()
    req(projections)
    brushedPs <- suppressMessages(plotly::event_data("plotly_selected", source = "subset"))
    # inpClusters <- (input$clusters)
    myshowCells <- (input$showCells)
    geneNames <- input$geneIds
    geneNames2 <- input$geneIds2
    dimY <- input$dimension_y
    dimX <- input$dimension_x
    scEx_log <- scEx_log()
    # moreOptions <- input$moreOptions
    retVal <- selectedCellNames()
    grpN <- grpNameDebounced()
    grpSelected <- make.names(input$groupNames, unique = TRUE)
    grpNs <- groupNames$namesDF
    myns <- ns("cellSelection")
    inputGroupName = isolate(input$groupName)
    # in case there is a name given for a group we only print the rownames belonging to this name
    if (inputGroupName != ""){
      if(grpN %in% names(grpNs) ){
        retVal = rownames(grpNs)[as.logical(grpNs[,grpN])]
      } else {
        return("")
      }
    }
    if (!myshowCells) {
      return("")
    }
    if (is.null(projections)) {
      return("")
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/clustercellSelection.RData", list = c(ls()))
    }
    # cp = load(file=paste0("~/SCHNAPPsDebug/clustercellSelection.RData"))
    
    # inpClusters <- levels(projections$dbCluster)
    # featureData <- rowData(scEx_log)
    # subsetData <- subset(projections, dbCluster %in% inpClusters)
    # geneid <- geneName2Index(geneNames, featureData)
    # subsetData <- updateProjectionsWithUmiCount(
    #   dimX = dimX, dimY = dimY,
    #   geneNames = geneNames,
    #   geneNames2 = geneNames2,
    #   scEx = scEx_log[, rownames(projections)], projections = projections
    # )
    #
    # cells.names <- rownames(projections)[subset(brushedPs, curveNumber == 0)$pointNumber + 1]
    # cells.names <- cells.names[!is.na(cells.names)]
    retVal <- paste(retVal, collapse = ", ")
    
    exportTestValues(ClusterCellSelection = {
      retVal
    })
    return(retVal)
  })
  
  # clusterServer - return ----
  return(reactive({
    returnValues
  }))
}


tableSelectionServer <- function(input, output, session,
                                 dataTab, caption = "Table") {
  if (DEBUG) cat(file = stderr(), paste("tableSelectionServer", session$ns("test"), "\n"))
  ns <- session$ns
  # modSelectedRows <- c()
  # colOrder <- list()
  # searchStr <- ""
  # colState <- list()
  
  # can be removed. changed to the code below because it is more consitant with other assignments
  # assign(ns("colState"), list(), envir = .schnappsEnv)
  # assign(ns("pageLength"), 15, envir = .schnappsEnv)
  # assign(ns("colOrder"), list(), envir = .schnappsEnv)
  # assign(ns("modSelectedRows"), c(), envir = .schnappsEnv)
  # assign(ns("currentStart"), 0, envir = .schnappsEnv)
  
  .schnappsEnv[[ns("colState")]] = list()
  .schnappsEnv[[ns("pageLength")]] = 15
  .schnappsEnv[[ns("colOrder")]] = list()
  .schnappsEnv[[ns("modSelectedRows")]] = c()
  .schnappsEnv[[ns("currentStart")]] = 0
  
  
  observe(label = "cellNameTablesave2Hist", {
    clicked <- input$save2HistTabUi
    myns <- session$ns("cellNameTable")
    if (DEBUG) cat(file = stderr(), "observe input$save2HistTabUi \n")
    start.time <- base::Sys.time()
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        printTimeEnd(start.time, "cellNameTable")
        removeNotification(id = "save2Hist")
      }
    )
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("save2Hist", id = "save2Hist", duration = NULL)
    }
    if (is.null(clicked)) {
      if (DEBUG) cat(file = stderr(), "observe input$save2HistTabUi clicked=NULL\n")
      return()
    }
    if (clicked < 1) {
      if (DEBUG) cat(file = stderr(), paste("observe input$save2HistTabUi clicked < 1:", clicked,"\n"))
      return()
    }
    req(.schnappsEnv[[ns("historyTable")]])
    add2history(
      type = "renderDT", input = isolate( reactiveValuesToList(input)), comment = "Table",
      tableData = .schnappsEnv[[ns("historyTable")]]
    )
  })
  
  
  output$rowSelection <- renderText({
    if (DEBUG) cat(file = stderr(), "rowSelection\n")
    start.time <- Sys.time()
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        printTimeEnd(start.time, "rowSelection")
        removeNotification(id = "rowSelection")
      }
    )
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("rowSelection", id = "rowSelection", duration = NULL)
    }
    
    ns <- session$ns
    nsStr <- ns("-")
    dataTables <- dataTab()
    selectedRows <- input$cellNameTable_rows_selected
    scEx <- scEx()
    inputData = inputData()
    
    # update if expanded and not showing
    # input$refreshtable
    
    # we only need this for the removed genes table, so to not use too much memory we introduce this if statement
    # inputData <- NULL
    # if (nsStr == "gsRMGenesMod--") {
    #   inputData <- rowData(inputData()$scEx)
    # } else {
    #   inputData <- dataTables
    # }
    
    if (is.null(inputData)) {
      return(NULL)
    }
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("rowSelection", id = "rowSelection", duration = NULL)
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = paste0("~/SCHNAPPsDebug/rowSelection-", ns("bkup"), ".RData", collapse = "."),
           list = c(ls())
      )
    }
    # load(file=paste0("~/SCHNAPPsDebug/cellSelection-coE_topExpGenes-bkup.RData"))
    # deepDebug()
    # in case there is a table with multiple same row ids (see crPrioGenesTable) the gene names has "_#_" appended plus a number
    # remove this here
    if (length(selectedRows) > 0) {
      retVal <- rownames(dataTables[selectedRows, ,drop=F])
      retVal <- retVal[!is.na(retVal)]
      retVal <- sub("(.?)_##_(.*)", "\\1,\\2", retVal)
      retVal <- unlist(strsplit(retVal, ","))
      retVal <- sub("(.*)_#_.*", "\\1", retVal)
      retVal <- unique(retVal)
      # this removes everything other than row or col names
      # with just scEx we will cannot display the genes in the removed table
      # TODO need to check what we are comparing here. Just added "rowData(scEx)$symbol" because genes were not showing in some tables
      # also check that removed genes / cells table are working
      # retVal <- retVal[retVal %in% c(rowData(scEx)$symbol, inputData$symbol, colnames(scEx))]
      retVal <- paste0(retVal, collapse = ", ")
    } else {
      retVal <- NULL
    }
    .schnappsEnv[[ns("modSelectedRows")]] = retVal
    printTimeEnd(start.time, "tableSelectionServer-rowSelection")
    exportTestValues(tableSelectionServercellSelection = {
      retVal
    })
    if (DEBUG) cat(file = stderr(), paste("rowSelection retval: ",retVal, "\n"))
    return(retVal)
  })
  
  proxy <- DT::dataTableProxy("cellNameTable")
  
  observeEvent(input$selectAll, {
    # deepDebug()
    if (DEBUG) cat(file = stderr(), paste("observe input$selectAll",ns("test"),"\n"))
    ipSelect <- input$selectAll
    # prox <- proxy
    allrows <- input$cellNameTable_rows_all
    
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = paste0("~/SCHNAPPsDebug/inputselectAll.RData", collapse = "."),
           list = c(ls(), ls(envir = globalenv()))
      )
    }
    # load(file=paste0("~/SCHNAPPsDebug/inputselectAll.RData", collapse = "."))
    if (ipSelect) {
      proxy %>% DT::selectRows(selected = allrows)
    } else {
      proxy %>% DT::selectRows(selected = NULL)
    }
  })
  
  # observe: cellNameTable_rows_selected ----
  # observe(label = "ob44", {
  #   if (DEBUG) cat(file = stderr(), "observe input$cellNameTable_rows_selected\n")
  #   assign(ns("modSelectedRows"), input$cellNameTable_rows_selected, envir = .schnappsEnv)
  # })
  
  # observe: cellNameTable_state ----
  observe(label = "ob45", {
    if (DEBUG) cat(file = stderr(), "observe input$cellNameTable_state\n")
    # colOrder <<- input$cellNameTable_state$order
    # colState <<- input$cellNameTable_state
    if(!is.null(input$cellNameTable_state)){
      deepDebug()
      .schnappsEnv[[ns("colState")]] = input$cellNameTable_state
      .schnappsEnv[[ns("pageLength")]] = input$cellNameTable_state$pageLength
      .schnappsEnv[[ns("colOrder")]] = input$cellNameTable_state$order
      .schnappsEnv[[ns("currentStart")]] = input$cellNameTable_state$start
      tmp <- input$cellNameTable_state$search
    }
  })
  
  # renderDT cellNameTable ----
  output$cellNameTable <- DT::renderDT({
    myns <- session$ns("cellNameTable")
    if (DEBUG) cat(file = stderr(), "output$cellNameTable\n")
    start.time <- base::Sys.time()
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        printTimeEnd(start.time, "cellNameTable")
        removeNotification(id = "cellNameTable")
      }
    )
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("cellNameTable", id = "cellNameTable", duration = NULL)
    }
    
    dataTables <- dataTab()
    ns <- session$ns
    nsStr <- ns("-")
    reorderCells <- input$reorderCells
    selectedRows <- isolate(input$cellNameTable_rows_selected)
    showAllCells <- input$showAllCells
    showRowNames <- input$showRowNames
    
    # searchStr <-
    if (is.null(dataTables)) {
      .schnappsEnv[[ns("historyTable")]] <- NULL
      return(NULL)
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = paste0("~/SCHNAPPsDebug/cellNameTable-", nsStr, ".RData", collapse = "."),
           # list = c(ls(),ls(envir = .schnappsEnv))
           list = c(ls())
      )
    }
    #cp =  load(file=paste0("~/SCHNAPPsDebug/cellNameTable-sCA_dgeTable--.RData", collapse = "."))
    # cp =  load(file="~/SCHNAPPsDebug/cellNameTable-gQC_projCombTableMod--.RData")
    
    maxCol <- min(20, ncol(dataTables))
    if (dim(dataTables)[1] > 1) {
      # if the number of columns is big this can take really long
      # but this only happens if we are plotting the count matrix, so we assume that all columns are numeric
      if( ncol(dataTables) > 1000) {
        numericCols <- colnames(dataTables)
      } else {
        numericCols <- colnames(dataTables %>% select_if(is.numeric)) 
      }
      nonNumericCols <- which(!colnames(dataTables) %in% numericCols) # to keep non numeric columns...
      numericCols <- which(colnames(dataTables) %in% numericCols)
      if (reorderCells && length(selectedRows) > 0) {
        csums <- colSums(dataTables[selectedRows, numericCols,drop=F])
        cols2disp <- numericCols[order(csums, decreasing = TRUE)]
      } else {
        cols2disp <- numericCols
      }
      cols2disp <- c(nonNumericCols, cols2disp)[1:maxCol]
      if (showAllCells) cols2disp <- colnames(dataTables)
      dataTables <- as.data.frame(dataTables[, cols2disp,drop=F])
      colState <- .schnappsEnv[[ns("colState")]]
      if (length(colState) == 0) {
        colState <- list(
          orderClasses = TRUE,
          # setting to false, see https://stackoverflow.com/questions/46287971/column-headers-of-shiny-data-table-gets-shifted
          autoWidth = FALSE,
          scrollX = TRUE,
          # search = list(search = searchStr),
          stateSave = TRUE,
          order = .schnappsEnv[[ns("colOrder")]]
        )
        searchColList <- list()
      } else {
        searchColList <- list()
        for (cIdx in 1:length(colState$columns)) {
          searchColList[[cIdx]] <- colState$columns[[cIdx]]$search
        }
      }
      # if (DEBUG) cat(file = stderr(), paste(colState$search,"\n"))
      # deepDebug()
      dtout <- DT::datatable(dataTables,
                             rownames = showRowNames,
                             filter = "top",
                             selection = list(mode = "multiple"
                                              # , selected = get(ns("modSelectedRows"), envir = .schnappsEnv)
                             ),
                             caption = caption,
                             options = list(
                               orderClasses = TRUE,
                               autoWidth = TRUE,
                               scrollX = TRUE,
                               pageLength = .schnappsEnv[[ns("pageLength")]],
                               search = colState$search,
                               searchCols = searchColList,
                               stateSave = TRUE,
                               displayStart = .schnappsEnv[[ns("currentStart")]],
                               order = .schnappsEnv[[ns("colOrder")]]
                             )
      )
      .schnappsEnv[[ns("historyTable")]] <- dtout
      return(
        dtout
      )
    } else {
      return(warning("test"))
    }
  })
  
  output$download_cellNameTable <- downloadHandler(
    filename = function() {
      paste("cellNameTable", "table.csv", sep = "_")
    },
    content = function(file) {
      if (DEBUG) cat(file = stderr(), "download_cellNameTable started.\n")
      start.time <- base::Sys.time()
      on.exit({
        printTimeEnd(start.time, "download_cellNameTable")
        if (!is.null(getDefaultReactiveDomain())) {
          removeNotification(id = "download_cellNameTable")
        }
      })
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("download_cellNameTable", id = "download_cellNameTable", duration = NULL)
      }
      
      dataTables <- dataTab()
      write.csv(dataTables, file)
    }
  )
  
  return(reactive({
    if (DEBUG) cat(file = stderr(), "rowSelection.return\n")
    start.time <- Sys.time()
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        printTimeEnd(start.time, "rowSelection.return")
        removeNotification(id = "rowSelection.return")
      }
    )
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("rowSelection.return", id = "rowSelection.return", duration = NULL)
    }
    
    ns <- session$ns
    nsStr <- ns("-")
    dataTables <- dataTab()
    selectedRows <- input$cellNameTable_rows_selected
    inputData = inputData()
    
    # update if expanded and not showing
    # input$refreshtable
    
    # we only need this for the removed genes table, so to not use too much memory we introduce this if statement
    # inputData <- NULL
    # if (nsStr == "gsRMGenesMod--") {
    #   inputData <- rowData(inputData()$scEx)
    # } else {
    #   inputData <- dataTables
    # }
    
    if (is.null(inputData)) {
      return(NULL)
    }
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("rowSelection.return", id = "rowSelection.return", duration = NULL)
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = paste0("~/SCHNAPPsDebug/rowSelection-", ns("bkup"), ".return.RData", 
                         collapse = "."),list = c(ls())
      )
    }
    # load(file=paste0("~/SCHNAPPsDebug/cellSelection-coE_topExpGenes-bkup.RData"))
    # deepDebug()
    # in case there is a table with multiple same row ids (see crPrioGenesTable) the gene names has "___" appended plus a number
    # remove this here
    if (length(selectedRows) > 0) {
      retVal <- rownames(dataTables[selectedRows, ,drop=F])
      retVal <- retVal[!is.na(retVal)]
      retVal <- sub("(.?)_##_(.*)", "\\1,\\2", retVal)
      retVal <- unlist(strsplit(retVal, ","))
      retVal <- sub("(.*)_#_.*", "\\1", retVal)
      retVal <- unique(retVal)
      # this removes everything other than row or col names
      # with just scEx we will cannot display the genes in the removed table
      # TODO need to check what we are comparing here. Just added "rowData(scEx)$symbol" because genes were not showing in some tables
      # also check that removed genes / cells table are working
      # retVal <- retVal[retVal %in% c(rowData(scEx)$symbol, inputData$symbol, colnames(scEx))]
      retVal <- paste0(retVal, collapse = ", ")
    } else {
      retVal <- NULL
    }
    
    printTimeEnd(start.time, "tableSelectionServer-rowSelection.return")
    exportTestValues(tableSelectionServercellSelection = {
      retVal
    })
    if (DEBUG) cat(file = stderr(), paste("rowSelection retval: ",retVal, "\n"))
    return(retVal)
    
  }))
  
}




pHeatMapModule <- function(input, output, session,
                           pheatmapList # list with arguments for pheatmap
) {
  require(glue)
  require(evaluate)
  if (DEBUG) cat(file = stderr(), paste("pHeatMapModule", session$ns("test"), "\n"))
  ns <- session$ns
  outfilePH <- NULL
  heatmap_id = NULL
  
  # deepDebug()
  
  ht_obj = reactiveVal(NULL)
  ht_pos_obj = reactiveVal(NULL)
  myretVal = reactiveVal(NULL)
  
  observe(label = ns("parameterChanged"),{
    if (DEBUG) cat(file = stderr(), paste0("observe: ",ns("parameterChange"),"\n"))
    .schnappsEnv[[ns('sortingCols')]]  = input$sortingCols
    .schnappsEnv[[ns('normRow')]]  = input$normRow
    .schnappsEnv[[ns('heatMapGrpName')]]  = input$heatMapGrpName
    .schnappsEnv[[ns('ColNames')]]  = input$ColNames
    .schnappsEnv[[ns('colPal')]]  = input$colPal
    .schnappsEnv[[ns('orderNames')]]  = input$orderNames
    .schnappsEnv[[ns('heatmapCellGrp')]]  = input$heatmapCellGrp
    .schnappsEnv[[ns('sortingRows')]]  = input$sortingRows
    
    .schnappsEnv$defaultValues[[ns('sortingCols')]]  = input$sortingCols
    .schnappsEnv$defaultValues[[ns('normRow')]]  = input$normRow
    .schnappsEnv$defaultValues[[ns('heatMapGrpName')]]  = input$heatMapGrpName
    .schnappsEnv$defaultValues[[ns('ColNames')]]  = input$ColNames
    .schnappsEnv$defaultValues[[ns('colPal')]]  = input$colPal
    .schnappsEnv$defaultValues[[ns('orderNames')]]  = input$orderNames
    .schnappsEnv$defaultValues[[ns('heatmapCellGrp')]]  = input$heatmapCellGrp
    .schnappsEnv$defaultValues[[ns('sortingRows')]]  = input$sortingRows
  })
  
  addOptions <- reactive(
    {
      ph = pheatmapList()
      req(ph)
      list(
        sortingCols = input$sortingCols,
        normRow = input$normRow,
        heatMapGrpName = input$heatMapGrpName,
        ColNames = input$ColNames,
        colPal = input$colPal,
        orderColNames = input$orderNames,
        heatmapCellGrp = input$heatmapCellGrp,
        sortingRows = input$sortingRows,
        heatmapMinMaxValue = ifelse({
          # deepDebug();
          "mat" %in% names(pheatmapList())},
          c(min(pheatmapList()$mat), max(pheatmapList()$mat)),
          c(0,1)
        )
      )}) %>% debounce(1000)
  
  
  
  heatmapMinMaxValueDeb <- reactive(
    input$heatmapMinMaxValue
  ) %>% debounce(1000)
  
  heatmapCellGrpDeb <- reactive(
    input$heatmapCellGrp
  ) %>% debounce(1000)
  
  # pHeatMapModule - pHeatMapPlot ----
  # output$pHeatMapPlot <- renderImage(deleteFile = T,
  observe({
    if (DEBUG) cat(file = stderr(), "pHeatMapPlot started.\n")
    # deepDebug()
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "pHeatMapPlot")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "pHeatMapPlot")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("pHeatMapPlot", id = "pHeatMapPlot", duration = NULL)
    }
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "pHeatMapPlotWARNING")
    }
    
    ns <- session$ns
    heatmapData <- pheatmapList()
    addColNames <- addOptions()$ColNames
    orderColNames <- addOptions()$orderColNames
    # moreOptions <- input$moreOptions
    sortingCols <- addOptions()$sortingCols
    sortingRows <- addOptions()$sortingRows
    scale <- addOptions()$normRow
    pc = projectionColors %>% reactiveValuesToList()
    myns <- ns("pHeatMap")
    save2History <- input$save2History
    # pWidth <- input$heatmapWidth
    heatmapCellGrp <- heatmapCellGrpDeb()
    # pHeight <- input$heatmapHeight
    colPal <- addOptions()$colPal
    minMaxVal <- heatmapMinMaxValueDeb()
    # maxVal <- input$heatmapMaxValue
    # sampCol <- projectionColors$sampleNames
    # ccols <- projectionColors$dbCluster
    # force redraw
    input$pHeatMapPlot__shinyjquiBookmarkState__resizable$width
    input$pHeatMapPlot__shinyjquiBookmarkState__resizable$height
    proje <- projections()
    if (DEBUG) cat(file = stderr(), "output$pHeatMapModule:pHeatMapPlot\n")
    if (.schnappsEnv$DEBUGSAVE) {
      cat(file = stderr(), "output$pHeatMapModule:pHeatMapPlot saving\n")
      save(file = "~/SCHNAPPsDebug/pHeatMapPlotModule.RData", 
           list = c(ls()))
      cat(file = stderr(), "output$pHeatMapModule:pHeatMapPlot saving done\n")
    }
    # cp = load(file = "~/SCHNAPPsDebug/pHeatMapPlotModule.RData")
    outfile <- paste0(tempdir(), "/heatmap", ns("debug"), base::sample(1:10000, 1), ".png")
    outfile <- normalizePath(outfile, mustWork = FALSE)
    if (!"annotation_colors" %in% names(heatmapData)) {
      heatmapData$annotation_colors = list()
    }
    # browser()
    colorList = lapply(addColNames, FUN = function(x){
      # cat(file = stderr(), x)
      if(x %in% names(pc)){
        heatmapData$annotation_colors[[x]] = pc[[x]]
      }
    })
    names(colorList) = addColNames
    heatmapData$annotation_colors = colorList
    # if (!"sampleNames" %in% names(heatmapData$annotation_colors)) {
    #   heatmapData$annotation_colors$sampleNames = sampCol
    # }
    # if (!"dbCluster" %in% names(heatmapData$annotation_colors)) {
    #   heatmapData$annotation_colors$dbCluster = ccols
    # }
    updateSliderInput(session = session, inputId = "heatmapCellGrp",
                      max = min(200,ncol(heatmapData$mat)-1))
    
    retVal <- heatmapModuleFunction_m(
      heatmapData = heatmapData,
      addColNames = addColNames,
      orderColNames = orderColNames,
      sortingCols = sortingCols,
      sortingRows = sortingRows,
      scale = scale,
      colPal = colPal,
      minMaxVal = minMaxVal,
      proje = proje,
      heatmapCellGrp = heatmapCellGrp,
      outfile = outfile
    )
    myretVal(retVal)
    
    af = heatmapModuleFunction
    # remove env because it is too big
    environment(af) = new.env(parent = emptyenv())
    .schnappsEnv[[ns("historyHeatmap")]] <- list(plotFunc = af,
                                                 heatmapData = heatmapData,
                                                 addColNames = addColNames,
                                                 orderColNames = orderColNames,
                                                 sortingCols = sortingCols,
                                                 sortingRows = sortingRows,
                                                 scale = scale,
                                                 colPal = colPal,
                                                 minMaxVal = minMaxVal,
                                                 proje = proje,
                                                 heatmapCellGrp = heatmapCellGrp,
                                                 outfile = outfile
                                                 
    )
    # }
    # return(retVal)
    # deepDebug()
    
    output$pHeatMapPlot = renderPlot({
      cat(file = stderr(), "output$pHeatMapModule:rendering\n")
      if (is.null(retVal)){
        cat(file = stderr(), "output$pHeatMapModule:renderingNULL\n")
        return(NULL)
      }
      # deepDebug()
      if (.schnappsEnv$DEBUGSAVE) {
        save(file = "~/SCHNAPPsDebug/pHeatMapModule.RData",
             list = ls()
        )
      }
      # cp = load(file="~/SCHNAPPsDebug/pHeatMapModule.RData")
      
      ht = ComplexHeatmap::draw(retVal)
      ht_pos = htPositionsOnDevice(ht, include_annotation = FALSE, calibrate = FALSE)
      # cat(file = stderr(), glue_collapse(list_components(), sep="\n"))
      
      # fill reactive values:
      ht_obj(ht)
      ht_pos_obj(ht_pos)
      # deepDebug()
      # heatmap_id = shiny_env$current_heatmap_id
    })
  }) 
  
  observeEvent(label = "ob54",
               eventExpr = input$heatMapGrpNameButton,{
                 if (DEBUG) cat(file = stderr(), "observe input$heatMapGrpNameButton \n")
                 start.time <- base::Sys.time()
                 on.exit({
                   printTimeEnd(start.time, "observe input$heatMapGrpNameButton")
                   if (!is.null(getDefaultReactiveDomain())) {
                     removeNotification(id = "input-heatMapGrpNameButton")
                   }
                 })
                 if (!is.null(getDefaultReactiveDomain())) {
                   showNotification("observe: heatMapGrpNameButton", 
                                    id = "input-heatMapGrpNameButton", duration = NULL)
                 }
                 
                 heatmap_brush =  isolate(input$heatmap_brush)
                 newPrj = isolate(input$heatMapGrpName)
                 htobj = ht_obj()
                 htpos_obj = ht_pos_obj()
                 projections <- projections()
                 newPrjs <- projectionsTable$newProjections
                 acn = allCellNames()
                 htDat = myretVal()
                 # heatmap_id2 = shiny_env$current_heatmap_id
                 if (is.null(projections)) return(NULL)
                 if (is.null(htDat)) return(NULL)
                 if (is.null(htobj)) return(NULL)
                 # deepDebug()
                 if (.schnappsEnv$DEBUGSAVE) {
                   save(file = "~/SCHNAPPsDebug/heatMapGrpNameButton.RData",
                        list = ls()
                   )
                 }
                 # cp = load(file="~/SCHNAPPsDebug/heatMapGrpNameButton.RData")
                 if (is.null(heatmap_brush)) {
                   showNotification(
                     "No cells selected",
                     type = "error",
                     duration = NULL
                   )
                   return(NULL)
                 }
                 if (str_length(newPrj)<1) {
                   showNotification(
                     "New column name not set",
                     type = "error",
                     duration = NULL
                   )
                   return(NULL)
                 }
                 if (newPrj %in% colnames(projections)) {
                   showNotification(
                     "New column name already used",
                     type = "error",
                     duration = NULL
                   )
                   return(NULL)
                 }
                 # deepDebug()
                 newPrj = make.names(newPrj)
                 pos = getPositionFromBrush(brush = isolate(heatmap_brush), 
                                            ratio = 1)
                 
                 selection = tryCatch(
                   selectArea(ht_list = htobj, pos1 = pos[[1]], pos2 = pos[[2]], 
                              mark = T, ht_pos = htpos_obj, 
                              verbose = T, calibrate = FALSE), 
                   error = function(e) {
                     cat(file = stderr(), paste("inputData: NULL", e,"\n"))
                     return(NULL)}
                 )
                 if (is.null(selection)) {
                   save(file = "~/SCHNAPPsDebug/pHeatMapAreaNULL2.RData", list = c( ls()  ))
                   return(NULL)
                 }
                 # cp = load("~/SCHNAPPsDebug/pHeatMapAreaNULL2.RData")
                 
                 addPrj = rep(FALSE, nrow(projections))
                 names(addPrj) = rownames(projections)
                 mat=NULL
                 # the matrix can be in different places
                 if ("ht_list" %in% slotNames(htDat)){ 
                   mat = htDat@ht_list[[1]]@matrix
                 } else{ 
                   mat = htDat@matrix
                 }
                 addPrj[colnames(mat)[unlist(selection$column_index)]] = TRUE
                 # 
                 if (ncol(newPrjs) == 0) {
                   newPrjs = data.frame(row.names = acn)
                 }
                 newPrjs[,newPrj] = FALSE
                 newPrjs[names(addPrj),newPrj] <- addPrj
                 newPrjs[,newPrj] = as.factor(as.character(newPrjs[,newPrj]))
                 # } else {
                 #   # newPrjs <- cbind(newPrjs[rownames(projections), , drop = FALSE], projections[, addPrj, drop = FALSE])
                 #   # deepDebug()
                 #   newPrjs <- dplyr::left_join(
                 #     tibble::rownames_to_column(newPrjs),
                 #     tibble::rownames_to_column(projections[, addPrj, drop = FALSE]),
                 #     by='rowname')
                 #   rownames(newPrjs) = newPrjs[,1]
                 #   newPrjs = newPrjs[,-1]
                 # }
                 # colnames(newPrjs)[ncol(newPrjs)] <- newPrj
                 projectionsTable$newProjections <- newPrjs
               }) 
  
  renderGeneName = reactiveVal()
  # observer click ----
  observe( {
    if (DEBUG) cat(file = stderr(), "observe input$heatmap_click \n")
    heatmap_click = input$heatmap_click
    htobj = ht_obj()
    htpos_obj = ht_pos_obj()
    projections <- projections()
    newPrjs <- isolate(projectionsTable$newProjections)
    acn = allCellNames()
    htDat = myretVal()
    scEx_log = scEx_log()
    orderNames = isolate(addOptions()$orderColNames)
    if(is.null(input$heatmap_click)) return(NULL)
    if(is.null(scEx_log)) return(NULL)
    pos = getPositionFromClick(click = input$heatmap_click, 
                               ratio = 1)
    # save(file = "~/SCHNAPPsDebug/pHeatMapClick.RData", list = c( ls()  ))
    # cp = load("~/SCHNAPPsDebug/pHeatMapClick.RData")
    
    if(is.null(pos)) return(NULL)
    if(is.null(htobj)) return(NULL)
    if (is.null(projections)) return(NULL)
    if (is.null(htDat)) return(NULL)
    # deepDebug()
    selection = tryCatch(
      selectPosition(ht_list = htobj, pos = pos, mark = T, ht_pos = htpos_obj, 
                     verbose = F, calibrate = FALSE), 
      error = function(e) {
        cat(file = stderr(), paste("inputData: NULL", e,"\n"))
        return(NULL)}
    )
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/heatMapGrpNameClickButton.RData",
           list = ls()
      )
    }
    # cp = load(file="~/SCHNAPPsDebug/heatMapGrpNameClickButton.RData")
    # cat(file=stderr(), paste("\n\n",current.viewport()$name,"\n\n"))
    if (is.null(selection)) {
      # save(file = "~/SCHNAPPsDebug/pHeatMapClickNULL.RData", list = c( ls()  ))
      return(NULL)
      # cp = load("~/SCHNAPPsDebug/pHeatMapClickNULL.RData")
    }
    output$pHeatMapPlotSelection = renderPrint({
      print(selection)
    })
    
    if (is.null(heatmap_click)) {
      showNotification(
        "No cells selected",
        type = "error",
        duration = NULL
      )
      return(NULL)
    }
    
    # deepDebug()
    
    if (ncol(newPrjs) == 0) {
      newPrjs = data.frame(row.names = acn)
    }
    if ("ht_list" %in% slotNames(htDat)) {
      htMat = htDat@ht_list[[1]]@matrix
    } else if("matrix" %in% slotNames(htDat)) {
      htMat = htDat@matrix
    }
    
    # unlist(selection$row_index)
    geneName = rowData(scEx_log)[rownames(htMat)[selection$row_index],"symbol"]
    renderGeneName(geneName)
    
    
    if (!addOptions()$sortingCols == "gene (click)") return(NULL)
    newPrj = make.names(geneName)
    newPrjs[colnames(htMat),newPrj] <- htMat[selection$row_index,]
    projectionsTable$newProjections <- newPrjs
    
    updateSelectInput(inputId = "orderNames", 
                      selected = c(newPrj, orderNames))
    
    # clusterServer - return ----
    return(reactive({
      returnValues
    }))
    
    
    
  })
  
  # clusterServer - output$heatmapSelectedGenes ----
  output$heatmapSelectedGenes <- renderText({
    if (DEBUG)
      cat(file = stderr(), "heatmapSelectedGenes started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "heatmapSelectedGenes")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "heatmapSelectedGenes")
      }
      
    })
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("heatmapSelectedGenes", id = "heatmapSelectedGenes", duration = NULL)
    }
    # browser()
    geneName = renderGeneName()
    req(geneName)
    
    retVal <- paste(geneName, collapse = ", ")
    
    # exportTestValues(ClusterCellSelection = {
    #   retVal
    # })
    return(retVal)
  })
  
  
  # observe brush -----
  observe( {
    heatmap_brush = input$heatmap_brush
    htobj = ht_obj()
    htpos_obj = ht_pos_obj()
    scEx_log = scEx_log()
    htDat = myretVal()
    
    if(is.null(input$heatmap_brush)) return(NULL)
    pos = getPositionFromBrush(brush = input$heatmap_brush, 
                               ratio = 1)
    if(is.null(scEx_log)) return(NULL)
    if(is.null(pos)) return(NULL)
    if(is.null(htobj)) return(NULL)
    if(is.null(htpos_obj)) return(NULL)
    # save(file = "~/SCHNAPPsDebug/pHeatMapClick.RData", list = c( ls()  ))
    # cp = load("~/SCHNAPPsDebug/pHeatMapClick.RData")
    
    selection = tryCatch(
      selectArea(ht_list = htobj, pos1 = pos[[1]], pos2 = pos[[2]], 
                 mark = T, ht_pos = htpos_obj, 
                 verbose = T, calibrate = FALSE), 
      error = function(e) {
        cat(file = stderr(), paste("inputData: NULL", e,"\n"))
        return(NULL)}
    )
    if (is.null(selection)) {
      save(file = "~/SCHNAPPsDebug/pHeatMapAreaNULL.RData", list = c( ls()  ))
      return(NULL)
    }
    # cp = load("~/SCHNAPPsDebug/pHeatMapAreaNULL.RData")
    
    if (!is.null(selection)){
      output$pHeatMapPlotSelection = renderPrint({
        print(selection)
      })
      
    }
    
    if ("ht_list" %in% slotNames(htDat)) {
      htMat = htDat@ht_list[[1]]@matrix
    } else if("matrix" %in% slotNames(htDat)) {
      htMat = htDat@matrix
    } else {
      save(file = "~/SCHNAPPsDebug/pHeatMaphtDat.RData", list = c( ls()  ))
      return(NULL)
    }
    if(!is.null(selection)){
      r_index = selection$row_index[[1]]
      geneName = rowData(scEx_log)[rownames(htMat)[r_index],"symbol"]
      
      # clusterServer - output$heatmapSelectedGenes ----
      output$heatmapSelectedGenes <- renderText({
        retVal <- paste(geneName, collapse = ", ")
        return(retVal)
      })
    }
    
  }) 
  
  # observe save 2 history ----
  observe(label = "save2histHM", {
    clicked <- input$save2HistHM
    if (DEBUG) cat(file = stderr(), "observe input$save2HistHM \n")
    myns <- ns("pHeatMap")
    # deepDebug()
    req(.schnappsEnv[[ns("historyHeatmap")]])
    start.time <- base::Sys.time()
    if (DEBUG) cat(file = stderr(), "cluster: save2HistHM\n")
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "save2HistHM")
      }
    )
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("save2HistHM", id = "save2HistHM", duration = NULL)
    }
    
    if (is.null(clicked)) {
      return()
    }
    if (clicked < 1) {
      return()
    }
    add2history(
      type = "save", input = isolate( reactiveValuesToList(input)),
      comment = paste("# ",myns, "\n",
                      "fun = plotData$plotData$plotFunc\n", 
                      "environment(fun) = environment()\n",
                      "plotData$plotData$outfile=NULL\n",
                      "print(do.call(\"fun\",plotData$plotData[2:length(plotData$plotData)]))\n"
      ),
      plotData = .schnappsEnv[[ns("historyHeatmap")]]
    )
  })
  
  # pHeatMapModule - updateInput ----
  # updateInput <-
  # this is calling projections during loading of data
  observe(label = "ob46", {
    if (DEBUG) cat(file = stderr(), "observer: updateInput started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "updateInput")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "updateInput")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("updateInput", id = "updateInput", duration = NULL)
    }
    
    proje <- projections()
    
    # Can use character(0) to remove all choices
    if (is.null(proje)) {
      return(NULL)
    }
    
    # Can also set the label and select items
    updateSelectInput(session, "ColNames",
                      choices = colnames(proje),
                      selected = .schnappsEnv[[ns('ColNames')]]
    )
    updateSelectInput(session, "orderNames",
                      choices = colnames(proje),
                      selected = .schnappsEnv[[ns('orderNames')]]
    )
    
    # updateSelectInput(session, "dimension_y",
    #                   choices = colnames(proje),
    #                   selected = colnames(proje)[2]
    # )
  }) 
  
  observe(label = "ob46a", {
    if (DEBUG) cat(file = stderr(), "observer: updateInput started 46a.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "updateInput")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "updateInput")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("updateInput", id = "updateInput", duration = NULL)
    }
    # deepDebug()
    # browser()
    heatmapData <- pheatmapList()
    if(is.null(heatmapData)) return(NULL)
    # to remove the warning message when nothing is loaded
    if (!"mat" %in% names(heatmapData)){
      heatmapData$mat = 0
      min = 0
      max = 200
    }else{
      min = signif(min(heatmapData$mat), digits = 4)
      max = signif(max(heatmapData$mat), digits = 4)
      if (input$normRow == "row_order"){
        min = 0
        max = ncol(heatmapData$mat)
      }
      if (input$normRow == "col_order"){
        min = 0
        max = nrow(heatmapData$mat)
      }
    }
    if (DEBUG) cat(file = stderr(), paste("minmax",min,max,"\n"))
    updateSliderInput(session,
                      inputId = "heatmapMinMaxValue",
                      min = min,
                      max = max
                      # ,
                      # value = c(min, max)
    )
    # updateNumericInput(session, inputId = "heatmapMaxValue",
    #                    min = min(heatmapData$mat),
    #                    max = max(heatmapData$mat)
    # )
  }) 
  
  
  
  
  # pHeatMapModule - additionalOptions ----
  # output$additionalOptions <- renderUI({
  #   if (DEBUG) cat(file = stderr(), "additionalOptions started.\n")
  #   start.time <- base::Sys.time()
  #   on.exit({
  #     printTimeEnd(start.time, "additionalOptions")
  #     if (!is.null(getDefaultReactiveDomain())) {
  #       removeNotification(id = "additionalOptions")
  #     }
  #   })
  #   if (!is.null(getDefaultReactiveDomain())) {
  #     showNotification("additionalOptions", id = "additionalOptions", duration = NULL)
  #   }
  #
  #   ns <- session$ns
  #   moreOptions <- (input$moreOptions)
  #   groupNs <- groupNames$namesDF
  #   proje <- projections()
  #   if (!moreOptions | is.null(proje)) {
  #     return("")
  #   }
  #
  #
  #   if (.schnappsEnv$DEBUGSAVE) {
  #     save(file = "~/SCHNAPPsDebug/heatMapadditionalOptions.RData", list = c(ls()))
  #   }
  #   # load(file="~/SCHNAPPsDebug/heatMapadditionalOptions.RData")
  #
  #   tagList(
  #    )
  # })
  #
  # pHeatMapModule - download_pHeatMapUI ----
  output$download_pHeatMapUI <- downloadHandler(
    filename = function() {
      paste("pHeatMap", "data.zip", sep = "_")
    },
    content = function(file) {
      if (DEBUG) cat(file = stderr(), "download_pHeatMapUI started.\n")
      start.time <- base::Sys.time()
      on.exit({
        printTimeEnd(start.time, "download_pHeatMapUI")
        if (!is.null(getDefaultReactiveDomain())) {
          removeNotification(id = "download_pHeatMapUI")
        }
      })
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("download_pHeatMapUI", id = "download_pHeatMapUI", duration = NULL)
      }
      
      heatmapData <- pheatmapList()
      addColNames <- addOptions()$ColNames
      orderColNames <- addOptions()$orderColNames
      # moreOptions <- input$moreOptions
      groupNs <- groupNames$namesDF
      scale <- addOptions()$normRow
      proje <- projections()
      if (.schnappsEnv$DEBUGSAVE) {
        save(file = "~/SCHNAPPsDebug/download_pHeatMapUI.RData", list = c(
          "outfilePH", ls(),
          ls(.schnappsEnv)
        ))
      }
      # load("~/SCHNAPPsDebug/download_pHeatMapUI.RData")
      dfilename <- paste0(.schnappsEnv$reportTempDir, "/sessionData.RData")
      base::save(file = dfilename, list =
                   c("heatmapData", "addColNames", "orderColNames", "proje", "groupNs", "scale")
      )
      
      
      # 2do: not sure why I cannot find the original file...
      # maybe there is an intermediate session created?
      outfile <- paste0(tempdir(), "/heatmap", ns("debug"), base::sample(1:10000, 1), ".png")
      outfile <- normalizePath(outfile, mustWork = FALSE)
      heatmapData$filename <- outfile
      
      # if (length(addColNames) > 0 & moreOptions) {
      if (length(addColNames) > 0) {
        heatmapData$annotation_col <- proje[rownames(heatmapData$annotation_col), addColNames, drop = FALSE]
      }
      # if (sum(orderColNames %in% colnames(proje)) > 0 & moreOptions) {
      if (sum(orderColNames %in% colnames(proje)) > 0) {
        heatmapData$cluster_cols <- FALSE
        colN <- rownames(psychTools::dfOrder(proje, orderColNames))
        colN <- colN[colN %in% colnames(heatmapData$mat)]
        heatmapData$mat <- heatmapData$mat[, colN, drop = FALSE]
      }
      do.call(pheatmap, heatmapData)
      
      
      zippedReportFiles <- c(
        dfilename,
        outfile
      )
      zip(file, zippedReportFiles, flags = "-9Xj")
    }
  )
}
#


cellSelectionModule <- function(input, output, session) {
  if (DEBUG) cat(file = stderr(), paste("cellSelectionModule", session$ns("my-test"), "\n"))
  ns <- session$ns
  
  
  # init valuse for Env to remember selections
  # and observe/save if changed
  
  # browser()
  .schnappsEnv[[ns("Mod_clusterPP")]] = "dbCluster"
  observe(label = "Mod_clusterPP", 
          {
            if (DEBUG) cat(file = stderr(), paste0("observe: Mod_clusterPP\n"))
            # assign(ns("Mod_clusterPP"), input$Mod_clusterPP, envir = .schnappsEnv)
            # browser()
            # initialization by history
            if(input$Mod_clusterPP=="") {
              if("defaultValues" %in% names(.schnappsEnv))
                if(ns("Mod_clusterPP") %in% names(.schnappsEnv$defaultValues))
                  .schnappsEnv[[ns("Mod_clusterPP")]] = .schnappsEnv$defaultValues[[ns("Mod_clusterPP")]]
              return()
            }
            
            .schnappsEnv[[ns("Mod_clusterPP")]] = input$Mod_clusterPP
            .schnappsEnv$defaultValues[[ns("Mod_clusterPP")]] = input$Mod_clusterPP
          })
  
  
  .schnappsEnv[[ns("Mod_PPGrp")]] = "1"
  observe(label = "Mod_PPGrp",
          {
            if (DEBUG) cat(file = stderr(), paste0("observe: Mod_PPGrp",input$Mod_PPGrp,"\n"))
            # assign(ns("Mod_PPGrp"), input$Mod_PPGrp, envir = .schnappsEnv)
            .schnappsEnv[[ns("Mod_PPGrp")]] = input$Mod_PPGrp
            .schnappsEnv$defaultValues[[ns("Mod_PPGrp")]] = input$Mod_PPGrp
            
          })
  
  
  # Mod_updateInputPPt if projections changed ====
  observe(label = "Mod_clusterPP", 
          {
            if (DEBUG) cat(file = stderr(), paste("Mod_updateInputPPt started.", 
                                                  .schnappsEnv[[ns("Mod_clusterPP")]],"\n"))
            start.time <- base::Sys.time()
            on.exit({
              printTimeEnd(start.time, "Mod_updateInputPPt")
              if (!is.null(getDefaultReactiveDomain())) {
                removeNotification(id = "Mod_updateInputPPt")
              }
            })
            if (!is.null(getDefaultReactiveDomain())) {
              showNotification("Mod_updateInputPPt", id = "Mod_updateInputPPt", duration = NULL)
            }
            projections <- projections()
            projFactors <- projFactors()
            # Can use character(0) to remove all choices
            if (is.null(projections)) {
              return(NULL)
            }
            # save(file = "~/SCHNAPPsDebug/Mod_updateInputPPt", list = c(ls(), ls(envir = globalenv())))
            # load(file = "~/SCHNAPPsDebug/Mod_updateInputPPt")
            
            # coln <- colnames(projections)
            # choices <- c()
            # for (cn in coln) {
            #   if (length(levels(as.factor(projections[, cn]))) < 50) {
            #     choices <- c(choices, cn)
            #   }
            # }
            # 
            # browser()
            
            updateSelectInput(
              session,
              "Mod_clusterPP",
              choices = unique(projFactors, .schnappsEnv[[ns("Mod_clusterPP")]]),
              selected = .schnappsEnv[[ns("Mod_clusterPP")]]
            )
          })
  
  observe({
    projections <- projections()
    
    if (DEBUG) cat(file = stderr(), paste("observeEvent:: input$Mod_clusterPP: \n"))
    # Can use character(0) to remove all choices
    if (is.null(projections)) {
      return(NULL)
    }
    if (!input$Mod_clusterPP %in% colnames(projections)) {
      return(NULL)
    }
    # deepDebug()
    choicesVal <- levels(projections[, input$Mod_clusterPP])
    updateSelectInput(
      session,
      "Mod_PPGrp",
      choices = choicesVal,
      selected = .schnappsEnv$defaultValues[[ns("Mod_PPGrp")]]
    )
  })
  
  
  returnValues <- reactiveValues(
    # return cell names
    # return description of selection
    cellNames = reactive(
      {
        prjNames <- input$Mod_clusterPP
        prjVals <- input$Mod_PPGrp
        projections <- projections()
        req(projections)
        req(prjNames)
        req(prjVals)
        # browser() 
        # debugControl("cellSelectionModule", list = c(ls()))
        # cp = load(file = "~/SCHNAPPsDebug/cellSelectionModule.RData")
        if (length(prjVals) == 0 | str_length(prjNames) == 0 |
            !prjNames %in% colnames(projections)) {
          return(NULL)
        } else {
          retVal <- rownames(projections[projections[, prjNames] %in% prjVals, ])
          retVal
        }
      }),
    selectionDescription = reactive(
      {
        prjNames <- input$Mod_clusterPP
        prjVals <- input$Mod_PPGrp
        
        if (.schnappsEnv$DEBUGSAVE) {
          save(file = "~/SCHNAPPsDebug/selectionDescription.RData", list = c(ls()))
        }
        # cp = load(file = "~/SCHNAPPsDebug/selectionDescription.RData")
        retVal <- paste("projections: ", prjNames, "with levels:", paste(prjVals, collapse = ", "))
        retVal
      }),
    ProjectionUsed = reactive(
      {
        prjNames <- input$Mod_clusterPP
        prjNames
      }),
    ProjectionValsUsed = reactive(
      {
        if (DEBUG) cat(file = stderr(), "ProjectionValsUsed: ",.schnappsEnv[[ns("Mod_PPGrp")]],"\n")
        prjGrp <- input$Mod_PPGrp
        prjGrp
      })
  )
  
  if (DEBUG) cat(file = stderr(), paste("cellSelectionModule", session$ns("my-test:end"), "\n"))
  return(reactive({
    returnValues
  }))
}


if (DEBUG) {
  cat(file = stderr(), "done loading Module server.\n")
}
