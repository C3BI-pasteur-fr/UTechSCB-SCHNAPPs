
# just in case?
# if ( exists(envir = .schnappsEnv, x=".SCHNAPPs_LiteData")){
#   
#   return (NULL)
#   
# }

# inputData ----
inputData <- reactive({NULL})

if (!exists('AllowClustering')) {
  AllowClustering = FALSE
  if (DEBUG) cat(file = stderr(), "ui-lite: AllowClustering not defined\n")
}



# scEx_log ----
scEx_log <- reactive({
  cat(file = stderr(), green("lite: scEx_log\n"))
  get(".SCHNAPPs_LiteData",envir = .schnappsEnv)$scEx_log
})

# scEx ----
scEx <- reactive({
  cat(file = stderr(), green("lite: scEx\n"))
  get(".SCHNAPPs_LiteData",envir = .schnappsEnv)$scEx
})

# pca ----
if (!AllowClustering)
  pcaReact <- reactive({
    cat(file = stderr(), green("lite: pca\n"))
    get(".SCHNAPPs_LiteData",envir = .schnappsEnv)$pca
  })

# dbCluster ----
if (!AllowClustering)
  dbCluster <- reactive({
    cat(file = stderr(), green("lite: dbCluster\n"))
    get(".SCHNAPPs_LiteData",envir = .schnappsEnv)$projections$dbCluster
  })

# projections ----
# if (!AllowClustering)
  projections <- reactive({
    if (DEBUG) {
      cat(file = stderr(), green("projections lite started.\n"))
    }
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "projections")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "projections")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("projections", id = "projections", duration = NULL)
    }
    # scEx is the fundamental variable with the raw data, which is available after loading
    # data. Here we ensure that everything is loaded and all varialbles are set by waiting
    # input data being loaded
    scEx <- scEx()
    
    # pca is already fixed in the calculation so we don't need to recalculate it.
    # pcaReact <- NULL
    
    # in schnapps-lite we are only interested in adding the sessionProjections
    # to the already calculated other projections
    prjs <- sessionProjections$prjs
    newPrjs <- projectionsTable$newProjections
    if (!exists("scEx") |
        is.null(scEx)) {
      if (DEBUG) {
        cat(file = stderr(), "sampleInfo: NULL\n")
      }
      return(NULL)
    }
    # if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/projections.RData", list = c(ls()))
    # }
    # cp = load(file="~/SCHNAPPsDebug/projections.RData"); DEBUGSAVE=FALSE
    
      
    
    # browser()
    
    # todo colData() now returns a s4 object of class DataFrame
    # not sure what else is effected...
    pd <- as.data.frame(colData(scEx))
    if (ncol(pd) < 2) {
      cat(file = stderr(), "phenoData for scEx has less than 2 columns\n")
      return(NULL)
    }
    projections <- pd
    
    # Commenting out during dev. here we are usually calculating the user defined new projections like tsne, umap
    # etc which are already calculated.
    # 
    # withProgress(message = "Performing projections", value = 0, {
    #   n <- length(.schnappsEnv$projectionFunctions)
    #   iter <- 1
    #   for (proj in .schnappsEnv$projectionFunctions) {
    #     start.time1 <- Sys.time()
    #     incProgress(1 / n, detail = paste("Creating ", proj[1]))
    #     if (DEBUG) {
    #       cat(file = stderr(), paste("calculation projection:  ", proj[1], "\n"))
    #     }
    #     if (DEBUG) cat(file = stderr(), paste("projection: ", proj[2], "\n"))
    #     assign("tmp", eval(parse(text = paste0(proj[2], "()"))))
    #     if (.schnappsEnv$DEBUGSAVE) {
    #       save(
    #         file = paste0("~/SCHNAPPsDebug/projections.", iter, ".RData"),
    #         list = c("tmp")
    #       )
    #       iter <- iter + 1
    #     }
    #     # load(file="~/SCHNAPPsDebug/projections.1.RData")
    #     # browser()
    #     # TODO here, dbCluster is probably overwritten and appended a ".1"
    #     if (is(tmp, "data.frame")) {
    #       cn <- make.names(c(colnames(projections), colnames(tmp)), unique = TRUE)
    #     } else {
    #       cn <- make.names(c(colnames(projections), make.names(proj[1])), unique = TRUE)
    #     }
    #     if (length(tmp) == 0) {
    #       next()
    #     }
    #     if (ncol(projections) == 0) {
    #       # never happening because we set pca first
    #       projections <- data.frame(tmp = tmp)
    #     } else {
    #       if (nrow(projections) == length(tmp)) {
    #         projections <- cbind(projections, tmp)
    #       } else {
    #         if (!is.null(nrow(tmp))) {
    #           if (nrow(projections) == nrow(tmp)) {
    #             projections <- cbind(projections, tmp)
    #           }
    #         } else {
    #           save(file = "~/SCHNAPPsDebug/projectionsError.RData", list = c(ls()))
    #           stop("error: ", proj[1], "didn't produce a result, please send file ~/SCHNAPPsDebug/projectionsError.RData to bernd")
    #         }
    #       }
    #       # else {
    #       #   stop("error: ", proj[1], "didn't produce a result")
    #       # }
    #     }
    #     if (!length(colnames(projections)) == length(cn)) {
    #       save(file = "~/SCHNAPPsDebug/projectionsError2.RData", list = c(ls()))
    #       stop("error: ", proj[1], "didn't produce a result, please send file ~/SCHNAPPsDebug/projectionsError2.RData to bernd")
    #     }
    #     colnames(projections) <- cn
    #     if (DEBUG) cat(file = stderr(), paste("colnames ", paste0(colnames(projections), collapse = " "), "\n"))
    #     if (DEBUG) cat(file = stderr(), paste("observe this: ", proj[2], "\n"))
    #     # observe(proj[2], quoted = TRUE)
    #   }
    # })
    # add a column for gene specific information that will be filled/updated on demand
    # projections$UmiCountPerGenes <- 0
    # projections$UmiCountPerGenes2 <- 0
    # for (pdIdx in colnames(pd)) {
    #   if (!pdIdx %in% colnames(projections)) {
    #     projections[, pdIdx] <- pd[, pdIdx]
    #   }
    # }
    
    if (ncol(prjs) > 0 & nrow(prjs) == nrow(projections)) {
      projections <- cbind(projections, prjs)
    } else if (ncol(prjs) > 0) {
      commIds = intersect(rownames(prjs), rownames(projections))
      missing = rownames(projections)[!rownames(projections) %in% rownames(prjs)]
      if (length(missing)>0) {
        prjs[missing,] = NA
      }
      projections <- cbind(projections, prjs[rownames(projections),])
    }
    # remove columns with only one unique value
    rmC <- c()
    for (cIdx in 1:ncol(projections)) {
      # ignore sampleNames
      if (colnames(projections)[cIdx] == "sampleNames") next()
      if (length(unique(projections[, cIdx])) == 1) rmC <- c(rmC, cIdx)
    }
    if (length(rmC) > 0) projections <- projections[, -rmC]
    
    if (ncol(newPrjs) > 0) {
      projections <- cbind(projections, newPrjs[rownames(projections), , drop = FALSE])
    }
    # in case (no Normalization) no clusters or sample names have been assigned
    if (!"dbCluster" %in% colnames(projections)) {
      projections$dbCluster <- 0
    }
    if (!"sampleNames" %in% colnames(projections)) {
      projections$sampleNames <- "1"
    }
    # TODO figure out how to limit this.
    # add2history(type = "save", input = isolate( reactiveValuesToList(input)), comment = "projections", projections = projections)
    
    
    exportTestValues(projections = {
      projections
    })
    return(projections)
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
  # if (clicked < 1) {
  #   return(list(
  #     src = "",
  #     contentType = "image/png",
  #     width = 10,
  #     height = 10,
  #     alt = "Scater plot will be here when 'apply changes' is checked"
  #   ))
  # }
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
    save(file = "~/SCHNAPPsDebug/scater.Rmd", list = c(ls()))
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
  
  outfile <- paste0(getwd(), "/scaterPlot.png")
  if (file.exists(normalizePath(outfile, mustWork = FALSE)) & clicked == 0){
    return(list(
      src = normalizePath(outfile, mustWork = FALSE),
      contentType = "image/png",
      width = width,
      height = height,
      alt = "Scater plot should be here"
    ))
  }
  cat(file = stderr(), paste("getwd:", getwd(), "didn't contain png file.\n"))
  # outfile <- paste0("~/SCHNAPPsDebug",'/scaterPlot.png')
  if (DEBUG) cat(file = stderr(), paste("output file: ", outfile, "\n"))
  if (DEBUG) cat(file = stderr(), paste("output file normalized: ", normalizePath(outfile, mustWork = FALSE), "\n"))
  n <- min(nrow(scaterReads), 50)
  
  rownames(scaterReads) <- rowData(scaterReads)$symbol
  # p1 <- scater::plotHighestExprs(scaterReads, colour_cells_by = "sampleNames", n = n)
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
  # browser()
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

# groupNames ----
# store cell groups that are defined on the fly using the modular 2D plot
groupNames <- reactiveValues(
  namesDF = if(!is.null(.schnappsEnv$.SCHNAPPs_LiteData[["namesDF"]])) {
    .schnappsEnv$.SCHNAPPs_LiteData[["namesDF"]]
  } else {
    data.frame()
  }
)


if (DEBUG) cat(file = stderr(), "end: reactives-lite\n")

