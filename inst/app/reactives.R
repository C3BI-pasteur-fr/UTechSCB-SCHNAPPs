suppressMessages(require(tibble))
suppressMessages(require(scran))
suppressMessages(require(irlba))
suppressMessages(require(BiocSingular))

# reactive values  ------------------------------------------------------------------
inputFileStats <- reactiveValues(stats = NULL)

# store cell groups that are defined on the fly using the modular 2D plot
groupNames <- reactiveValues(namesDF = data.frame())

# colors for samples
sampleCols <- reactiveValues(colPal = allowedColors)

# colors for clusters
clusterCols <- reactiveValues(colPal = allowedColors)

# Here, we store projections that are created during the session. These can be selections of cells or other values that
# are not possible to precalculate.
sessionProjections <- reactiveValues(prjs = data.frame())

DEBUGSAVE <- get(".SCHNAPPs_DEBUGSAVE", envir = .schnappsEnv)

# Input file either rdata file or csv file

inputFile <- reactiveValues(
  inFile = "",
  annFile = ""
)
if ("crayon" %in% rownames(installed.packages()) == FALSE) {
  green <- function(x) {
    x
  }
} else {
  require(crayon)
}


# inputDataFunc ----
# loads singleCellExperiment
#   only counts, rowData, and colData are used. Everything else needs to be recomputed
inputDataFunc <- function(inFile) {
  if (DEBUG) {
    cat(file = stderr(), "inputDataFunc started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "inputDataFunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "inputDataFunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("inputDataFunc", id = "inputDataFunc", duration = NULL)
  }

  stats <- tibble(.rows = length(inFile$datapath))
  stats$names <- inFile$name
  stats$nFeatures <- 0
  stats$nCells <- 0

  #
  cat(file = stderr(), paste("reading", inFile$name[1], "\n"))
  fp <- inFile$datapath[1]
  # fp ="scEx.RData"
  # fp ="../SCHNAPPsData/patty1A.v2.RData"
  # a bit of cleanup
  for (v in c("scEx", "scEx_log", "featureData")) {
    if (exists(v)) rm(v)
  }
  fpLs <- tryCatch(load(fp), error = function(e) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("Error reading input file!!",
        id = "inputDataFuncERROR",
        type = "error", duration = NULL
      )
    }
    if (DEBUG) {
      cat(file = stderr(), "\n\nERROR reading data.\n\n\n")
    }
    NULL
  })

  # in case we are loading from a report
  scExFound <- FALSE
  if ("scEx" %in% fpLs) {
    scExFound <- TRUE
    varName <- "scEx"
  } else {
    scExFound <- FALSE
    for (varName in fpLs) {
      # if ("SingleCellExperiment" %in% class(get(varName))) {
      if (is(get(varName), "SingleCellExperiment")) {
        scEx <- get(varName)
        scExFound <- TRUE
        break()
      }
    }
  }
  if (!scExFound) {
    return(NULL)
  }
  cat(file = stderr(), paste("file ", inFile$name[1], "contains variable", varName, " as SingleCellExperiment.\n"))
  # save(file = "~/SCHNAPPsDebug/inputProblem.RData", list = ls())
  # load("~/SCHNAPPsDebug/inputProblem.RData")
  fdAll <- rowData(scEx)
  pdAll <- colData(scEx)
  exAll <- assays(scEx)[["counts"]]
  stats[1, "nFeatures"] <- nrow(fdAll)
  stats[1, "nCells"] <- nrow(pdAll)

  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/readInp1.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/SCHNAPPsDebug/readInp1.RData')

  # read multiple files
  if (length(inFile$datapath) > 1) {
    for (fpIdx in 2:length(inFile$datapath)) {
      inFile$datapath[fpIdx] <- "~/Downloads/paper1.RData"
      cat(file = stderr(), paste("reading", inFile$name[fpIdx], "\n"))
      fp <- inFile$datapath[fpIdx]
      fpLs <- load(fp)
      scExFound <- FALSE
      for (varName in fpLs) {
        # if ("SingleCellExperiment" %in% class(get(varName))) {
        if (is(get(varName), "SingleCellExperiment")) {
          scEx <- get(varName)
          scExFound <- TRUE
          break()
        }
      }
      cat(file = stderr(), paste("file ", inFile$datapath[fpIdx], "contains variable", varName, "\n"))
      if (!scExFound) {
        next()
      }
      rListAll <- rownames(fdAll)[!rownames(fdAll) %in% rownames(rowData(scEx))]
      rListNew <- rownames(rowData(scEx))[!rownames(rowData(scEx)) %in% rownames(fdAll)]

      # fdIdx <- intersect(rownames(fdAll), rownames(rowData(scEx)))
      # if (length(fdIdx) != nrow(fd)) {
      #   cat(file = stderr(), "Houston, there is a problem with the features\n")
      # }
      # fd <- featuredata[fdIdx, ]
      rowdataNew <- rowData(scEx)[rListNew, ]
      fdAll[rListNew, ] <- NA
      for (cIdx in colnames(rowdataNew)) {
        if (!cIdx %in% colnames(fdAll)) {
          fdAll <- cbind(fdAll, data.frame(row.names = rownames(fdAll), rep(NA, nrow(fdAll))))
          colnames(fdAll)[length(colnames(fdAll))] <- cIdx
        }
        fdAll[rListNew, cIdx] <- rowdataNew[, cIdx]
      }
      # fdAll[rListNew,"id"] = rListNew

      # append new genes
      tmpMat <- matrix(nrow = length(rListNew), ncol = ncol(exAll), data = 0)
      rownames(tmpMat) <- rListNew
      exAll <- rbind(exAll, tmpMat)



      # fdAll <- fdAll[fdIdx, ]
      pd1 <- colData(scEx)
      rownames(pd1) <- paste0(rownames(pd1), "-", fpIdx)
      if (!"counts" %in% names(assays(scEx))) {
        cat(file = stderr(), paste("file ", inFile$datapath[fpIdx], "with variable", varName, "didn't contain counts slot\n"))
        next()
      }
      ex1 <- assays(scEx)[["counts"]][fdIdx, ]
      pdAllCols2add <- colnames(pdAll)[!colnames(pdAll) %in% colnames(pd1)]
      pd1Cols2add <- colnames(pd1)[!colnames(pd1) %in% colnames(pdAll)]
      for (pd1C in pd1Cols2add) {
        pdAll[, pd1C] <- NA
      }
      for (pd1C in pdAllCols2add) {
        pd1[, pd1C] <- NA
      }
      if (sum(rownames(pdAll) %in% rownames(pd1)) > 0) {
        cat(green("Houston, there are cells with the same name\n"))
        save(file = "~/SCHNAPPsDebug/inputProblem.RData", list = c("pdAll", "pd1", "exAll", "ex1", "fpIdx", "scEx", "stats", "fdAll"))
        # load("~/SCHNAPPsDebug/inputProblem.RData")
        cat(file = stderr(), "Houston, there are cells with the same name\n")
        rownames(pd1) <- paste0(rownames(pd1), "_", fpIdx)
        pdAll <- rbind(pdAll, pd1)
        colnames(ex1) <- rownames(pd1)
      } else {
        pdAll <- rbind(pdAll, pd1)
      }
      # stats[fpIdx, "nFeatures"] <- nrow(fd)
      stats[fpIdx, "nCells"] <- nrow(pd1)

      nrow(ex1)
      nrow(exAll)
      exAll <- Matrix::cbind2(exAll[rownames(exAll), ], ex1[rownames(exAll), ])
    }
  }
  exAll <- as(exAll, "dgTMatrix")


  if ("sampleNames" %in% colnames(pdAll)) {
    if (!is(pdAll$sampleNames, "factor")) {
      pdAll$sampleNames <- factor(pdAll$sampleNames)
    }
    sampNames <- levels(pdAll$sampleNames)
    isolate({
      # sampleCols$colPal <- colorRampPalette(brewer.pal(
      #   n = 6, name =
      #     "PRGn"
      # ))(length(sampNames))
      sampleCols$colPal <- allowedColors[seq_along(sampNames)]
      names(sampleCols$colPal) <- sampNames
    })
  } else {
    showNotification(
      "scEx - colData doesn't contain sampleNames",
      duration = NULL,
      type = "error"
    )
    pdAll$sampleNames <- 1
  }


  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/readInp.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/SCHNAPPsDebug/readInp.RData')

  scEx <- SingleCellExperiment(
    assay = list(counts = exAll),
    colData = pdAll,
    rowData = fdAll
  )

  dataTables <- list()
  featuredata <- rowData(scEx)
  # handle different extreme cases for the symbol column (already encountered)
  if (is.factor(featuredata$symbol)) {
    if (levels(featuredata$symbol) == "NA") {
      featuredata$symbol <- make.unique(rownames(featuredata))
      rowData(scEx) <- featuredata
    }
  }
  if ("symbol" %in% colnames(featuredata)) {
    featuredata$symbol <- make.unique(featuredata$symbol)
    rowData(scEx) <- featuredata
  }

  # dataTables$featuredataOrg <- rowData(scEx)
  dataTables$scEx <- scEx
  dataTables$featuredata <- featuredata

  if (is.null(scEx$barcode)) {
    showNotification("scEx doesn't contain barcode column", type = "error")
    return(NULL)
  }
  # some checks

  if (sum(is.infinite(assays(scEx)[["counts"]])) > 0) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("scEx contains infinite values",
        type = "error"
      )
    }
    return(NULL)
  }


  if (sum(c("id", "symbol") %in% colnames(rowData(scEx))) < 2) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "scEx - rowData doesn't contain id and/or symbol columns",
        duration = NULL,
        type = "error"
      )
    }
  }

  if (!sum(c("symbol", "Description") %in% colnames(featuredata)) == 2) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "featuredata - one of is missing: symbol,  Description)",
        duration = NULL,
        type = "error"
      )
    }
    if (!"Description" %in% colnames(featuredata)) {
      featuredata$"Description" <- "not given"
    }
    featuredata$symbol <- featuredata$symbol
    dataTables$featuredata <- featuredata
  }
  # if (is.null(rowData(dataTables$scEx)$symbol)){
  #
  # }

  inputFileStats$stats <- stats
  return(dataTables)
}

readMM <- function(inFile) {
  data <- Matrix::readMM(file = inFile$datapath)

  return(NULL)
}

# inFile$datapath="~/Downloads/GSE122084_RAW/GSM3855868_Salmonella_exposed_cells.txt"
# load("data/scEx.RData")
# scMat = as.matrix(assays(scEx)[[1]])
# write.file.csv(scMat, row.names=TRUE, file="data/scEx.csv" )

readCSV <- function(inFile) {
  # check.names = T will change the rownames. Since this is not enforced for the singleExperiment we shouldn't do it here either.
  con <- file(inFile$datapath, "r")
  first_line <- readLines(con, n = 1)
  close(con)
  commaCount <- length(gregexpr(",", first_line, perl = T)[[1]])
  tabCount <- length(gregexpr("\t", first_line, perl = T)[[1]])
  spaceCount <- length(gregexpr(" ", first_line, perl = T)[[1]])
  sep <- ","
  if (spaceCount > commaCount) sep <- " "
  if (commaCount > tabCount) sep <- ","
  if (tabCount > commaCount) sep <- "\t"
  data <- read.table(file = inFile$datapath, check.names = FALSE, header = TRUE, sep = sep, stringsAsFactors = F)
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/readCSV.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/SCHNAPPsDebug/readCSV.RData')
  if (colnames(data)[1] %in% c("", "rownames", "ROWNAMES", "genes")) {
    rownames(data) <- make.unique(data[, 1])
    data <- data[, -1]
  }
  exAll <- as(as.matrix(data), "dgTMatrix")
  rownames(exAll) <- rownames(data)
  colnames(exAll) <- colnames(data)
  scExnew <- SingleCellExperiment(
    assay = list(counts = exAll),
    colData = list(
      sampleNames = as.factor(rep(tools::file_path_sans_ext(inFile$name), ncol(data))),
      barcode = colnames(data)
    ),
    rowData = list(
      id = rownames(data),
      Description = rownames(data),
      symbol = rownames(data)
    )
  )
  dataTables <- list()
  dataTables$scEx <- scExnew
  dataTables$featuredata <- rowData(scExnew)

  stats <- tibble(.rows = length(inFile$datapath))
  stats$names <- inFile$name
  stats$nFeatures <- 0
  stats$nCells <- 0
  stats[1, "nFeatures"] <- nrow(data)
  stats[1, "nCells"] <- ncol(data)
  inputFileStats$stats <- stats

  return(dataTables)
}

#
# inFile$datapath="data/scEx.csv"
# load("data/scEx.RData")
# write.file.csv(colData(scEx), row.names=TRUE, file="data/scExCells.csv" )
# write.file.csv(rowData(scEx), row.names=TRUE, file="data/scExGenes.csv" )
# annFile$datapath="data/scExGenes.csv"
#' appendAnnotation
#'
#' append annotation to singleCellExperiment object
#' uses colData
appendAnnotation <- function(scEx, annFile) {
  rDat <- rowData(scEx)
  cDat <- colData(scEx)

  for (fpIdx in 1:length(annFile$datapath)) {
    data <- read.table(file = annFile$datapath[fpIdx], check.names = FALSE, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/appendAnnotation.RData", list = c(ls(), ls(envir = globalenv())))
    }
    # load(file = "~/SCHNAPPsDebug/appendAnnotation.RData")
    if (colnames(data)[1] %in% c("", "rownames", "ROWNAMES")) {
      rownames(data) <- data[, 1]
      data <- data[, -1]
    }
    # feature data
    if (all(rownames(data) %in% rownames(rDat))) {
      # if any factor is already set we need to avoid having different levles
      if (any(colnames(rDat) %in% colnames(data))) {
        commonCols <- colnames(rDat) %in% colnames(data)
        for (cCol in colnames(rDat)[commonCols]) {
          if (is(rDat[, cCol], "factor")) {
            rDat[, cCol] <- factor(data[, cCol])
          }
        }
      }
      rDat[rownames(data), colnames(data)] <- data
    }
    # symbol used as row name
    if (any(rownames(data) %in% rDat$symbol)) {
      rowColfound <- FALSE
      nrComm <- c()
      for (cIdx in colnames(data)) {
        ci <- which(colnames(data) == cIdx)
        nrComm[ci] <- sum(rownames(rDat) %in% data[, cIdx])
      }
      if (max(nrComm) > (nrow(data) / 2)) {
        rowColfound <- which(nrComm == max(nrComm))
      }
      if (rowColfound == FALSE) {
        cat(file = stderr(), paste("couldn't find column with rownames when symbols are used as rownames for ", annFile$name, "\n"))
        break()
      } else {
        whichRowNames <- which(data[, cIdx] %in% rownames(rDat))
        rownames(data)[whichRowNames] <- as.character(data[whichRowNames, rowColfound])
        data <- data[, -rowColfound]
      }
      # if any factor is already set we need to avoid having different levles
      if (any(colnames(rDat) %in% colnames(data))) {
        commonCols <- colnames(rDat) %in% colnames(data)
        for (cCol in colnames(rDat)[commonCols]) {
          if (is(rDat[, cCol], "factor")) {
            rDat[, cCol] <- factor(data[, cCol])
          }
        }
      }
      rDat[rownames(data), colnames(data)] <- data
    }

    # colData
    if (all(rownames(data) %in% rownames(cDat))) {
      # if any factor is already set we need to avoid having different levles
      if (any(colnames(cDat) %in% colnames(data))) {
        commonCols <- colnames(cDat) %in% colnames(data)
        for (cCol in colnames(cDat)[commonCols]) {
          if (is(cDat[, cCol], "factor")) {
            cDat[, cCol] <- factor(data[, cCol])
          }
        }
      }
      # a vector with less than 20 levels (unique values) will be converted in a factor
      for (cn in colnames(data)) {
        if (is(data[, cn], "factor")) next()
        lv <- levels(factor(data[, cn]))
        if (length(lv) <= 20) {
          data[, cn] <- factor(data[, cn])
        }
      }
      cDat[rownames(data), colnames(data)] <- data
    }
  }
  cDat$sampleNames <- factor(cDat$sampleNames)
  colData(scEx) <- cDat
  rowData(scEx) <- rDat
  return(scEx)
}

# inputData ----
# load RData file with singlecellExperiment object
# internal, should not be used by plug-ins
inputData <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "inputData started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "inputData")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "inputData")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("inputData", id = "inputData", duration = NULL)
  }

  inFile <- input$file1
  annFile <- input$annoFile
  sampleCells <- input$sampleInput
  subsampleNum <- input$subsampleNum

  if (is.null(inFile)) {
    if (DEBUG) {
      cat(file = stderr(), "inputData: NULL\n")
    }
    return(NULL)
  }
  if (DEBUG) cat(file = stderr(), paste("inFile.", inFile$datapath[1], "\n"))
  if (!file.exists(inFile$datapath[1])) {
    if (DEBUG) {
      cat(file = stderr(), "inputData: ", inFile$datapath[1], " doesn't exist\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/inputData.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/SCHNAPPsDebug/inputData.RData')

  # inFile$datapath = "data/sessionData.RData"
  # annFile$datapath = "data/selectedCells.csv"
  # TODO either multiple files with rdata/rds or single file of csv/other
  fpExtension <- tools::file_ext(inFile$datapath[1])
  if (toupper(fpExtension) %in% c("RDATA", "RDS")) {
    retVal <- inputDataFunc(inFile)
  } else {
    retVal <- readCSV(inFile)
  }

  if (is.null(retVal)) {
    return(NULL)
  }

  if (is.null(retVal[["scEx"]])) {
    return(NULL)
  }

  if (!is.null(annFile)) {
    if (file.exists(annFile$datapath[1])) {
      retVal$scEx <- appendAnnotation(scEx = retVal$scEx, annFile = annFile)
    }
  }

  sampNames <- levels(colData(retVal$scEx)$sampleNames)
  isolate({
    sampleCols$colPal <- allowedColors[seq_along(sampNames)]
    names(sampleCols$colPal) <- sampNames
  })
  inputFile$inFile <- paste(inFile$name, collapse = ", ")
  inputFile$annFile <- paste(annFile$name, collapse = ", ")

  if (sampleCells) {
    samp <- base::sample(
      x = ncol(assays(retVal$scEx)[["counts"]]),
      size = min(subsampleNum, ncol(assays(retVal$scEx)[["counts"]])),
      replace = FALSE
    )
    retVal$scEx <- retVal$scEx[, samp]
  }
  exportTestValues(inputData = {
    list(
      assays(retVal$scEx)[["counts"]],
      rowData(retVal$scEx),
      colData(retVal$scEx)
    )
  })
  return(retVal)
})

medianENSGfunc <- function(scEx) {
  if (DEBUG) {
    cat(file = stderr(), "medianENSGfunc started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "medianENSGfunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "medianENSGfunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("medianENSGfunc", id = "medianENSGfunc", duration = NULL)
  }

  geneC <- Matrix::colSums(scEx > 0, na.rm = TRUE)
  return(median(t(geneC)))
}

# medianENSG ----
medianENSG <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "medianENSG started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "medianENSG")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "medianENSG")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("medianENSG", id = "medianENSG", duration = NULL)
  }

  scEx_log <- scEx_log()
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "medianENSG:NULL\n")
    }
    return(0)
  }
  scEx_log <- assays(scEx_log)[[1]]
  if (ncol(scEx_log) <= 1 | nrow(scEx_log) < 1) {
    return(0)
  }
  retVal <- medianENSGfunc(scEx_log)

  exportTestValues(medianENSG = {
    retVal
  })
  return(retVal)
})

medianUMIfunc <- function(scEx) {
  if (DEBUG) {
    cat(file = stderr(), "medianUMIfunc started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "medianUMIfunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "medianUMIfunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("medianUMIfunc", id = "medianUMIfunc", duration = NULL)
  }

  umiC <- Matrix::colSums(scEx, na.rm = TRUE)

  return(median(t(umiC)))
}

# medianUMI ----
medianUMI <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "medianUMI started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "medianUMI")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "medianUMI")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("medianUMI", id = "medianUMI", duration = NULL)
  }

  scEx <- scEx()
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "medianUMI:NULL\n")
    }
    return(0)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/medianUMI.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/SCHNAPPsDebug/medianUMI.RData')
  scEx <- assays(scEx)[["counts"]]
  retVal <- medianUMIfunc(scEx)

  exportTestValues(medianUMI = {
    retVal
  })
  return(retVal)
})




# for now we don't have a way to specifically select cells
# we could cluster numbers or the like
# internal, should not be used by plug-ins
useCellsFunc <-
  function(dataTables,
             geneNames,
             rmCells,
             rmPattern,
             keepCells,
             cellKeepOnly,
             geneNamesNEG) {
    if (DEBUG) {
      cat(file = stderr(), "useCellsFunc started.\n")
    }
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "useCellsFunc")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "useCellsFunc")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("useCellsFunc", id = "useCellsFunc", duration = NULL)
    }

    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/useCellsFunc.RData", list = c(ls()))
    }
    # load(file='~/SCHNAPPsDebug/useCellsFunc.Rdata')
    goodCols <- rep(TRUE, dim(dataTables$scEx)[2])
    scEx <- assays(dataTables$scEx)[[1]]
    #### start: cells with genes expressed
    # take only cells where these genes are expressed with at least one read
    genesin <- toupper(geneNames)
    genesin <- gsub(" ", "", genesin, fixed = TRUE)
    genesin <- strsplit(genesin, ",")
    genesin <- genesin[[1]]

    cellKeep <- toupper(keepCells)
    cellKeep <- gsub(" ", "", cellKeep, fixed = TRUE)
    cellKeep <- strsplit(cellKeep, ",")
    cellKeep <- cellKeep[[1]]

    cellKeepOnly <- toupper(cellKeepOnly)
    cellKeepOnly <- gsub(" ", "", cellKeepOnly, fixed = TRUE)
    cellKeepOnly <- strsplit(cellKeepOnly, ",")
    cellKeepOnly <- cellKeepOnly[[1]]

    # specifically remove cells
    if (nchar(rmCells) > 0) {
      cellsRM <- toupper(rmCells)
      cellsRM <- gsub(" ", "", cellsRM, fixed = TRUE)
      cellsRM <- strsplit(cellsRM, ",")
      cellsRM <- cellsRM[[1]]
      goodCols[which(toupper(colnames(dataTables$scEx)) %in% cellsRM)] <-
        FALSE
    }

    # remove cells by pattern
    if (nchar(rmPattern) > 0) {
      goodCols[grepl(rmPattern, colnames(dataTables$scEx), ignore.case = TRUE)] <- FALSE
    }

    if (!length(cellKeep) == 0) {
      ids <- which(toupper(colnames(dataTables$scEx)) %in% cellKeep)
      goodCols[ids] <- TRUE
    }

    # genes that have to be expressed at least in one of them.
    selCols <- rep(FALSE, length(goodCols))
    if (!length(genesin) == 0) {
      ids <- which(toupper(dataTables$featuredata$symbol) %in% genesin)
      if (length(ids) == 1) {
        selCols <- scEx[ids, ] > 0
      } else if (length(ids) == 0) {
        showNotification(
          "not enough cells, check gene names for min coverage",
          type = "warning",
          duration = NULL
        )
        return(NULL)
      } else {
        selCols <- Matrix::colSums(scEx[ids, ]) > 0
      }
      goodCols <- goodCols & selCols
    }

    # genes that should NOT be expressed.
    genesin <- toupper(geneNamesNEG)
    genesin <- gsub(" ", "", genesin, fixed = TRUE)
    genesin <- strsplit(genesin, ",")
    genesin <- genesin[[1]]

    selCols <- rep(FALSE, length(goodCols))
    if (!length(genesin) == 0) {
      ids <- which(toupper(dataTables$featuredata$symbol) %in% genesin)
      if (length(ids) == 1) {
        selCols <- scEx[ids, ] > 0
      } else if (length(ids) == 0) {
        showNotification(
          "not enough cells for NonExpGenes, check gene names for min coverage",
          type = "error",
          duration = NULL
        )
        return(NULL)
      } else {
        selCols <- Matrix::colSums(scEx[ids, ]) > 0
      }
      # now the cells that should be removed are set to T
      # we inverse this and then keep only the ones that are not found
      selCols <- !selCols
      goodCols <- goodCols & selCols
    }

    if (!length(cellKeepOnly) == 0) {
      goodCols[c(1:length(goodCols))] <- FALSE
      ids <-
        which(toupper(colnames(dataTables$scEx)) %in% cellKeepOnly)
      goodCols[ids] <- TRUE
    }

    #### end: cells with genes expressed
    return(goodCols)
  }

# useCells ----
# works on cells only
# internal, should not be used by plug-ins
useCells <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "useCells started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "useCells")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "useCells")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("useCells", id = "useCells", duration = NULL)
  }

  dataTables <- inputData()
  cellSelectionValues <- cellSelectionValues()
  geneNames <- cellSelectionValues$minExpGenes
  geneNamesNEG <- cellSelectionValues$minNonExpGenes
  rmCells <- cellSelectionValues$cellsFiltersOut
  rmPattern <- cellSelectionValues$cellPatternRM
  keepCells <- cellSelectionValues$cellKeep
  cellKeepOnly <- cellSelectionValues$cellKeepOnly
  # useGenes = isolate(useGenes())
  if (!exists("dataTables") || is.null(dataTables)) {
    if (DEBUG) {
      cat(file = stderr(), "useCells:NULL\n")
    }
    return(NULL)
  }

  retVal <- useCellsFunc(
    dataTables,
    geneNames,
    rmCells,
    rmPattern,
    keepCells,
    cellKeepOnly,
    geneNamesNEG
  )

  exportTestValues(useCells = {
    retVal
  })
  return(retVal)
})


useGenesFunc <-
  function(dataTables,
             ipIDs,
             # regular expression of genes to be removed
             geneListSelection,
             genesKeep,
             geneLists) {
    if (DEBUG) {
      cat(file = stderr(), "useGenesFunc started.\n")
    }
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "useGenesFunc")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "useGenesFunc")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("useGenesFunc", id = "useGenesFunc", duration = NULL)
    }

    gList <-
      geneLists # global variable, assigning it locally ensures that it will be saved
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/useGenesFunc.Rdata", list = c(ls(), ls(envir = globalenv())))
    }
    # load(file='~/SCHNAPPsDebug/useGenesFunc.Rdata')
    # regular expression with gene names to be removed
    if (nchar(ipIDs) > 0) {
      keepIDs <- !grepl(ipIDs, dataTables$featuredata$symbol, ignore.case = TRUE)
    } else {
      keepIDs <- rep(TRUE, dim(dataTables$scEx)[1])
    }
    # explicit list of genes to keep
    genesKeep <- toupper(genesKeep)
    genesKeep <- gsub(" ", "", genesKeep, fixed = TRUE)
    genesKeep <- strsplit(genesKeep, ",")
    genesKeep <- genesKeep[[1]]
    keepGeneIds <-
      which(toupper(dataTables$featuredata$symbol) %in% genesKeep)

    # dataTables$featuredata$symbol[keepIDs]
    # gene groups to be included
    # work on the geneList tree
    if (!is.null(geneListSelection)) {
      selectedgeneList <- get_selected(geneListSelection)
      if (length(selectedgeneList) > 0) {
        selGenes <- c()
        for (sIdx in 1:length(selectedgeneList)) {
          print(sIdx)
          att <- attr(selectedgeneList[[sIdx]], "ancestry")
          if (length(att) > 0) {
            selGenes <- c(selGenes, gList[[att]][[selectedgeneList[[sIdx]]]])
          }
        }
        selGenes <- unique(selGenes)
        keepIDs <-
          (rownames(dataTables$scEx) %in% selGenes) & keepIDs
      }
    }

    keepIDs[keepGeneIds] <- TRUE
    return(keepIDs)
  }

# selected genes table ----
gsSelectedGenesTable <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "output$selectedGenesTable\n")
  }
  dataTables <- inputData()
  useGenes <- useGenes()
  useCells <- useCells()
  minGenes <- input$minGenesGS

  if (is.null(dataTables) | is.null(useGenes) | is.null(useCells)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(
      file = "~/SCHNAPPsDebug/selectedGenesTable.RData",
      list = c("normaliztionParameters", ls(), ls(envir = globalenv()))
    )
  }
  # load("~/SCHNAPPsDebug/selectedGenesTable.RData")

  scEx <- assays(dataTables$scEx)[[1]]
  fd <- rowData(dataTables$scEx)
  dt <- fd[useGenes, ]
  dt$rowSums <- Matrix::rowSums(scEx[useGenes, useCells])
  dt$rowSamples <- Matrix::rowSums(scEx[useGenes, useCells] > 0)
  # get the order of the frist two columns correct
  firstCol <- which(colnames(dt) == "symbol")
  firstCol <- c(firstCol, which(colnames(dt) == "Description"))
  # those we created so we know they are there
  firstCol <- firstCol <- c(firstCol, which(colnames(dt) %in% c("rowSums", "rowSamples")))
  colOrder <- c(firstCol, (1:ncol(dt))[-firstCol])
  dt <- dt[, colOrder]
  dt <- dt[dt$rowSums >= minGenes, ]
  exportTestValues(selectedGenesTable = {
    as.data.frame(dt)
  })
  # DT::datatable(as.data.frame(dt),
  #               options = list(scrollX = TRUE))
  rownames(dt) <- dt$symbol
  as.data.frame(dt)
})


# removed genes table ----
gsRMGenesTable <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "output$removedGenesTable\n")
  }
  dataTables <- inputData()
  useGenes <- useGenes()
  useCells <- useCells()
  minGenes <- input$minGenesGS

  if (is.null(dataTables) | is.null(useGenes) | is.null(useCells)) {
    return(NULL)
  }
  useGenes <- !useGenes

  if (.schnappsEnv$DEBUGSAVE) {
    save(
      file = "~/SCHNAPPsDebug/removedGenesTable.RData",
      list = c("normaliztionParameters", ls(), ls(envir = globalenv()))
    )
  }
  # load("~/SCHNAPPsDebug/removedGenesTable.RData")
  scEx <- assays(dataTables$scEx)[[1]]
  fd <- rowData(dataTables$scEx)
  dt <- fd[useGenes, c("symbol", "Description")]
  dt$rowSums <- Matrix::rowSums(scEx[useGenes, useCells])
  dt$rowSamples <- Matrix::rowSums(scEx[useGenes, useCells] > 0)

  # dt <- dt[dt$rowSums < minGenes, ]
  exportTestValues(removedGenesTable = {
    as.data.frame(dt)
  })
  if (nrow(dt) == 0) {
    return(NULL)
  }
  rownames(dt) <- dt$symbol
  as.data.frame(dt)
  # DT::datatable(as.data.frame(dt))
})

# before gene filtering -----
beforeFilterCounts <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "beforeFilterCounts started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "beforeFilterCounts")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "beforeFilterCounts")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("beforeFilterCounts", id = "beforeFilterCounts", duration = NULL)
  }

  dataTables <- inputData()
  geneSelectionValues <- geneSelectionValues()
  # ipIDs <- input$selectIds # regular expression of genes to be removed
  ipIDs <- geneSelectionValues$selectIds
  if (!exists("dataTables") |
    is.null(dataTables) |
    length(dataTables$featuredata$symbol) == 0) {
    if (DEBUG) {
      cat(file = stderr(), "beforeFilterCounts: NULL\n")
    }
    return(NULL)
  }

  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/beforeFilterCounts.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/beforeFilterCounts.RData")

  geneIDs <- NULL
  if (nchar(ipIDs) > 0) {
    geneIDs <- grepl(ipIDs, dataTables$featuredata$symbol, ignore.case = TRUE)
  }
  if (is.null(geneIDs)) {
    return(rep(0, nrow(dataTables$featuredata)))
  }
  retVal <- Matrix::colSums(assays(dataTables$scEx)[["counts"]][geneIDs, ])

  exportTestValues(beforeFilterCounts = {
    retVal
  })
  return(retVal)
})

# useGenes ----
# collects information from all places where genes being removed or specified
useGenes <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "useGenes started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "useGenes")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "useGenes")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("useGenes", id = "useGenes", duration = NULL)
  }

  dataTables <- inputData()
  geneSelectionValues <- geneSelectionValues()
  # ipIDs <-
  #   input$selectIds # regular expression of genes to be removed
  ipIDs <- geneSelectionValues$selectIds
  # genesKeep <- input$genesKeep
  # geneListSelection <- input$geneListSelection
  genesKeep <- geneSelectionValues$genesKeep
  geneListSelection <- geneSelectionValues$geneListSelection

  if (!exists("dataTables") |
    is.null(dataTables) |
    length(dataTables$featuredata$symbol) == 0) {
    if (DEBUG) {
      cat(file = stderr(), "useGenes: NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("which genes to use", id = "useGenes", duration = NULL)
  }

  retVal <-
    useGenesFunc(dataTables, ipIDs, geneListSelection, genesKeep, geneLists)

  exportTestValues(useGenes = {
    retVal
  })
  return(retVal)
})


# will be called recursively to ensure that nothing changes when cells/genes are changing.
scExFunc <-
  function(scExOrg,
             useCells,
             useGenes,
             minGene,
             minG,
             maxG) {
    if (DEBUG) {
      cat(file = stderr(), "scExFunc started.\n")
    }
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "scExFunc")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "scExFunc")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("scExFunc", id = "scExFunc", duration = NULL)
    }
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scExFunc1")
      removeNotification(id = "scExFunc2")
    }

    # change names to be hopefully a bit more clear
    changed <- FALSE # trace if something changed
    keepGenes <- useGenes
    keepCells <- useCells
    scEx <- assays(scExOrg)[[1]]

    # overall gene expression Min
    if (!is.null(minGene)) {
      selGenes <- Matrix::rowSums(scEx[, keepCells]) >= minGene
      selGenes <- keepGenes & selGenes
      if (!all(selGenes == keepGenes)) {
        keepGenes <- selGenes
        changed <- TRUE
      }
    }

    # min reads per cell
    if (!is.null(minG)) {
      selCols <- Matrix::colSums(scEx[keepGenes, ], na.rm = FALSE) > minG
      selCols[is.na(selCols)] <- FALSE
      selCols <- keepCells & selCols
      if (!all(selCols == keepCells)) {
        keepCells <- selCols
        changed <- TRUE
      }
    }

    # max reads per cell
    if (!is.null(maxG)) {
      selCols <- Matrix::colSums(scEx[keepGenes, ], na.rm = FALSE) <= maxG
      selCols[is.na(selCols)] <- FALSE
      selCols <- selCols & keepCells
      if (!all(selCols == keepCells)) {
        changed <- TRUE
        keepCells <- selCols
      }
    }

    if (sum(keepCells) == 0) {
      showNotification(
        "not enough cells left",
        type = "warning",
        id = "scExFunc1",
        duration = NULL
      )
      return(NULL)
    }
    if (sum(keepGenes) == 0) {
      showNotification(
        "not enough genes left",
        type = "warning",
        id = "scExFunc1",
        duration = NULL
      )
      return(NULL)
    }
    # cells with same expression pattern
    # save(file = "~/SCHNAPPsDebug/scEx3.RData", list = c(ls(), ls(envir = globalenv())))
    # load(file = "~/SCHNAPPsDebug/scEx3.RData")
    # which(duplicated(as.matrix(assays(scExOrg[keepGenes, keepCells])[[1]])))
    # which(duplicated(as.matrix(assays(scEx)[[1]])))
    #
    # sum(duplicated(t(as.matrix(assays(scExOrg[keepGenes, keepCells])[[1]]))))
    # if something changed, check that it doesn't change again
    scExNew <- scExOrg[keepGenes, keepCells]
    if (changed) {
      scExNew <-
        scExFunc(scExOrg[keepGenes, keepCells], useCells[keepCells], useGenes[keepGenes], minGene, minG, maxG)
      if (is.null(scExNew)) {
        return(NULL)
      }
    }

    pD <- colData(scExNew)
    for (colN in colnames(pD)) {
      if (colN == "barcode") {
        next()
      }
      if (is(pD[, colN], "character")) {
        pD[, colN] <- factor(as.character(pD[, colN]))
      }
    }
    colData(scExNew) <- pD

    return(scExNew)
  }

# scEx ----
# apply filters that depend on genes & cells
# it is here that useCells and useGenes are combined and applied to select for
scEx <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "scEx started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scEx")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scEx")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scEx", id = "scEx", duration = NULL)
  }

  dataTables <- inputData()
  useCells <- useCells()
  useGenes <- useGenes()
  cellSelectionValues <- cellSelectionValues()
  minGene <- cellSelectionValues$minGenesGS # min number of reads per gene
  minG <- cellSelectionValues$minGenes # min number of reads per cell
  maxG <- cellSelectionValues$maxGenes # max number of reads per cell

  if (!exists("dataTables") |
    is.null(dataTables) | is.null(useGenes) | is.null(useCells)) {
    if (DEBUG) {
      cat(file = stderr(), "scEx: NULL\n")
    }
    return(NULL)
  }

  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scEx.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/scEx.RData")

  # TODO:??? should be keep the cell names from the projections?
  # here we are just reinitializing.
  retVal <- scExFunc(
    scExOrg = dataTables$scEx,
    useCells = useCells,
    useGenes = useGenes,
    minGene = minGene,
    minG = minG,
    maxG = maxG
  )


  exportTestValues(scEx = {
    list(rowData(retVal), colData(retVal))
  })
  return(retVal)
})

# rawNormalization ----
rawNormalization <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "rawNormalization started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "rawNormalization")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "rawNormalization")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("rawNormalization", id = "rawNormalization", duration = NULL)
  }

  scEx <- scEx()
  names(assays(scEx)) <- "logcounts"

  exportTestValues(rawNormalization = {
    str(scEx)
  })
  return(scEx)
})

# scEx_log ----
scEx_log <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "scEx_log started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scEx_log")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scEx_log")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scEx_log", id = "scEx_log", duration = NULL)
  }

  scEx <- scEx()
  normMethod <- input$normalizationRadioButton
  disableNorm <- input$disablescEx_log
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "scEx_log:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scEx_log.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/scEx_log.RData")

  if (disableNorm) {
    return(NULL)
  }
  scEx_log <- do.call(normMethod, args = list())

  exportTestValues(scEx_log = {
    assays(scEx_log)["logcounts"]
  })
  return(scEx_log)
})


scEx_log_sha <- reactive({
  scEx_log <- scEx_log()
  require(digest)
  if (is.null(scEx_log)) {
    return(NULL)
  }
  return(sha1(as.matrix(assays(scEx_log)[[1]])))
})
# scExLogMatrixDisplay ----
# scExLog matrix with symbol as first column
# TODO
# we should probably just rename the rows and then have an option to tableSelectionServer that shows (or not) rownames
scExLogMatrixDisplay <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "scExLogMatrixDisplay started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scExLogMatrixDisplay")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scExLogMatrixDisplay")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scExLogMatrixDisplay",
      id = "scExLogMatrixDisplay",
      duration = NULL
    )
  }

  # dataTables = inputData()
  # useCells = useCells()
  # useGenes = useGenes()
  scEx <- scEx()
  scEx_log <- scEx_log()
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "scExLogMatrixDisplay:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scExLogMatrixDisplay.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/scExLogMatrixDisplay.RData")

  # TODO
  if (dim(scEx)[1] > 20000) {

  }
  retVal <-
    data.frame(
      symbol = make.names(rowData(scEx)$symbol, unique = TRUE),
      stringsAsFactors = FALSE
    )
  if (!is.null(scEx_log)) {
    retVal <- cbind(
      retVal,
      as.matrix(assays(scEx_log)[[1]])
    )
  }
  rownames(retVal) <-
    make.names(rowData(scEx)$symbol, unique = TRUE)

  return(retVal)
})

pcaFunc <- function(scEx_log, rank, center, scale, pcaGenes, featureData, pcaN) {
  if (DEBUG) {
    cat(file = stderr(), "pcaFunc started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "pcaFunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "pcaFunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("pcaFunc", id = "pcaFunc", duration = NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "pcawarning")
  }

  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/pcaFunc.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/pcaFunc.RData")
  genesin <- geneName2Index(pcaGenes, featureData)
  if (is.null(genesin) || length(genesin) == 0) {
    genesin <- rownames(scEx_log)
  }
  if (length(genesin) < rank) {
    rank <- length(genesin)
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        paste("Problem with PCA: Too many PCs requested, setting to ", rank),
        type = "warning",
        id = "pcawarning",
        duration = NULL
      )
    }
  }

  withWarnings <- function(expr) {
    wHandler <- function(w) {
      if (DEBUG) {
        cat(file = stderr(), paste("runPCA created a warning:", w, "\n"))
      }
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification(
          "Problem with PCA?",
          type = "warning",
          id = "pcawarning",
          duration = NULL
        )
      }
      invokeRestart("muffleWarning")
    }
    eHandler <- function(e) {
      cat(file = stderr(), paste("error in PCA:", e))
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification(
          paste("Problem with PCA, probably not enough cells?", e),
          type = "warning",
          id = "pcawarning",
          duration = NULL
        )
      }
      cat(file = stderr(), "PCA FAILED!!!\n")
      return(NULL)
    }
    val <- withCallingHandlers(expr, warning = wHandler, error = eHandler)
    return(val)
  }

  scaterPCA <- withWarnings({
    # not sure, but this works on another with dgTMatrix
    if (is(assays(scEx_log)[["logcounts"]], "dgTMatrix")) {
      assays(scEx_log)[["logcounts"]] <-
        as(assays(scEx_log)[["logcounts"]], "dgCMatrix")
    }
    BiocSingular::runPCA(
      # t(assays(scEx_log)[["logcounts"]]),
      scEx_log[genesin, ],
      ncomponents = rank,
      ntop = pcaN,
      exprs_values = "logcounts",
      # rank = rank,
      #  center = center,
      scale = scale
      # ,
      # method = "irlba",
      # BPPARAM = bpparam(),
      # BPPARAM = SnowParam(workers = 3, type = "SOCK),
      # BPPARAM = MulticoreParam(
      #   workers = ifelse(detectCores()>1, detectCores()-1, 1))
      # BSPARAM = IrlbaParam()
    )
  })

  if (is.null(scaterPCA)) {
    return(NULL)
  }
  # pca = reducedDim(scaterPCA, "PCA")
  # attr(pca,"percentVar")
  #
  # rownames(scaterPCA$x) = colnames(scEx_log)
  return(list(
    # x =  scaterPCA$x,
    x = SingleCellExperiment::reducedDim(scaterPCA, "PCA"),
    # var_pcs = scaterPCA$sdev
    var_pcs = attr(
      SingleCellExperiment::reducedDim(scaterPCA, "PCA"),
      "percentVar"
    )
  ))
}

# pca ----
pca <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "pca started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "pca")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "pca")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("pca", id = "pca", duration = NULL)
  }
  rank <- input$pcaRank
  pcaN <- input$pcaN
  center <- input$pcaCenter
  scale <- input$pcaScale
  pcaGenes <- input$genes4PCA
  scEx_log <- scEx_log()

  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "pca:NULL\n")
    }
    return(NULL)
  }

  if (is.null(rank)) rank <- 10
  if (is.null(center)) centr <- TRUE
  if (is.null(scale)) scale <- FALSE
  featureData <- rowData(scEx_log)
  retVal <- pcaFunc(scEx_log, rank, center, scale, pcaGenes, featureData, pcaN)

  printTimeEnd(start.time, "pca")
  exportTestValues(pca = {
    retVal
  })
  return(retVal)
})

# scranCluster -----
scranCluster <- function(pca,
                         scEx,
                         scEx_log,
                         seed,
                         clusterSource,
                         geneSelectionClustering = "",
                         minClusterSize = 2,
                         clusterMethod = "igraph",
                         featureData,
                         useRanks) {
  if (DEBUG) {
    cat(file = stderr(), "scranCluster started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scranCluster")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scranCluster")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scranCluster", id = "scranCluster", duration = NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scranCluster.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/scranCluster.RData")

  set.seed(seed)
  geneid <- geneName2Index(geneSelectionClustering, featureData)

  params <- list(
    # min.size: An integer scalar specifying the minimum size of each cluster.
    min.size = minClusterSize,
    method = clusterMethod,
    use.ranks = useRanks,
    # d: An integer scalar specifying the number of principal components to retain.
    # If NULL and use.ranks=TRUE, this defaults to 50. If use.rank=FALSE, the number of PCs
    # is chosen by denoisePCA.
    # If NA, no dimensionality reduction is performed and the gene expression values
    # (or their rank equivalents) are directly used in clustering.
    d = ncol(pca$x),
    # When use.ranks=TRUE, the function will also filter out genes with average counts
    # (as defined by calcAverage) below min.mean. This removes low-abundance genes with many
    # tied ranks, especially due to zeros, which may reduce the precision of the clustering.
    # We suggest setting min.mean to 1 for read count data and 0.1 for UMI data.
    min.mean = 0.1,
    BPPARAM = MulticoreParam(
      workers = ifelse(detectCores() > 1, detectCores() - 1, 1)
    ),
    block.BPPARAM = MulticoreParam(
      workers = ifelse(detectCores() > 1, detectCores() - 1, 1)
    )
  )
  switch(clusterSource,
    "PCA" = {
      reducedDims(scEx) <- SimpleList(PCA = pca$x)
      params$assay.type <- "counts"
      params$x <- scEx
    },
    "logcounts" = {
      params$x <- scEx_log
      reducedDims(scEx_log) <- SimpleList(PCA = pca$x)
      params$assay.type <- "logcounts"
      if (length(geneid) > 0) {
        params$subset.row <- geneid
      }
      use.ranks <- NA
    },
    "counts" = {
      params$x <- scEx
      reducedDims(scEx) <- SimpleList(PCA = pca$x)
      params$assay.type <- "counts"
      if (length(geneid) > 0) {
        params$subset.row <- geneid
      }
      use.ranks <- NA
    }
  )
  require(scran)
  retVal <- tryCatch({
    suppressMessages(do.call("quickCluster", params))
  },
  error = function(e) {
    cat(file = stderr(), paste("\nProblem with clustering:\n\n", as.character(e), "\n\n"))
    if (grepl("rank variances of zero", as.character(e))) {
      if (useRanks) {
        if (!is.null(getDefaultReactiveDomain())) {
          showNotification("quickClusterError", id = "quickClusterError", type = "error", duration = NULL)
        }
        cat(file = stderr(), "\nNot using ranks due to identical ranks\n")
        params$use.ranks <- F
        return(do.call("quickCluster", params))
      }
      cat(file = stderr(), "\nRemove duplicated cells\n")
    }
    return(NULL)
  },
  warning = function(e) {
    if (DEBUG) cat(file = stderr(), paste("\nclustering produced Warning:\n", e, "\n"))
    return(suppressMessages(do.call("quickCluster", params)))
  }
  )
  if ("barcode" %in% colnames(colData(scEx_log))) {
    barCode <- colData(scEx_log)$barcode
  } else {
    barCode <- rownames(colData(scEx_log))
  }

  retVal <- data.frame(
    Barcode = barCode,
    Cluster = retVal
  )
  rownames(retVal) <- rownames(colData(scEx_log))
  return(retVal)
}




scran_Cluster <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "scran_Cluster started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scran_Cluster")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scran_Cluster")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scran_Cluster", id = "scran_Cluster", duration = NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "dbClusterError")
  }

  pca <- pca()
  scEx <- scEx()
  scEx_log <- scEx_log()
  seed <- input$seed
  # kNr <- input$kNr
  useRanks <- input$useRanks

  clusterSource <- clusterMethodReact$clusterSource
  geneSelectionClustering <- input$geneSelectionClustering
  minClusterSize <- input$minClusterSize
  clusterMethod <- clusterMethodReact$clusterMethod

  if (is.null(pca) | is.null(scEx_log) | is.na(minClusterSize)) {
    if (DEBUG) {
      cat(file = stderr(), "scran_Cluster:NULL\n")
    }
    return(NULL)
  }
  # if (.schnappsEnv$DEBUGSAVE) {
  #   save(file = "~/SCHNAPPsDebug/scran_Cluster.RData", list = c(ls(), ls(envir = globalenv())))
  # }
  # load(file="~/SCHNAPPsDebug/scran_Cluster.RData")

  featureData <- rowData(scEx_log)

  if (is.null(seed)) {
    seed <- 1
  }
  retVal <- scranCluster(
    pca,
    scEx,
    scEx_log,
    seed,
    clusterSource,
    geneSelectionClustering,
    minClusterSize,
    clusterMethod,
    featureData,
    useRanks
  )
  if (is.null(retVal)) {
    showNotification(
      paste("error: clustering didn't produce a result"),
      type = "error",
      id = "dbClusterError",
      duration = NULL
    )
    exit()
  }

  exportTestValues(scran_Cluster = {
    retVal
  })
  return(retVal)
})

# dbCluster ----
dbCluster <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "dbCluster started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "dbCluster")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "dbCluster")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("dbCluster", id = "dbCluster", duration = NULL)
  }

  # kNr <- input$kNr
  clustering <- scran_Cluster()

  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/dbCluster.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/dbCluster.RData")

  if (is.null(clustering)) {
    if (DEBUG) {
      cat(file = stderr(), "dbCluster: NULL\n")
    }
    return(NULL)
  }

  dbCluster <- clustering$Cluster

  # cluster colors
  inCols <- list()
  lev <- levels(dbCluster)

  inCols <- allowedColors[1:length(lev)]
  names(inCols) <- lev

  clusterCols$colPal <- unlist(inCols)

  exportTestValues(dbCluster = {
    dbCluster
  })
  return(dbCluster)
})


clusterMethodReact <- reactiveValues(
  clusterMethod = "igraph",
  clusterSource = "counts"
)

# collect copied/renamed projections
projectionsTable <- reactiveValues(
  newProjections = data.frame()
)

# projections ----
#' projections
#' each column is of length of number of cells
#' if factor than it is categorical and can be cluster number of sample etc
#' if numeric can be projection
#' projections is a reactive and cannot be used in reports. Reports have to organize
#' themselves as it is done here with tsne.data.
projections <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "projections started.\n")
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
  pca <- pca()
  prjs <- sessionProjections$prjs
  newPrjs <- projectionsTable$newProjections
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/projections.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/projections.RData"); DEBUGSAVE=FALSE
  if (!exists("scEx") |
    is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "sampleInfo: NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/projections.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/projections.RData"); DEBUGSAVE=FALSE



  # todo colData() now returns a s4 object of class DataFrame
  # not sure what else is effected...
  pd <- as.data.frame(colData(scEx))
  if (ncol(pd) < 2) {
    cat(file = stderr(), "phenoData for scEx has less than 2 columns\n")
    return(NULL)
  }
  projections <- pd

    if (!is.null(pca)) {
    projections <- cbind(projections, pca)
  }

  withProgress(message = "Performing projections", value = 0, {
    n <- length(.schnappsEnv$projectionFunctions)
    iter <- 1
    for (proj in .schnappsEnv$projectionFunctions) {
      start.time1 <- Sys.time()
      incProgress(1 / n, detail = paste("Creating ", proj[1]))
      if (DEBUG) {
        cat(file = stderr(), paste("calculation projection:  ", proj[1], "\n"))
      }
      if (DEBUG) cat(file = stderr(), paste("projection: ", proj[2], "\n"))
      assign("tmp", eval(parse(text = paste0(proj[2], "()"))))
      if (.schnappsEnv$DEBUGSAVE) {
        save(
          file = paste0("~/SCHNAPPsDebug/projections.", iter, ".RData"),
          list = c("tmp")
        )
        iter <- iter + 1
      }
      # load(file="~/SCHNAPPsDebug/projections.1.RData")
      # browser()
      # TODO here, dbCluster is probably overwritten and appended a ".1"
      if (is(tmp, "data.frame")) {
        cn <- make.names(c(colnames(projections), colnames(tmp)))
      } else {
        cn <- make.names(c(colnames(projections), make.names(proj[1])))
      }
      if (length(tmp) == 0) {
        next()
      }
      if (ncol(projections) == 0) {
        # never happening because we set pca first
        projections <- data.frame(tmp = tmp)
      } else {
        if (nrow(projections) == length(tmp)) {
          projections <- cbind(projections, tmp)
        } else if (nrow(projections) == nrow(tmp)) {
          projections <- cbind(projections, tmp)
        } else {
          stop("error: ", proj[1], "didn't produce a result")
        }
      }

      colnames(projections) <- cn
      if (DEBUG) cat(file = stderr(), paste("colnames ", paste0(colnames(projections), collapse = " "), "\n"))
      if (DEBUG) cat(file = stderr(), paste("observe this: ", proj[2], "\n"))
      # observe(proj[2], quoted = TRUE)
    }
  })
  # add a column for gene specific information that will be filled/updated on demand
  # projections$UmiCountPerGenes <- 0
  # projections$UmiCountPerGenes2 <- 0
  for (pdIdx in colnames(pd)) {
    if (!pdIdx %in% colnames(projections)) {
      projections[, pdIdx] <- pd[, pdIdx]
    }
  }

  if (ncol(prjs) > 0 & nrow(prjs) == nrow(projections)) {
    projections <- cbind(projections, prjs)
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

  exportTestValues(projections = {
    projections
  })
  return(projections)
})

# initializeGroupNames ----
# TODO shouldn't this be an observer???
initializeGroupNames <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "initializeGroupNames started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "initializeGroupNames")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "initializeGroupNames")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("initializeGroupNames",
      id = "initializeGroupNames",
      duration = NULL
    )
  }

  scEx <- scEx()
  if (is.null(scEx)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/initializeGroupNames.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/initializeGroupNames.RData")
  isolate({
    df <-
      data.frame(
        all = rep(TRUE, dim(scEx)[2]),
        none = rep(FALSE, dim(scEx)[2])
      )
    rownames(df) <- colnames(scEx)
    groupNames[["namesDF"]] <- df
  })
})

observe(initializeGroupNames())
# sample --------
sample <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "sample started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sample")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sample")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sample", id = "sample", duration = NULL)
  }

  scEx <- scEx()
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "sample: NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sample.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/sample.RData")

  pd <- colData(scEx)
  retVal <- NULL
  for (pdColName in colnames(pd)) {
    # in case dbCluster is already in the colData this would create dbCluster.1 later on
    if (pdColName == "dbCluster") next()
    if (length(levels(factor(pd[, pdColName]))) < 30) {
      if (is.null(retVal)) {
        retVal <- data.frame(pd[, pdColName])
        colnames(retVal) <- pdColName
      }
      retVal[, pdColName] <- factor(as.character(pd[, pdColName]))
    }
  }
  rownames(retVal) <- rownames(pd)
  exportTestValues(sample = {
    retVal
  })
  return(retVal)
})

# geneCount --------
geneCount <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "geneCount started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "geneCount")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "geneCount")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("geneCount", id = "geneCount", duration = NULL)
  }

  scEx_log <- scEx_log()

  if (is.null(scEx_log)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/geneCount.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/geneCount.RData")

  retVal <- Matrix::colSums(assays(scEx_log)[["logcounts"]] > 0)

  exportTestValues(geneCount = {
    retVal
  })
  return(retVal)
})

# beforeFilterPrj ----
beforeFilterPrj <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "beforeFilterPrj started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "beforeFilterPrj")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "beforeFilterPrj")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("beforeFilterPrj", id = "beforeFilterPrj", duration = NULL)
  }

  scEx <- scEx()
  bfc <- beforeFilterCounts()

  if (is.null(scEx) | is.null(bfc)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/beforeFilterPrj.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/beforeFilterPrj.RData")

  cn <- colnames(scEx)
  retVal <- bfc[cn]

  exportTestValues(beforeFilterPrj = {
    retVal
  })
  return(retVal)
})

# umiCount ----
umiCount <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "umiCount started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "umiCount")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "umiCount")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("umiCount", id = "umiCount", duration = NULL)
  }

  scEx <- scEx()

  if (is.null(scEx)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/umiCount.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/umiCount.RData")

  retVal <- Matrix::colSums(assays(scEx)[["counts"]])

  exportTestValues(umiCount = {
    retVal
  })
  return(retVal)
})

# sampleInfoFunc ----
sampleInfoFunc <- function(scEx) {
  # gsub(".*-(.*)", "\\1", scEx$barcode)
  colData(scEx)$sampleNames
}


# sampleInfo -------
# sample information
sampleInfo <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "sampleInfo started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sampleInfo")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sampleInfo")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sampleInfo", id = "sampleInfo", duration = NULL)
  }

  scEx <- scEx()
  if (!exists("scEx")) {
    if (DEBUG) {
      cat(file = stderr(), "sampleInfo: NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sampleInfo.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/sampleInfo.RData")

  retVal <- sampleInfoFunc(scEx)

  exportTestValues(sampleInfo = {
    retVal
  })
  return(retVal)
})


# table of input cells with sample information
# TODO: used in tableSeletionServer table; should be divided into function and reactive
inputSample <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "inputSample started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "inputSample")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "inputSample")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("inputSample", id = "inputSample", duration = NULL)
  }

  dataTables <- inputData()

  if (is.null(dataTables)) {
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("inputSample", id = "inputSample", duration = NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/inputSample.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/SCHNAPPsDebug/inputSample.RData')

  # TODO should come from sampleInfo
  sampInf <- gsub(".*-(.*)", "\\1", dataTables$scEx$barcode)
  cellIds <- data.frame(
    cellName = colnames(dataTables$scEx),
    sample = sampInf,
    ngenes = Matrix::colSums(assays(dataTables$scEx)[[1]])
  )

  if (DEBUG) {
    cat(file = stderr(), "inputSample: done\n")
  }
  if (dim(cellIds)[1] > 1) {
    return(cellIds)
  } else {
    return(NULL)
  }
})


getMemoryUsed <- reactive({
  suppressMessages(require(pryr))
  if (DEBUG) {
    cat(file = stderr(), "getMemoryUsed started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "getMemoryUsed")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "getMemoryUsed")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("getMemoryUsed", id = "getMemoryUsed", duration = NULL)
  }
  # number of times memory was calculated
  # not used anymore
  # umu <- updateMemUse$update
  paste(utils:::format.object_size(mem_used(), "auto"))
})

# used in coExpression, subclusterAnalysis, moduleServer, generalQC, DataExploration
# TODO change to scEx_log everywhere and remove
# log2cpm <- reactive({
#   if (DEBUG) {
#     cat(file = stderr(), "log2cpm\n")
#   }
#   scEx_log <- scEx_log()
#   if (is.null(scEx_log)) {
#     if (DEBUG) {
#       cat(file = stderr(), "log2cpm: NULL\n")
#     }
#     return(NULL)
#   }
#   log2cpm <- as.data.frame(as.matrix(assays(scEx_log)[[1]]))
#
#   return(log2cpm)
# })


reportFunction <- function(tmpPrjFile) {
  return(NULL)
}

# reacativeReport ----
reacativeReport <- function() {
  if (DEBUG) {
    cat(file = stderr(), "reacativeReport started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "reacativeReport")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "reacativeReport")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("reacativeReport", id = "reacativeReport", duration = NULL)
  }

  scEx <- scEx()
  projections <- projections()
  scEx_log <- scEx_log()
  inputNames <- names(input)
  # browser()

  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "output$report:NULL\n")
    }
    return(NULL)
  }

  tmpPrjFile <-
    tempfile(
      pattern = "file",
      tmpdir = .schnappsEnv$reportTempDir,
      fileext = ".RData"
    )

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
  pca <- pca()
  tsne <- tsne()

  scEx <- consolidateScEx(scEx, projections, scEx_log, pca, tsne)
  reportTempDir <- get("reportTempDir", envir = .schnappsEnv)
  base::save(
    file = tmpPrjFile,
    list = c(
      "reportTempDir",
      "projections",
      "scEx_log",
      "scEx",
      "report.env",
      ".schnappsEnv"
    )
  )
  userDataEnv <-
    as.environment(as.list(session$userData, all.names = TRUE))
  # browser()
  if (.schnappsEnv$DEBUGSAVE) {
    save(
      file = "~/SCHNAPPsDebug/tempReport.1.RData",
      list = c("session", "report.env", "file", ls(), ls(envir = globalenv()))
    )
  }
  # load('~/SCHNAPPsDebug/tempReport.1.RData')

  outZipFile <- paste0(.schnappsEnv$reportTempDir, "/report.zip")

  reactiveFiles <- ""

  # fixed files -----------
  tmpFile <-
    tempfile(
      pattern = "file",
      tmpdir = .schnappsEnv$reportTempDir,
      fileext = ".RData"
    )
  file.copy(paste0(packagePath, "/geneLists.RData"),
    tmpFile,
    overwrite = TRUE
  )
  file.copy(paste0(packagePath, "/Readme.txt"),
    .schnappsEnv$reportTempDir,
    overwrite = TRUE
  )
  reactiveFiles <-
    paste0(reactiveFiles,
      "#geneLists.RData\nload(file=\"",
      tmpFile,
      "\")\n",
      collapse = "\n"
    )

  # global variables
  tmpFile <-
    tempfile(
      pattern = "file",
      tmpdir = .schnappsEnv$reportTempDir,
      fileext = ".RData"
    )
  save(file = tmpFile, list = c("allowedColors"))
  reactiveFiles <-
    paste0(reactiveFiles,
      "#allowedColors",
      "\nload(file=\"",
      tmpFile,
      "\")\n",
      collapse = "\n"
    )

  # we load this twice since we need .schnappsEnv here for the first time...
  reactiveFiles <-
    paste0(
      reactiveFiles,
      "#load internal data\nload(file=\"",
      tmpPrjFile,
      "\")\nfor (n in names(report.env)) {assign(n,report.env[[n]])}\n",
      collapse = "\n"
    )
  # Projections -----
  # projections can contain mannually annotated groups of cells and different normalizations.
  # to reduce complexity we are going to save those in a separate RData file

  # return(tmpPrjFile)
  # the reactive.R can hold functions that can be used in the report to reduce the possibility of code replication
  # we copy them to the temp directory and load them in the markdown
  # localContributionDir <- .SCHNAPPs_locContributionDir
  uiFiles <-
    dir(
      path = c(
        paste0(packagePath, "/contributions"),
        localContributionDir
      ),
      pattern = "reactives.R",
      full.names = TRUE,
      recursive = TRUE
    )
  for (fp in c(
    paste0(packagePath, "/serverFunctions.R"),
    paste0(packagePath, "/reactives.R"),
    uiFiles
  )) {
    if (DEBUG) {
      cat(file = stderr(), paste("loading: ", fp, "\n"))
    }
    tmpFile <-
      tempfile(
        pattern = "file",
        tmpdir = .schnappsEnv$reportTempDir,
        fileext = ".R"
      )
    file.copy(fp, tmpFile, overwrite = TRUE)
    reactiveFiles <-
      paste0(
        reactiveFiles,
        "# load ",
        fp,
        "\nsource(\"",
        tmpFile,
        "\", local = TRUE)\n",
        collapse = "\n"
      )
  }
  # otherwise reactive might overwrite projections...
  reactiveFiles <-
    paste0(
      reactiveFiles,
      "#load internal data\nload(file=\"",
      tmpPrjFile,
      "\")\nfor (n in names(report.env)) {assign(n,report.env[[n]])}\n",
      collapse = "\n"
    )
  # encapsulte the load files in an R block
  LoadReactiveFiles <-
    paste0(
      "\n\n```{r load-reactives, include=FALSE}\n",
      reactiveFiles,
      "\n```\n\n"
    )

  # handle plugin reports
  # load contribution reports
  # parse all report.Rmd files under contributions to include in application
  uiFiles <-
    dir(
      path = c(
        paste0(packagePath, "/contributions"),
        localContributionDir
      ),
      pattern = "report.Rmd",
      full.names = TRUE,
      recursive = TRUE
    )
  pluginReportsString <- ""
  fpRidx <- 1
  for (fp in uiFiles) {
    if (DEBUG) {
      cat(file = stderr(), paste("loading: ", fp, "\n"))
    }
    tmpFile <-
      tempfile(
        pattern = "file",
        tmpdir = .schnappsEnv$reportTempDir,
        fileext = ".Rmd"
      )
    file.copy(fp, tmpFile, overwrite = TRUE)
    pluginReportsString <- paste0(
      pluginReportsString,
      "\n\n```{r child-report-",
      fpRidx,
      ", child = '",
      tmpFile,
      "', eval=TRUE}\n```\n\n"
    )
    fpRidx <- fpRidx + 1
  }

  # Copy the report file to a temporary directory before processing it, in
  # case we don't have write permissions to the current working dir (which
  # can happen when deployed).
  tempReport <- file.path(.schnappsEnv$reportTempDir, "report.Rmd")

  # tempServerFunctions <- file.path(.schnappsEnv$reportTempDir, "serverFunctions.R")
  # file.copy("serverFunctions.R", tempServerFunctions, overwrite = TRUE)

  # create a new list of all parameters that can be passed to the markdown doc.

  params <- list( # tempServerFunctions = tempServerFunctions,
    # tempprivatePlotFunctions = tempprivatePlotFunctions,
    calledFromShiny = TRUE # this is to notify the markdown that we are running the script from shiny. used for debugging/development
    # save the outputfile name for others to use to save
    # params$outputFile <- file$datapath[1])
  )

  # for (idx in 1:length(names(input))) {
  #   params[[inputNames[idx]]] <- input[[inputNames[idx]]]
  # }
  params[["reportTempDir"]] <- .schnappsEnv$reportTempDir

  file.copy(paste0(packagePath, "/report.Rmd"), tempReport, overwrite = TRUE)

  # read the template and replace parameters placeholder with list
  # of paramters
  x <- readLines(tempReport)
  # x <- readLines("report.Rmd")
  paramString <-
    paste0("  ", names(params), ": NA", collapse = "\n")
  y <- gsub("#__PARAMPLACEHOLDER__", paramString, x)
  y <- gsub("__CHILDREPORTS__", pluginReportsString, y)
  y <- gsub("__LOAD_REACTIVES__", LoadReactiveFiles, y)
  # cat(y, file="tempReport.Rmd", sep="\n")
  cat(y, file = tempReport, sep = "\n")

  if (DEBUG) {
    cat(file = stderr(), "output$report:scEx:\n")
  }
  if (DEBUG) {
    cat(file = stderr(), paste("\n", tempReport, "\n"))
  }
  # Knit the document, passing in the `params` list, and eval it in a
  # child of the global environment (this isolates the code in the document
  # from the code in this app)
  if (DEBUG) {
    file.copy(tempReport, "~/SCHNAPPsDebug/tempReport.Rmd")
  }
  myparams <-
    params # needed for saving as params is already taken by knitr
  # if (.schnappsEnv$DEBUGSAVE)
  # save(file = "~/SCHNAPPsDebug/tempReport.RData", list = c("session", "myparams", ls(), "zippedReportFiles"))
  # load(file = '~/SCHNAPPsDebug/tempReport.RData')
  if (DEBUG) cat(file = stderr(), paste("workdir: ", getwd()))
  suppressMessages(require(callr))
  # if (.schnappsEnv$DEBUGSAVE)
  # file.copy(tempReport, "~/SCHNAPPsDebug/tmpReport.Rmd", overwrite = TRUE)

  # tempReport = "~/SCHNAPPsDebug/tmpReport.Rmd"
  # file.copy("contributions/gQC_generalQC//report.Rmd",
  #           '/var/folders/tf/jwlc7r3d48z7pkq0w38_v7t40000gp/T//RtmpTx4l4G/file1a6e471a698.Rmd', overwrite = TRUE)
  tryCatch(
    callr::r(
      function(input, output_file, params, envir)
        rmarkdown::render(
          input = input,
          output_file = output_file,
          params = params,
          envir = envir
        ),
      args = list(
        input = tempReport,
        output_file = "report.html",
        params = params,
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

  tDir <- paste0(.schnappsEnv$reportTempDir, "/")
  base::file.copy(tmpPrjFile, paste0(.schnappsEnv$reportTempDir, "/sessionData.RData"))
  write.csv(as.matrix(assays(scEx_log)[[1]]),
    file = paste0(.schnappsEnv$reportTempDir, "/normalizedCounts.csv")
  )
  base::save(
    file = paste0(.schnappsEnv$reportTempDir, "/inputUsed.RData"),
    list = c("scEx", "projections")
  )
  zippedReportFiles <- c(paste0(tDir, zippedReportFiles))
  zip(outZipFile, zippedReportFiles, flags = "-9Xj")
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

consolidateScEx <-
  function(scEx, projections, scEx_log, pca, tsne) {
    reducedDims(scEx) <- SimpleList(PCA = pca$x, TSNE = tsne)
    assays(scEx)[["logcounts"]] <- assays(scEx_log)[[1]]
    colData(scEx)[["before.Filter"]] <- projections$before.filter
    colData(scEx)[["dbCluster"]] <- projections$dbCluster
    colData(scEx)[["UmiCountPerGenes"]] <- projections$UmiCountPerGenes
    colData(scEx)[["UmiCountPerGenes2"]] <- projections$UmiCountPerGenes2

    return(scEx)
  }


if (DEBUG) {
  cat(file = stderr(), "\n\ndone loading reactives.R.\n\n\n")
}
