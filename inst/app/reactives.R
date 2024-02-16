# inst/app/reactives.R

suppressMessages(require(tibble))
suppressMessages(require(Seurat))
suppressMessages(require(scran))
suppressMessages(require(irlba))
suppressMessages(require(BiocSingular))
suppressMessages(require(dplyr))
suppressMessages(require(BPCells))
suppressMessages(require(ggalluvial))
library(SingleCellExperiment)
# require(tidySingleCellExperiment)
# base::source(paste0(packagePath, "/outputs.R"), local = TRUE)
if ("crayon" %in% rownames(installed.packages()) == FALSE) {
  green <- function(x) {
    x
  }
} else {
  require(crayon)
}

# reactive values  ------------------------------------------------------------------
inputFileStats <- reactiveValues(stats = NULL)

# store groups of cells that are defined on the fly using the modular 2D plot
groupNames <- reactiveValues(namesDF = data.frame()) 

# colors for samples
sampleCols <- reactiveValues(colPal = allowedColors)

# colors for clusters
clusterCols <- reactiveValues(colPal = allowedColors)



# 
allCellNames <- reactiveVal(
  label = "CellNames",
  value = ""
)
if (exists(".SCHNAPPs_DEBUGSAVE", envir = .schnappsEnv)) {
  DEBUGSAVE <- get(".SCHNAPPs_DEBUGSAVE", envir = .schnappsEnv)
}
# Input file either rdata file or csv file

inputFile <- reactiveValues(
  inFile = "",
  annFile = ""
)

# add comment to history ----

commentModal <- function(failed = FALSE) {
  modalDialog(
    # TODO
    #  mce not working, maybe this helps eventually: https://github.com/twbs/bootstrap/issues/549
    # if ("shinyMCE" %in% rownames(installed.packages())) {
    #   shinyMCE::tinyMCE(
    #     "Comment4history",
    #     "Please describe your work. This will be included in the history"
    #   )
    # } else {
    sc_textInput("Comment4history", "Please describe your work. This will be included in the history")
    # }
    ,
    footer = tagList(
      modalButton("Cancel"),
      actionButton("commentok", "OK")
    )
  )
}

# Show modal when button is clicked.
observeEvent(input$comment2History, {
  deepDebug()
  showModal(commentModal())
})
# Show modal when button is clicked.
# # TODO modal to confirm
observeEvent(input$Quit, {
  deepDebug()
  showModal(modalDialog(
    title="Quit SCHNAPPs",
    footer = tagList(actionButton("confirmQuit", "yes, quit"),
                     modalButton("Cancel")
    )
  ))
  
})

observeEvent(input$confirmQuit, {
  if(DEBUG) {cat(file = stderr(), "\nquit the app\n\n")}
  add2history(type = "save", input=isolate( reactiveValuesToList(input)), 
              comment = "manual quit", sessionInfo = sessionInfo())
  
  stopApp()
  removeModal()
})





observeEvent(input$openBrowser, {
  deepDebug()
  require(rmarkdown)
  # browser()
  knitr::opts_chunk$set(
    message = FALSE,
    warning = FALSE,
    echo = FALSE,
    include = TRUE
  )
  inputList = reactiveValuesToList(input)
  scEx_log = scEx_log()
  scEx = scEx()
  projections = projections()
  pca = pcaReact()
  if(is.null(scEx_log)) scEx_log = "NULL"
  tryCatch(
    rmarkdown::render("test.Rmd",
                      output_file = paste0("test.report.html"),
                      output_format = "html_document",
                      params = list(
                        scEx = scEx,
                        input = inputList,
                        scEx_log = scEx_log,
                        projections = projections(),
                        projectionColors = isolate( projectionColors) %>% reactiveValuesToList(),
                        pa = pca,
                        projections
                      )
    ),
    error = function(e) {cat(file = stderr(),paste("Error\n",e,"\n")); NULL}
  )
  
  
}
)
# When OK button is pressed, attempt to load the data set. If successful,
# remove the modal. If not show another modal, but this time with a failure
# message.
observeEvent(input$commentok, {
  # cat(file = stderr(), paste0("commentok: \n"))
  deepDebug()
  comment <- input$Comment4history
  add2history(type = "text", comment = "",  input=isolate( reactiveValuesToList(input)), text2add = comment)
  removeModal()
})
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
  if(is.null(fpLs)) {return(NULL)}
  
  # in case we are loading from a report
  scExFound <- FALSE
  allScEx_log <- NULL
  if ("scEx" %in% fpLs) { # we prefer the scEx object ;)
    scExFound <- TRUE
    varName <- "scEx"
    #when loading history files
    if (is(scEx, "list")){
      scEx = scEx[[1]]
    }
  } else {
    scExFound <- FALSE
    for (varName in fpLs) {
      # if ("SingleCellExperiment" %in% class(get(varName))) {
      if (is(get(varName), "SingleCellExperiment")) {
        scEx <- get(varName)
        scExFound <- TRUE
        break()
      } else{
        ## list as from history file
        weiter = tryCatch({
          is(get(varName)[[1]], "SingleCellExperiment")} , error = function(e) {
            return(FALSE)
          })
        if (weiter) {
          scEx <- get(varName)[[1]]
          if (!"counts" == names(assays(scEx))[1]) {
            assays(scEx)["counts"] = assays(scEx)[1]
          }
          scExFound <- TRUE
          break()
          
          
        }
      }
    }
  }
  if (!scExFound) {
    return(NULL)
  }
  cat(file = stderr(), green(paste("file ", inFile$name[1], "contains variable", varName, " as SingleCellExperiment.\n")))
  
  # reducedDims
  if(!isEmpty(reducedDims(scEx))){
    colData(scEx)
    for(rdN in 1:length(reducedDims(scEx))){
      rDim = reducedDims(scEx)[[rdN]]
      colnames(rDim) = make.unique(c(colnames(colData(scEx)),colnames(rDim)))[-c(1:ncol(colData(scEx)))]
      colData(scEx) <- cbind(colData(scEx),rDim)
    }
    reducedDims(scEx) <- NULL
  }
  
  # save(file = "~/SCHNAPPsDebug/inputProblem.RData", list = ls())
  # load("~/SCHNAPPsDebug/inputProblem.RData")
  # fdAll <- rowData(scEx)
  # pdAll <- colData(scEx)
  # exAll <- assays(scEx)[["counts"]]
  stats[1, "nFeatures"] <- nrow(rowData(scEx))
  stats[1, "nCells"] <- nrow(colData(scEx))
  allScEx <- scEx # save to allScEx to be able to combine if multiple SCE objects are given
  # In case we are loading precalculated normalizations.
  if ("logcounts" %in% names(assays(scEx))) {
    allScEx_log <- scEx
    if(! "counts" %in% names(assays(scEx)))
      assays(scEx)[["counts"]] = assays(scEx)[["logcounts"]]
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/readInp1.RData", list = c(ls()))
  }
  # cp = load(file='~/SCHNAPPsDebug/readInp1.RData')
  
  # read multiple files [2:1] => c(2,1) and not empty list
  # inFile$datapath[2] = "/Volumes/@Single_cell/SingleCellCourse_2021/schnapps/Tcells_d2_10X.schnapps.RData"
  if (length(inFile$datapath) > 1) {
    for (fpIdx in 2:length(inFile$datapath)) {
      # inFile$datapath[fpIdx] <- "~/Downloads/paper1.RData"
      cat(file = stderr(), paste("reading", inFile$name[fpIdx], "\n"))
      fp <- inFile$datapath[fpIdx]
      fpLs <- load(fp)
      scExFound <- FALSE
      varName = "not found"
      if ("scEx" %in% fpLs) {
        scExFound <- TRUE
        varName <- "scEx" # we don't have to set scEx
        # In case we are loading precalculated normalizations.
      } else {
        for (varName in fpLs) {
          # if ("SingleCellExperiment" %in% class(get(varName))) {
          if (is(get(varName), "SingleCellExperiment")) {
            scEx <- get(varName)
            scExFound <- TRUE
            break()
          }
        }
      }
      # reducedDims
      if(!isEmpty(reducedDims(scEx))){
        colData(scEx)
        for(rdN in 1:length(reducedDims(scEx))){
          rDim = reducedDims(scEx)[[rdN]]
          colnames(rDim) = make.unique(c(colnames(colData(scEx)),colnames(rDim)))[-c(1:ncol(colData(scEx)))]
          colData(scEx) <- cbind(colData(scEx),rDim)
        }
        reducedDims(scEx) <- NULL
      }
      
      cat(file = stderr(), green(paste("file ", inFile$datapath[fpIdx], "contains variable", varName, "\n")))
      if (!scExFound) {
        next()
      }
      if ("counts" %in% names(assays(scEx))) {
        if (is.null(allScEx)) {
          allScEx <- scEx
        } else {
          
          
          newColnames <- make.unique(c(colnames(allScEx), colnames(scEx)))
          colnames(scEx) <- newColnames[(ncol(allScEx) + 1):length(newColnames)]
          genesUnion <- intersect(rownames(scEx), rownames(allScEx))
          
          allScEx <- addColData(allScEx, scEx)
          scEx <- addColData(scEx, allScEx)
          
          # if (!is.null(allScEx_log)) scEx <- addColData(scEx, allScEx_log)
          # colData(allScEx_log)
          
          # in case one of the loaded files doesn't have all the assays. (since cbind cannot handle this)
          if(!all(union(names(assays(allScEx)), names(assays(scEx))) == 
                  intersect(names(assays(allScEx)), names(assays(scEx))))){
            cTemp <- counts(allScEx)
            assays(allScEx) <- S4Vectors::SimpleList()
            counts(allScEx) <- cTemp
            cTemp <- counts(scEx)
            assays(scEx) <- S4Vectors::SimpleList()
            counts(scEx) <- cTemp
          }
          allScEx <- SingleCellExperiment::cbind(allScEx[genesUnion, ], scEx[genesUnion, ])
        }
      } else {
        cat(file = stderr(), (paste("!!!!!file ", inFile$datapath[fpIdx], "contains variable", varName, " but no counts assay!!!!\n")))
        next()
      }
      cat(file = stderr(), paste("nCells: ", ncol(allScEx),"\n"))
      
      # In case we are loading precalculated normalizations we need to also keep the logcounts.
      if ("logcounts" %in% names(assays(scEx))) {
        if (is.null(allScEx_log)) {
          allScEx_log <- scEx
        } else {
          newColnames <- make.unique(c(colnames(allScEx_log), colnames(scEx)))
          colnames(scEx) <- newColnames[(ncol(allScEx_log) + 1):length(newColnames)]
          genesUnion <- intersect(rownames(scEx), rownames(allScEx_log))
          
          allScEx_log <- addColData(allScEx_log, scEx)
          scEx <- addColData(scEx, allScEx_log)
          # colData(allScEx_log)
          allScEx_log <- SingleCellExperiment::cbind(allScEx_log[genesUnion, ], scEx[genesUnion, ])
        }
      }
      stats[fpIdx, "nFeatures"] <- nrow(scEx)
      stats[fpIdx, "nCells"] <- ncol(scEx)
    }
  }
  # save(file = "~/SCHNAPPsDebug/inputProblem.RData", list = c("allScEx", "allScEx_log", "sampleCols"))
  # load(file = "~/SCHNAPPsDebug/inputProblem.RData")
  if ("sampleNames" %in% colnames(colData(allScEx))) {
    if (!is(colData(allScEx)$sampleNames, "factor")) {
      colData(allScEx)$sampleNames <- factor(colData(allScEx)$sampleNames)
    }
    sampNames <- levels(colData(allScEx)$sampleNames)
    # isolate({
    #   projectionColors$sampleNames <- allowedColors[seq_along(sampNames)]
    #   names(projectionColors$sampleNames) <- sampNames
    # })
  } else {
    showNotification(
      "scEx - colData doesn't contain sampleNames",
      duration = NULL,
      type = "error"
    )
    colData(allScEx)$sampleNames <- 1
    colData(allScEx)$sampleNames = as.factor(colData(allScEx)$sampleNames)
  }
  
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/readInp.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/readInp.RData')
  
  
  #' save orginial data
  dataTables <- list()
  featuredata <- rowData(allScEx)
  # handle different extreme cases for the symbol column (already encountered)
  if (is.factor(featuredata$symbol)) {
    if (levels(featuredata$symbol) == "NA") {
      featuredata$symbol <- make.unique(rownames(featuredata))
      rowData(allScEx) <- featuredata
    }
  }
  if ("symbol" %in% colnames(featuredata)) {
    featuredata$symbol <- make.unique(featuredata$symbol)
    rowData(allScEx) <- featuredata
  }
  
  # dataTables$featuredataOrg <- rowData(allScEx)
  dataTables$scEx <- allScEx
  dataTables$featuredata <- featuredata
  if (!is.null(allScEx_log)) {
    assays(allScEx_log)[["counts"]] <- NULL
  }
  dataTables$scEx_log <- allScEx_log
  
  if (is.null(allScEx$barcode)) {
    showNotification("scEx doesn't contain barcode column", type = "error")
    return(NULL)
  }
  # some checks
  
  if (sum(is.infinite(assays(allScEx)[["counts"]])) > 0) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("scEx contains infinite values",
                       type = "error"
      )
    }
    return(NULL)
  }
  
  
  if (sum(c("id", "symbol") %in% colnames(rowData(allScEx))) < 2) {
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
  # deepDebug()
  
  isolate(allCellNames(rownames(colData(dataTables$scEx))))
  inputFileStats$stats <- stats
  return(dataTables)
} # end of inputDataFunc

readMM <- function(inFile) {
  data <- Matrix::readMM(file = inFile$datapath)
  
  return(NULL)
}

# inFile$datapath="~/Downloads/GSE122084_RAW/GSM3855868_Salmonella_exposed_cells.txt"
# load("data/scEx.RData")
# scMat = as.matrix(assays(scEx)[[1]])
# write.file.csv(scMat, row.names=TRUE, file="data/scEx.csv" )


## readGMT ----
#' Read GMT File
#'
#' Read GMT file and return a list containing the name, description, and filtered genes for each element.
#'
#' @param inFile A list containing the path of the file to be read (\code{inFile$datapath}).
#' @param scEx A SingleCellExperiment object.
#' @return A list containing the name, description, and filtered genes for each element.
#' @export
#'
#' @examples
#' inFile <- list()
#' inFile$datapath <- 'h.all.v2022.1.Hs.symbols.gmt'
#' scEx <- data.frame(symbol = c("A1BG", "A2M", "A4GALT"))
#' readGMT(inFile, scEx)

readGMT <- function(inFile, scEx){
  # inFile = list()
  # inFile$datapath='h.all.v2022.1.Hs.symbols.gmt'
  
  # Open the file specified in inFile$datapath for reading
  con <- file(inFile$datapath, 'r')
  
  # Read the first line of the file
  first_line <- readLines(con, n = 1)
  
  # Close the file
  close(con)
  
  # Count the number of occurrences of commas, tabs, and spaces in the first line
  commaCount <- length(gregexpr(',', first_line, perl = T)[[1]])
  tabCount <- length(gregexpr('\t', first_line, perl = T)[[1]])
  spaceCount <- length(gregexpr(' ', first_line, perl = T)[[1]])
  
  # Set the separator based on the counts
  sep <- ','
  if (spaceCount > commaCount) sep <- ' '
  if (commaCount > tabCount) sep <- ','
  if (tabCount > commaCount) sep <- '\t'
  
  # Open the file again
  con <- file(inFile$datapath,'r')
  
  # Read all the lines from the file
  dat <- readLines(con = con)
  
  # Close the file
  close(con)
  
  # Split each line into a list of elements using the separator
  dat <- strsplit(dat, sep)
  
  # If the DEBUGSAVE flag is set, save the variables to a file
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = '~/SCHNAPPsDebug/readGMT.RData', list = c(ls()))
  }
  
  # Define a function to process each element in the list
  FUN = function(x){
    name = x[1]
    desc = x[2]
    genes = x[3:length(x)]

    # Find the genes that are not in the rowData of scEx
    noGenes = genes[which(!genes %in% rowData(scEx)$symbol)]
    
    # If the DEBUG flag is set, print a message listing the genes not found
    if(.schnappsEnv$DEBUG){
      cat(file = stderr(), paste('GMT reader: ', name, 'genes not found:', paste(noGenes, collapse = ','), '\n'))
    }
    
    # Filter the genesthat are found in the rowData of scEx
    genes = genes[which(genes %in% rowData(scEx)$symbol)]
    
    # If no genes are found, print a message and return NULL
    if (length(genes) == 0) {
      if(.schnappsEnv$DEBUG){
        cat(file = stderr(), paste('GMT reader: ', name, 'no genes found:\n'))
      }
      return(NULL)
    }
    
    # Return a list containing the name, description, and filtered genes
    return(list(name = name, desc = desc, genes = genes))
    
  }
  
  # Apply the FUN function to each element in dat and store the results in retVal
  retVal = lapply(dat,FUN = FUN)
  
  # Filter out elements in retVal that have length 0
  retVal = retVal[lapply(retVal,length)>0]
  
  # Set the names of retVal based on the name attribute of each element
  names(retVal) = lapply(retVal, FUN=function(x)x$name) %>% unlist
  
  # Return retVal
  retVal
}

## readCSV ----
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
    save(file = "~/SCHNAPPsDebug/readCSV.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/readCSV.RData')
  if (colnames(data)[1] %in% c("", "rownames", "ROWNAMES", "genes")) {
    rownames(data) <- make.unique(data[, 1])
    data <- data[, -1]
  }
  exAll <- as(as.matrix(data), "TsparseMatrix")
  rownames(exAll) <- rownames(data)
  colnames(exAll) <- colnames(data)
  if (all(stringr::str_detect(colnames(data),"-"))) {
    sampleNames = unlist(lapply(colnames(data), function(x) stringr::str_split(x,"-")[[1]][2]))
  } else {
    sampleNames = as.factor(rep(tools::file_path_sans_ext(inFile$name) %>% make.names(), ncol(data)))
  }
  scExnew <- SingleCellExperiment(
    assay = list(counts = exAll),
    colData = list(
      sampleNames = sampleNames,
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
  isolate(allCellNames(rownames(colData(dataTables$scEx))))
  return(dataTables)
}

#
# inFile$datapath="data/scEx.csv"
# load("data/scEx.RData")
# write.file.csv(colData(scEx), row.names=TRUE, file="data/scExCells.csv" )
# write.file.csv(rowData(scEx), row.names=TRUE, file="data/scExGenes.csv" )
# annFile$datapath="data/scExGenes.csv"
# appendAnnotation ----
#
# append annotation to singleCellExperiment object
# uses colData
appendAnnotation <- function(scEx, annFile) {
  if (DEBUG) {
    cat(file = stderr(), "appendAnnotation\n")
  }
  rDat <- rowData(scEx)
  cDat <- colData(scEx)
  success = FALSE
  for (fpIdx in 1:length(annFile$datapath)) {
    data <- read.table(file = annFile$datapath[fpIdx], check.names = FALSE, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/appendAnnotation.RData", list = c(ls()))
    }
    # cp = load(file = "~/SCHNAPPsDebug/appendAnnotation.RData")
    if (colnames(data)[1] %in% c("", "rownames", "ROWNAMES", "CELL_ID")) {
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
      success = TRUE
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
      success = TRUE
    }
    
    # colData
    if (any(rownames(data) %in% rownames(cDat))) {
      # if any factor is already set we need to avoid having different levels
      if (any(colnames(cDat) %in% colnames(data))) {
        commonCols <- colnames(cDat) %in% colnames(data)
        for (cCol in colnames(cDat)[commonCols]) {
          if (is(cDat[, cCol], "factor")) {
            cDat[, cCol] <- factor(data[, cCol])
            levels(cDat[, cCol]) = make.names(levels(cDat[, cCol]))
          }
        }
      }
      # a vector with less than 20 levels (unique values) will be converted in a factor
      for (cn in colnames(data)) {
        if (is(data[, cn], "factor")) next()
        lv <- levels(factor(data[, cn]))
        if (length(lv) <= 20) {
          data[, cn] <- factor(data[, cn])
          levels(data[,cn]) = make.names(levels(data[,cn]))
        }
      }
      # only take cells that are also in the data set
      data = data[rownames(data) %in% colnames(scEx),]
      data = data[rownames(cDat),]
      cDat = cbind(cDat, data)
      success = TRUE
    }
  }
  cDat$sampleNames <- factor(cDat$sampleNames)
  colData(scEx) <- cDat
  rowData(scEx) <- rDat
  if(!success){
    showNotification(
      "loading annotations didn't result in any added annotations",
      type = "error",
      duration = NULL
    )
    
  }
  if (DEBUG) {
    cat(file = stderr(), "appendAnnotation done\n")
  }
  return(scEx)
}


### dataFile reactive ----
dataFile <- reactive({
  deepDebug()
  if(!exists("restoreHistory",envir = .schnappsEnv)){
    .schnappsEnv$restoreHistory = FALSE
  }
  if( .schnappsEnv$restoreHistory & !is.null(.schnappsEnv$inputFile)){
    iFile = .schnappsEnv$inputFile
  }else{
    iFile = input$file1
  }
  return(iFile)
})

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
  inFile <- dataFile()
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
    save(file = "~/SCHNAPPsDebug/inputData.RData", list = c(ls()))
  }
  # cp =load(file='~/SCHNAPPsDebug/inputData.RData')
  
  # inFile$datapath = "data/sessionData.RData"
  # annFile$datapath = "data/selectedCells.csv"
  # TODO either multiple files with rdata/rds or single file of csv/other
  fpExtension <- tools::file_ext(inFile$datapath[1])
  if (toupper(fpExtension) %in% c("RDATA", "RDS")) {
    retVal <- tryCatch({inputDataFunc(inFile)}, 
                       error = function(e) {
                         cat(file = stderr(), paste("inputData: NULL", e,"\n"))
                         return(NULL)})
  } else {
    retVal <- tryCatch({readCSV(inFile)},
                       error = function(e){
                         cat(file = stderr(), paste("inputData: NULL", e,"\n"))
                         return(NULL)}
    )
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
  # deepDebug()
  # browser()
  
  isolate({
    if (is.null(names(projectionColors$sampleNames)) | !all(names(projectionColors$sampleNames) %in% sampNames)){
      projectionColors$sampleNames <- rev(allowedColors[rep(1:length(allowedColors),ceiling(length(sampNames) / length(allowedColors)))[1:length(sampNames)]])
      # projectionColors$sampleNames <- rev(allowedColors)[seq_along(sampNames)]
      names(projectionColors$sampleNames) <- sampNames
      add2history(type = "save", input=isolate( reactiveValuesToList(input)), comment = "scol", scol = projectionColors$sampleNames)
    }
  })
  inputFile$inFile <- paste(inFile$name, collapse = ", ")
  inputFile$annFile <- paste(annFile$name, collapse = ", ")
  
  
  # subsample this number of cells per sample
  # no replacement. it is possible that one sample is more represented if the required number of cells is larger than there are cells for this sample.
  # save(file = "~/SCHNAPPsDebug/inputData.RData", list = c(ls()), compress = F)
  # load(file='~/SCHNAPPsDebug/inputData.RData')
  if (sampleCells) {
    tb = table(colData(retVal$scEx)$sampleNames)
    smpIdx = c()
    for (nm in names(tb)) {
      nidx = which(colData(retVal$scEx)$sampleNames == nm)
      if (length(nidx) <= subsampleNum){ 
        smpIdx = c(smpIdx, nidx)
      } else { 
        samp <- base::sample(
          x = length(nidx),
          size = subsampleNum,
          replace = FALSE
        )
        smpIdx = c(smpIdx, nidx[samp])
      }
    }
    # samp <- base::sample(
    #   x = ncol(assays(retVal$scEx)[["counts"]]),
    #   size = min(subsampleNum, ncol(assays(retVal$scEx)[["counts"]])),
    #   replace = FALSE
    # )
    retVal$scEx <- retVal$scEx[, smpIdx]
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
  gmtData()
  scEx <- scEx()
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "medianUMI:NULL\n")
    }
    return(0)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/medianUMI.RData", list = c(ls()))
  }
  #cp =  load(file='~/SCHNAPPsDebug/medianUMI.RData')
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
    # cp = load(file='~/SCHNAPPsDebug/useCellsFunc.Rdata')
    req(dataTables$scEx)
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
  if (!exists("dataTables") || is.null(dataTables) || ! "scEx" %in% names(dataTables)) {
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
  
  setRedGreenButton(
    vars = list(
      c("minExpGenes", isolate(input$minExpGenes)),
      c("minGenes", isolate(input$minGenes)),
      c("maxGenes", isolate(input$maxGenes)),
      c("cellPatternRM", isolate(input$cellPatternRM)),
      c("cellKeep", isolate(input$cellKeep)),
      c("cellKeepOnly", isolate(input$cellKeepOnly)),
      c("cellsFiltersOut", isolate(input$cellsFiltersOut)),
      c("minNonExpGenes", isolate(input$minNonExpGenes))
    ),
    button = "updateCellSelectionParameters"
  )
  
  exportTestValues(useCells = {
    retVal
  })
  return(retVal)
})

# useGenesFunc ----
useGenesFunc <-
  function(dataTables,
           ipIDs,
           # regular expression of genes to be removed
           geneListSelection,
           genesKeep,
           geneLists) {
    if (DEBUG) {
      cat(file = stderr(), "useGenesFunc started.\n")
      cat(file = stderr(), paste("useGenesFunc ipIDs:", ipIDs,"\n"))
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
      save(file = "~/SCHNAPPsDebug/useGenesFunc.Rdata", list = c(ls()))
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
    if (DEBUG) cat(file = stderr(), paste("useGenesFunc sum ret", sum(keepIDs),"\n"))
    return(keepIDs)
  }

# HVAinfoTable ----
HVAinfoTable <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "output$HVAinfoTable\n")
  }
  scEx_log <- scEx_log()
  scEx <- scEx()
  pca = pcaReact()
  hvgSelection <- isolate(input$hvgSelection)
  
  if (is.null(pca)) {
    return(NULL)
  }
  if (is.null(scEx_log)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/HVAinfoTable.RData",
         list = c( ls())
    )
  }
  # cp = load("~/SCHNAPPsDebug/HVAinfoTable.RData")
  
  hvgenes = rownames(pca$rotation)
  pcaN = length(hvgenes)
  # if(hvgSelection %in% c("vst", "mvp", "disp")){
  tryCatch(  switch(hvgSelection,
                    "vst" = {
                      pbmc <-
                        FindVariableFeatures(assays(scEx)[[1]], selection.method = "vst")
                      return(pbmc[order(pbmc$vst.variance.standardized, decreasing = T)[1:pcaN],])
                    },
                    "mvp" = {
                      pbmc <-
                        FindVariableFeatures(assays(scEx_log)[[1]], selection.method = "mean.var.plot")
                      pbmc <-
                        FindVariableFeatures(assays(scEx_log)[[1]], selection.method = "mean.var.plot")
                      pbmc$mvp.dispersion.scaled[which(is.na(pbmc$mvp.dispersion.scaled))] = max(pbmc$mvp.dispersion.scaled, na.rm = T) + 1
                      return(pbmc[order(pbmc$mvp.dispersion.scaled, decreasing = T)[1:pcaN],])
                    },
                    "disp" = {
                      pbmc <-
                        Seurat::FindVariableFeatures(assays(scEx_log)[[1]],selection.method = "dispersion", nfeatures = pcaN)
                      # NA values seem to have the highest variance. so we make sure to include them
                      pbmc$mvp.dispersion.scaled[which(is.na(pbmc$mvp.dispersion.scaled))] = max(pbmc$mvp.dispersion.scaled, na.rm = T) + 1
                      return(pbmc[order(pbmc$mvp.dispersion.scaled, decreasing = T)[1:pcaN],])
                      
                    },
                    "getTopHVGs" = {
                      # apply(as.matrix(assays(scEx_log)[[1]]), 1, max)
                      scEx = logNormCounts(scEx)
                      dec <- scran::modelGeneVar(scEx, assay.type = "logcounts")
                      # geneshvg <- getTopHVGs(dec, prop=0.1)
                      return(as.data.frame(dec[hvgenes,]))
                    }
  ),
  error = function(e) {
    cat(file = stderr(), paste("\n\n!!!Error during HVAinfoTable:\n", e, "\n\n"))
    return(NULL)
  }
  )
  
  #   seurDat <- tryCatch( 
  #     {
  #       # seurDat <- CreateAssayObject(
  #       #   data = as.matrix(assays(scEx_log)[[1]])
  #       # )
  #       # seurDat[["SCT"]] = 
  #       # vf = FindVariableFeatures(assays(scEx_log)[[1]], selection.method = selection.method)
  #       vf = vf[hvgenes,]
  #       return(vf@meta.features[order(vf@meta.features$vst.variance.standardized),])
  #     },
  #     error = function(e) {
  #       cat(file = stderr(), paste("\n\n!!!Error during HVAinfoTable:\n", e, "\n\n"))
  #       return(NULL)
  #     }
  #   )
  #   
  # }else {  #"getTopHVGs"
  #   # apply(as.matrix(assays(scEx_log)[[1]]), 1, max)
  #   scEx = logNormCounts(scEx)
  #   dec <- scran::modelGeneVar(scEx, assay.type = "logcounts")
  #   geneshvg <- getTopHVGs(dec, prop=0.1)
  #   return(as.data.frame(dec[geneshvg,]))
  # }
  
  cat(file = stderr(), "ERROR output$HVAinfoTable This should not happen!\n")
  return(NULL)
})


# PCAloadingsTable ----
PCAloadingsTable <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "output$PCAloadingsTable\n")
  }
  pca <- pcaReact()
  
  if (is.null(pca)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/PCAloadingsTable.RData",
         list = c( ls())
    )
  }
  # load("~/SCHNAPPsDebug/PCAloadingsTable.RData")
  as.data.frame(pca$rotation)
})



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
    save(file = "~/SCHNAPPsDebug/selectedGenesTable.RData",
         list = c("normaliztionParameters", ls())
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
    save(file = "~/SCHNAPPsDebug/removedGenesTable.RData",
         list = c("normaliztionParameters", ls())
    )
  }
  # load("~/SCHNAPPsDebug/removedGenesTable.RData")
  
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
      is.null(dataTables) ) {
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
  
  setRedGreenButton(
    vars = list(
      c("selectIds", isolate(input$selectIds)),
      c("geneListSelection", isolate(input$geneListSelection)),
      c("minGenesGS", isolate(input$minGenesGS)),
      c("genesKeep", isolate(input$genesKeep))
    ),
    button = "updateGeneSelectionParameters"
  )
  
  
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
    # deepDebug()
    # change names to be hopefully a bit more clear
    changed <- FALSE # trace if something changed
    keepGenes <- useGenes
    keepCells <- useCells
    scEx <- assays(scExOrg)[[1]]
    
    # overall gene expression Min (minGenesGS)
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
        id = "scExFunc2",
        duration = NULL
      )
      return(NULL)
    }
    # cells with same expression pattern
    # save(file = "~/SCHNAPPsDebug/scEx3.RData", list = c(ls()))
    # load(file = "~/SCHNAPPsDebug/scEx3.RData")
    # which(duplicated(as.matrix(assays(scExOrg[keepGenes, keepCells])[[1]])))
    # which(duplicated(as.matrix(assays(scEx)[[1]])))
    #
    # sum(duplicated(t(as.matrix(assays(scExOrg[keepGenes, keepCells])[[1]]))))
    # if something changed, check that it doesn't change again
    # TODO something is fishy here
    scExNew <- scExOrg[keepGenes, keepCells]
    if (changed) {
      if (DEBUG) cat(file = stderr(), "--- scExFunc changed\n")
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
    
    # remove potentially existing projections (PC, tTSNE, ...)
    knownProj = c("UMI.count", "Feature.count", "before.filter", "tsne1", "tsne2", "tsne3", "dbCluster", "all", "none" )
    rmCols = which(colnames(pD) %in% knownProj)
    rmCols = c(rmCols, which(base::grepl("^PC_[[:digit:]]+$", colnames(pD))))
    if(length(rmCols) > 0) {
      pD = pD[,-rmCols]
    }
    colData(scExNew) <- pD
    
    return(scExNew)
  }

# scEx ----
# apply filters that depend on genes & cells
# it is here that useCells and useGenes are combined and applied to select for
scEx <- reactive({
  deepDebug()
  if (DEBUG) {
    cat(file = stderr(), "scEx started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scEx")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scEx")
      # removeNotification(id = "scExError")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scEx", id = "scEx", duration = NULL)
  }
  
  dataTables <- inputData()
  useCells <- useCells()
  useGenes <- useGenes()
  cellSelectionValues <- cellSelectionValues()
  geneSelectionValues <- geneSelectionValues()
  minGene <- geneSelectionValues$minGenesGS # min number of reads per gene
  minG <- cellSelectionValues$minGenes # min number of reads per cell
  maxG <- cellSelectionValues$maxGenes # max number of reads per cell
  
  # browser()
  if (!exists("dataTables") |
      is.null(dataTables) | is.null(useGenes) | is.null(useCells)) {
    if (DEBUG) {
      cat(file = stderr(), "scEx: NULL\n")
    }
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scEx.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/scEx.RData")
  
  # TODO:??? should be keep the cell names from the projections?
  # here we are just reinitializing.
  retVal <- tryCatch({scExFunc(
    scExOrg = dataTables$scEx,
    useCells = useCells,
    useGenes = useGenes,
    minGene = minGene,
    minG = minG,
    maxG = maxG
  )}, error = function(e) {return(NULL)})
  scEx = retVal
  # ensure that we take the counts slot
  # when reloading additional slots might be loaded as well.
  # deepDebug()
  if(!is.null(scEx)){
    assays(scEx) = list(counts = assays(scEx)[["counts"]])
  }
  
  if (min(dim(scEx)) <100) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("Less than 100 cells or genes", id = "scExError",type = "error", duration = NULL)
    }
    return(NULL)
  }
  
  add2history(type = "save", input = isolate( reactiveValuesToList(input)), comment = "scEx", scEx = retVal)
  exportTestValues(scEx = {
    list(rowData(retVal), colData(retVal))
  })
  return(retVal)
})

## scEx_Hash ----
scEx_Hash <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "scEx_Hash started.\n")
  }
  scEx <- scEx()
  require(digest)
  if (is.null(scEx)) {
    return(NULL)
  }
  hash = sha1(as.matrix(assays(scEx)[[1]]))
  if (DEBUG) {
    cat(file = stderr(), "scEx_Hash ended\n")
  }
  return(hash)
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
  
  shinyjs::addClass("updateNormalization", "green")
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
  # deepDebug()
  scEx <- scEx()
  dataTables <- inputData()
  whichscLog <- input$whichscLog # from the input page
  # update if button is clicked is haneled in individual normalizations reactives/observers
  # update <- input$updateNormalization
  # don't update if parameters are changed
  normMethod <- isolate(input$normalizationRadioButton)
  clicked <- input$updateNormalization
  
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "scEx_log:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scEx_log.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/scEx_log.RData")
  
  if (whichscLog == "disablescEx_log") {
    return(NULL)
  }
  
  if (whichscLog == "useLog" & !is.null(dataTables$scEx_log)) {
    scEx_log <- dataTables$scEx_log[rownames(scEx), colnames(scEx)]
  } else {
    scEx_log <- do.call(normMethod, args = list())
  }
  if (is.null(scEx_log)) {
    # problem with normalization
    return(NULL)
  }
  .schnappsEnv$calculated_normalizationRadioButton <- normMethod
  
  add2history(type = "save", input = isolate( reactiveValuesToList(input)), comment = "scEx_log", scEx_log = scEx_log)
  
  exportTestValues(scEx_log = {
    assays(scEx_log)["logcounts"]
  })
  return(scEx_log)
})

## scEx_log_Hash ----
scEx_log_Hash <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "scEx_log_Hash started.\n")
  }
  scEx_log <- scEx_log()
  require(digest)
  if (is.null(scEx_log)) {
    return(NULL)
  }
  hash = sha1(as.matrix(assays(scEx_log)[[1]]))
  if (DEBUG) {
    cat(file = stderr(), "scEx_log_Hash ended\n")
  }
  return(hash)
})
# scEx_log_sha <- reactive({
#   scEx_log <- scEx_log()
#   require(digest)
#   if (is.null(scEx_log)) {
#     return(NULL)
#   }
#   return(sha1(as.matrix(assays(scEx_log)[[1]])))
# })
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
  if (is.null(scEx) | is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "scExLogMatrixDisplay:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scExLogMatrixDisplay.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/scExLogMatrixDisplay.RData")
  
  # TODO
  if (dim(scEx)[1] > 20000) {
    
  }
  retVal <- NULL
  if (!is.null(scEx_log)) {
    retVal <-
      data.frame(
        symbol = make.names(rowData(scEx_log)$symbol, unique = TRUE),
        stringsAsFactors = FALSE
      )
    retVal <- cbind(
      retVal,
      as.matrix(assays(scEx_log)[[1]])
    )
    # rownames(retVal) <-
    #   make.names(rowData(scEx)$symbol, unique = TRUE)
  }
  rownames(retVal) <-
    retVal$symbol
  
  return(retVal)
})

# pcaFunc ----
#' Principal Component Analysis (PCA) Function
#'
#' This function performs Principal Component Analysis (PCA) on single-cell RNA-seq data.
#'
#' @param scEx A SingleCellExperiment object containing the raw expression data.
#' @param scEx_log A SingleCellExperiment object containing the log-transformed expression data.
#' @param rank The number of principal components to compute.
#' @param center Logical, indicating whether to center the data.
#' @param scale Logical, indicating whether to scale the data.
#' @param useSeuratPCA Logical, indicating whether to use Seurat's PCA method. If FALSE, scran's PCA method will be used.
#' @param pcaGenes A character vector of gene names to include in the analysis.
#' @param rmGenes A character vector of gene names to exclude from the analysis.
#' @param featureData A DataFrame containing gene-related metadata.
#' @param pcaN The maximum number of genes to use for PCA.
#' @param maxGenes The maximum number of genes to consider in variable gene selection.
#' @param hvgSelection The method for selecting highly variable genes ("vst", "mvp", "disp", or "getTopHVGs").
#' @param inputNormalization The method of input data normalization ("DE_something", "DE_seuratSCTnorm", etc.).
#'
#' @return A list containing the PCA results, including 'x' (PCA coordinates), 'var_pcs' (standard deviations), and 'rotation' (loadings).
#'
#' @details This function performs PCA on single-cell RNA-seq data using either Seurat's PCA method or scran's PCA method, depending on the value of `useSeuratPCA`. It allows you to specify various parameters for PCA analysis and variable gene selection.
#'
#' @seealso \code{\link{SingleCellExperiment}}, \code{\link{FindVariableFeatures}}, \code{\link{ScaleData}}, \code{\link{RunPCA}}, \code{\link{runPCA}}
#'
#' @examples
#' \dontrun{
#' # Load required libraries and create mock data
#' library(SingleCellExperiment)
#' scEx <- SingleCellExperiment(assays = list(logcounts = matrix(rnorm(100), nrow = 10)))
#' scEx_log <- scEx
#' rank <- 3
#' center <- FALSE
#' scale <- FALSE
#' useSeuratPCA <- TRUE
#' pcaGenes <- colnames(assays(scEx)[["logcounts"]])
#' rmGenes <- c()
#' featureData <- rowData(scEx)
#' pcaN <- 10
#' maxGenes <- 1000
#' hvgSelection <- "vst"
#' inputNormalization <- "DE_something"
#'
#' # Perform PCA analysis
#' result <- pcaFunc(
#'   scEx,
#'   scEx_log,
#'   rank,
#'   center,
#'   scale,
#'   useSeuratPCA,
#'   pcaGenes,
#'   rmGenes,
#'   featureData,
#'   pcaN,
#'   maxGenes,
#'   hvgSelection,
#'   inputNormalization
#' )
#'
#' # View PCA results
#' str(result)
#' }
#'
#' @export

pcaFunc <- function(scEx, scEx_log, 
                    rank, center, scale,
                    cale, useSeuratPCA, 
                    pcaGenes, rmGenes=c(), 
                    featureData, pcaN, 
                    maxGenes = 1000, hvgSelection, 
                    inputNormalization="DE_something") {
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
    removeNotification(id = "pcawarning2")
    removeNotification(id = "pcawarning3")
  }
  # browser()
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/pcaFunc.RData", list = c(ls()))
  }
  # cp =load(file="~/SCHNAPPsDebug/pcaFunc.RData")
  
  # if there are no genes provided take all genes
  genesin <- geneName2Index(pcaGenes, featureData)
  genesout <- geneName2Index(rmGenes, featureData)
  # if (is.null(genesin) || length(genesin) == 0) {
  #   genesin <- rownames(scEx_log)
  # }
  geneshvg = c()
  if (inputNormalization == "DE_seuratSCTnorm") {
    geneshvg = rownames(assays(scEx_log)[[1]])
    
  } else {
    tryCatch({
      
      switch(hvgSelection,
             "vst" = {
               pbmc <-
                 FindVariableFeatures(assays(scEx)[[1]], selection.method = "vst")
               
               geneshvg = rownames(pbmc)[order(pbmc$vst.variance.standardized, decreasing = T)[1:pcaN]]
             },
             "mvp" = {
               pbmc <-
                 FindVariableFeatures(assays(scEx_log)[[1]], selection.method = "mean.var.plot")
               pbmc$mvp.dispersion.scaled[which(is.na(pbmc$mvp.dispersion.scaled))] = max(pbmc$mvp.dispersion.scaled, na.rm = T) + 1
               geneshvg = rownames(pbmc)[order(pbmc$mvp.dispersion.scaled, decreasing = T)[1:pcaN]]
             },
             "disp" = {
               pbmc <-
                 FindVariableFeatures(assays(scEx_log)[[1]], selection.method = "dispersion", nfeatures = pcaN)
               # NA values seem to have the highest variance. so we make sure to include them
               pbmc$mvp.dispersion.scaled[which(is.na(pbmc$mvp.dispersion.scaled))] = max(pbmc$mvp.dispersion.scaled, na.rm = T) + 1
               geneshvg = rownames(pbmc)[order(pbmc$mvp.dispersion.scaled, decreasing = T)[1:pcaN]]
             },
             "getTopHVGs" = {
               # apply(as.matrix(assays(scEx_log)[[1]]), 1, max)
               scEx = logNormCounts(scEx)
               dec <- scran::modelGeneVar(scEx, assay.type = "logcounts")
               geneshvg <- getTopHVGs(dec, prop=0.1)
             }
      )
    }, error = function(e){
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification(
          paste("Problem with PCA: Find variables not working ", e),
          type = "error",
          id = "pcawarning",
          duration = NULL
        )
      }
    })
  }
  # We insist that the genes we select manually are included and then take 
  
  genesin = c(genesin, geneshvg, rownames(scEx_log))[1:pcaN]
  
  if (length(genesin) < rank) {
    rank <- length(genesin)
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        paste("Problem with PCA: Too many PCs requested, setting to ", rank),
        type = "error",
        id = "pcawarning2",
        duration = NULL
      )
    }
  }
  if (rank < 3 ) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        paste("Problem with PCA: not enough genes, aboarding!"),
        type = "error",
        id = "pcawarning3",
        duration = NULL
      )
    }
    return(NULL)
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
    # not sure, but this works on another with TsparseMatrix
    if (!is(assays(scEx_log)[["logcounts"]], "CsparseMatrix")) {
      assays(scEx_log)[["logcounts"]] <-
        as(assays(scEx_log)[["logcounts"]], "CsparseMatrix")
    }
    # x <- assays(scEx_log)[["logcounts"]]
    genesin = genesin[genesin %in% rownames(scEx_log)]
    # x <- x[genesin, , drop = FALSE]
    mi = BPCells::write_matrix_memory(assays(scEx_log)[["logcounts"]][genesin,], compress = F)
    if (scale | center) {
      set.seed(1)
      # parallel
      # plan("multiprocess", workers = 4)
      mi <- Seurat::ScaleData(mi, do.scale = scale, do.center = center, verbose = T)
      genesin = rownames(mi) # ScaleData can remove genes.
    }
    if (useSeuratPCA){
      # Seurat:
      reductObj = RunPCA(mi, npcs = rank, verbose = FALSE, assay = "RNA")
      pca = list()
      pca$rotation = Loadings(reductObj)
      pca$x = Embeddings(reductObj)
      pca$sdev = Stdev(reductObj)
    } else {
      # # scran:
      x <- assays(scEx_log)[["logcounts"]][genesin,]
      x <- t(x)
      pca <- runPCA(x, rank=rank, get.rotation=TRUE)
    }

    pca
  })
  
  if (is.null(scaterPCA)) {
    return(NULL)
  }
  # pca = reducedDim(scaterPCA, "PCA")
  # attr(pca,"percentVar")
  #
  # rownames(scaterPCA$x) = colnames(scEx_log)
  return(list(
    x =  scaterPCA$x,
    # x = SingleCellExperiment::reducedDim(scaterPCA, "PCA"),
    var_pcs = scaterPCA$sdev,
    rotation = scaterPCA$rotation
    # var_pcs = attr(
    #   SingleCellExperiment::reducedDim(scaterPCA, "PCA"),
    #   "percentVar"
    # )
  ))
}
# pca ----
runPCAclicked <- reactive({
  runme = input$updatePCAParameters
  if(runme>0) return(1)
  return(0)
})
pcaReact <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "pca started.\n")
  }
  nWorkers = nbrOfWorkers()
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
  # browser()
  scEx <- scEx()
  scEx_log <- scEx_log()
  # only redo calculations if button is pressed.
  input$updatePCAParameters
  rank <- isolate(input$pcaRank)
  pcaN <- isolate(input$pcaN)
  center <- isolate(input$pcaCenter)
  scale <- isolate(input$pcaScale)
  pcaGenes <- isolate(input$genes4PCA)
  rmGenes <- isolate(input$genesRMPCA)
  hvgSelection <- isolate(input$hvgSelection)
  useSeuratPCA <- isolate(input$useSeuratPCA)
  inputNormalization <- isolate(input$normalizationRadioButton)
  
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "pcaReact: scEx_log:NULL\n")
    }
    return(NULL)
  }
  # deepDebug()
  if (is.null(rank)) rank <- 10
  if (is.null(center)) centr <- TRUE
  if (is.null(scale)) scale <- FALSE
  featureData <- rowData(scEx_log)
  retVal <- pcaFunc(scEx = scEx, scEx_log = scEx_log, 
                    rank = rank, center = center, 
                    scale = scale, useSeuratPCA = useSeuratPCA, 
                    pcaGenes = pcaGenes, rmGenes = rmGenes, featureData = featureData,
                    pcaN =  pcaN, 
                    hvgSelection = hvgSelection,
                    inputNormalization = inputNormalization)
  
  setRedGreenButton(
    vars = list(
      c("pcaRank", rank),
      c("pcaN", pcaN),
      c("pcaCenter", center),
      c("pcaScale", scale),
      c("hvgSelection", hvgSelection),
      c("genes4PCA", pcaGenes)
    ),
    button = "updatePCAParameters"
  )
  
  printTimeEnd(start.time, "pca")
  exportTestValues(pca = {
    retVal
  })
  return(retVal)
})
# %>%
#   bindCache(scEx(),
#             scEx_log(),
#             runPCAclicked(),
#             isolate(input$pcaRank),
#             isolate(input$pcaN),
#             isolate(input$pcaCenter),
#             isolate(input$pcaScale),
#             isolate(input$genes4PCA),
#             isolate(input$genesRMPCA),
#             isolate(input$hvgSelection),
#             isolate(input$useSeuratPCA),
#             isolate(input$normalizationRadioButton)
#   )

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
    save(file = "~/SCHNAPPsDebug/scranCluster.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/scranCluster.RData")
  
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
    min.mean = 0.1
    # ,
    # BPPARAM = MulticoreParam(
    #   workers = ifelse(detectCores() > 1, detectCores() - 1, 1)
    # ),
    # block.BPPARAM = MulticoreParam(
    #   workers = ifelse(detectCores() > 1, detectCores() - 1, 1)
    # )
  )
  switch(clusterSource,
         # "PCA" = {
         #   params$x <- scEx
         #   reducedDims(scEx) <- SimpleList(PCA = pca$x)
         #   params$assay.type <- "counts"
         # },
         "logcounts" = {
           reducedDims(scEx_log) <- SimpleList(PCA = pca$x)
           params$x <- scEx_log
           params$assay.type <- "logcounts"
           if (length(geneid) > 0) {
             params$subset.row <- geneid
           }
           use.ranks <- NA
         },
         "counts" = {
           reducedDims(scEx) <- SimpleList(PCA = pca$x)
           params$x <- scEx
           params$assay.type <- "counts"
           if (length(geneid) > 0) {
             params$subset.row <- geneid
           }
           use.ranks <- NA
         }
  )
  require(scran)
  retVal <- tryCatch(
    {
      # parallel (BSPARAM)
      suppressMessages(BiocGenerics::do.call("quickCluster", params))
    },
    error = function(e) {
      cat(file = stderr(), paste("\nProblem with clustering:\n\n", as.character(e), "\n\n"))
      # if (grepl("rank variances of zero", as.character(e))) {
      if (useRanks) {
        if (!is.null(getDefaultReactiveDomain())) {
          showNotification(paste("Not using ranks due to identical ranks\n"), id = "quickClusterError", type = "error", duration = NULL)
        }
        cat(file = stderr(), paste("\nNot using ranks due to identical ranks\n", e))
        params$use.ranks <- F
        save(file = "~/SCHNAPPsDebug/scranClusterError1.RData", list = c(ls()))
        # cp = load(file="~/SCHNAPPsDebug/scranClusterError1.RData")
        return(do.call("quickCluster", params))
      }
      cat(file = stderr(), "\nRemove duplicated cells\n")
      # }
      return(NULL)
    },
    warning = function(e) {
      if (DEBUG) cat(file = stderr(), paste("\nclustering produced Warning:\n", e, "\n"))
      save(file = "~/SCHNAPPsDebug/scranClusterError.RData", list = c(ls()))
      # cp = load(file="~/SCHNAPPsDebug/scranClusterError.RData")
      return(suppressMessages(do.call("quickCluster", params)))
    }
  )
  if ("barcode" %in% colnames(colData(scEx_log))) {
    barCode <- colData(scEx_log)$barcode
  } else {
    barCode <- rownames(colData(scEx_log))
  }
  if (is.null(retVal)){
    retVal = rep(0,length(barCode))
  }
  retVal <- data.frame(
    Barcode = barCode,
    Cluster = retVal
  )
  rownames(retVal) <- rownames(colData(scEx_log))
  return(retVal)
}

clusterBootstrapReactive <- reactive({
  
})
# memoise functionality 
#   scranCluster_m <- memoise::memoise(scranCluster,cache=do.call(cachem::cache_disk,.schnappsEnv$cacheDir))
scranCluster_m <- scranCluster


scran_Cluster <- function(){
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
  
  # react to the following changes
  # input$updateClusteringParameters
  scEx <- scEx()
  scEx_log <- scEx_log()
  pca <- pcaReact()
  
  # ignore these changes
  seed <- isolate(input$seed)
  useRanks <- isolate(input$useRanks)
  clusterSource <- isolate(clusterMethodReact$clusterSource)
  geneSelectionClustering <- isolate(input$geneSelectionClustering)
  minClusterSize <- isolate(input$minClusterSize)
  clusterMethod <- isolate(clusterMethodReact$clusterMethod)
  tabsetCluster = isolate(input$tabsetCluster)
  
  if (is.null(pca) | is.null(scEx_log) | is.na(minClusterSize) | tabsetCluster != "scran_Cluster") {
    if (DEBUG) {
      cat(file = stderr(), "scran_Cluster:NULL\n")
    }
    return(NULL)
  }
  # if (.schnappsEnv$DEBUGSAVE) {
  #   save(file = "~/SCHNAPPsDebug/scran_Cluster.RData", list = c(ls()))
  # }
  # load(file="~/SCHNAPPsDebug/scran_Cluster.RData")
  
  featureData <- rowData(scEx_log)
  
  if (is.null(seed)) {
    seed <- 1
  }
  retVal <- scranCluster_m(
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
    save(file = "~/SCHNAPPsDebug/scran_ClusterError.RData", list = c(ls()))
    return(NULL)
    # stop("error: clustering didn't produce a result")
  }
  setRedGreenButton(
    vars = list(
      c("seed", isolate(input$seed)),
      c("useRanks", isolate(input$useRanks)),
      c("clusterSource", isolate(clusterMethodReact$clusterSource)),
      c("geneSelectionClustering", isolate(input$geneSelectionClustering)),
      c("minClusterSize", isolate(input$minClusterSize)),
      c("tabsetCluster", isolate(input$tabsetCluster)),
      c("clusterMethod", isolate(clusterMethodReact$clusterMethod))
    ),
    button = "updateClusteringParameters"
  )
  
  exportTestValues(scran_Cluster = {
    retVal
  })
  return(retVal)
}

BPCellsLog <- reactive({
  scEx_log = scEx_log()
  req(scEx_log)
  # browser()
  mi = NULL
  mi = tryCatch({
    if(!is(assay(scEx_log, "logcounts"),"CsparseMatrix")){
      m = as(assay(scEx_log, "logcounts"),"CsparseMatrix")
    }else{
      m = assay(scEx_log, "logcounts")
    }
    BPCells::write_matrix_memory(m, compress = F)
  }, error = function(e) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("Problem BPCellsLog", type = "warning", duration = NULL)
    }
    return(NULL)
  }
  )
  return(mi)
})


BPCellsCounts <- reactive({
  scEx = scEx()
  req(scEx)
  if(!is(assay(scEx, "counts"),"CsparseMatrix")){
    m = as(assay(scEx, "counts"),"CsparseMatrix")
  }else{
    m = assay(scEx, "counts")
  }
  mi = BPCells::write_matrix_memory(m, compress = F)
  return(mi)
})


# Seurat clustering ----
runSeuratClustering <- function(scEx, meta.data, dims, pca, k.param, resolution) {
  # creates object @assays$RNA@data and @assays$RNA@counts
  seurDat <- NULL
  if(is(scEx,"SingleCellExperiment")){
    seurDat <- CreateSeuratObject(
      counts = as(assays(scEx)[[1]], "CsparseMatrix"),
      meta.data = meta.data
    )
  }else{
    # if(is(scEx,"BPCells")){
    if(is(scEx,"UnpackedMatrixMem_double")){
      seurDat <- CreateSeuratObject(
        counts = scEx,
        meta.data = meta.data
      )
    } else{
      cat(file = stderr(), "Neither SingleCellExperiment nor UnpackedMatrixMem_double. \n")
      stop()
    }
  }
  # bowser()
  # we remove e.g. "genes" from total seq (CD3-TotalSeqB)
  useGenes = which(rownames(seurDat[["RNA"]]$counts) %in% rownames(scEx))
  # seurDat[["RNA"]]$counts = as(assays(scEx)[[1]], "CsparseMatrix")[useGenes,]
  seurDat = seurDat[useGenes,]
  dims = min(dims,ncol(pca$x)) 
  
  seurDat[["pca"]] = CreateDimReducObject(embeddings = pca$x[colnames(seurDat),], 
                                          loadings = pca$rotation, 
                                          stdev = pca$var_pcs, 
                                          key = "PC_", 
                                          assay = "RNA")
  
  
  seurDat = FindNeighbors(seurDat, dims = 1:dims, k.param = k.param)
  # parallel
  # plan("multiprocess", workers = 4)
  seurDat <- FindClusters(seurDat, resolution = resolution, method = "igraph", algorithm=1 )
  retVal = data.frame(Barcode = colnames(seurDat),
                      Cluster = Idents(seurDat))
}

# cacheDir is not known before and messes up things
# if(!is.null(.schnappsEnv$cacheDir)){
#   cat(file = stderr(), unlist(.schnappsEnv$cacheDir))
#   # heatmapModuleFunction_m = memoise::memoise(heatmapModuleFunction,cache=do.call(cachem::cache_disk,.schnappsEnv$cacheDir)) is called too often
#   heatmapModuleFunction_m = heatmapModuleFunction
#   runSeuratClustering_m <- memoise::memoise(runSeuratClustering,cache=do.call(cachem::cache_disk,.schnappsEnv$cacheDir))
#   panelPlotFunc_m = memoise::memoise(panelPlotFunc,cache=do.call(cachem::cache_disk,.schnappsEnv$cacheDir))
#   scranCluster_m <- memoise::memoise(scranCluster,cache=do.call(cachem::cache_disk,.schnappsEnv$cacheDir))

#   
#   }
# } else {
runSeuratClustering_m = runSeuratClustering

# 
seurat_Clustering <- function() {
  if (DEBUG) {
    cat(file = stderr(), "seurat_Clustering started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "seurat_Clustering")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "seurat_Clustering")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("seurat_Clustering", id = "seurat_Clustering", duration = NULL)
  }
  # browser()
  scEx = scEx() # need to be run when updated
  scExMat <- BPCellsCounts()
  scEx_log = scEx_log()
  pca = pcaReact()
  tabsetCluster = isolate(input$tabsetCluster)
  dims = isolate(input$seurClustDims)
  k.param = isolate(input$seurClustk.param)
  resolution = isolate(input$seurClustresolution)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/seurat_Clustering.RData", list = c(ls()))
  }
  # cp =load(file="~/SCHNAPPsDebug/seurat_Clustering.RData")
  if(any(c(is.null(k.param), is.null(tabsetCluster), is.null(dims), is.null(resolution)))) {
    return(NULL)
  }
  if (tabsetCluster != "seurat_Clustering" | is.null(pca)){
    return(NULL)
  }
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "seurat_Clustering: scEx_log:NULL\n")
    }
    return(NULL)
  }
  cellMeta <- colData(scEx)
  rData <- rowData(scEx)
  meta.data <- cellMeta[, "sampleNames", drop = FALSE]
  # browser()
  retVal = NULL
  retVal <- tryCatch({
    runSeuratClustering_m(scExMat, meta.data, dims, pca, k.param, resolution)
  },
  error = function(e) {
    cat(file = stderr(), paste("\n\n!!!Error during Seurat clustering:\ncheck console", e, "\n\n"))
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("ERROR in seurat_Clustering", id = "seurat_ClusteringERROR", duration = NULL, type = "error")
    }
    return(NULL)
  }
  )
  # runSeuratClustering(scEx, meta.data, dims, pca, k.param, resolution)
  
  if(is.null(retVal)){
    return(NULL)
  }
  
  setRedGreenButton(
    vars = list(
      c("seurClustDims", isolate(input$seurClustDims)),
      c("seurClustk.param", isolate(input$seurClustk.param)),
      c("seurClustresolution", isolate(input$seurClustresolution)),
      c("tabsetCluster", isolate(input$tabsetCluster))
    ),
    button = "updateClusteringParameters"
  )
  
  return(retVal)
} 

# snnGraph ----
snnGraph <- function(){
  if (DEBUG) {
    cat(file = stderr(), "snnGraph started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "snnGraph")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "snnGraph")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("snnGraph", id = "snnGraph", duration = NULL)
  }
  
  scEx = scEx() # need to be run when updated
  scEx_log = scEx_log()
  pca = pcaReact()
  type = isolate(input$snnType)
  # clusterSource = isolate(input$snnClusterSource)
  tabsetCluster = isolate(input$tabsetCluster)
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "snnGraph scEx_log:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/snnGraph_Clustering.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/snnGraph_Clustering.RData")
  
  if (tabsetCluster != "snnGraph" | is.null(pca)){
    return(NULL)
  }
  
  # switch(clusterSource,
  #        # "PCA" = {
  #        #   params$x <- scEx
  #        #   reducedDims(scEx) <- SimpleList(PCA = pca$x)
  #        #   params$assay.type <- "counts"
  #        # },
  #        "logcounts" = {
  reducedDims(scEx_log) <- SimpleList(PCA = pca$x)
  sce = scEx_log
  assay.type <- "logcounts"
  #        },
  #        "counts" = {
  #          reducedDims(scEx) <- SimpleList(PCA = pca$x)
  #          sce = scEx
  #          assay.type <- "counts"
  #        }
  # )
  require(scran)
  retVal <- tryCatch(
    {
      g1 <- buildSNNGraph(sce, use.dimred="PCA", type = type)
      cluster <- factor(igraph::cluster_walktrap(g1)$membership)
    },
    error = function(e) {
      cat(file = stderr(), paste("\nProblem with SNNGRAPH:\n\n", as.character(e), "\n\n"))
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("snnClusterError", id = "snnClusterError", type = "error", duration = NULL)
      }
      
      return(NULL)
    },
    warning = function(e) {
      if (DEBUG) cat(file = stderr(), paste("\nSNN clustering produced Warning:\n", e, "\n"))
    }
  )
  if ("barcode" %in% colnames(colData(scEx_log))) {
    barCode <- colData(sce)$barcode
  } else {
    barCode <- rownames(colData(sce))
  }
  
  if(is.null(retVal)) {return(NULL)}
  
  
  retVal <- data.frame(
    Barcode = barCode,
    Cluster = retVal
  )
  rownames(retVal) <- rownames(colData(scEx_log))
  
  setRedGreenButton(
    vars = list(
      # c("snnClusterSource", isolate(input$snnClusterSource)),
      c("snnType", isolate(input$snnType)),
      c("tabsetCluster", isolate(input$tabsetCluster))
    ),
    button = "updateClusteringParameters"
  )
  return(retVal)
}



# simlrFunc ----
simlrFunc  <- function(){
  if (DEBUG) {
    cat(file = stderr(), "simlrFunc started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "simlrFunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "simlrFunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("simlrFunc", id = "simlrFunc", duration = NULL)
  }
  
  # scEx = scEx() # need to be run when updated
  scEx_log = scEx_log()
  # pca = pcaReact()
  nClust = isolate(input$simlr_nClust)
  maxClust = isolate(input$simlr_maxClust)
  # clusterSource = isolate(input$snnClusterSource)
  tabsetCluster = isolate(input$tabsetCluster)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/simlrFunc_Clustering.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/simlrFunc_Clustering.RData")
  
  if (tabsetCluster != "simlrFunc" | is.null(scEx_log)){
    return(NULL)
  }
  
  if(!"SIMLR" %in% rownames(installed.packages())) {
    return(NULL)
  }
  require(SIMLR)
  # why????
  if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
      Sys.info()["sysname"] == "Darwin" && getRversion() == "4.0.2") {
    parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  }
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "simlrFunc scEx_log:NULL\n")
    }
    return(NULL)
  }
  
  withWarnings <- function(expr) {
    wHandler <- function(w) {
      if (DEBUG) {
        cat(file = stderr(), paste("simlrFunc created a warning:", w, "\n"))
      }
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification(
          "Problem with simlrFunc?",
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
  
  retVal <- withWarnings({
    if (nClust == 0) {
      estimates = SIMLR_Estimate_Number_of_Clusters(X = as.matrix(assays(scEx_log)[[1]]),
                                                    NUMC=2:maxClust,
                                                    cores.ratio = 1)
      nClust = max(which.min(estimates$K1), which.min(estimates$K2))
    }
    sim = SIMLR(X = as.matrix(assays(scEx_log)[[1]]), 
                c = nClust, 
                cores.ratio = 1)
    
    cluster <- factor(sim$y$cluster)
    cluster
  })
  
  
  if (is.null(retVal)) {return(NULL)}
  
  
  retVal <- data.frame(
    Cluster = retVal
  )
  # rownames(retVal) = barCode
  rownames(retVal) <- rownames(colData(scEx_log))
  
  setRedGreenButton(
    vars = list(
      # c("snnClusterSource", isolate(input$snnClusterSource)),
      c("snnType", isolate(input$snnType)),
      c("tabsetCluster", isolate(input$tabsetCluster)),
      c("simlr_nClust", isolate(input$simlr_nClust)),
      c("simlr_maxClust", isolate(input$simlr_maxClust))
    ),
    button = "updateClusteringParameters"
  )
  return(retVal)
}




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
  
  clicked = input$updateClusteringParameters
  scEx = scEx() # need to be run when updated
  scEx_log = scEx_log()
  pca = pcaReact()
  tabsetCluster = isolate(input$tabsetCluster)
  
  
  # kNr <- input$kNr
  clustering <- do.call(tabsetCluster, args = list())
  
  if (.schnappsEnv$DEBUGSAVE) {
    prjCols = isolate( reactiveValuesToList(projectionColors))
    save(file = "~/SCHNAPPsDebug/dbCluster.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/dbCluster.RData")
  
  if (is.null(clustering)) {
    if (DEBUG) {
      cat(file = stderr(), "dbCluster: clustering NULL\n")
    }
    return(NULL)
  }
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "dbCluster scEx_log:NULL\n")
    }
    return(NULL)
  }
  dbCluster <- clustering$Cluster
  
  # cluster colors
  inCols <- list()
  lev <- levels(dbCluster)
  
  # repeat allowedColors to match cluster colors with repetion
  inCols <- allowedColors[rep(1:length(allowedColors),ceiling(length(lev) / length(allowedColors)))[1:length(lev)]]
  # inCols <- allowedColors[1:length(lev)]
  names(inCols) <- lev
  # browser()
  isolate(
    if(is.null(names(projectionColors$dbCluster)) | !all(levels(dbCluster) %in% names(projectionColors$dbCluster))){
      projectionColors$dbCluster <- unlist(inCols)
      add2history(type = "save", input=isolate( reactiveValuesToList(input)), comment = "projectionColors", projectionColors = projectionColors)
    }
  )
  # setRedGreenButton(
  #   vars = list(
  #     c("sampleNamecol", isolate(projectionColors$sampleNames)),
  #     c("clusterCols", isolate(projectionColors$dbCluster))
  #   ),
  #   button = "updateColors"
  # )
  exportTestValues(dbCluster = {
    dbCluster
  })
  return(dbCluster)
})

observe({
  pc = projectionColors
  cat(file = stderr(), paste("colornames", paste(names(pc), collapse = ", "), "\n"))
  cat(file = stderr(), paste("dbCluster colnames", paste(names(pc$dbCluster), collapse = ", "), "\n"))
})

observe({
  pc = projections()
  cat(file = stderr(), paste("projection names", paste(names(pc), collapse = ", "), "\n"))
  cat(file = stderr(), paste("projection dbCluster colnames", paste(levels(pc$dbCluster), collapse = ", "), "\n"))
})

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
  # .schnappsEnv$DEBUGSAVE = T
  
  scEx <- scEx()
  printTimeEnd(start.time, "projections scEx")
  pca <- pcaReact()
  printTimeEnd(start.time, "projections pca")
  # manually specified groups of cells (see 2D plot in moduleServer.R)
  prjs <- sessionProjections$prjs
  # derived/modified projections from projections tab
  newPrjs <- projectionsTable$newProjections
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/projections.bf.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/projections.bf.RData"); DEBUGSAVE=FALSE
  if (!exists("scEx") |
      is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "sampleInfo: NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/projections.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/projections.RData"); DEBUGSAVE=FALSE
  
  # save session specific projections for restore
  if (exists("historyPath", envir = .schnappsEnv)) {
    # if this variable is not set we are not saving
    # browser()
    tfile <- paste0(.schnappsEnv$historyPath, "/userProjections.RData")
    save(file = tfile, list = c("prjs", "newPrjs"))
  }
  printTimeEnd(start.time, "projections save1")
  
  # deepDebug()
  # browser()
  
  
  # todo colData() now returns a s4 object of class DataFrame
  # not sure what else is effected...
  pd <- as.data.frame(colData(scEx))
  if (ncol(pd) < 2) {
    cat(file = stderr(), "phenoData for scEx has less than 2 columns\n")
    return(NULL)
  }
  projections <- pd
  
  if (!is.null(pca)) {
    pcax = pca$x
    comColNames = colnames(projections) %in% colnames(pcax)
    if(any(comColNames)){
      colnames(projections)[comColNames] = paste0(colnames(projections)[comColNames], ".old")
    }
    if (!all(c(rownames(pcax) %in% rownames(projections) , 
               rownames(projections) %in% rownames(pcax)))){
      save(file = "~/SCHNAPPsDebug/projectionsError.RData", list = c(ls()))
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("projFactors", id = "projFactors", duration = NULL, type = "error")
      }
      return(NULL)
    }
    projections <- cbind(projections, pcax)
  }
  printTimeEnd(start.time, "projections pca2")
  
  withProgress(message = "Performing projections", value = 0, {
    n <- length(.schnappsEnv$projectionFunctions)
    iter <- 1
    for (proj in .schnappsEnv$projectionFunctions) {
      start.time1 <- Sys.time()
      if (!is.null(getDefaultReactiveDomain())) 
        incProgress(1 / n, detail = paste("Creating ", proj[1]))
      if (DEBUG) {
        cat(file = stderr(), paste("calculation projection (", iter,"):  ", proj[1], "\n"))
      }
      if (DEBUG) cat(file = stderr(), paste("projection: ", proj[2], "\n"))
      assign("tmp", eval(parse(text = paste0(proj[2], "()"))))
      #browser()
      if (.schnappsEnv$DEBUGSAVE) {
        save(file = paste0("~/SCHNAPPsDebug/projections.", iter, ".RData"),
             list = c("tmp")
        )
        iter <- iter + 1
      }
      # cp = load(file="~/SCHNAPPsDebug/projections.5.RData")
      # deepDebug()
      # TODO here, dbCluster is probably overwritten and appended a ".1"
      if (is(tmp, "data.frame")) {
        cn <- make.names(c(colnames(projections), colnames(tmp)), unique = TRUE)
      } else {
        cn <- make.names(c(colnames(projections), make.names(proj[1])), unique = TRUE)
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
        } else {
          if (!is.null(nrow(tmp))) {
            if (nrow(projections) == nrow(tmp)) {
              projections <- cbind(projections, tmp)
            }
          } else {
            save(file = "~/SCHNAPPsDebug/projectionsError.RData", list = c(ls()))
            stop("error: ", proj[1], "didn't produce a result, please send file ~/SCHNAPPsDebug/projectionsError.RData to bernd")
          }
        }
        # else {
        #   stop("error: ", proj[1], "didn't produce a result")
        # }
      }
      if (!length(colnames(projections)) == length(cn)) {
        save(file = "~/SCHNAPPsDebug/projectionsError2.RData", list = c(ls()))
        stop("error: ", proj[1], "didn't produce a result, please send file ~/SCHNAPPsDebug/projectionsError2.RData to bernd")
      }
      colnames(projections) <- cn
      if (DEBUG) cat(file = stderr(), paste("colnames ", paste0(colnames(projections), collapse = " "), "\n"))
      if (DEBUG) cat(file = stderr(), paste("observe this: ", proj[2], "\n"))
      # observe(proj[2], quoted = TRUE)
    }
  })
  printTimeEnd(start.time, "projections after for")
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/projections.af.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/projections.af.RData"); DEBUGSAVE=FALSE
  
  # add a column for gene specific information that will be filled/updated on demand
  # projections$UmiCountPerGenes <- 0
  # projections$UmiCountPerGenes2 <- 0
  for (pdIdx in colnames(pd)) {
    if (!pdIdx %in% colnames(projections)) {
      projections[, pdIdx] <- pd[, pdIdx]
    }
  }
  printTimeEnd(start.time, "projections gene specific")
  
  if (ncol(prjs) > 0 ) {
    errorPrjs = FALSE
    if(is.null(rownames(prjs))) {
      errorPrjs = TRUE
    } else {
      if (!all(rownames(projections) %in% rownames(prjs))) errorPrjs = TRUE
    }
    if(errorPrjs){
      cat(file = stderr(), paste("\n\n\nprjs ERROR\n\n\n\n"))
      save(file = "~/SCHNAPPsDebug/prjs.ERROR.RData", list = c(ls()))
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("prjs ERROR", id = "prjs", duration = NULL, type = "error")
      }
      stop("\n\n\nERROR: prjs didn't work correctly, please send file ~/SCHNAPPsDebug/prjs.ERROR.RData to bernd\n\n\n")
    }else{
      projections <- cbind(projections, prjs[rownames(projections),,drop=FALSE])
    }
  }
  printTimeEnd(start.time, "projections error?")
  
  # # remove columns with only one unique value
  # rmC <- c()
  # for (cIdx in 1:ncol(projections)) {
  #   # ignore sampleNames
  #   if (colnames(projections)[cIdx] == "sampleNames") next()
  #   if (length(unique(projections[, cIdx])) == 1) rmC <- c(rmC, cIdx)
  # }
  # if (length(rmC) > 0) projections <- projections[, -rmC]
  
  if (ncol(newPrjs) > 0) {
    projections <- cbind(projections, newPrjs[rownames(projections), , drop = FALSE])
  }
  # in case (no Normalization) no clusters or sample names have been assigned
  if (!"dbCluster" %in% colnames(projections)) {
    projections$dbCluster <- 0
    projections$dbCluster = as.factor(projections$dbCluster)
  }
  if (!"sampleNames" %in% colnames(projections)) {
    projections$sampleNames <- "1"
  }
  # TODO
  # this takes too long
  # UNCOMMENT
  printTimeEnd(start.time, "projections add history")
  add2history(type = "save", input=isolate( reactiveValuesToList(input)), comment = "projections", projections = projections)
  
  cat(file = stderr(), paste("\n\nscLog: ",isolate(input$whichscLog),"\n\n"))
  # add2history(type = "save", comment = "projections", projections = projections)
  
  exportTestValues(projections = {
    projections
  })
  # browser()
  # .schnappsEnv$DEBUGSAVE = F
  return(projections)
})

## projections_Hash ----
projections_Hash <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "projections_Hash started.\n")
  }
  projections <- projections()
  require(digest)
  if (is.null(projections)) {
    return(NULL)
  }
  hash = sha1(projections)
  if (DEBUG) {
    cat(file = stderr(), "projections_Hash started.\n")
  }
  return(hash)
})

# projFactors ----
# which projections are factors?
projFactors <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "projFactors started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "projFactors")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "projFactors")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("projFactors", id = "projFactors", duration = NULL)
  }
  
  projections <- projections()
  deepDebug()
  if (is.null(projections)) {
    choices <- c("no valid columns")
    return(choices)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/projFactors.RData", list = c(ls()))
  }
  # cp=load(file="~/SCHNAPPsDebug/projFactors.RData")
  
  coln <- colnames(projections)
  choices <- c()
  for (cn in coln) {
    if (length(unique(projections[, cn])) < 50) {
      choices <- c(choices, cn)
    }
  }
  if (is.null(choices)) {
    choices <- c("no valid columns")
  }
  if (length(choices) == 0) {
    choices <- c("no valid columns")
  }
  return(choices)
})


# initializeGroupNames ----
# TODO shouldn't this be an observer??? or just a function???
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
    if (DEBUG) cat(file = stderr(), "     scEx NULL.\n")
    return(NULL)
  }
  isolate({
    grpNs <- groupNames$namesDF
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/initializeGroupNames.RData", list = c(ls()))
    }
    # cp = load(file="~/SCHNAPPsDebug/initializeGroupNames.RData")
    # TODO ??? if cells have been removed it is possible that other cells that were excluded previously show up
    # this will invalidate all previous selections.
    # browser()
    if (rlang::is_empty(grpNs) | !all(colnames(scEx) %in% rownames(grpNs))) {
      df <-
        data.frame(
          all = rep("TRUE", dim(scEx)[2]),
          none = rep("FALSE", dim(scEx)[2])
        )
      rownames(df) <- colnames(scEx)
      groupNames$namesDF <- df
    } else {
      groupNames$namesDF = groupNames$namesDF[colnames(scEx),]
    }
  })
})

# since initializeGroupNames depends on scEx only this will be set when the org data is changed.
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
    save(file = "~/SCHNAPPsDebug/sample.RData", list = c(ls()))
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
    save(file = "~/SCHNAPPsDebug/geneCount.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/geneCount.RData")
  
  retVal <- Matrix::colSums(assays(scEx_log)[["logcounts"]] > 0)
  
  exportTestValues(geneCount = {
    retVal
  })
  return(retVal)
})

# beforeFilterCounts -----
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
  ipIDs = input$beforeFilterRegEx
  scEx = scEx()
  
  # ipIDs <- input$selectIds # regular expression of genes to be removed
  # ipIDs <- geneSelectionValues$selectIds
  if (!exists("dataTables") |
      is.null(dataTables) |
      length(dataTables$featuredata$symbol) == 0) {
    if (DEBUG) {
      cat(file = stderr(), "beforeFilterCounts: NULL\n")
    }
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/beforeFilterCounts.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/beforeFilterCounts.RData")
  
  geneIDs <- NULL
  if (nchar(ipIDs) > 0) {
    geneIDs <- grepl(ipIDs, dataTables$featuredata$symbol, ignore.case = TRUE)
  }
  if (is.null(geneIDs)) {
    retVal = rep(0, ncol(dataTables$scEx))
    names(retVal) = colnames(dataTables$scEx)
    return(retVal[colnames(scEx)])
  }
  retVal <- Matrix::colSums(assays(dataTables$scEx)[["counts"]][geneIDs, ])
  names(retVal) = colnames(dataTables$scEx)
  
  exportTestValues(beforeFilterCounts = {
    retVal[colnames(scEx)]
  })
  return(retVal[colnames(scEx)])
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
    save(file = "~/SCHNAPPsDebug/beforeFilterPrj.RData", list = c(ls()))
  }
  # cp =load(file="~/SCHNAPPsDebug/beforeFilterPrj.RData")
  if (is(bfc, "numeric")) {
    bfc = data.frame(before.filter = bfc)
    # rownames(bfc) = colnames(scEx)
  }
  cn <- colnames(scEx)
  retVal <- bfc[cn,]
  
  exportTestValues(beforeFilterPrj = {
    retVal
  })
  return(retVal)
})


# featureCount ----
featureCount <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "featureCount started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "featureCount")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "featureCount")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("featureCount", id = "featureCount", duration = NULL)
  }
  
  scEx <- scEx()
  
  if (is.null(scEx)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/featureCount.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/featureCount.RData")
  
  retVal <- Matrix::colSums(assays(scEx)[["counts"]]>0)
  
  exportTestValues(featureCount = {
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
    save(file = "~/SCHNAPPsDebug/umiCount.RData", list = c(ls()))
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
  factor(colData(scEx)$sampleNames)
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
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sampleInfo.RData", list = c(ls()))
  }
  # cp=load(file="~/SCHNAPPsDebug/sampleInfo.RData")
  if (!exists("scEx")) {
    if (DEBUG) {
      cat(file = stderr(), "sampleInfo: NULL\n")
    }
    return(NULL)
  }
  
  retVal <- sampleInfoFunc(scEx)
  
  exportTestValues(sampleInfo = {
    retVal
  })
  return(retVal)
})


# table of input cells with sample information ----
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
  
  scEx <- scEx()
  
  if (is.null(scEx)) {
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("inputSample", id = "inputSample", duration = NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/inputSample.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/inputSample.RData')
  
  # TODO should come from sampleInfo
  sampInf <- gsub(".*-(.*)", "\\1", scEx$barcode)
  cellIds <- data.frame(
    cellName = colnames(scEx),
    sample = sampInf,
    # number of genes per cell
    ngenes = Matrix::colSums(assays(scEx)[[1]]>0),
    nUMI = Matrix::colSums(assays(scEx)[[1]])
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

inputSampleOrg <- reactive({
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
    save(file = "~/SCHNAPPsDebug/inputSample.RData", list = c(ls()))
  }
  # load(file='~/SCHNAPPsDebug/inputSample.RData')
  
  # TODO should come from sampleInfo
  sampInf <- gsub(".*-(.*)", "\\1", dataTables$scEx$barcode)
  cellIds <- data.frame(
    cellName = colnames(dataTables$scEx),
    sample = sampInf,
    # number of genes per cell
    ngenes = Matrix::colSums(assays(dataTables$scEx)[[1]]>0),
    nUMI = Matrix::colSums(assays(dataTables$scEx)[[1]])
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
  #
  if ("pryr" %in% rownames(installed.packages())) {
    suppressMessages(require(pryr))
  } else {
    mem_used <- function() {
      showNotification("Please install pryr", id = "noPryr", type = "error", duration = NULL)
      return(0)
    }
  }
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
  if (!"callr" %in% rownames(installed.packages())) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("please install 'callr' to enable reports",
                       id = "noCallR", duration = NULL, type = "error"
      )
    }
    return(NULL)
  }
  scEx <- scEx()
  projections <- projections()
  scEx_log <- scEx_log()
  inputNames <- names(input)
  # deepDebug()
  
  if (is.null(scEx) | is.null(scEx_log)) {
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
  
  report.env <- getReactEnv(DEBUG = DEBUG)
  pca <- pcaReact()
  tsne <- tsne()
  
  scEx <- consolidateScEx(scEx, projections, scEx_log, pca, tsne)
  reportTempDir <- get("reportTempDir", envir = .schnappsEnv)
  base::save(file = tmpPrjFile,
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
  # deepDebug()
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/tempReport.1.RData",
         list = c("session", "report.env", "file", ls())
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
      function(input, output_file, params, envir) {
        rmarkdown::render(
          input = input,
          output_file = output_file,
          params = params,
          envir = envir
        )
      },
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
  
  tDir <- paste0(.schnappsEnv$reportTempDir, .Platform$file.sep)
  base::file.copy(tmpPrjFile, paste0(.schnappsEnv$reportTempDir, "/sessionData.RData"))
  write.csv(as.matrix(assays(scEx_log)[[1]]),
            file = paste0(.schnappsEnv$reportTempDir, "/normalizedCounts.csv")
  )
  base::save(file = paste0(.schnappsEnv$reportTempDir, "/inputUsed.RData"),
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

# gmtData ----
# load RData file with singlecellExperiment object
# internal, should not be used by plug-ins

## combine user defined and file defined GMT lists ----
observe({
  gmtfDat = gmtFileData()
  gmtuDat = gmtUserData() 
  
  if(is.null(gmtfDat) & is.null(gmtuDat)){
    return(NULL)
  }
  retVal = append(gmtfDat, gmtuDat) %>% compact
  gmtData(retVal)
  gmt = gmtData()
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/gmtFileDataObs.RData", list = c(ls()))
  }
  # cp =load(file='~/SCHNAPPsDebug/gmtFileDataObs.RData')
  
})


## load GMT file ----
gmtFileData <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "gmtData started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gmtData")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gmtData")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gmtData", id = "gmtData", duration = NULL)
  }
  
  gmtFile <- input$geneSetFile
  scEx = scEx() # needed to validate input
  
  if (is.null(gmtFile) |is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "gmtData: NULL\n")
    }
    return(NULL)
  }
  
  if (DEBUG) cat(file = stderr(), paste("gmtFile.", gmtFile$datapath[1], "\n"))
  if (!file.exists(gmtFile$datapath[1])) {
    if (DEBUG) {
      cat(file = stderr(), "gmtData: ", gmtFile$datapath[1], " doesn't exist\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/gmtFileData.RData", list = c(ls()))
  }
  # cp =load(file='~/SCHNAPPsDebug/gmtFileData.RData')
  
  # gmtFile$datapath = "data/sessionData.RData"
  # annFile$datapath = "data/selectedCells.csv"
  # TODO either multiple files with rdata/rds or single file of csv/other
  fpExtension <- tools::file_ext(gmtFile$datapath[1])
  retVal <- tryCatch({readGMT(gmtFile, scEx)},
                     error = function(e){
                       cat(file = stderr(), paste("gmtData: NULL", e,"\n"))
                       return(NULL)}
  )
  
  
  if (is.null(retVal)) {
    return(NULL)
  }
  
  
  exportTestValues(gmtData = {
    list(
      retVal
    )
  })
  return(retVal)
})


if (DEBUG) {
  cat(file = stderr(), "\n\ndone loading reactives.R.\n\n\n")
}
