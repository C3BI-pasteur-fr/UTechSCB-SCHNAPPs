
# normalization parameters

# choice for the radio buttion
myNormalizationChoices <- list(
  scEx_log = "DE_logNormalization",
  gene_norm = "DE_logGeneNormalization"
)

# value should be of class shiny.tag
# will be displayed via renderUI
myNormalizationParameters <- list(
  DE_logNormalization = h5("no Parameters implemented"),
  DE_logGeneNormalization = textInput("DE_geneIds_norm", "comma separated list of genes used for normalization", value = "")
)

# DE_logGeneNormalization ----
#' DE_logGeneNormalization
#' Normalization just as scEx_log only that sum of counts of genes is used if available.
DE_logGeneNormalization <- reactive({
  if (DEBUG) cat(file = stderr(), "DE_logGeneNormalization started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_logGeneNormalization")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_logGeneNormalization")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_logGeneNormalization", id = "DE_logGeneNormalization", duration = NULL)
  }
  
  scEx <- scEx()
  inputGenes <- input$DE_geneIds_norm
  
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "DE_logGeneNormalization:NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_logGeneNormalization.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/DE_logGeneNormalization.RData")
  
  # TODO ?? define scaling factor somewhere else???
  sfactor = max(max(assays(scEx)[["counts"]]),1000)
  retVal <- DE_logNormalizationGenefunc(scEx, inputGenes, scalingFactor = sfactor)
  
  exportTestValues(DE_logGeneNormalization = {assays(retVal)[["logcounts"]]})  
  return(retVal)
})

#' DE_logNormalizationGenefunc
#' actual computation of the normalization as it is done in seurat
DE_logNormalizationGenefunc <- function(scEx, inputGenes, scalingFactor = 10000) {
  if (DEBUG) cat(file = stderr(), "DE_logNormalizationGenefunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_logNormalizationGenefunc")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_logNormalizationGenefunc")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_logNormalizationGenefunc", id = "DE_logNormalizationGenefunc", duration = NULL)
  }
  
  # use_genes <- sort(unique(1 + slot(as(assays(scEx)[[1]], "dgTMatrix"), 
  #                                   "i")))
  # 
  # bc_sums <- Matrix::colSums(assays(scEx)[[1]])
  # median_sum <- median(bc_sums)
  genesin <- geneName2Index(inputGenes, rowData(scEx))
  if (is_empty(genesin)) {
    genesin = rownames(scEx)
  }
  A <- as(assays(scEx)[[1]], "dgTMatrix")
  assays(scEx)[[1]] = as(assays(scEx)[[1]], "dgTMatrix")
  nenner <- Matrix::colSums(assays(scEx)[[1]][genesin, ,drop=FALSE])
  nenner[nenner==0] = 1
  # A@x <- A@x / Matrix::colSums(A)[assays(scEx)[[1]]@j + 1L]
  A@x <- A@x / (nenner[A@j + 1L])
  scEx_bcnorm <- SingleCellExperiment(assay = list(logcounts = as(A,"dgTMatrix")),
                                      colData = colData(scEx),
                                      rowData = rowData(scEx))
  # assays(scEx_bcnorm[genesin,])[[1]]
  
  
  x <- uniqTsparse(assays(scEx_bcnorm)[[1]])
  slot(x, "x") <- log(1 + slot(x, "x"), base = 2) * scalingFactor
  assays(scEx_bcnorm)[[1]] <- x
  return(scEx_bcnorm)
}


# DE_logNormalization ----
#' DE_logNormalization
#' reactive for normalizing data according to seurat
DE_logNormalization <- reactive({
  if (DEBUG) cat(file = stderr(), "DE_logNormalization started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_logNormalization")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_logNormalization")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_logNormalization", id = "DE_logNormalization", duration = NULL)
  }
  
  scEx <- scEx()
  
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "DE_logNormalization:NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_logNormalization.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/DE_logNormalization.RData")
  
  # TODO ?? define scaling factor somewhere else???
  sfactor = max(max(assays(scEx)[["counts"]]),1000)
  retVal <- DE_logNormalizationfunc(scEx, scalingFactor = sfactor)
  
  exportTestValues(DE_logNormalization = {assays(retVal)[["logcounts"]]})  
  return(retVal)
})

#' DE_logNormalizationfunc
#' actual computation of the normalization as it is done in seurat
DE_logNormalizationfunc <- function(scEx, scalingFactor = 10000) {
  if (DEBUG) cat(file = stderr(), "DE_logNormalizationfunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_logNormalizationfunc")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_logNormalizationfunc")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_logNormalizationfunc", id = "DE_logNormalizationfunc", duration = NULL)
  }
  
  # use_genes <- sort(unique(1 + slot(as(assays(scEx)[[1]], "dgTMatrix"), 
  #                                   "i")))
  # 
  # bc_sums <- Matrix::colSums(assays(scEx)[[1]])
  # median_sum <- median(bc_sums)
  A <- as(assays(scEx)[[1]], "dgCMatrix")
  assays(scEx)[[1]] = as(assays(scEx)[[1]], "dgTMatrix")
  A@x <- A@x / Matrix::colSums(A)[assays(scEx)[[1]]@j + 1L]
  scEx_bcnorm <- SingleCellExperiment(assay = list(logcounts = as(A,"dgTMatrix")),
                                      colData = colData(scEx),
                                      rowData = rowData(scEx))
  
  x <- uniqTsparse(assays(scEx_bcnorm)[[1]])
  slot(x, "x") <- log(1 + slot(x, "x"), base = 2) * scalingFactor
  assays(scEx_bcnorm)[[1]] <- x
  return(scEx_bcnorm)
}
