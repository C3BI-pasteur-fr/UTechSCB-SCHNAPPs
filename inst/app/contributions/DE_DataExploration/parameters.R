suppressMessages(require(Matrix))
# normalization parameters

# choice for the radio buttion
myNormalizationChoices <- list(
  scEx_log = "DE_logNormalization",
  gene_norm = "DE_logGeneNormalization",
  SeuratStandard = "DE_seuratStandard",
  SeuratSCtransform = "DE_seuratSCtransform",
  SeuratRefBased = "DE_seuratRefBased"
)

# value should be of class shiny.tag
# will be displayed via renderUI
myNormalizationParameters <- list(
  DE_logNormalization = h5("no Parameters implemented"),
  DE_logGeneNormalization = textInput(inputId = "DE_geneIds_norm", label = "comma separated list of genes used for normalization", value = ""),
  DE_seuratStandard = tagList(
    numericInput("DE_seuratStandard_dims",
      label = "Which dimensions to use from the CCA to specify the neighbor search space",
      min = 2, max = 20000, step = 1, value = 30
    ),
    numericInput("DE_seuratStandard_anchorF",
      label = "select the provided number of features to be used in anchor finding",
      min = 60, max = 30000, step = 10,
      value = 2000
    ),
    numericInput("DE_seuratStandard_kF",
      label = "How many neighbors (k) to use when filtering anchors",
      min = 2, max = 3000, step = 1,
      value = 200
    ),
    numericInput("DE_seuratStandard_k.weight",
      label = "Number of neighbors to consider when weighting",
      min = 2, max = 3000, step = 1,
      value = 100
    )
  ),
  DE_seuratSCtransform = tagList(
    numericInput("DE_seuratSCtransform_nfeatures",
      label = "Number of features to retain/use",
      min = 200, max = 20000, step = 10, value = 3000
    ),
    numericInput("DE_seuratSCtransform_k.filter",
      label = "How many neighbors (k) to use when filtering anchors, should be smaller than the lowest number of cells per sample",
      min = 60, max = 30000, step = 10,
      value = 200
    ),
    numericInput("DE_seuratSCtransform_scaleFactor",
                 label = "Scaling to use for transformed data",
                 min = 1, max = 30000, step = 10,
                 value = 1000
    ),
    textInput("DE_seuratSCtransformm_keepfeatures", "comma separated list of genes keep", value = "")
    
  ),
  DE_seuratRefBased = tagList(
    numericInput("DE_seuratRefBased_nfeatures",
      label = "Number of features to retain/use",
      min = 200, max = 20000, step = 10, value = 3000
    ),
    numericInput("DE_seuratRefBased_k.filter",
      label = "How many neighbors (k) to use when filtering anchors, should be smaller than the lowest number of cells per sample",
      min = 60, max = 30000, step = 10,
      value = 200
    ),
    numericInput("DE_seuratRefBased_scaleFactor",
                 label = "Scaling to use for transformed data",
                 min = 1, max = 30000, step = 10,
                 value = 1000
    ),
    textInput("DE_seuratRefBased_keepfeatures", "comma separated list of genes keep", value = "")
  )
)

# DE_seuratRefBasedFunc ----
DE_seuratRefBasedFunc <- function(scEx, nfeatures = 3000, k.filter = 100, 
                                  scalingFactor = 1000, keep.features = "") {
  require(Seurat)
  cellMeta <- colData(scEx)
  # split in different samples
  meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
  seurDat <- CreateSeuratObject(
    counts = assays(scEx)[[1]],
    meta.data = meta.data
  )
  seur.list <- SplitObject(seurDat, split.by = "sampleNames")
  for (i in 1:length(seur.list)) {
    seur.list[[i]] <- SCTransform(seur.list[[i]], verbose = TRUE)
  }
  integrated <- tryCatch(
    {
      # save(file = "~/SCHNAPPsDebug/DE_seuratRefBased.RData", list = c(ls()))
      # load(file = "~/SCHNAPPsDebug/DE_seuratRefBased.RData")
      features <- SelectIntegrationFeatures(object.list = seur.list, nfeatures = nfeatures)
     
      keep.features = keep.features[keep.features %in% rownames(scEx)]
      features = unique(c(features, keep.features))

      seur.list <- PrepSCTIntegration(
        object.list = seur.list, anchor.features = features,
        verbose = TRUE
      )
      # take the sample with the highest number of cells as reference
      reference_dataset <- order(unlist(lapply(seur.list, FUN = function(x) {ncol(x)})), decreasing =T)[1]
      
 
      anchors <- FindIntegrationAnchors(
        object.list = seur.list, normalization.method = "SCT",
        anchor.features = features, verbose = TRUE, k.filter = k.filter,
        reference = reference_dataset
      )
      integrated <- IntegrateData(
        anchorset = anchors, normalization.method = "SCT",
        verbose = TRUE
      )
      # integrated <- NormalizeData(integrated, verbose = TRUE)

      # integrated <- RunPCA(integrated, verbose = TRUE)
      # integrated <- RunUMAP(integrated, dims = 1:30)
      # plots <- DimPlot(integrated, group.by = c("sampleNames"), combine = TRUE, dims = )
      # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") +
      #                   guides(color = guide_legend(nrow = 3,
      #                                               byrow = TRUE, override.aes = list(size = 3))))
      # CombinePlots(plots)

      # FeaturePlot(integrated, c("CCR7", "S100A4", "GZMB", "GZMK", "GZMH"))
    },
    error = function(e) {
      cat(file = stderr(), paste("\n\n!!!Error during Seurat normalization:\n", e, "\n\n"))
      return(NULL)
    }
  )
  if (is.null(integrated)) {
    return(NULL)
  }
  A <- integrated@assays$integrated@data * scalingFactor
  scEx_bcnorm <- SingleCellExperiment(
    assay = list(logcounts = as(A, "dgTMatrix")),
    colData = colData(scEx)[colnames(A), , drop = FALSE],
    rowData = rowData(scEx)[rownames(A), , drop = FALSE]
  )
  return(scEx_bcnorm)
}

DE_seuratRefBasedButton <- reactiveVal(
  label = "DE_seuratRefBasedButton",
  value = ""
)

# DE_seuratRefBased ----
DE_seuratRefBased <- reactive({
  if (DEBUG) cat(file = stderr(), "DE_seuratRefBased started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_seuratRefBased")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_seuratRefBased")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_seuratRefBased", id = "DE_seuratRefBased", duration = NULL)
  }

  scEx <- scEx()
  runThis <- DE_seuratRefBasedButton()
  nfeatures <- isolate(input$DE_seuratRefBased_nfeatures)
  k.filter <- isolate(input$DE_seuratRefBased_k.filter)
  scalingFactor <- isolate(input$DE_seuratRefBased_scaleFactor)
  geneNames <-  isolate(input$DE_seuratRefBased_keepfeatures)

  if (is.null(scEx) | runThis == "") {
    if (DEBUG) {
      cat(file = stderr(), "DE_seuratRefBased:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_seuratRefBased.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/DE_seuratRefBased.RData")
  .schnappsEnv$normalizationFactor = scalingFactor
  featureData <- rowData(scEx)
  geneid <- geneName2Index(geneNames, featureData)
  
  
  # # TODO ?? define scaling factor somewhere else???
  # sfactor = max(max(assays(scEx)[["counts"]]),1000)
  retVal <- DE_seuratRefBasedFunc(scEx = scEx, nfeatures = nfeatures, k.filter = k.filter, 
                                  scalingFactor = scalingFactor, keep.features = geneid)
  
  if (is.null(retVal)) {
    showNotification("An error occurred during Seurat normalization, please check console", id = "DE_seuratError", duration = NULL, type = "error")
  }

  .schnappsEnv$calculated_DE_seuratRefBased_nfeatures <- nfeatures
  .schnappsEnv$calculated_DE_seuratRefBased_k.filter <- k.filter
  .schnappsEnv$calculated_DE_seuratRefBased_scaleFactor <- scalingFactor

  addClass("updateNormalization", "green")

  exportTestValues(DE_seuratSCtransform = {
    assays(retVal)[["logcounts"]]
  })
  return(retVal)
})

# DE_seuratSCtransformFunc =======
DE_seuratSCtransformFunc <- function(scEx, nfeatures = 3000, k.filter = 100, 
                                     scalingFactor = 1000, keep.features = "") {
  require(Seurat)
  cellMeta <- colData(scEx)
  # split in different samples
  meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
  integrated <- tryCatch(
    {
      seurDat <- CreateSeuratObject(
        counts = assays(scEx)[[1]],
        meta.data = meta.data
      )
      seur.list <- SplitObject(seurDat, split.by = "sampleNames")
      for (i in 1:length(seur.list)) {
        seur.list[[i]] <- SCTransform(seur.list[[i]], verbose = TRUE)
      }

      features <- SelectIntegrationFeatures(object.list = seur.list, nfeatures = nfeatures)
      keep.features = keep.features[keep.features %in% rownames(scEx)]
      features = unique(c(features, keep.features))
      
      seur.list <- PrepSCTIntegration(
        object.list = seur.list, anchor.features = features,
        verbose = TRUE
      )
      
      anchors <- FindIntegrationAnchors(
        object.list = seur.list, normalization.method = "SCT",
        anchor.features = features, verbose = TRUE, k.filter = k.filter
      )
      keep.features = keep.features[keep.features %in% rownames(scEx)]
      anchors = unique(c(anchors, keep.features))
      integrated <- IntegrateData(
        anchorset = anchors, normalization.method = "SCT",
        verbose = TRUE
      )
    },
    error = function(e) {
      cat(file = stderr(), paste("\n\n!!!Error during Seurat normalization:\n", e, "\n\n"))
      return(NULL)
    }
  )
  # integrated <- NormalizeData(integrated, verbose = TRUE)

  # integrated <- RunPCA(integrated, verbose = TRUE)
  # integrated <- RunUMAP(integrated, dims = 1:30)
  # plots <- DimPlot(integrated, group.by = c("sampleNames"), combine = TRUE, dims = )
  # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") +
  #                   guides(color = guide_legend(nrow = 3,
  #                                               byrow = TRUE, override.aes = list(size = 3))))
  # CombinePlots(plots)

  # FeaturePlot(integrated, c("CCR7", "S100A4", "GZMB", "GZMK", "GZMH"))

  if (is.null(integrated)) {
    return(NULL)
  }
  A <- integrated@assays$integrated@data * scalingFactor
  scEx_bcnorm <- SingleCellExperiment(
    assay = list(logcounts = as(A, "dgTMatrix")),
    colData = colData(scEx)[colnames(A), , drop = FALSE],
    rowData = rowData(scEx)[rownames(A), , drop = FALSE]
  )
  return(scEx_bcnorm)
}

DE_seuratSCtransformButton <- reactiveVal(
  label = "DE_seuratSCtransformButton",
  value = ""
)

# DE_seuratSCtransform ----
DE_seuratSCtransform <- reactive({
  if (DEBUG) cat(file = stderr(), "DE_seuratSCtransform started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_seuratSCtransform")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_seuratSCtransform")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_seuratSCtransform", id = "DE_seuratSCtransform", duration = NULL)
  }

  scEx <- scEx()
  runThis <- DE_seuratSCtransformButton()
  nfeatures <- isolate(input$DE_seuratSCtransform_nfeatures)
  k.filter <- isolate(input$DE_seuratSCtransform_k.filter)
  scalingFactor <- isolate(input$DE_seuratSCtransform_scaleFactor)
  geneNames <-  isolate(input$DE_seuratSCtransformm_keepfeatures)

  if (is.null(scEx) | runThis == "") {
    if (DEBUG) {
      cat(file = stderr(), "DE_seuratSCtransform:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_seuratSCtransform.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/DE_seuratSCtransform.RData")
  .schnappsEnv$normalizationFactor <- scalingFactor
  featureData <- rowData(scEx)
  geneid <- geneName2Index(geneNames, featureData)
  
  # # TODO ?? define scaling factor somewhere else???
  # sfactor = max(max(assays(scEx)[["counts"]]),1000)
  retVal <- DE_seuratSCtransformFunc(scEx = scEx, nfeatures = nfeatures, k.filter = k.filter, 
                                     scalingFactor = scalingFactor, keep.features = geneid)
  
  if (is.null(retVal)) {
    if (DEBUG) green(cat(file = stderr(), "An error occurred during Seurat normalization, please check console\n"))
    showNotification("An error occurred during Seurat normalization, please check console", id = "DE_seuratError", duration = NULL, type = "error")
  }

  .schnappsEnv$calculated_DE_seuratSCtransform_nfeatures <- nfeatures
  .schnappsEnv$calculated_DE_seuratSCtransform_k.filter <- k.filter
  .schnappsEnv$calculated_DE_seuratSCtransform_scaleFactor <- scalingFactor

  addClass("updateNormalization", "green")
  exportTestValues(DE_seuratSCtransform = {
    assays(retVal)[["logcounts"]]
  })
  return(retVal)
})



DE_seuratStandardfunc <- function(scEx, dims = 10, anchorsF = 2000, kF = 200, k.weight = 100) {
  require(Seurat)
  cellMeta <- colData(scEx)
  # split in different samples
  meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
  # creates object @assays$RNA@data and @assays$RNA@counts
  integrated <- NULL
  seurDat <- tryCatch(
    {
      seurDat <- CreateSeuratObject(
        counts = assays(scEx)[[1]],
        meta.data = meta.data
      )
      seur.list <- SplitObject(seurDat, split.by = "sampleNames")
      for (i in 1:length(seur.list)) {
        seur.list[[i]] <- NormalizeData(seur.list[[i]], verbose = FALSE)
        seur.list[[i]] <- FindVariableFeatures(seur.list[[i]],
          selection.method = "vst",
          nfeatures = 2000, verbose = FALSE
        )
      }
      anchors <- FindIntegrationAnchors(
        object.list = seur.list, dims = 1:dims, anchor.features = anchorsF,
        k.filter = kF
      )
      integrated <- IntegrateData(anchorset = anchors, dims = 1:dims, k.weight = k.weight, verbose = TRUE)
      DefaultAssay(integrated) <- "integrated"

      # Run the standard workflow for visualization and clustering
      integrated <- ScaleData(integrated, verbose = TRUE)
      # integrated@assays


      # NormalizeData(seurDat, normalization.method = "LogNormalize", scale.factor = 10000)
    },
    error = function(e) {
      cat(file = stderr(), paste("\n\n!!!Error during Seurat normalization:\n", e, "\n\n"))
      return(NULL)
    }
  )
  if (is.null(seurDat)) {
    return(NULL)
  }
  A <- seurDat@assays$integrated@data
  scEx_bcnorm <- SingleCellExperiment(
    assay = list(logcounts = as(A, "dgTMatrix")),
    colData = colData(scEx)[colnames(A), , drop = FALSE],
    rowData = rowData(scEx)[rownames(A), , drop = FALSE]
  )
  return(scEx_bcnorm)
}

# DE_seuratStandard ----
DE_seuratStandardButton <- reactiveVal(
  label = "DE_seuratStandardButton",
  value = ""
)
DE_seuratStandard <- reactive({
  if (DEBUG) cat(file = stderr(), "DE_seuratStandard started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_seuratStandard")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_seuratStandard")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_seuratStandard", id = "DE_seuratStandard", duration = NULL)
  }

  scEx <- scEx()
  runThis <- DE_seuratStandardButton()

  sDims <- isolate(input$DE_seuratStandard_dims)
  anchorF <- isolate(input$DE_seuratStandard_anchorF)
  kF <- isolate(input$DE_seuratStandard_kF)
  k.weight <- isolate(input$DE_seuratStandard_k.weight)

  if (is.null(scEx) | runThis == "") {
    if (DEBUG) {
      cat(file = stderr(), "DE_seuratStandard:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_seuratStandard.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/DE_seuratStandard.RData")



  # # TODO ?? define scaling factor somewhere else???
  # sfactor = max(max(assays(scEx)[["counts"]]),1000)
  retVal <- DE_seuratStandardfunc(scEx = scEx, dims = sDims, anchorsF = anchorF, kF = kF, k.weight = k.weight)

  if (is.null(retVal)) {
    showNotification("An error occurred during Seurat normalization, please check console", id = "DE_seuratError", duration = NULL, type = "error")
  }

  .schnappsEnv$calculated_DE_seuratStandard_dims <- sDims
  .schnappsEnv$calculated_DE_seuratStandard_anchorF <- anchorF
  .schnappsEnv$calculated_DE_seuratStandard_kF <- kF
  .schnappsEnv$calculated_DE_seuratStandard_k.weight <- k.weight

  addClass("updateNormalization", "green")
  exportTestValues(DE_seuratStandard = {
    assays(retVal)[["logcounts"]]
  })
  return(retVal)
})


# observe input$DE_geneIds_norm ----


# DE_logGeneNormalization ----
#' DE_logGeneNormalization
#' Normalization just as scEx_log only that sum of counts of genes is used if available.
#'

DE_logGeneNormalizationButton <- reactiveVal(
  value = NULL,
  label = "pressed"
)


DE_logGeneNormalization <- reactive(label = "rlogGene", {
  if (DEBUG) cat(file = stderr(), "DE_logGeneNormalization started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_logGeneNormalization")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_logGeneNormalization")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_logGeneNormalization", id = "DE_logGeneNormalization", duration = NULL)
  }

  scEx <- scEx()
  DE_logGeneNormalizationButton()
  # input$updateNormalizationButton
  # inputGenes <- isolate(input$DE_geneIds_norm)
  inputGenes <- isolate(input$DE_geneIds_norm)

  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "DE_logGeneNormalization:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_logGeneNormalization.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/DE_logGeneNormalization.RData")

  # TODO ?? define scaling factor somewhere else???
  sfactor <- max(max(assays(scEx)[["counts"]]), 1000)
  retVal <- DE_logNormalizationGenefunc(scEx, inputGenes, scalingFactor = sfactor)

  if (is.null(retVal)) {
    # print error message
  }
  if (DEBUG) {
    cat(file = stderr(), "assign: .schnappsEnv$calculated_DE_geneIds_norm\n")
  }
  .schnappsEnv$calculated_DE_geneIds_norm <- inputGenes

  # turn normalization button green
  addClass("updateNormalization", "green")

  .schnappsEnv$normalizationFactor <- sfactor
  exportTestValues(DE_logGeneNormalization = {
    assays(retVal)[["logcounts"]]
  })
  return(retVal)
})

#' DE_logNormalizationGenefunc
#' actual computation of the normalization as it is done in seurat
DE_logNormalizationGenefunc <- function(scEx, inputGenes, scalingFactor = 10000) {
  if (DEBUG) cat(file = stderr(), "DE_logNormalizationGenefunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_logNormalizationGenefunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_logNormalizationGenefunc")
    }
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
    genesin <- rownames(scEx)
  }
  A <- as(assays(scEx)[[1]], "dgTMatrix")
  assays(scEx)[[1]] <- as(assays(scEx)[[1]], "dgTMatrix")
  nenner <- Matrix::colSums(assays(scEx)[[1]][genesin, , drop = FALSE])
  nenner[nenner == 0] <- 1
  # A@x <- A@x / Matrix::colSums(A)[assays(scEx)[[1]]@j + 1L]
  A@x <- A@x / (nenner[A@j + 1L])
  scEx_bcnorm <- SingleCellExperiment(
    assay = list(logcounts = as(A, "dgTMatrix")),
    colData = colData(scEx),
    rowData = rowData(scEx)
  )
  # assays(scEx_bcnorm[genesin,])[[1]]


  x <- uniqTsparse(assays(scEx_bcnorm)[[1]])
  slot(x, "x") <- log(1 + slot(x, "x"), base = 2) * scalingFactor
  assays(scEx_bcnorm)[[1]] <- x
  return(scEx_bcnorm)
}


# DE_logNormalization ----
#' DE_logNormalization
#' reactive for normalizing data according to seurat
DE_logNormalizationButton <- reactiveVal(
  value = NULL,
  label = "pressed"
)

DE_logNormalization <- reactive(label = "rlogNorm", {
  if (DEBUG) cat(file = stderr(), "DE_logNormalization started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_logNormalization")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_logNormalization")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_logNormalization", id = "DE_logNormalization", duration = NULL)
  }

  scEx <- scEx()
  DE_logNormalizationButton()

  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "DE_logNormalization:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_logNormalization.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/DE_logNormalization.RData")

  # TODO ?? define scaling factor somewhere else???
  sfactor <- max(max(assays(scEx)[["counts"]]), 1000)
  retVal <- DE_logNormalizationfunc(scEx, scalingFactor = sfactor)

  # turn normalization button green
  addClass("updateNormalization", "green")

  .schnappsEnv$normalizationFactor <- sfactor
  exportTestValues(DE_logNormalization = {
    assays(retVal)[["logcounts"]]
  })
  return(retVal)
})

#' DE_logNormalizationfunc
#' actual computation of the normalization as it is done in seurat
DE_logNormalizationfunc <- function(scEx, scalingFactor = 10000) {
  if (DEBUG) cat(file = stderr(), "DE_logNormalizationfunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_logNormalizationfunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_logNormalizationfunc")
    }
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
  assays(scEx)[[1]] <- as(assays(scEx)[[1]], "dgTMatrix")
  A@x <- A@x / Matrix::colSums(A)[assays(scEx)[[1]]@j + 1L]
  scEx_bcnorm <- SingleCellExperiment(
    assay = list(logcounts = as(A, "dgTMatrix")),
    colData = colData(scEx),
    rowData = rowData(scEx)
  )

  x <- uniqTsparse(assays(scEx_bcnorm)[[1]])
  slot(x, "x") <- log(1 + slot(x, "x"), base = 2) * scalingFactor
  assays(scEx_bcnorm)[[1]] <- x
  return(scEx_bcnorm)
}
