suppressMessages(require(Matrix))
suppressMessages(require(Seurat))
# normalization parameters

# choice for the radio buttion
myNormalizationChoices <- list(
  scEx_log = "DE_logNormalization",
  scaterNorm = "DE_scaterNormalization",
  gene_norm = "DE_logGeneNormalization",
  "seurat::SCT normalization" = "DE_seuratSCTnorm",
  "seurat::(NormalizeData + FindVariableFeatures + ScaleData)" = "DE_seuratLogNorm",
  "seurat::SCTransform" = "DE_seuratStandard",
  "seurat::integrate and SCtransform" = "DE_seuratSCtransform",
  "seurat::SeuratRefBased" = "DE_seuratRefBased"
)

# value should be of class shiny.tag
# will be displayed via renderUI
myNormalizationParameters <- list(
  DE_logNormalization = tagList(
    numericInput("DE_logNormalization_sf",
                 label = "scale by (0 => minvalue)",
                 min = 0, max = 200000, step = 1, value = 0
    )
    ),
  DE_scaterNormalization = h5("no Parameters implemented"),
  DE_seuratSCTnorm = tagList(
    numericInput("DE_seuratSCTnorm_nHVG",
                 label = "Number of variable genes",
                 min = 10, max = 200000, step = 1, value = 3000
    ),
    selectInput(
      "DE_seuratSCTnorm_var2reg",
      label = "Vars to regress out (only factors are allowed)",
      multiple = TRUE,
      choices = unique(c(defaultValue("DE_seuratSCTnorm_var2reg", ""), ""))
      ,
      selected = defaultValue("DE_seuratSCTnorm_var2reg", "")
      # , 
      # options = list(maxItems = 20)
    )
  ),
  DE_seuratLogNorm = tagList(
    numericInput("DE_seuratLogNorm_nHVG",
                 label = "Number of variable genes",
                 min = 10, max = 200000, step = 1, value = 3000
    ),
    selectInput(
      "DE_seuratLogNorm_var2reg",
      label = "Vars to regress out (only factors are allowed)",
      multiple = TRUE,
      choices = unique(c(defaultValue("DE_seuratLogNorm_var2reg", ""), ""))
      ,
      selected = defaultValue("DE_seuratLogNorm_var2reg", "")
      # , 
      # options = list(maxItems = 20)
    )
  ),
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
    ),
    selectInput(
      "DE_seuratStandard_splitby",
      label = "apply scTransform to individual entries (levels of factor)",
      multiple = FALSE,
      choices = unique(c(defaultValue("DE_seuratStandard_splitby", ""), ""))
      , selected = defaultValue("DE_seuratStandard_splitby", "")
      # , 
      # options = list(maxItems = 20)
    )
  ),
  DE_seuratSCtransform = tagList(
    numericInput("DE_seuratSCtransform_nhvg",
                 label = "number of highly variable genes",
                 min = 200, max = 400000, step = 10, value = 3000
    ),
    selectInput(
      "DE_seuratSCtransform_vars2regress",
      label = "Vars to regress out (only factors are allowed)",
      multiple = TRUE,
      choices = unique(c(defaultValue("DE_seuratSCtransform_vars2regress", ""), ""))
      ,
      selected = defaultValue("DE_seuratSCtransform_vars2regress", "")
      # , 
      # options = list(maxItems = 20)
    ),
    numericInput("DE_seuratSCtransform_dimsMin",
                 label = "minimum dimension (PCA) to use",
                 min = 1, max = 100, step = 1, value = 1
    ),
    numericInput("DE_seuratSCtransform_dimsMax",
                 label = "maximum dimension (PCA) to use",
                 min = 3, max = 100, step = 1, value = 30
    ),
    numericInput("DE_seuratSCtransform_nfeatures",
                 label = "number of genes to keep for the integration step",
                 min = 200, max = 400000, step = 10, value = 3000
    ),
    numericInput("DE_seuratSCtransform_k.anchor",
                 label = "Number of anchors to use",
                 min = 1, max = 100, step = 1,
                 value = 5
    ),
    numericInput("DE_seuratSCtransform_k.filter",
                 label = "How many neighbors (k) to use when filtering anchors, should be smaller than the lowest number of cells per sample",
                 min = 1, max = 1000, step = 10,
                 value = 200
    ),
    numericInput("DE_seuratSCtransform_k.score",
                 label = "k score",
                 min = 1, max = 1000, step = 1,
                 value = 30
    ),
    # numericInput("DE_seuratSCtransform_scalingFactor",
    #              label = "Scaling to use for transformed data",
    #              min = 1, max = 30000, step = 10,
    #              value = 1000
    # ),
    selectInput(
      "DE_seuratSCtransformm_split.by",
      label = "apply scTransform to individual entries (levels of factor)",
      multiple = FALSE,
      choices = unique(c(defaultValue("DE_seuratSCtransformm_split.by", ""), ""))
      ,
      selected = defaultValue("DE_seuratSCtransformm_split.by", "")
      # , 
      # options = list(maxItems = 20)
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
    # numericInput("DE_seuratRefBased_scaleFactor",
    #              label = "Scaling to use for transformed data",
    #              min = 1, max = 30000, step = 10,
    #              value = 1000
    # ),
    selectInput(
      "DE_seuratRefBased_splitby",
      label = "apply scTransform to individual entries (levels of factor)",
      multiple = FALSE,
      choices = unique(c(defaultValue("DE_seuratRefBased_splitby", ""), ""))
      , selected = defaultValue("DE_seuratRefBased_splitby", "")
      # , 
      # options = list(maxItems = 20)
    ),
    textInput("DE_seuratRefBased_keepfeatures", "comma separated list of genes keep", value = "")
  )
)

# DE_seuratRefBasedFunc ----
DE_seuratRefBasedFunc <- function(scEx, nfeatures = 3000, k.filter = 100, 
                                  keep.features = "",
                                  splitby = NULL) {
  require(Seurat)
  cellMeta <- colData(scEx)
  # split in different samples
  A <- tryCatch(
    {
      # save(file = "~/SCHNAPPsDebug/DE_seuratRefBased.RData", list = c(ls()))
      # load(file = "~/SCHNAPPsDebug/DE_seuratRefBased.RData")
      if(is.null(splitby)) {
        splitby = "ident"
        meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
      } else if(splitby == "" ) {
        splitby = "ident"
        meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
      } else {
        meta.data <- as.data.frame(cellMeta[, splitby, drop = FALSE])
        limitCells = meta.data[,1] %in% levels(meta.data[,1])[table(meta.data[,1]) > 30]
        # we cannot remove cell here because this would change scEX and projections
        # not sure that NA would be a good solution
        # so we are asking to remove the cells manually
        if (sum(limitCells) < ncol(scEx)) {
          errStr = paste("please remove the following cells:\n", 
                         paste(colnames(scEx)[!limitCells],
                               collapse = ", "))
          cat (file = stderr(), errStr)
          if (!is.null(getDefaultReactiveDomain())) {
            showNotification(errStr, id = "DE_seuratError", duration = NULL, type = "error")  
          }
          return(NULL)
        }
        # 
        # scEx = scEx[, limitCells]
        # meta.data = meta.data[limitCells,, drop = FALSE]
      }
      seurDat <- CreateSeuratObject(
        counts = assays(scEx)[[1]],
        meta.data = meta.data
      )
      seur.list <- SplitObject(seurDat, split.by = splitby)
      # for (i in 1:length(seur.list)) {
      #   seur.list[[i]] <- SCTransform(seur.list[[i]], verbose = TRUE)
      # }
      seur.list <- lapply(seur.list, FUN = function(x) SCTransform(object = x, 
                                                                   # variable.features.n = nhvg, 
                                                                   # vars.to.regress=vars2regress, 
                                                                   verbose = DEBUG)
      )
      
      features <- SelectIntegrationFeatures(object.list = seur.list, nfeatures = nfeatures)
      
      if (length(keep.features)>1) {
        keep.features = keep.features[keep.features %in% rownames(scEx)]
      }
      
      features = unique(c(features, keep.features))
      
      seur.list <- PrepSCTIntegration(
        object.list = seur.list, anchor.features = features,
        verbose = DEBUG
      )
      
      # take the sample with the highest number of cells as reference
      reference_dataset <- order(unlist(lapply(seur.list, FUN = function(x) {ncol(x)})), decreasing =T)[1]
      
      if (length(seur.list) > 1) {
        anchors <- FindIntegrationAnchors(
          object.list = seur.list, 
          normalization.method = "SCT",
          anchor.features = features, 
          verbose = DEBUG, 
          # k.anchor = k.anchor,
          k.filter = k.filter, 
          reference = reference_dataset
          # k.score=k.score,
          # dims=dimsMin:dimsMax
        )
        # keep.features = keep.features[keep.features %in% rownames(scEx)]
        # anchors = unique(c(anchors, keep.features))
        integrated <- IntegrateData(
          anchorset = anchors, normalization.method = "SCT",
          verbose = DEBUG,
          k.weight = min(100, min(unlist(lapply(seur.list, ncol))))
        )
        integrated@assays$integrated@scale.data
      } else {
        seur.list[[1]]@assays$SCT@scale.data 
      }
      
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
      cat(file = stderr(), paste("\n\n!!!Error during Seurat normalization:\nNeed multiple samples", e, "\n\n"))
      return(NULL)
    }
  )
  if (is.null(integrated)) {
    return(NULL)
  }
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
  # scalingFactor <- isolate(input$DE_seuratRefBased_scaleFactor)
  geneNames <-  isolate(input$DE_seuratRefBased_keepfeatures)
  splitby <- isolate(input$DE_seuratRefBased_splitby)
  
  if (is.null(scEx) | runThis == "") {
    if (DEBUG) {
      cat(file = stderr(), "DE_seuratRefBased:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_seuratRefBased.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/DE_seuratRefBased.RData")
  # .schnappsEnv$normalizationFactor = scalingFactor
  featureData <- rowData(scEx)
  geneid <- geneName2Index(geneNames, featureData)
  
  
  # # TODO ?? define scaling factor somewhere else???
  # sfactor = max(max(assays(scEx)[["counts"]]),1000)
  retVal <- DE_seuratRefBasedFunc(scEx = scEx, nfeatures = nfeatures, k.filter = k.filter, 
                                  # scalingFactor = scalingFactor,
                                  keep.features = geneid,
                                  splitby = splitby)
  
  if (is.null(retVal)) {
    showNotification("An error occurred during Seurat normalization, please check console", id = "DE_seuratError", duration = NULL, type = "error")
  }
  
  .schnappsEnv$calculated_DE_seuratRefBased_nfeatures <- nfeatures
  .schnappsEnv$calculated_DE_seuratRefBased_k.filter <- k.filter
  # .schnappsEnv$calculated_DE_seuratRefBased_scaleFactor <- scalingFactor
  
  addClass("updateNormalization", "green")
  
  exportTestValues(DE_seuratSCtransform = {
    assays(retVal)[["logcounts"]]
  })
  return(retVal)
})

# DE_seuratSCtransformFunc =======
DE_seuratSCtransformFunc <- function(scEx, 
                                     nhvg = 3000, 
                                     nfeatures = 3000,
                                     vars2regress = "", 
                                     dimsMin = 1, 
                                     dimsMax = 30, 
                                     k.anchor = 5,
                                     k.filter = 200, 
                                     k.score = 30,
                                     splitby = "sampleNames",
                                     # scalingFactor = 1000, 
                                     keep.features = "") {
  require(Seurat)
  cellMeta <- colData(scEx)
  # split in different samples
  # meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
  
  A <- tryCatch(
    {
      if(is.null(splitby)) {
        splitby = "ident"
        meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
      } else if(splitby == "" ) {
        splitby = "ident"
        meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
      } else {
        meta.data <- as.data.frame(cellMeta[, splitby, drop = FALSE])
        # only groups of cells with more than 20 cells
        limitCells = meta.data[,1] %in% levels(meta.data[,1])[table(meta.data[,1]) > 30]
        # we cannot remove cell here because this would change scEX and projections
        # not sure that NA would be a good solution
        # so we are asking to remove the cells manually
        if (sum(limitCells) < ncol(scEx)) {
          errStr = paste("please remove the following cells:\n", 
                         paste(colnames(scEx)[!limitCells],
                               collapse = ", "))
          cat (file = stderr(), errStr)
          if (!is.null(getDefaultReactiveDomain())) {
            showNotification(errStr, id = "DE_seuratError", duration = NULL, type = "error")  
          }
          return(NULL)
        }
        # scEx = scEx[, limitCells]
        # meta.data = meta.data[limitCells,, drop = FALSE]
      }
      
      seurDat <- CreateSeuratObject(
        counts = assays(scEx)[[1]],
        meta.data = meta.data
      )
      
      seur.list <- SplitObject(seurDat, split.by = splitby)
      seur.list <- lapply(seur.list, FUN = function(x) SCTransform(object = x, 
                                                                   variable.features.n = nhvg, 
                                                                   vars.to.regress=vars2regress, verbose = DEBUG)
      )
      #       for (i in 1:length(seur.list)) {
      #   seur.list[[i]] <- SCTransform(object = seur.list[[i]], variable.features.n = nhvg, vars.to.regress=vars2regress, verbose = TRUE)
      # }
      sel.features <- SelectIntegrationFeatures(object.list = seur.list, nfeatures = nfeatures)
      
      
      if (length(keep.features)>1) {
        keep.features = keep.features[keep.features %in% rownames(scEx)]
      }
      sel.features = unique(c(sel.features, keep.features))
      sel.features = sel.features[ sel.features %in% rownames(seur.list[[1]])]
      seur.list <- PrepSCTIntegration(
        object.list = seur.list, anchor.features = sel.features,
        verbose = DEBUG
      )
      
      if (length(seur.list) > 1) {
        anchors <- FindIntegrationAnchors(
          object.list = seur.list, 
          normalization.method = "SCT",
          anchor.features = sel.features, 
          verbose = DEBUG, 
          k.anchor = k.anchor,
          k.filter = k.filter, 
          k.score=k.score,
          dims=dimsMin:dimsMax
        )
        # keep.features = keep.features[keep.features %in% rownames(scEx)]
        # anchors = unique(c(anchors, keep.features))
        integrated <- IntegrateData(
          anchorset = anchors, normalization.method = "SCT",
          verbose = DEBUG,
          k.weight = min(100, min(unlist(lapply(seur.list, ncol))))
        )
        integrated@assays$integrated@scale.data
      } else {
        seur.list[[1]]@assays$SCT@scale.data 
      }
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
  
  if (is.null(A)) {
    return(NULL)
  }
  
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
# Seurat SCT and integration
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
  nhvg <- isolate(input$DE_seuratSCtransform_nhvg)
  nfeatures <- isolate(input$DE_seuratSCtransform_nfeatures)
  vars2regress <- isolate(input$DE_seuratSCtransform_vars2regress)
  dimsMin <- isolate(input$DE_seuratSCtransform_dimsMin)
  dimsMax <- isolate(input$DE_seuratSCtransform_dimsMax)
  k.anchor <- isolate(input$DE_seuratSCtransform_k.anchor )
  k.filter <- isolate(input$DE_seuratSCtransform_k.filter)
  k.score <- isolate(input$DE_seuratSCtransform_k.score)
  # scalingFactor <- isolate(input$DE_seuratSCtransform_scalingFactor)
  geneNames <-  isolate(input$DE_seuratSCtransformm_keepfeatures)
  split.by <-  isolate(input$DE_seuratSCtransformm_split.by)
  
  if (is.null(scEx) | runThis == "") {
    if (DEBUG) {
      cat(file = stderr(), "DE_seuratSCtransform:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_seuratSCtransform.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/DE_seuratSCtransform.RData")
  # .schnappsEnv$normalizationFactor <- scalingFactor
  featureData <- rowData(scEx)
  geneid <- geneName2Index(geneNames, featureData)
  
  # # TODO ?? define scaling factor somewhere else???
  # sfactor = max(max(assays(scEx)[["counts"]]),1000)
  retVal <- DE_seuratSCtransformFunc(scEx = scEx, 
                                     nhvg = nhvg, 
                                     vars2regress = vars2regress, 
                                     nfeatures = nfeatures, 
                                     dimsMin = dimsMin, 
                                     dimsMax = dimsMax, 
                                     k.anchor = k.anchor,
                                     k.filter = k.filter, 
                                     k.score = k.score,
                                     splitby = split.by,
                                     # scalingFactor = scalingFactor, 
                                     keep.features = geneid)
  
  if (is.null(retVal)) {
    if (DEBUG) green(cat(file = stderr(), "An error occurred during Seurat normalization, please check console\n"))
    showNotification("An error occurred during Seurat normalization, please check console", id = "DE_seuratError", duration = NULL, type = "error")
  }
  
  .schnappsEnv$calculated_DE_seuratSCtransform_nfeatures <- nfeatures
  .schnappsEnv$calculated_DE_seuratSCtransform_k.filter <- k.filter
  # .schnappsEnv$calculated_DE_seuratSCtransform_scaleFactor <- scalingFactor
  
  addClass("updateNormalization", "green")
  exportTestValues(DE_seuratSCtransform = {
    assays(retVal)[["logcounts"]]
  })
  return(retVal)
})


# DE_seuratStandardfunc ----
DE_seuratStandardfunc <- function(scEx, dims = 10, anchorsF = 2000, kF = 200, k.weight = 100, 
                                  splitby = "sampleNames") {
  require(Seurat)
  cellMeta <- colData(scEx)
  # split in different samples
  # creates object @assays$RNA@data and @assays$RNA@counts
  integrated <- NULL
  A <- tryCatch(
    {
      if (is.null(splitby)) {
        splitby = "ident"
        meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
      } else if(splitby == "" ) {
        splitby = "ident"
        meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
      } else {
        meta.data <- as.data.frame(cellMeta[, splitby, drop = FALSE])
        limitCells = meta.data[,1] %in% levels(meta.data[,1])[table(meta.data[,1]) > 30]
        # we cannot remove cell here because this would change scEX and projections
        # not sure that NA would be a good solution
        # so we are asking to remove the cells manually
        if (sum(limitCells) < ncol(scEx)) {
          errStr = paste("please remove the following cells:\n", 
                         paste(colnames(scEx)[!limitCells],
                               collapse = ", "))
          cat (file = stderr(), errStr)
          if (!is.null(getDefaultReactiveDomain())) {
            showNotification(errStr, id = "DE_seuratError", duration = NULL, type = "error")  
          }
          return(NULL)
        }
        # 
        # scEx = scEx[, limitCells]
        # meta.data = meta.data[limitCells,, drop = FALSE]
      }
      seurDat <- CreateSeuratObject(
        counts = assays(scEx)[[1]],
        meta.data = meta.data
      )
      seur.list <- SplitObject(seurDat, split.by = splitby)
      seur.list <- lapply(seur.list, FUN = function(x) {
        x <- NormalizeData(x, verbose = DEBUG)
        x <- FindVariableFeatures(x,
                                  selection.method = "vst",
                                  nfeatures = 2000, verbose = DEBUG
        )
      }
      )
      # for (i in 1:length(seur.list)) {
      #   seur.list[[i]] <- NormalizeData(seur.list[[i]], verbose = FALSE)
      #   seur.list[[i]] <- FindVariableFeatures(seur.list[[i]],
      #     selection.method = "vst",
      #     nfeatures = 2000, verbose = FALSE
      #   )
      # }
      if (length(seur.list) > 1){
        anchors <- FindIntegrationAnchors(
          object.list = seur.list, dims = 1:dims, anchor.features = anchorsF,
          k.filter = kF
        )
        integrated <- IntegrateData(anchorset = anchors, dims = 1:dims, k.weight = k.weight, verbose = DEBUG)
        DefaultAssay(integrated) <- "integrated"
        
        # Run the standard workflow for visualization and clustering
        integrated <- ScaleData(integrated, verbose = DEBUG)
        integrated@assays$integrated@data
      } else {
        seur.list[[1]]@assays$RNA@data 
      }
      # integrated@assays
      # NormalizeData(seurDat, normalization.method = "LogNormalize", scale.factor = 10000)
    },
    error = function(e) {
      cat(file = stderr(), paste("\n\n!!!Error during Seurat normalization:\n", e, "\n\n"))
      return(NULL)
    }
  )
  if (is.null(A)) {
    return(NULL)
  }
  # A <- seurDat@assays$integrated@data
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
  splitby <- isolate(input$DE_seuratStandard_splitby)
  
  if (is.null(scEx) | runThis == "") {
    if (DEBUG) {
      cat(file = stderr(), "DE_seuratStandard:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_seuratStandard.RData", list = c(ls()))
  }
  # cp =load(file="~/SCHNAPPsDebug/DE_seuratStandard.RData")
  
  
  
  # # TODO ?? define scaling factor somewhere else???
  # sfactor = max(max(assays(scEx)[["counts"]]),1000)
  retVal <- DE_seuratStandardfunc(scEx = scEx, dims = sDims, 
                                  anchorsF = anchorF, kF = kF, 
                                  k.weight = k.weight, splitby = splitby)
  
  if (is.null(retVal)) {
    showNotification("An error occurred during Seurat normalization, please check console", 
                     id = "DE_seuratError", duration = NULL, type = "error")
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

# Seurat DE_seuratSCTnorm ----

DE_seuratSCTnormButton  <- reactiveVal(
  label = "DE_seuratSCTnormButton",
  value = ""
)


DE_seuratSCTnorm <- reactive({
  if (DEBUG) cat(file = stderr(), "DE_seuratSCTnorm started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_seuratSCTnorm")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_seuratSCTnorm")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_seuratSCTnorm", id = "DE_seuratSCTnorm", duration = NULL)
  }
  
  scEx <- scEx()
  runThis <- DE_seuratSCTnormButton()
  nHVG =  isolate(input$DE_seuratSCTnorm_nHVG)
  var2reg =  isolate(input$DE_seuratSCTnorm_var2reg)

  # if (length(var2reg)<1)
  if (is.null(scEx) | runThis == "") {
    if (DEBUG) {
      cat(file = stderr(), "DE_seuratSCTnorm:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_seuratSCTnorm.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/DE_seuratSCTnorm.RData")
  
  # make sure only factors are used.
  var2reg = var2reg[var2reg %in% names(Filter(is.factor, colData(scEx)))]
  
  # # TODO ?? define scaling factor somewhere else???
  # sfactor = max(max(assays(scEx)[["counts"]]),1000)
  retVal <- DE_seuratSCTnormfunc(scEx = scEx, nHVG, var2reg)
  
  if (is.null(retVal)) {
    showNotification("An error occurred during Seurat normalization, please check console", id = "DE_seuratError", duration = NULL, type = "error")
  }
  
  addClass("updateNormalization", "green")
  exportTestValues(DE_seuratSCTnorm = {
    assays(retVal)[["logcounts"]]
  })
  return(retVal)
})

DE_seuratSCTnormfunc <- function(scEx, nHVG, var2reg) {
  require(Seurat)
  # save(file = "~/SCHNAPPsDebug/DE_seuratSCTnormf.RData", list = c(ls()))
  # cp = load(file="~/SCHNAPPsDebug/DE_seuratSCTnormf.RData")
  
  cellMeta <- colData(scEx)
  # split in different samples
  meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
  # creates object @assays$RNA@data and @assays$RNA@counts
  # integrated <- NULL
  
  choicesVal = names(Filter(is.factor, cellMeta))
  choicesVal = choicesVal[unlist(lapply(choicesVal, FUN = function(x) {length(levels(cellMeta[,x]))>1}))]
  var2reg= var2reg[var2reg %in% choicesVal]
  if (is.null(var2reg)) {
    var2reg = NULL
    meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
  } else if(var2reg == "" ) {
    var2reg = NULL
    meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
  } else {
    meta.data <- as.data.frame(cellMeta[, var2reg, drop = FALSE])
    limitCells = meta.data[,1] %in% levels(meta.data[,1])[table(meta.data[,1]) > 30]
    # we cannot remove cell here because this would change scEX and projections
    # not sure that NA would be a good solution
    # so we are asking to remove the cells manually
    if (sum(limitCells) < ncol(scEx)) {
      errStr = paste("please remove the following cells:\n", 
                     paste(colnames(scEx)[!limitCells],
                           collapse = ", "))
      cat (file = stderr(), errStr)
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification(errStr, id = "DE_seuratError", duration = NULL, type = "error")  
      }
      return(NULL)
    }
    # 
    # scEx = scEx[, limitCells]
    # meta.data = meta.data[limitCells,, drop = FALSE]
  }
  seurDat <- tryCatch(
    {
      seurDat <- CreateSeuratObject(
        counts = assays(scEx)[[1]],
        meta.data = as.data.frame(cellMeta)
      )
    },
    error = function(e) {
      cat(file = stderr(), paste("\n\n!!!Error during Seurat normalization:\n", e, "\n\n"))
      return(NULL)
    }
  )
  # UMI-based normalisation & logTransformation
  seurDat = Seurat::NormalizeData(seurDat)
  
  seurDat <- SCTransform(object = seurDat, 
                         variable.features.n = nHVG, 
                         vars.to.regress = var2reg, #No additional correction, eg. no cell cycle correction
                         verbose = DEBUG) 
  
  
  A <- seurDat[["SCT"]]@scale.data
  scEx_bcnorm <- SingleCellExperiment(
    assay = list(logcounts = as(A, "dgTMatrix")),
    colData = colData(scEx)[colnames(A), , drop = FALSE],
    rowData = rowData(scEx)[rownames(A), , drop = FALSE]
  )
  
  return(scEx_bcnorm)
}



# Seurat LogNorm ----

DE_seuratLogNormButton  <- reactiveVal(
  label = "DE_seuratLogNormButton",
  value = ""
)


DE_seuratLogNorm <- reactive({
  if (DEBUG) cat(file = stderr(), "DE_seuratLogNorm started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_seuratLogNorm")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_seuratLogNorm")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_seuratLogNorm", id = "DE_seuratLogNorm", duration = NULL)
  }
  
  scEx <- scEx()
  runThis <- DE_seuratLogNormButton()
  nHVG = 3000
  var2reg = ""
  if (DEBUG) cat(file = stderr(), "DE_seuratLogNorm started2.\n")
  
  nHVG = isolate(input$DE_seuratLogNorm_nHVG)
  if (DEBUG) cat(file = stderr(), "DE_seuratLogNorm started4.\n")
  
  var2reg = isolate(input$DE_seuratLogNorm_var2reg)
  if (DEBUG) cat(file = stderr(), "DE_seuratLogNorm started5.\n")
  
  # projections = projections()
  
  if (DEBUG) cat(file = stderr(), "DE_seuratLogNorm started3.\n")
  
  # if (length(var2reg)<1)
  if (is.null(scEx) | runThis == "") {
    if (DEBUG) {
      cat(file = stderr(), "DE_seuratStandard:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_seuratStandard.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/DE_seuratStandard.RData")
  
  # make sure only factors are used.
  var2reg = var2reg[var2reg %in% names(Filter(is.factor, colData(scEx)))]
  
  # # TODO ?? define scaling factor somewhere else???
  # sfactor = max(max(assays(scEx)[["counts"]]),1000)
  retVal <- DE_seuratLogNormfunc(scEx = scEx, nHVG, var2reg)
  
  if (is.null(retVal)) {
    showNotification("An error occurred during Seurat normalization, please check console", id = "DE_seuratError", duration = NULL, type = "error")
  }
  
  addClass("updateNormalization", "green")
  exportTestValues(DE_seuratLogNorm = {
    assays(retVal)[["logcounts"]]
  })
  return(retVal)
})

DE_seuratLogNormfunc <- function(scEx, nHVG, var2reg) {
  require(Seurat)
  # save(file = "~/SCHNAPPsDebug/DE_seuratStandardf.RData", list = c(ls()))
  # cp = load(file="~/SCHNAPPsDebug/DE_seuratStandardf.RData")
  
  cellMeta <- colData(scEx)
  # split in different samples
  meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
  # creates object @assays$RNA@data and @assays$RNA@counts
  # integrated <- NULL
  
  choicesVal = names(Filter(is.factor, cellMeta))
  choicesVal = choicesVal[unlist(lapply(choicesVal, FUN = function(x) {length(levels(cellMeta[,x]))>1}))]
  var2reg= var2reg[var2reg %in% choicesVal]
  if (is.null(var2reg)){
    var2reg = NULL
    meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
    
  } else if(var2reg == "" ) {
    var2reg = NULL
    meta.data <- as.data.frame(cellMeta[, "sampleNames", drop = FALSE])
  } else {
    meta.data <- as.data.frame(cellMeta[, var2reg, drop = FALSE])
    limitCells = meta.data[,1] %in% levels(meta.data[,1])[table(meta.data[,1]) > 30]
    # we cannot remove cell here because this would change scEX and projections
    # not sure that NA would be a good solution
    # so we are asking to remove the cells manually
    if (sum(limitCells) < ncol(scEx)) {
      errStr = paste("please remove the following cells:\n", 
                     paste(colnames(scEx)[!limitCells],
                           collapse = ", "))
      cat (file = stderr(), errStr)
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification(errStr, id = "DE_seuratError", duration = NULL, type = "error")  
      }
      return(NULL)
    }
    # 
    # scEx = scEx[, limitCells]
    # meta.data = meta.data[limitCells,, drop = FALSE]
  }
  seurDat <- tryCatch(
    {
      seurDat <- CreateSeuratObject(
        counts = assays(scEx)[[1]],
        meta.data = as.data.frame(cellMeta)
      )
    },
    error = function(e) {
      cat(file = stderr(), paste("\n\n!!!Error during Seurat normalization:\n", e, "\n\n"))
      return(NULL)
    }
  )
  # UMI-based normalisation & logTransformation
  seurDat = Seurat::NormalizeData(seurDat)
  # Finding variable genes
  seurDat = Seurat::FindVariableFeatures(object = seurDat,
                                         selection.method = "vst",
                                         nfeatures = nHVG)
  if (length(var2reg) == 0)
    var2reg = NULL
  # # scaling the data (only variable genes)
  seurDat = Seurat::ScaleData(object = seurDat,
                              vars.to.regress = var2reg)
  
  A <- seurDat[["RNA"]]@scale.data
  scEx_bcnorm <- SingleCellExperiment(
    assay = list(logcounts = as(A, "dgTMatrix")),
    colData = colData(scEx)[colnames(A), , drop = FALSE],
    rowData = rowData(scEx)[rownames(A), , drop = FALSE]
  )
  
  return(scEx_bcnorm)
}

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
  retVal <- DE_logNormalizationGenefunc(scEx, inputGenes)
  
  if (is.null(retVal)) {
    # print error message
  }
  if (DEBUG) {
    cat(file = stderr(), "assign: .schnappsEnv$calculated_DE_geneIds_norm\n")
  }
  .schnappsEnv$calculated_DE_geneIds_norm <- inputGenes
  
  # turn normalization button green
  addClass("updateNormalization", "green")
  
  #used for DGE for some of the methods
  .schnappsEnv$normalizationFactor <- sfactor
  exportTestValues(DE_logGeneNormalization = {
    assays(retVal)[["logcounts"]]
  })
  return(retVal)
})

#' DE_logNormalizationGenefunc
#' actual computation of the normalization as it is done in seurat
DE_logNormalizationGenefunc <- function(scEx, inputGenes) {
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
  if (rlang::is_empty(genesin)) {
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
  slot(x, "x") <- log(1 + slot(x, "x"), base = 2) 
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

# DE_scaterNormalization ----
DE_scaterNormalization <- reactive(label = "scaterNorm", {
  if (DEBUG) cat(file = stderr(), "DE_scaterNormalization started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_scaterNormalization")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_scaterNormalization")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_scaterNormalization", id = "DE_scaterNormalization", duration = NULL)
  }
  
  scEx <- scEx()
  DE_logNormalizationButton()
  
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "DE_scaterNormalization:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_scaterNormalization.RData", list = c(ls()))
  }
  # cp=load(file="~/SCHNAPPsDebug/DE_scaterNormalization.RData")
  
  # TODO ?? define scaling factor somewhere else???
  sfactor <- max(max(assays(scEx)[["counts"]]), 1000)
  retVal <- DE_scaterNormalizationfunc(scEx)
  
  # turn normalization button green
  addClass("updateNormalization", "green")
  
  .schnappsEnv$normalizationFactor <- sfactor
  exportTestValues(DE_scaterNormalization = {
    assays(retVal)[["logcounts"]]
  })
  return(retVal)
  
})

DE_scaterNormalizationfunc <- function(scEx) {
  require(scater)
  require(scran)
  if (DEBUG) cat(file = stderr(), "DE_scaterNormalizationfunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "DE_scaterNormalizationfunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_scaterNormalizationfunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_scaterNormalizationfunc", id = "DE_scaterNormalizationfunc", duration = NULL)
  }
  # sampinfo = colData(scEx)$sampleNames
  # genes2use = rownames(scEx)
  # ta = table(sampinfo)
  # ta = ta[ta>0]
  # mtab = min(ta)
  # stp = round((mtab - 21) /10)
  # scaterReads <- scran::computeSumFactors(scEx, sizes = seq(21, mtab, stp), 
  #                                         clusters = sampinfo, 
  #                                         subset.row = genes2use)
  # # scaterReads <- scran::computeSumFactors(scEx, clusters = sampinfo)
  # # dt = data.frame(x=librarySizeFactors(scaterReads), y=sizeFactors(scaterReads))
  # # p = plot_ly(dt, x=~x,y=~y)  %>% add_markers()
  # # p <- layout(p, xaxis = list(type = "log"),
  # #             yaxis = list(type = "log"))
  # # 
  # # p
  # scaterReads <- normalize(scaterReads)
  # assays(scaterReads)['counts'] = NULL
  # return(scaterReads)
  clusters <- quickCluster(scEx)
  scEx <- computeSumFactors(scEx, clusters=clusters)
  
  sce <- scuttle::logNormCounts(scEx)
  assays(sce)['counts'] = NULL
  return(sce)
}



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
  sfactor = isolate(input$DE_logNormalization_sf)
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "DE_logNormalization:NULL\n")
    }
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/DE_logNormalization.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/DE_logNormalization.RData")
  
  # TODO ?? define scaling factor somewhere else???
  # sfactor <- max(max(assays(scEx)[["counts"]]), 1000)
  retVal <- DE_logNormalizationfunc(scEx, sfactor)
  
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
DE_logNormalizationfunc <- function(scEx, sfactor) {
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
  if (sfactor <=0){
    sfactor = min(A@x[A@x>0])
  }
  A@x = A@x / sfactor
  scEx_bcnorm <- SingleCellExperiment(
    assay = list(logcounts = as(A, "dgTMatrix")),
    colData = colData(scEx),
    rowData = rowData(scEx)
  )
  
  x <- uniqTsparse(assays(scEx_bcnorm)[[1]])
  slot(x, "x") <- log(1 + slot(x, "x"), base = 2) 
  assays(scEx_bcnorm)[[1]] <- x
  return(scEx_bcnorm)
}
