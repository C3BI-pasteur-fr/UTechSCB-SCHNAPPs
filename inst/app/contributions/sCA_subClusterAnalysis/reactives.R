#' sCA_selectedDge
#' stores table with differentially expressed genes
#' is used to export file and write csv file
sCA_selectedDge <- reactiveValues(
  sCA_dgeTable = data.frame()
)

#' sCA_getCells
#' get list of cells from selctions
sCA_getCells <- function(projections, cl1, db1, db2, db1x, db1y) {
  if (DEBUG) cat(file = stderr(), "sCA_getCells started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_getCells")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_getCells")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_getCells", id = "sCA_getCells", duration = NULL)
  }
  
  # save(file = "~/SCHNAPPsDebug/sCA_getCells.RData", list = c( ls()))
  # cp =load("~/SCHNAPPsDebug/sCA_getCells.RData")
  # dbCluster = projections$dbCluster
  subsetData <- projections[cl1,]
  # db1x = db1$mapping$x
  # db1y = db1$mapping$y
  db1$mapping$x = db1x
  db1$mapping$y = db1y
  if (is(subsetData[,db1x], "logical")) {
    subsetData[,db1x] = as.numeric(subsetData[,db1x]) + 1
  }
  if (is(subsetData[, db1y], "logical")) {
    subsetData[, db1y] <- as.numeric(subsetData[, db1y]) + 1
  }
  
  
  # factors and brushedPoints don't work together.
  # so we change a factor into a numeric
  # TODO WHY is discrete_limits set/misused?????
  # ggplot is not displaying levels for which there are no values
  # thus, the numbering might be off 
  if (is(subsetData[, db1x], "factor")) {
    subsetData[, db1x] <- as.numeric(subsetData[, db1x])
    db1$domain$discrete_limits <- NULL
  }
  if (is(subsetData[, db1y], "factor")) {
    subsetData[, db1y] <- as.numeric(subsetData[, db1y])
    db1$domain$discrete_limits <- NULL
  }
  db1$domain$discrete_limits <- NULL
  cells.1 <- rownames(shiny::brushedPoints(df = subsetData, brush = db1))
  
  cells.2 = c()
  if (!is.null(db2)) {
    db2$domain$discrete_limits <- NULL
    db2x = db1x
    db2y = db1y
    db2$mapping$x = db1x
    db2$mapping$y = db1y
    
    # factors and brushedPoints don't work together.
    # so we change a factor into a numeric
    if (is(subsetData[, db2x], "factor")) {
      subsetData[, db2x] <- as.numeric(subsetData[, db2x])
      db2$domain$discrete_limits <- NULL
    }
    if (is(subsetData[, db2y], "factor")) {
      subsetData[, db2y] <- as.numeric(subsetData[, db2y])
      db2$domain$discrete_limits <- NULL
    }
    cells.2 <- rownames(shiny::brushedPoints(df = subsetData, brush = db2))
  } else {
    cells.2 <- rownames(subsetData)[!rownames(subsetData) %in% cells.1]
  }
  
  # cells.2 <- rownames(shiny::brushedPoints(df = subsetData, brush = db2))
  retVal <- list(c1 = cells.1, c2 = cells.2)
  # save(file = "~/SCHNAPPsDebug/sCA_getCells.RData", list = c( ls()))
  # cp = load("~/SCHNAPPsDebug/sCA_getCells.RData")
  return(retVal)
}

#' define different methods for calculating diff. expressed genes
#' first entry is displayed in Radio box, second is function to be called.
#' to set sCA_dgeRadioButton
myDiffExpFunctions <- list(
  c("Chi-square test of an estimated binomial distribution", "sCA_dge_CellViewfunc"),
  c("t-test", "sCA_dge_ttest"),
  c("DESeq2", "sCA_dge_deseq2"),
  c("seurat:wilcox", "sCA_dge_s_wilcox"),
  c("seurat:bimod", "sCA_dge_s_bimod"),
  c("seurat:t-test", "sCA_dge_s_t"),
  c("seurat:LR", "sCA_dge_s_LR"),
  c("seurat:neg-binomial", "sCA_dge_s_negbinom"),
  c("seurat:Poisson", "sCA_dge_s_poisson")
)

if("scDEA" %in% rownames(installed.packages())){
  myDiffExpFunctions[[length(myDiffExpFunctions)+1]] = c("scDEA - !! this might take hours !!","sCA_scDEA")
} else {
  cat(file = stderr(), "Please install scDEA:
      devtools::install_github('nghiavtr/BPSC')
      BiocManager::install('DEsingle')
      devtools::install_github('nghiavtr/BPSC')
      BiocManager::install('DESeq2')
      BiocManager::install('edgeR')
      BiocManager::install('MAST')
      BiocManager::install('monocle')
      BiocManager::install('limma')
      BiocManager::install('Seurat')
      devtools::install_github('statOmics/zingeR')
      BiocManager::install('SingleCellExperiment')
      BiocManager::install('scater')
      devtools::install_github('Zhangxf-ccnu/scDEA')
  ")
}




# ! TODO
# some dge functions are based on counts most on log-counts
# this is hard-coded here, but should be a parameter somehow like myDiffExpFunctions
# such it can be used by contributed functions

#' Seurat FindMarkers
#'
#' cellMeta = colData(scEx)
sCA_seuratFindMarkers <- function(scEx, scEx_logMat, cells.1, cells.2, test="wilcox", normFact = 1){
  if (DEBUG) cat(file = stderr(), "sCA_seuratFindMarkers started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_seuratFindMarkers")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_seuratFindMarkers")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_seuratFindMarkers", id = "sCA_seuratFindMarkers", duration = NULL)
  }
  if (!"DESeq2" %in% rownames(installed.packages())) {
    warning("Please install DESeq2 - learn more at https://bioconductor.org/packages/release/bioc/html/DESeq2.html")
    showNotification("Please install DESeq2", id = "sCA_dge_deseq2NOTFOUND", duration = NULL, type = "error")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_seuratFindMarkers.RData", list = c(ls()))
  }
  #cp = load(file='~/SCHNAPPsDebug/sCA_seuratFindMarkers.RData')
  cellMeta <- colData(scEx)
  rData <- rowData(scEx)
  meta.data <- cellMeta[, "sampleNames", drop = FALSE]
  # creates object @assays$RNA@data and @assays$RNA@counts
  # browser()
  if(!is.null(scEx_logMat)){
    seurDat <- tryCatch.W.E(Seurat::CreateSeuratObject(
      counts = scEx_logMat, 
      meta.data = as.data.frame(meta.data)
    ))
  }else{
    seurDat <- tryCatch.W.E(Seurat::CreateSeuratObject(
      counts = assays(scEx)[[1]],
      meta.data = as.data.frame(meta.data)
    ))
  }
  
  if (is(seurDat$value, "Seurat")) {
    seurDat = seurDat$value
  } else {
    cat(file = stderr(), "something went wrong with seurat CreateSeuratObject\n")
    cat(file = stderr(), as.character(seurDat$value))
    cat(file = stderr(), "\n")
    return(NULL)
  }
  
  
  # change "_" to "-"
  rowNameLookup = data.frame(org = rownames(scEx), seurat = stringr::str_replace_all(rownames(scEx),"_","-"))
  rownames(rowNameLookup) = rowNameLookup$org
  # rownames(scEx) = stringr::str_replace_all(rownames(scEx),"_","-")
  # we remove e.g. "genes" from total seq (CD3-TotalSeqB)
  useGenes = rowNameLookup[which(rowNameLookup$seurat %in% Features(seurDat)),"seurat"]
  rownames(scEx) = rowNameLookup[rownames(scEx),"seurat"]
  
  
  
  # seurDat@assays$RNA$counts = as(assays(scEx)[[1]], "CsparseMatrix")[useGenes,]
  # seurDat@assays$RNA$scale.data = as.matrix(seurDat@assays$RNA$counts)
  
  # not sure we need the normalization factor
  # markers <- Seurat::FindMarkers(seurDat@assays$RNA@data/normFact, 
  require(Seurat)
  markers <- tryCatch.W.E(
    Seurat::FindMarkers(seurDat, 
                        ident.1 = cells.1,
                        ident.2 = cells.2,
                        min.pct = 0.00000000001,
                        test.use = test,
                        logfc.threshold = 0.000001,
                        slot = Layers(seurDat)[1]
                        # test.use = "wilcox" # p_val  avg_logFC pct.1 pct.2    p_val_adj
                        # test.use = "bimod"  # p_val  avg_logFC pct.1 pct.2    p_val_adj
                        # test.use = "roc"    # myAUC   avg_diff power pct.1 pct.2
                        # test.use = "t"      # p_val avg_logFC pct.1 pct.2 p_val_adj
                        # test.use = "negbinom" # needs UMI; p_val  avg_logFC pct.1 pct.2    p_val_adj
                        # test.use = "poisson"  # needs UMI; p_val  avg_logFC pct.1 pct.2    p_val_adj
                        # test.use = "LR"     # p_val  avg_logFC pct.1 pct.2 p_val_adj
                        # test.use = "MAST" # not working: Assay in position 1, with name et is unlogged. Set `check_sanity = FALSE` to override and then proceed with caution.
                        # test.use = "DESeq2" # needs UMI # done separately because the estimating process isn't working with 0s
    )
  )
  
  if (is(markers$value , "data.frame")) {
    markers = markers$value
    rownames(rowNameLookup) = rowNameLookup$seurat
    rownames(markers) = rowNameLookup[rownames(markers),"org"]
  }else {
    cat(file = stderr(), "something went wrong with seurat find markers\n")
    cat(file = stderr(), as.character(markers$value))
    cat(file = stderr(), "\n")
    return(NULL)
  }
  if(is.null(markers)) {
    return(NULL)
  }
  if(length(markers)<1) {
    return(NULL)
  }
  
  if (nrow(markers) > 0) {
    markers$symbol <- rData[rownames(markers), "symbol"]
  }
  if ("Description" %in% colnames(rData)) {
    markers$Description <- rData[rownames(markers), "Description"]
  }
  return(markers)
}
#   sCA_seuratFindMarkers_m = memoise::memoise(sCA_seuratFindMarkers,cache=do.call(cachem::cache_disk,.schnappsEnv$cacheDir))
sCA_seuratFindMarkers_m = sCA_seuratFindMarkers


sCA_dge_s_wilcox <- function(scEx_log, scEx_logMat, cells.1, cells.2){
  normFact = .schnappsEnv$normalizationFactor
  sCA_seuratFindMarkers_m(scEx_log, scEx_logMat, cells.1, cells.2, test="wilcox", normFact) 
}
#Differential expression with BPCells currently only supports the 'wilcox' method. Please rerun with test.use = 'wilcox'
sCA_dge_s_bimod <- function(scEx_log, scEx_logMat, cells.1, cells.2){
  normFact = .schnappsEnv$normalizationFactor
  sCA_seuratFindMarkers_m(scEx_log, scEx_logMat = NULL, cells.1, cells.2, test="bimod") 
}
sCA_dge_s_t <- function(scEx_log, scEx_logMat, cells.1, cells.2){
  if (.schnappsEnv$DEBUG) cat(file = stderr(), paste0("sCA_dge_s_t: \n"))
  normFact = .schnappsEnv$normalizationFactor
  sCA_seuratFindMarkers_m(scEx_log, scEx_logMat = NULL, cells.1, cells.2, test="t", normFact) 
}
sCA_dge_s_LR<- function(scEx_log, scEx_logMat, cells.1, cells.2){
  normFact = .schnappsEnv$normalizationFactor
  sCA_seuratFindMarkers_m(scEx_log, scEx_logMat = NULL, cells.1, cells.2, test="LR", normFact)
}
sCA_dge_s_negbinom <- function(scEx_log, scEx_logMat, cells.1, cells.2){
  sCA_seuratFindMarkers_m(scEx_log, scEx_logMat = NULL, cells.1, cells.2, test="negbinom", normFact = 1)
}
sCA_dge_s_poisson <- function(scEx_log, scEx_logMat, cells.1, cells.2){
  sCA_seuratFindMarkers_m(scEx_log, scEx_logMat = NULL, cells.1, cells.2, test="poisson", normFact = 1) 
}

#' scDEA
#'  



sCA_scDEA <- function(scEx_log, scEx_logMat, cells.1, cells.2){
  withWarnings <- function(expr) {
    wHandler <- function(w) {
      if (DEBUG) {
        cat(file = stderr(), paste("runscDEA created a warning:", w, "\n"))
      }
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification(
          "Problem with scDEA?",
          type = "warning",
          id = "scDEAwarning",
          duration = NULL
        )
      }
      invokeRestart("muffleWarning")
    }
    eHandler <- function(e) {
      cat(file = stderr(), paste("error in scDEA:", e))
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification(
          paste("Problem with scDEA, probably not enough cells?", e),
          type = "warning",
          id = "scDEAwarning",
          duration = NULL
        )
      }
      cat(file = stderr(), "scDEA FAILED!!!\n")
      return(NULL)
    }
    val <- withCallingHandlers(expr, warning = wHandler, error = eHandler)
    return(val)
  }
  
  if (DEBUG) cat(file = stderr(), "sCA_scDEA started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_scDEA")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_scDEA")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_scDEA - this might take several hours", id = "sCA_scDEA", duration = NULL)
  }
  if (!"scDEA" %in% rownames(installed.packages())) {
    warning("Please install scDEA - learn more at https://github.com/Zhangxf-ccnu/scDEA")
    showNotification("Please install scDEA", id = "sCA_dge_deseq2NOTFOUND", duration = NULL, type = "error")
    return(NULL)
  }
  require(scDEA)
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_scDEA.RData", list = c(ls()))
  }
  # cp = load(file='~/SCHNAPPsDebug/sCA_scDEA.RData')
  
  dups = c(cells.1, cells.2)[duplicated(c(cells.1, cells.2))]
  if(length(dups) > 0){
    cells.1 = cells.1[!cells.1 %in% dups]
    cells.2 = cells.2[!cells.2 %in% dups]
    warning("Duplicated cells, using only unique ones")
    showNotification("Duplicated cells, using only unique ones", id = "sCA_dge_deseq2NOTUnique", duration = NULL, type = "warning")
  }
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  group.info$wellKey <- rownames(x = group.info)
  # TODO how to handle data / transformation ?
  # it uses non-transformed org data. What happens if we loaded normalized data?
  # can back tronsform the data in counts?
  if (is.null(.schnappsEnv$normalizationFactor)) {
    .schnappsEnv$normalizationFactor = 1
  }
  counts = assays(scEx_log)[[1]][,rownames(group.info)]
  
  Pvals <- withWarnings(
    suppressMessages(scDEA::scDEA_individual_methods(raw.count = as.matrix(counts), 
                                                     cell.label = group.info$group,
                                                     BPSC = isolate(input$scDEA_BPSC), 
                                                     DEsingle = isolate(input$scDEA_BPSC), 
                                                     DESeq2 = isolate(input$scDEA_DESeq2), 
                                                     edgeR = isolate(input$scDEA_edgeR), 
                                                     MAST = isolate(input$scDEA_MAST), 
                                                     monocle = isolate(input$scDEA_monocle), 
                                                     scDD = isolate(input$scDEA_scDD),
                                                     Ttest = isolate(input$scDEA_Ttest),
                                                     
                                                     
                                                     
                                                     Wilcoxon = isolate(input$scDEA_Wilcoxon), 
                                                     limma = isolate(input$scDEA_limma), 
                                                     Seurat = isolate(input$scDEA_Seurat),
                                                     zingeR.edgeR = isolate(input$scDEA_zingeR.edgeR),
                                                     
                                                     BPSC.parallel = isolate(input$scDEA_parallel),
                                                     DEsingle.parallel = isolate(input$scDEA_parallel),
                                                     DESeq2.parallel = isolate(input$scDEA_parallel),
                                                     MAST.parallel = isolate(input$scDEA_parallel),
                                                     monocle.cores = ifelse(isolate(input$scDEA_parallel),
                                                                            max(floor(parallel::detectCores()/2),1), #assure min 1 core used and not all, i.e. about 50% for memory reasons
                                                                            1)
                                                     
    ))
  )
  if(is.null(Pvals)) return(NULL)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_scDEA.res.RData", list = c(ls()))
  }
  # cp = load(file='~/SCHNAPPsDebug/sCA_scDEA.res.RData')
  allOnes = which(apply(Pvals,1,FUN=function(x)all(x==1)))
  if(length(allOnes)==0) allOnes = -(1:nrow(Pvals))
  PvalsO = Pvals[-allOnes,]
  combination.Pvals <- lancaster.combination(PvalsO, weight = TRUE, trimmed = 0.2)
  
  adjusted.Pvals <- scDEA.p.adjust(combination.Pvals, adjusted.method = "bonferroni")
  res = as.data.frame(PvalsO)
  res$p_val_adj = adjusted.Pvals
  res$p_val = combination.Pvals
  res$avg_log2FC = log2(apply(counts[-allOnes,cells.1],1,mean) / apply(counts[-allOnes,cells.2],1,mean))
  newMax = max(abs(res$avg_log2FC[!is.infinite(res$avg_log2FC)])) * 1.5
  res$avg_log2FC[res$avg_log2FC==-Inf] = -newMax
  res$avg_log2FC[res$avg_log2FC==Inf] = newMax
  
  
  return(res)
  
}

#' sCA_dge_deseq2
#' calculate dge using DESeq2
runDESEQ2 <- function(data.use, group.info) {
  if (DEBUG) cat(file = stderr(), "DESeq2 setup.\n")
  
  # TODO DESeqParallel 
  
  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data.use,
    colData = group.info,
    design = ~group
  )
  if (DEBUG) cat(file = stderr(), "DESeq2 estimate size factors\n")
  dds1 <- DESeq2::estimateSizeFactors(object = dds1, type = "poscounts")
  if (DEBUG) cat(file = stderr(), "DESeq2 estimateDispersions\n")
  dds1 <- DESeq2::estimateDispersions(object = dds1, fitType = "local")
  if (DEBUG) cat(file = stderr(), "DESeq2 nbinomWaldTest\n")
  dds1 <- DESeq2::nbinomWaldTest(object = dds1)
  if (DEBUG) cat(file = stderr(), "DESeq2 results\n")
  res <- DESeq2::results(
    object = dds1,
    contrast = c("group", "Group1", "Group2"),
    alpha = 0.05, parallel = T
  )
  if (DEBUG) cat(file = stderr(), "DESeq2 done\n")
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
  res$padj[is.na(res$padj)] <- 1
  res$pvalue[is.na(res$pvalue)] <- 1
  return(res)
}
#   runDESEQ2_m <- memoise::memoise(runDESEQ2,cache=do.call(cachem::cache_disk,.schnappsEnv$cacheDir))
sCA_seuratFindMarkers_m = sCA_seuratFindMarkers
runDESEQ2_m <- runDESEQ2


sCA_dge_deseq2 <- function(scEx_log, scEx_logMat, cells.1, cells.2) {
  if (DEBUG) cat(file = stderr(), "sCA_dge_deseq2 started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dge_deseq2")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_dge_deseq2")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dge_deseq2", id = "sCA_dge_deseq2", duration = NULL)
  }
  if (!"DESeq2" %in% rownames(installed.packages())) {
    warning("Please install DESeq2 - learn more at https://bioconductor.org/packages/release/bioc/html/DESeq2.html")
    showNotification("Please install DESeq2", id = "sCA_dge_deseq2NOTFOUND", duration = NULL, type = "error")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dge_deseq2.RData", list = c(ls()))
  }
  # cp = load(file='~/SCHNAPPsDebug/sCA_dge_deseq2.RData')
  
  dups = c(cells.1, cells.2)[duplicated(c(cells.1, cells.2))]
  if(length(dups) > 0){
    cells.1 = cells.1[!cells.1 %in% dups]
    cells.2 = cells.2[!cells.2 %in% dups]
    warning("Duplicated cells, using only unique ones")
    showNotification("Duplicated cells, using only unique ones", id = "sCA_dge_deseq2NOTUnique", duration = NULL, type = "warning")
  }
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  group.info$wellKey <- rownames(x = group.info)
  # TODO how to handle data / transformation ?
  # it uses non-transformed org data. What happens if we loaded normalized data?
  # can back tronsform the data in counts?
  if (is.null(.schnappsEnv$normalizationFactor)) {
    .schnappsEnv$normalizationFactor = 1
  }
  data.use = assays(scEx_log)[[1]][,rownames(group.info)]
  res = runDESEQ2_m(data.use, group.info) 
  featureData <- rowData(scEx_log)
  
  to.return <- data.frame(
    p_val = res$padj,
    avg_diff = res$log2FoldChange, row.names = rownames(res), symbol = rownames(res)
  )
  
  
  return(to.return)
}


#' sCA_dge_CellViewfunc
#' calculate differentically expressed genes given 2 sets of cells
sCA_dge_CellViewfunc <- function(scEx_log, scEx_logMat, cells.1, cells.2) {
  if (DEBUG) cat(file = stderr(), "sCA_dge_CellViewfunc started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dge_CellViewfunc")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_dge_CellViewfunc")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dge_CellViewfunc", id = "sCA_dge_CellViewfunc", duration = NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dge_CellViewfunc.RData", list = c(ls()))
  }
  # cp =load(file='~/SCHNAPPsDebug/sCA_dge_CellViewfunc.RData')
  
  featureData <- rowData(scEx_log)
  scEx_log <- as.matrix(assays(scEx_log)[[1]])
  subsetExpression <- scEx_log[complete.cases(scEx_log[, union(cells.1, cells.2)]), ]
  genes.use <- rownames(subsetExpression)
  # expMean exponential mean
  dat <- subsetExpression[genes.use, cells.1]
  data.1 <- apply(dat, 1, function(x) expMean(x, .schnappsEnv$normalizationFactor))
  dat <- subsetExpression[genes.use, cells.2]
  data.2 <- apply(dat, 1, function(x) expMean(x, .schnappsEnv$normalizationFactor))
  total.diff <- (data.1 - data.2)
  
  genes.diff <- names(which(abs(total.diff) > .2))
  genes.use <- ainb(genes.use, genes.diff)
  
  retVal <-
    DiffExpTest(subsetExpression, cells.1, cells.2, genes.use = genes.use)
  if(is.null(retVal)) return(NULL)
  if(nrow(retVal)==0) return(NULL)
  retVal[is.na(retVal[, "p_val"]), ] <- 1
  retVal[, "avg_diff"] <- total.diff[rownames(retVal)]
  retVal$symbol <-
    featureData[rownames(retVal), "symbol"]
  return(retVal)
}

#' sCA_dge_ttest
#' t-test on selected cells
sCA_dge_ttest <- function(scEx_log, scEx_logMat, cells.1, cells.2) {
  if (DEBUG) cat(file = stderr(), "sCA_dge_ttest started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dge_ttest")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_dge_ttest")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dge_ttest", id = "sCA_dge_ttest", duration = NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dge_ttest.RData", list = c(ls()))
  }
  #cp = load(file='~/SCHNAPPsDebug/sCA_dge_ttest.RData')
  # browser()
  featureData <- rowData(scEx_log)
  scEx_log <- as.matrix(assays(scEx_log)[[1]])
  subsetExpression <- scEx_log[complete.cases(scEx_log[, union(cells.1, cells.2)]), ]
  genes.use <- rownames(subsetExpression)
  
  if(!exists("t.test.cores", envir = .schnappsEnv)) .schnappsEnv$t.test.cores = 1
  cl <- makeCluster(.schnappsEnv$t.test.cores)
  
  p_val <- parApply(cl,subsetExpression, 1, cells.1=cells.1,cells.2=cells.2,function(x, cells.1, cells.2) t.test(x[cells.1], x[cells.2])$p.value)
  # p_val <- apply(subsetExpression, 1, function(x) t.test(x[cells.1], x[cells.2])$p.value)
  p_val[is.na(p_val)] <- 1
  dat <- subsetExpression[genes.use, cells.1]
  normFact = 1000
  if(!is.null(.schnappsEnv$normalizationFactor))
    if(.schnappsEnv$normalizationFactor == 0) normFact = 1000
  clusterExport(cl, "expMean")
  data.1 <- parApply(cl, dat, 1, normFact=normFact, function(x, normFact) expMean(x, normFactor = normFact))
  dat <- subsetExpression[genes.use, cells.2]
  data.2 <- parApply(cl, dat, 1, normFact=normFact, function(x, normFact) expMean(x, normFactor = normFact))
  avg_diff <- (data.1 - data.2)
  stopCluster(cl)
  
  retVal <- data.frame(p_val = p_val, avg_diff = avg_diff, symbol = featureData[names(p_val), "symbol"])
  return(retVal)
}

# sCA_dge reactive ----
#' manage calculation for differential expression analysis
sCA_dge <- reactive({
  if (DEBUG) cat(file = stderr(), "sCA_dge started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "sCA_dge")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "sCA_dge")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("sCA_dge", id = "sCA_dge", duration = NULL)
    removeNotification(id = "dgewarning")
  }
  
  selectedCells <- isolate(Liana_dataInput())
  
  scEx_log <- scEx_log()
  scEx <- scEx()
  scEx_logMat <- BPCellsLog()
  scExMat <- BPCellsCounts()
  
  projections <- projections()
  # cl1 <- input$sCA_dgeClustersSelection
  # selectedCells <- isolate(dePanelCellSelection())
  # cellNs <- isolate(selectedCells$cellNames())
  # sampdesc <- isolate(selectedCells$selectionDescription())
  
  clicked <- input$updateDGEParameters
  selectedCells <- isolate(sCA_dataInp())
  cellNs <- isolate(selectedCells$cellNames())
  sampdesc <- isolate(selectedCells$selectionDescription())
  prj <- isolate(selectedCells$ProjectionUsed())
  prjVals <- isolate(selectedCells$ProjectionValsUsed())
  db1x <- isolate(input$sCA_subscluster_x1)
  db1y <- isolate(input$sCA_subscluster_y1)
  db1 <- isolate(input$db1)
  db2 <- isolate(input$db2)
  method <- isolate(input$sCA_dgeRadioButton)
  
  if (is.null(scEx_log) | is.null(projections)  || is.null(db1)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sCA_dge.RData", list = c(ls(), ".schnappsEnv"))
  }
  # browser()
  if(method==""){
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("dge: NULL,did you select a method?", id = "dgewarning", 
                       duration = NULL, type = "error")
    }
    return(NULL)
  }
  # cp = load(file='~/SCHNAPPsDebug/sCA_dge.RData')
  # require(archivist)
  # library(tools)
  # lazyLoad = local({load("~/SCHNAPPsDebug/sCA_dge.RData"); environment()})
  # tools:::makeLazyLoadDB(lazyLoad, "Huge")
  # lazyLoad("Huge")
  # objNames <- ls()
  # deepDebug()
  methodIdx <- ceiling(which(unlist(.schnappsEnv$diffExpFunctions) == method) / 2)
  dgeFunc <- .schnappsEnv$diffExpFunctions[[methodIdx]][2]
  gCells <- sCA_getCells(projections, cl1 = cellNs, db1, db2, db1x, db1y)
  
  # in case we need counts and not normalized counts
  if (dgeFunc %in% c("sCA_dge_deseq2", "sCA_dge_s_negbinom", "sCA_dge_s_poisson", "sCA_scDEA")) {
    scEx_log = scEx
    scEx_logMat = scExMat
  }
  # retVal <- do.call(dgeFunc, args = list(
  #   scEx_log = scEx_log,
  #   cells.1 = gCells$c1, cells.2 = gCells$c2
  # ))
  
  withWarnings <- function(expr) {
    wHandler <- function(w) {
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification(as.character(w), id = "dgeFuncWarning", duration = NULL,type = "warning")
      }
      cat(file = stderr(), "something went wrong with dge\n")
      cat(file = stderr(), as.character(w))
      cat(file = stderr(), "\n")
      # at some point I introduced return NULL, but I don\t remember why. Now, there is warning issued
      # about 1GB memory being used but the function still returns correctly. Thus removing and this has
      # to be handled later
      # return(NULL)
      invokeRestart("muffleWarning")
    }
    eHandler <- function(e) {
      require(Seurat)
      if(is.null(e)) e = "NULL"
      if (!is.null(getDefaultReactiveDomain()) & !dgeFunc ==  "sCA_scDEA") {
        showNotification(e, id = "dgeFuncError", duration = NULL,type = "error")
      }
      cat(file = stderr(), "something went wrong with dge\n")
      cat(file = stderr(), as.character(e))
      cat(file = stderr(), "\n")
      return(NULL)
    }
    val <- withCallingHandlers(expr, warning = wHandler, error = eHandler)
    return(val)
  }
  
  # browser()
  
  retVal <- withWarnings({
    do.call(dgeFunc, args = list(
      scEx_log = scEx_log,
      scEx_logMat = scEx_logMat,
      cells.1 = gCells$c1, cells.2 = gCells$c2
    ))})
  # browser()
  if(is.null(retVal)) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("dge: NULL,did you install everything?", id = "dgewarning", duration = 10, type = "warning")
    }
    return(NULL)
  }
  
  if (nrow(retVal) == 0) {
    if (DEBUG) cat(file = stderr(), "dge: nothing found\n")
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("dge: nothing found", id = "dgewarning", duration = 10, type = "warning")
    }
  }
  
  # update reactiveValue
  sCA_selectedDge$sCA_dgeTable <- retVal
  
  setRedGreenButton(
    vars = list(
      # c("sCA_dataInpSelected_cells", isolate(sCA_dataInp()$selectedCells())),
      c("db1", isolate(input$db1)),
      c("db2", isolate(input$db2)),
      c("sCA_dgeRadioButton", isolate(input$sCA_dgeRadioButton)),
      c("sCA_dataInput-Mod_PPGrp", prjVals),
      c("sCA_dataInput-Mod_clusterPP", prj)
      
    ),
    button = "updateDGEParameters"
  )
  
  exportTestValues(sCA_dge = {
    retVal
  })
  return(retVal)
})

# using these global variables allows us to store a set values even when the projections are changing
.schnappsEnv$subClusterDim1 <- "PC1"
.schnappsEnv$subClusterDim2 <- "PC2"
.schnappsEnv$subClusterClusters <- NULL


#' subCluster2Dplot
#' plots cells in 2D for the subcluster anlaysis. The handling of the selection is done
#' outside this function
subCluster2Dplot <- function() {
  if (DEBUG) cat(file = stderr(), "subCluster2Dplot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "subCluster2Dplot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "subCluster2Dplot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("subCluster2Dplot", id = "subCluster2Dplot", duration = NULL)
  }
  
  renderPlot({
    if (DEBUG) cat(file = stderr(), "output$sCA_dge_plot2\n")
    
    projections <- projections()
    x1 <- input$sCA_subscluster_x1
    y1 <- input$sCA_subscluster_y1
    # c1 <- input$sCA_dgeClustersSelection
    # selectedCells <- isolate(dePanelCellSelection())
    # cellNs <- isolate(selectedCells$cellNames())
    # sampdesc <- isolate(selectedCells$selectionDescription())
    selectedCells <- sCA_dataInp()
    if(is.null(selectedCells)) return(NULL)
    cellNs <- selectedCells$cellNames()
    sampdesc <- selectedCells$selectionDescription()
    prjs <- selectedCells$ProjectionUsed()
    sampCol <- projectionColors$sampleNames
    ccols <- projectionColors$dbCluster
    
    if (is.null(projections)) {
      return(NULL)
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/sCA_dge_plot2.RData", list = c(ls()))
    }
    # cp = load(file="~/SCHNAPPsDebug/sCA_dge_plot2.RData")
    
    subsetData <- projections[cellNs,]
    #     xAxis <- list(
    #       title = x1,
    #       titlefont = f
    #     )
    #     yAxis <- list(
    #       title = y1,
    #       titlefont = f
    #     )
    #     
    #     p1 <- plotly::plot_ly(
    #       data = subsetData, source = "subset",
    #       key = rownames(subsetData)
    #     ) %>%
    #       add_trace(
    #         x = ~ get(x1),
    #         # x = ~ get(dimX),
    #         # y = ~ get(dimY),
    #         y = ~ get(y1),
    #         type = "scatter", mode = "markers",
    #         text = ~ paste(1:nrow(subsetData), " ", rownames(subsetData), "<br />", subsetData$exprs),
    #         # color = ~ get(dimCol),
    #         colors = colors,
    #         showlegend = TRUE
    #       ) %>%
    #       layout(
    #         xaxis = xAxis,
    #         yaxis = yAxis,
    #         # title = "gtitle",
    #         dragmode = "select"
    #       )
    # p1
    prj = subsetData[,prjs]
    mycolPal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(
      n = 12, name =
        "Paired"
    ))(length(levels(prj)))
    
    if (prjs == "sampleNames") {
      mycolPal <- sampCol
    }
    if (prjs == "dbCluster") {
      mycolPal <- ccols
    }
    
    p1 <- subsetData %>%
      ggplot(
        # this version causes a crash with unknown column selected
        aes(x = .data[[x1]], y = .data[[y1]]),
        # aes_string(x = x1, y = y1),
        colour = mycolPal[subsetData[,prjs]]
      ) +
      geom_point(colour = mycolPal[subsetData[,prjs]]) +
      geom_point(
        shape = 1,
        size = 4,
        colour = mycolPal[subsetData[,prjs]]
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          size = 12,
          vjust = 0.5
        ),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16),
        legend.position = "none"
      ) + ggtitle(sampdesc) 
    if (is.factor(subsetData[,x1])) {
      p1 <- p1 + scale_x_discrete(drop=FALSE) 
    }
    if (is.factor(subsetData[,y1])) {
      p1 <- p1 + scale_y_discrete(drop=FALSE) 
    }
    
    # + scale_y_discrete(drop=FALSE)
    p1
  })
}

# save to history violoin observer ----
observe({
  clicked  = input$save2HistVolc
  if (DEBUG) cat(file = stderr(), "observe input$save2HistVolc \n")
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
  add2history(type = "renderPlot", 
              input = isolate( reactiveValuesToList(input)), 
              comment = "volcano plot",  
              plotData = .schnappsEnv[["sCA_volcanoPlot"]])
  
})


