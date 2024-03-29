
# panelplotFactFunc ----
# scEx_log singlecell Experiment object
# projections as used in schnapps
# genesin gene names to be plotted
# dimx4, dimy4 dimensions to be plotted on
# sameScale True/False
# nCol number of columns for final plot
# sampdes header for plot
# cellNs cell names to be used

require(dplyr)
require(BiocParallel)
require(SingleCellExperiment)
require(BPCells)
require(cowplot)
require(ggpubr)

panelPlotFactFunc <- function(scEx_log, projections, factsin, dimx4, dimy4, dimCol, sameScale, nCol, sampdesc, cellNs,
                              lowCol = "blue", highCol = "red", midCol = "white",
                              midFunc = function(x){(max(x)-min(x))/2}, applyPvalue=FALSE,
                              projectionColors=NULL,
                              .schnappsEnv=.schnappsEnv) {
  prjCol = projectionColors[[dimCol]]
  if(exists(".schnappsEnv") & is.environment(.schnappsEnv)){
    if (.schnappsEnv$DEBUGSAVE) {
      # prjCol= reactiveValuesToList(projectionColors)
      save(file = "~/SCHNAPPsDebug/panelPlotFactFunc.RData", list = c(ls()))
    }
  }
  # cp=load(file='~/SCHNAPPsDebug/panelPlotFactFunc.RData')
  # save(file = "~/SCHNAPPsDebug/panelPlotFactFunc.RData", list = c(ls()))
  scEx_log = scEx_log[,cellNs]
  projections = projections[cellNs,]
  # browser()
  finalLevels = projections[,factsin] %>% unique()
  featureData <- rowData(scEx_log)
  
  par(mfrow = c(ceiling(length(finalLevels) / 4), 4), mai = c(0., .3, .3, .3))
  rbPal <- colorRampPalette(c("#018c0f", "red"))
  ylim <- c(min(projections[, dimy4]), max(projections[, dimy4]))
  if (dimy4 == "UMI.count") {
    ymax <- 0
    # for (i in 1:length(genesin)) {
    #   geneIdx <- which(toupper(featureData$symbol) == genesin[i])
    ymax <- max(ymax, max(Matrix::colSums(assays(scEx_log)[[1]][, , drop = FALSE])))
    # }
    ylim <- c(0, ymax)
    if (!sameScale) {
      ylim <- NULL
    }
  }
  plotList <- list()
  plotIdx <- 0
  printData = data.frame(level = character(), dimx = numeric(), dimy = numeric(), col = numeric)
  
  
  # projectionColors[[dimCol]]
  # color for each cell based on the expression
  # this will always show the max value per cell
  # this should apply when sameScale = FALSE
  # cuts is a factor with the ranges
  # cuts = cut(
  #   as.numeric(as.character(projections[,dimCol])),
  #   breaks = 10
  # )
  # Col <- rbPal(10)[
  #   as.numeric(
  #     cuts
  #   )
  # ]
  # names(Col) <- cuts
  # colF = Col[unique(names(Col))] %>% sort(decreasing=T)
  
  printData =  bplapply (seq(finalLevels), FUN = function(i) {
    prjIdx = projections[,factsin]==finalLevels[i]
    subsetTSNE <- projections[prjIdx,]
    # plotCol <- Col[prjIdx]
    data.frame(level = finalLevels[i], 
               dimx  = subsetTSNE[, dimx4, drop=T], 
               dimy  = subsetTSNE[, dimy4, drop=T],
               col   =  if(is.factor(subsetTSNE[, dimCol])){
                 prjCol[subsetTSNE[, dimCol, drop=T]]}
               else{
                 subsetTSNE[, dimCol, drop=T]
               }
    )
  }) %>% bind_rows()
  printData = printData[order(printData$col, decreasing = F),]
  # printData$gene = factor(printData$gene, levels = unique(genesin))
  if (is.factor(projections[,dimx4])) {
    retVal <-  ggboxplot(printData, x="dimx", y="dimy",add="jitter", 
                         color  = "col")  
      if(applyPvalue){
      # retVal <- ggplot(printData, aes(dimx, dimy)) + geom_boxplot(show.legend = FALSE)
      # compare_means(dimy ~ dimx,  data = printData, method = "anova")
      cm = compare_means(dimy ~ dimx,  data = printData, method = "t.test")
      maxcomp = min(3,nrow(cm))
      cm = cm[order(cm$p.adj)[1:maxcomp],]
      comparisons = apply(cbind(cm$group1,cm$group2) ,1, list) %>% unlist(recursive = F)
      retVal <-  retVal + stat_compare_means(method = "t.test", comparisons = comparisons)
    }
  } else {
    retVal <- ggplot(printData, aes(dimx, dimy, color = col)) + 
      geom_point(alpha = 0.4) 
  }
  if (is.factor(projections[,dimCol])) {
    retVal = retVal + 
      scale_color_identity(name = dimCol, 
                           labels = names(prjCol),
                           guide = "legend", aesthetics = "colour")
  }
  # if(!is.factor(projections[,dimCol])){
  #   retVal <- retVal + scale_color_continuous(name = dimCol)
  # } else {
    retVal <- retVal + labs(color=dimCol)
  # }
  if (sameScale){
    retVal = retVal + 
      facet_wrap(~level,ncol = nCol) 
  } else {
    retVal = retVal + facet_wrap(~level ,ncol = nCol,scales = "free")
  }
  retVal = retVal + xlab(dimx4) + 
    ylab(dimy4) 
  # dotcolors = printData$col
  # names(dotcolors) = cuts
  
  # retVal = retVal + scale_color_identity(name = dimy4, breaks = colF, labels = names(colF),guide = "legend", aesthetics = "colour") 
  
  # , labels = c("A", "B", "C")
  # +
  #   scale_colour_gradient2()
  
  # retVal =  retVal + scale_color_gradient2(low = lowCol, high = highCol, mid = midCol, midpoint = midFunc(colData))
  
  retVal <-
    ggpubr::annotate_figure(retVal,
                            top = text_grob(sampdesc)
    )
  # retVal
  return(retVal)
} 
#   panelPlotFunc_m = memoise::memoise(panelPlotFunc,cache=do.call(cachem::cache_disk,.schnappsEnv$cacheDir))
