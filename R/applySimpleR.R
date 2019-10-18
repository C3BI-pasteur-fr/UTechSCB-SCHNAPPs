require(SingleR)


#' attach data on cell types using SingleR
#'
#' @details using the SingleR functionaly to attach that information to a scEx (SingleCellExperiment object)
#'
#' 
#' @param scEx singleCellExperiment object
#' @param hpca.se e.g. HumanPrimaryCellAtlasData()
#' @param colName name for the new col
#' 
#' @importFrom SingleR SingleR pruneScores
#'
#' @export applySimpleR
#'
applySimpleR <- function(scEx, hpca.se, colName) {
  if (!"logcounts" %in% names(SummarizedExperiment::assays(scEx))){
    scEx <- scater::normalize(scEx)
  }
  
  common <- dplyr::intersect(rownames(scEx), rownames(hpca.se))
  hpca.se <- hpca.se[common,]
  scExC <- scEx[common,]
  pred.hpca <- SingleR::SingleR(test = scExC, ref = hpca.se, labels = hpca.se$label.main)
  to.remove <- SingleR::pruneScores(pred.hpca)
  new.pruned <- pred.hpca$labels
  new.pruned[pruneScores(pred.hpca, nmads=5)] <- NA
  cdata = SingleCellExperiment::colData(scEx)
  cdata[, colName] <- NA
  cdata[colnames(scExC), colName]  <- new.pruned
  SingleCellExperiment::colData(scEx) <- cdata
  return(scEx)
}
