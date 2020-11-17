#' checkbsTT
#' all variables have to be stored in .schnappsEnv and have to be called with "sbTT_" + variable names (see function definition at the end of this file.)
if ("shinyBS" %in% rownames(installed.packages())) {
  suppressMessages(require(shinyBS))

  # ui.R
  .schnappsEnv$sbTT_summaryStatsSideBar <- bsPopover("summaryStatsSideBar",
    title = "",
    "<h3>Data summary</h3> <ul><li>medium UMI: shows how many genes are  expressed in log2 space of normalized data</li> </ul> "
  )

  .schnappsEnv$sbTT_RDSsave <- bsPopover("RDSsave",
    title = "",
    "<h3>download current cell/gene configuration for reimport to this app</h3>"
  )
  # tabs.R
  # must not contain "\n", or other special characters.
  .schnappsEnv$sbTT_file1 <- bsPopover("file1", title = "", "input file with csv data or an RData/RDs file with SingleCellExperiment objects. If multiple files are given only common genes are used. Projections/colData is filled with NA values. csv files cannot be combined.")
  .schnappsEnv$sbTT_sampleInput <- bsPopover("sampleInput", title = "", "sub sample this amount of cells from each sample. In case logcount data is provided and used (see use scEx from loaded data) only the cells sampled will be used.")
  .schnappsEnv$sbTT_subsampleNum <- bsPopover("subsampleNum", title = "", "max number of cells")
  .schnappsEnv$sbTT_disablescEx_log <- bsPopover("disablescEx_log", title = "", "disable Normalization")
  .schnappsEnv$sbTT_beforeFilterRegEx <- bsPopover("beforeFilterRegEx", title = "", "regular expression to count genes/cell")
  .schnappsEnv$sbTT_selectIds <- bsPopover("selectIds", title = "", "regular expression for selection of genes to be removed")
  .schnappsEnv$sbTT_minGenesGS <- bsPopover("minGenesGS", title = "", "Min # of UMIs over all cells")
  .schnappsEnv$sbTT_genesKeep <- bsPopover("genesKeep", title = "", "genes to keep")
  .schnappsEnv$sbTT_minExpGenes <- bsPopover("minExpGenes", title = "", "<h3>Cells must have one or more</h3> <ul><li>These cells must have at least one of those genes expressed</li> </ul> ")
  .schnappsEnv$sbTT_minNonExpGenes <- bsPopover("minNonExpGenes", title = "", "<h3>Cells must have NOT one or more</h3> <ul><li>These cells must NOT have at least one of those genes expressed</li> </ul> ")
  .schnappsEnv$sbTT_minGenes <- bsPopover("minGenes", title = "", "Min # of UMIs")
  .schnappsEnv$sbTT_maxGenes <- bsPopover("maxGenes", title = "", "Max # of UMIs")
  .schnappsEnv$sbTT_cellSelectionComment <- bsPopover("cellSelectionComment", title = "", "Comment for selection of cells")
  .schnappsEnv$sbTT_cellPatternRM <- bsPopover("cellPatternRM", title = "", "cells to be filtered out by pattern")
  .schnappsEnv$sbTT_cellKeep <- bsPopover("cellKeep", title = "", "cells to keep")
  .schnappsEnv$sbTT_cellKeepOnly <- bsPopover("cellKeepOnly", title = "", "cells to keep (remove others)")
  .schnappsEnv$sbTT_cellsFiltersOut <- bsPopover("cellsFiltersOut", title = "", "Cells to be removed")
  .schnappsEnv$sbTT_pcaRank <- bsPopover("pcaRank", title = "", "of PCA to use for clustering and other calculations")
  .schnappsEnv$sbTT_pcaCenter <- bsPopover("pcaCenter", title = "", "center data")
  .schnappsEnv$sbTT_pcaScale <- bsPopover("pcaScale", title = "", "Whether data should be scaled to unit variance before the PCA is performed.")
  .schnappsEnv$sbTT_genes4PCA <- bsPopover("genes4PCA", title = "", "Genes to be used for PCA")
  .schnappsEnv$sbTT_pcaN <- bsPopover("pcaN", title = "Number of genes to use", "number of genes to use")
  .schnappsEnv$sbTT_pcaCenter <- bsPopover("pcaCenter", title = "", "center data")
  .schnappsEnv$sbTT_pcaScale <- bsPopover("pcaScale", title = "", "scale data")
  .schnappsEnv$sbTT_genes4PCA <- bsPopover("genes4PCA", title = "", "Genes to be used for PCA")
  .schnappsEnv$sbTT_clusterSource <- bsPopover("clusterSource", title = "", "use raw counts or normalized data?")
  .schnappsEnv$sbTT_minClusterSize <- bsPopover("minClusterSize", title = "", "minimum size of each cluster.")
  .schnappsEnv$sbTT_clusterMethod <- bsPopover("clusterMethod", title = "", "clustering method to use")
  .schnappsEnv$sbTT_useRanks <- bsPopover("useRanks", title = "", "whether to use ranks")
  .schnappsEnv$sbTT_geneSelectionClustering <- bsPopover("geneSelectionClustering", title = "", "Genes to be used for clustering")
  .schnappsEnv$sbTT_descriptionOfWork <- bsPopover("descriptionOfWork", title = "", "Please describe your work. This will be included in the report.")
  .schnappsEnv$sbTT_oldPrj <- bsPopover("oldPrj", title = "", "projections to copy + rename")
  .schnappsEnv$sbTT_newPrj <- bsPopover("newPrj", title = "", "new name of Projection")
  .schnappsEnv$sbTT_delPrj <- bsPopover("delPrj", title = "", "projections to delete")

  .schnappsEnv$sbTT_tabsetPCA <- bsPopover("tabsetPCA", title = "PCA", "performs principle component analysis using BiocSingular::runPCA", placement = "left")

  checkbsTT <- function(item) {
    if (exists(".schnappsEnv")) 
      if (is.environment(.schnappsEnv)) {
      if (exists(paste0("sbTT_", item), envir = .schnappsEnv)) {
        .schnappsEnv[[paste0("sbTT_", item)]]
      }
    }
  }
} else {
  checkbsTT <- function(...) {
    return("")
  }
}
# checkbsTT(item="pcaN")
#####
