#' checkbsTT
#' all variables have to be stored in .schnappsEnv and have to be called with "sbTT_" + variable names (see function definition at the end of this file.)
#' no "." in variable name
#' 
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
  .schnappsEnv$sbTT_whichscLog <- bsPopover("whichscLog", title = "comment on Normalization option", "A SingleCellExperiment object can contain already calculated normalized counts in an assay called logcounts. This can be used (2nd option). It is also possible to not calcualte any (1st option), or use the SCHNAPPs to populate this. There are no other assays being used by default in SCHNAPPs.")
  .schnappsEnv$sbTT_subsampleNum <- bsPopover("subsampleNum", title = "", "max number of cells per sample. Use this to either reduce memory footprint or get similar number of cells per sample.")
  .schnappsEnv$sbTT_disablescEx_log <- bsPopover("disablescEx_log", title = "", "disable Normalization")
  .schnappsEnv$sbTT_beforeFilterRegEx <- bsPopover("beforeFilterRegEx", title = "", "regular expression to count genes/cell")
  .schnappsEnv$sbTT_selectIds <- bsPopover("selectIds", title = "", "regular expression for selection of genes to be removed")
  .schnappsEnv$sbTT_minGenesGS <- bsPopover("minGenesGS", title = "", "Min # of UMIs over all cells")
  .schnappsEnv$sbTT_genesKeep <- bsPopover("genesKeep", title = "", "genes to keep")
  .schnappsEnv$sbTT_minExpGenes <- bsPopover("minExpGenes", title = "", "<h3>Cells must have one or more</h3> <ul><li>These cells must have at least one of those genes expressed</li> </ul> ")
  .schnappsEnv$sbTT_minNonExpGenes <- bsPopover("minNonExpGenes", title = "", "<h3>Cells must have NOT one or more</h3> <ul><li>These cells must NOT have at least one of those genes expressed</li> </ul> ")
  .schnappsEnv$sbTT_minGenes <- bsPopover("minGenes", title = "", "Min # of UMIs per cell")
  .schnappsEnv$sbTT_maxGenes <- bsPopover("maxGenes", title = "", "Max # of UMIs per cell")
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

  
  .schnappsEnv[['sbTT_coE_selected-geneIds']] <- bsPopover('coE_selected-geneIds', title = "", "sum of expression values if normalization is done, otherwise sum of UMI count.")
  .schnappsEnv[['sbTT_DE_panelplotPvalue']] <- bsPopover('DE_panelplotPvalue-geneIds', title = "", "Calculate t-test if x: factorial, y='UMI-count'. Max 5 are shown.")
  # should be read from contribution directory
  .schnappsEnv$sbTT_temporaCluster <-  bsPopover("temporaCluster", title = "", "name of cluster centroids.")
  .schnappsEnv$sbTT_temporaFactor <-  bsPopover("temporaFactor", title = "", "Factor used as time variable")
  .schnappsEnv$sbTT_temporaLevels <-  bsPopover("temporaLevels", title = "", "Ordered list of levels for time variable")
  .schnappsEnv$sbTT_temporaGMTFile <-  bsPopover("temporaGMTFile", title = "", "File with annotation, see http://download.baderlab.org/EM_Genesets/current_release/")
  .schnappsEnv$sbTT_temporaMinSz <-  bsPopover("temporaMinSz", title = "", "Minimum size of the genesets used in enrichment estimation, set to 5 genes by default.")
  .schnappsEnv$sbTT_temporaMaxSz <-  bsPopover("temporaMaxSz", title = "", "Maximum size of the genesets used in enrichment estimation, set to 200 genes by default.")
  .schnappsEnv$sbTT_temporaNPCs  <-  bsPopover("temporaNPCs", title = "", "Number of principal components to be used in building the network.")
  .schnappsEnv$sbTT_temporaDiff_thresh  <-  bsPopover("temporaDiff_thresh", title = "", "Percent of permissible difference between the temporal scores of two clusters to determine the direction of their connections. The temporal scores are calculated based on based on the clusters composition of cells from each timepoint. The directions of edges connecting pairs of clusters will only be determined for cluster pairs with difference in their time scores higher than the threshold. Other edges will remain undirected. Default at 0.01")
  .schnappsEnv$sbTT_temporaPval_thresh  <-  bsPopover("temporaPval_thresh", title = "", "P-value threshold to determine the significance of pathway enrichment over time. Set rel. high because filtering is done later. Default to 0.5.")
  
  
  .schnappsEnv$sbTT_gQC_tsneTheta <-  bsPopover("gQC_tsneTheta", title = "", "Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)")
  .schnappsEnv$sbTT_gQC_tsnePerplexity  <- bsPopover("gQC_tsnePerplexity", title = "", "text")
  .schnappsEnv$sbTT_gQC_um_n_neighbors  <- bsPopover("gQC_um_n_neighbors", title = "", "text")
  .schnappsEnv$sbTT_gQC_um_spread <- bsPopover("gQC_um_spread", title = "", "text")
  .schnappsEnv$sbTT_gQC_um_local_connectivity <- bsPopover("gQC_um_local_connectivity", title = "", "text")
  .schnappsEnv$sbTT_gQC_um_init <- bsPopover("gQC_um_init", title = "", "text")
  .schnappsEnv$sbTT_gQC_um_negative_sample_rate <- bsPopover("gQC_um_negative_sample_rate", title = "", "text")
  .schnappsEnv$sbTT_gQC_um_min_dist <- bsPopover("gQC_um_min_dist", title = "", "text")
  .schnappsEnv$sbTT_gQC_um_metric <- bsPopover("gQC_um_metric", title = "", "text")
  .schnappsEnv$sbTT_gQC_um_set_op_mix_ratio <- bsPopover("gQC_um_set_op_mix_ratio", title = "", "text")
  .schnappsEnv$sbTT_gQC_um_bandwidth  <- bsPopover("gQC_um_bandwidth", title = "", "text")
  
  
  
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
