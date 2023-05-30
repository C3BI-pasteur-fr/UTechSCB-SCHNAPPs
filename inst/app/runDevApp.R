#' this is used to run the app without installing it.
#'
#'

# git tag 1.6.134
# git --tags

library(doParallel)
library("future.callr")

registerDoParallel(cores=3)
suppressMessages(require(shinyjqui))

library(reactlog)
# if (!exists(".schnappsEnv")) {
.schnappsEnv <- new.env(parent=emptyenv())
# }
#
library(future)
# plan("multisession", workers = 3)
plan(callr, workers = 3
     )

library("BiocParallel")
register(MulticoreParam(9))
# register(SerialParam())

localContributionDir = "~/Rstudio/SCHNAPPsContributions/"
localContributionDir = ""
defaultValueSingleGene = "IL7R" # CD52
defaultValueMultiGenes = "IL7R, CCR7 IL7R, S100A4, CD8A, CD8A ,GNLY, NKG7,PPBP, FCER1A, MS4A7,CD14, LYZ,FCGR3A, MS4A7,MS4A"
defaultValueMultiGenes = "LINC00115, NOC2L, HES4, ISG15, TNFRSF18, CD52, SH3BGRL3"
# defaultValueMultiGenes = "prf1, Gzmb, IFNG, PDCD1, HAVCR2, LAG3, TSC22D3,ZFP36L2"
defaultValueRegExGene = "" # tip: '^CD7$|^KIT$; genes with min expression
DEBUG = T
DEBUGSAVE = F
historyPath = "/Volumes/LaCie2022/RStudio_history/celia/hist_2023-Jan-18.09.02"
historyPath = "/Volumes/LaCie2022/RStudio_history/celia/hist_2023-Feb-22.17.42/"
historyPath = "/Volumes/LaCie2022/RStudio_history/MPI/hist_2023-May-08.10.02/"
historyPath = "/Volumes/LaCie2022/RStudio_history/MPI/hist_grp3/"
historyPath = "/Volumes/LaCie2022/RStudio_history/MPI/hist_3 samples/"
# historyPath = "/Volumes/LaCie2022/RStudio_history/MPI/"
historyPath = "demoHistory/hist_2023-May-11.12.43/"
# historyPath = "demoHistory/hist_2023-May-11.12.43/"
historyPath = "/Volumes/CBUtechsZeus/bernd/celia/hist_2023-May-15.09.49/"
# historyPath = "demoHistory/MPI"
# historyPath = "/Volumes/LaCie2022/RStudio_history/marielle/hist_2022-Dec-15.18.15/"
# historyPath = NULL
assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueSingleGene", defaultValueSingleGene, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueMultiGenes", defaultValueMultiGenes, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueRegExGene", defaultValueRegExGene, envir = .schnappsEnv)
assign(".SCHNAPPs_DEBUG", DEBUG, envir = .schnappsEnv)
assign(".SCHNAPPs_DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
assign("localContributionDir", localContributionDir, envir = .schnappsEnv)
assign("defaultValueSingleGene", defaultValueSingleGene, envir = .schnappsEnv)
assign("defaultValueMultiGenes", defaultValueMultiGenes, envir = .schnappsEnv)
assign("defaultValueRegExGene", defaultValueRegExGene, envir = .schnappsEnv)
assign("DEBUG", DEBUG, envir = .schnappsEnv)
assign("DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
assign("DEBUGSAVE_cellSelectionModule", TRUE, envir = .schnappsEnv)
assign("historyPath", historyPath, envir = .schnappsEnv)
assign("allowFunctionality--shiny-tab-clusterParameters", T, envir = .schnappsEnv)
assign("enableTrajectories", T, envir = .schnappsEnv)

# ls(.schnappsEnv)


# # Scran parameters
defaultValues = list()


# Seurat parameters
defaultValues = list()

# defaultValues[["selectIds"]] = ""
# defaultValues[["pcaN"]] = 200
# defaultValues[["pcaScale"]] = TRUE
# defaultValues[["pcaRank"]]  = 15
# defaultValues[["sampleInput"]] =FALSE
# defaultValues[["hvgSelection"]] = "getTopHVGs"
# # defaultValues[["alluiv1"]] = "seurartCluster"
# defaultValues[["alluiv2"]] = "dbCluster"
# # defaultValues[["tabsetCluster"]] = "seurat_Clustering"
# defaultValues[["minGenesGS"]] = 1
# defaultValues[["minGenes"]] = 1
# defaultValues[["maxGenes"]] = 15000
# defaultValues[["seurClustDims"]] = 15
# defaultValues[["seurClustk.param"]] = 10
# defaultValues[["cellPatternRM"]] = ""
# defaultValues[["gQC_binSize"]] = 200
# defaultValues[["selectIds"]] = "^MT-|^RP|^MRP|MALAT1|B2M|EEF1A1"
# defaultValues[["selectIds"]] = ""
# defaultValues[["whichscLog"]] = "calcLog"
# defaultValues[["whichscLog"]] = "disablescEx_log"
# defaultValues[["normalizationRadioButton"]] = "rawNormalization"

# defaultValues[["gQC_um_n_neighbors"]] = 20
# defaultValues[["gQC_um_spread"]] = 6
# defaultValues[["gQC_um_local_connectivity"]] = 2
# defaultValues[["useSeuratPCA"]] = TRUE

# defaultValues[["DE_panelplotids"]] = c("CD8A", "CD4", "CD8B", "FCER1G", "CCR7", "GZMK", "FoxP3", "GZMK", "GZMB", "CCR7", "LEF1", "CCL5", "VIM", "CCL5", "TCF7", "NKG7", "LGALs1", "NKG7", "SELL", "CST7", "ANXA2", "CST7", "IL7R", "HLA-DRB1", "KLRG1", "CD27", "CTLA4A")
assign("defaultValues", defaultValues, envir = .schnappsEnv)

devscShinyApp = TRUE
packagePath <<- "inst/app/"
source(paste0(packagePath,  "/ui.R"))
source(paste0(packagePath,  "/server.R"))

options(shinyjqui.debug = TRUE)
options(shinyjquiui.debug = TRUE)

app <- shinyApp(ui = scShinyUI, server = scShinyServer, enableBookmarking = "server")
options(shiny.reactlog=TRUE)

#
#
options(keep.source=TRUE)
# p <- profvis::profvis({
# })
# htmlwidgets::saveWidget(p, '~/profvis1.html')

# schnapps(
# defaultValueMultiGenes = "IL7R, CCR7,CD14, LYZ ,IL7R, S100A4,MS4A1 ,CD8A,FCGR3A, MS4A7 ,GNLY, NKG7,FCER1A, CST3,PPBP",
# defaultValueSingleGene = "MS4A1", DEBUG=TRUE
# )
#
#
# sctkEx = SCtkExperiment(assays=list(counts=as.matrix(assays(scEx)[['counts']]),
#                                     logcounts = as.matrix(assays(scEx)[['logcounts']])),
#                         colData = colData(scEx),
#                         rowData = rowData(scEx))
# singleCellTK(inSCE = sctkEx)
#
#
#
# library(SCHNAPPs)
# schnapps(DEBUG = T, historyPath = "/Volumes/Oct2020/RStudio/history/celia/")
#
#
if (DEBUG) base::cat(file = stderr(), paste("\n\n\n", packagePath,"\n\n\n"))
# zz <- file("schnapps.output.txt", open = "wt")
# sink(zz, type = "message")

options("future.globals.maxSize")
options(future.globals.maxSize= 2024^3)
options(shinyjqui.debug = TRUE)
shiny::addResourcePath(
  prefix = "www",
  directoryPath = "./inst/www/"
)

chromeBrowser = function(url){
  system(paste("open -a 'Google Chrome' ", url))
}
edgeBrowser = function(url){
  system(paste("open -a '/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge' ", url))
}
# chromeBrowser = getOption("shiny.launch.browser", interactive())
source("R/DotPlotwithModuleScore.R")
# if(!is.null(historyPath)){
#   shinyOptions(cache = cachem::cache_disk(paste0(historyPath, "/app_cache/cache/")))
# }
# options(error = dump.frames)
# options(shiny.error=NULL)
# options(shiny.error=recover)
# options(shiny.jquery.version=3)
# options(shiny.trace=TRUE)
# options(shiny.devmode=F)
runApp(app, port=3838, launch.browser = chromeBrowser)
# 
# last.dump[length(last.dump)]
# debugger(last.dump)

