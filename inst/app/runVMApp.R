#' this is used to run the app without installing it.
#'
#'
library(reactlog)
# if (!exists(".schnappsEnv")) {
.schnappsEnv <- new.env(parent=emptyenv())
# }
# 
library(future)
library(future.callr)
# devtools::install_github("C3BI-pasteur-fr/UTechSCB-SCHNAPPs", dependencies = TRUE)


if(!exists("WORKERS")) WORKERS = parallel::detectCores()

plan('callr', workers = 2)
# plan("multicore", workers = WORKERS)
# plan(sequential)
library(doParallel)
registerDoParallel(cores=WORKERS)

library("BiocParallel")
register(safeBPParam(WORKERS))
# register(SerialParam())

localContributionDir = normalizePath("/home/schnapps/SCHNAPPsContributions/")
#localContributionDir = NULL
defaultValueSingleGene = "CD3g" # CD52
defaultValueMultiGenes = "cd3g, cd4, cd8b, ms4a1, TCF4, LILRA2, LYZ, cd79a, bcl11b, IL32, hbb, nkg7,MNDA" # itgae, cd69, itga1" # CD52, S100A9, S100A4
# defaultValueMultiGenes = "prf1, Gzmb, IFNG, PDCD1, HAVCR2, LAG3, TSC22D3,ZFP36L2"
defaultValueRegExGene = "" # tip: '^CD7$|^KIT$; genes with min expression
DEBUG = T
DEBUGSAVE = F
historyPath = normalizePath("/home/schnapps/history")
# historyPath = "/Volumes/LaCie2022//RStudio_history/"
#historyPath = NULL

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
assign("historyPath", historyPath, envir = .schnappsEnv)
ls(.schnappsEnv)


# # Scran parameters
defaultValues = list()


# Seurat parameters
# defaultValues = list()

# # defaultValues[["selectIds"]] = ""
# defaultValues[["pcaN"]] = 500
# defaultValues[["pcaScale"]] = TRUE
# defaultValues[["sampleInput"]] =FALSE
# defaultValues[["hvgSelection"]] = "vst"
# # defaultValues[["alluiv1"]] = "seurartCluster"
# defaultValues[["alluiv2"]] = "dbCluster"
# # defaultValues[["tabsetCluster"]] = "seurat_Clustering"
# defaultValues[["minGenesGS"]] = 100
# defaultValues[["minGenes"]] = 1
# defaultValues[["maxGenes"]] = 50000
# defaultValues[["seurClustDims"]] = 15
# defaultValues[["seurClustk.param"]] = 15
# defaultValues[["cellPatternRM"]] = ""
# defaultValues[["gQC_binSize"]] = 200
# defaultValues[["selectIds"]] = "^MT-|^RP|^MRP" #"^MT-|^RP|^MRP|MALAT1|B2M|EEF1A1"
# defaultValues[["selectIds"]] = ""
# defaultValues[["whichscLog"]] = "calcLog"
# defaultValues[["whichscLog"]] = "disablescEx_log"

# defaultValues[["gQC_um_n_neighbors"]] = 20 
# defaultValues[["gQC_um_spread"]] = 6 
# defaultValues[["gQC_um_local_connectivity"]] = 2
# defaultValues[["useSeuratPCA"]] = TRUE


# commented out for paper
# assign("defaultValues", defaultValues, envir = .schnappsEnv)

devscShinyApp = FALSE
packagePath <<- normalizePath("/home/schnapps/rstudio/schnappsGit/inst/app")
packagePath <<- paste0("inst",.Platform$file.sep, "app/")

source(normalizePath(paste0(packagePath,  "/ui.R")))
source(normalizePath(paste0(packagePath,  "/server.R")))

app <- shinyApp(ui = scShinyUI, server = scShinyServer, enableBookmarking = "server")
options(shiny.reactlog=FALSE)
shiny::addResourcePath(
  prefix = "www",
  directoryPath = "./inst/www/"
)
# options(keep.source=TRUE)
# p <- profvis::profvis({
runApp(app, host = "0.0.0.0", port = 6149, launch.browser = FALSE)
# })
# htmlwidgets::saveWidget(p, '~/profvis1.html')

# 
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
