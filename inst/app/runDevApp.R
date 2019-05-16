#' this is used to run the app without installing it.
#'

localContributionDir = "~/Rstudio/shHubgit/Dummy/"
defaultValueSingleGene = "CD52"
defaultValueMultiGenes = "CD52, S100A4, S100A9, S100A8"
defaultValueRegExGene = "" # tip: '^CD7$|^KIT$; genes with min expression
DEBUG = TRUE
DEBUGSAVE = FALSE
assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = globalenv())
assign(".SCHNAPPs_defaultValueSingleGene", defaultValueSingleGene, envir = globalenv())
assign(".SCHNAPPs_defaultValueMultiGenes", defaultValueMultiGenes, envir = globalenv())
assign(".SCHNAPPs_defaultValueRegExGene", defaultValueRegExGene, envir = globalenv())
assign(".SCHNAPPs_DEBUG", DEBUG, envir = globalenv())
assign(".SCHNAPPs_DEBUGSAVE", DEBUGSAVE, envir = globalenv())

devscShinyApp = TRUE
packagePath = "inst/app"
source(paste0(packagePath,  "/server.R"))
source(paste0(packagePath,  "/ui.R"))
library(shiny)
library(profvis)
  shinyApp(ui = scShinyUI, server = scShinyServer)

# sctkEx = SCtkExperiment(assays=list(counts=as.matrix(assays(scEx)[['counts']]), 
#                                     logcounts = as.matrix(assays(scEx)[['logcounts']])),
#                         colData = colData(scEx),
#                         rowData = rowData(scEx))
# singleCellTK(inSCE = sctkEx)
