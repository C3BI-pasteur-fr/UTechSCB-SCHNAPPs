#' this is used to run the app without installing it.
#'
#'

if (!exists(".schnappsEnv")) {
  .schnappsEnv <- new.env(parent=emptyenv())
}

localContributionDir = "~/Rstudio/shHubgit/Dummy/"
defaultValueSingleGene = "CD52"
defaultValueMultiGenes = "CD52, S100A4, S100A9, S100A8"
defaultValueRegExGene = "" # tip: '^CD7$|^KIT$; genes with min expression
DEBUG = FALSE
DEBUGSAVE = FALSE

assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueSingleGene", defaultValueSingleGene, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueMultiGenes", defaultValueMultiGenes, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueRegExGene", defaultValueRegExGene, envir = .schnappsEnv)
assign(".SCHNAPPs_DEBUG", DEBUG, envir = .schnappsEnv)
assign(".SCHNAPPs_DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)

devscShinyApp = TRUE
packagePath <<- "inst/app"
source(paste0(packagePath,  "/ui.R"))
source(paste0(packagePath,  "/server.R"))

shinyApp(ui = scShinyUI, server = scShinyServer)

# schnapps(
# defaultValueMultiGenes = "CD3e, CD3d, CD4, IL7R, CD8, CD8a,GNLY, NKG7, CD14, LYZ, MS4A1, CD79A, CD74, CST3, FCER1A",
# defaultValueSingleGene = "IL7R", DEBUG=TRUE
# )


# sctkEx = SCtkExperiment(assays=list(counts=as.matrix(assays(scEx)[['counts']]), 
#                                     logcounts = as.matrix(assays(scEx)[['logcounts']])),
#                         colData = colData(scEx),
#                         rowData = rowData(scEx))
# singleCellTK(inSCE = sctkEx)




