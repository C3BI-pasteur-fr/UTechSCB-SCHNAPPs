#' this is used to run the app without installing it.
#'
#'

# if (!exists(".schnappsEnv")) {
  .schnappsEnv <- new.env(parent=emptyenv())
# }

  localContributionDir = "~/Rstudio/scShinyHub-github/bjContributions/"
  # localContributionDir = ""
  defaultValueSingleGene = "LYZ"
defaultValueMultiGenes = "Gzmb, prf1, ifng, pdcd1, havcr2, lag3, tsc22d3, zfp36l2"
defaultValueRegExGene = "" # tip: '^CD7$|^KIT$; genes with min expression
DEBUG = TRUE
DEBUGSAVE = F
historyFile = "~/Rstudio/Schnapps/history.File.pdf"

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
assign("historyFile", historyFile, envir = .schnappsEnv)
ls(.schnappsEnv)

devscShinyApp = TRUE
packagePath <<- "inst/app"
source(paste0(packagePath,  "/ui.R"))
source(paste0(packagePath,  "/server.R"))

app <- shinyApp(ui = scShinyUI, server = scShinyServer)

runApp(app)

# schnapps(
# defaultValueMultiGenes = "IL7R, CCR7,CD14, LYZ ,IL7R, S100A4,MS4A1 ,CD8A,FCGR3A, MS4A7 ,GNLY, NKG7,FCER1A, CST3,PPBP",
# defaultValueSingleGene = "MS4A1", DEBUG=TRUE
# )


# sctkEx = SCtkExperiment(assays=list(counts=as.matrix(assays(scEx)[['counts']]), 
#                                     logcounts = as.matrix(assays(scEx)[['logcounts']])),
#                         colData = colData(scEx),
#                         rowData = rowData(scEx))
# singleCellTK(inSCE = sctkEx)




