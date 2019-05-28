require(shiny)
require(reactlog)
require(shinyTree)
require(tibble)
require(plotly)
require(shinythemes)
require(ggplot2)
require(DT)
require(pheatmap)
require(threejs)
require(sm)
require(RColorBrewer)
require(mclust)
require(reshape2)
require(ggplot2)
require(knitr)
require(kableExtra)
require(shinyWidgets)
require(scater)
require(shinyMCE)
require(kohonen)
require(Rsomoclu)
require(gtools)
require(SingleCellExperiment)
require(Matrix)
require(colourpicker)
require(shinytest)
require(scran)
require(callr)
require(debugme)
require(shinyBS)

localContributionDir = ""
defaultValueSingleGene = "CD52"
defaultValueMultiGenes = "CD52, S100A4, S100A9, S100A8"
defaultValueRegExGene = "" # tip: '^CD7$|^KIT$; genes with min expression
DEBUG = FALSE
DEBUGSAVE = FALSE
# I still don't understand how to pass a variable to a shinyApp without going through the globalenv.
assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = globalenv())
assign(".SCHNAPPs_defaultValueSingleGene", defaultValueSingleGene, envir = globalenv())
assign(".SCHNAPPs_defaultValueMultiGenes", defaultValueMultiGenes, envir = globalenv())
assign(".SCHNAPPs_defaultValueRegExGene", defaultValueRegExGene, envir = globalenv())
assign(".SCHNAPPs_DEBUG", DEBUG, envir = globalenv())
assign(".SCHNAPPs_DEBUGSAVE", DEBUGSAVE, envir = globalenv())

packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE)
packagePath <- paste0(packagePath,"/app/")
# cat (file = stderr(), paste("app.R\n"))
packagePath <- "/Users/bernd/Rstudio/Schnapps/inst/app"
source(paste0(packagePath,  "/server.R"))
source(paste0(packagePath,  "/ui.R"))

shinyApp(ui = scShinyUI, server = scShinyServer)
