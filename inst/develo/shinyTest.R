# shiny Testing

library(shinytest)

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
source(paste0(packagePath,  "/server.R"))
source(paste0(packagePath,  "/ui.R"))
app <- shinyApp(ui = scShinyUI, server = scShinyServer)


recordTest("/Users/bernd/Rstudio/Schnapps/inst/develo/", loadTimeout = 100000)

testApp("inst/develo", "mytest.R")

