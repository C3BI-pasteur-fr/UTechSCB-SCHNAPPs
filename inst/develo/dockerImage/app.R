options(shiny.sanitize.errors = FALSE)
library(SCHNAPPs)
library(dplyr)
library(shiny)
#
#schnappsLite(data = "epdc.rn-sham2.v2.lite.RData", DEBUG = T, historyPath = "history")
#schnappsLite(data = "rnAllEPDC.lite.RData", DEBUG = T, historyPath = "history")
defaultValueMultiGenes = "wt1,pdgfrb,agtr1a,col3a1,col1a1,postn,tbx18,scx,npr1,vegfa,npr2,notch2,tagln,acta2,ctgf, Rack1, Bmp4,  Eef2, col3a1,col1a1,tgfb3,postn"
defaultValueSingleGene = "wt1"
scShinyUI <- NULL
scShinyServer <- NULL
packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/")
.schnappsEnv <<- new.env(parent=emptyenv())

localContributionDir = "~/Rstudio/shHubgit/Dummy/"
defaultValueSingleGene = "CD3g"
defaultValueMultiGenes = "cd3g, cd4, cd8b, ms4a1, TCF4, LILRA2, LYZ, cd79a, bcl11b, IL32, hbb, nkg7,MNDA"
defaultValueRegExGene = ""
DEBUG = TRUE
DEBUGSAVE = FALSE
historyPath = NULL
defaultValues = list()
port = NULL
launch.browser = getOption("shiny.launch.browser", interactive())
  assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = .schnappsEnv)
  assign(".SCHNAPPs_defaultValueSingleGene", defaultValueSingleGene, envir = .schnappsEnv)
  assign(".SCHNAPPs_defaultValueMultiGenes", defaultValueMultiGenes, envir = .schnappsEnv)
  assign(".SCHNAPPs_defaultValueRegExGene", defaultValueRegExGene, envir = .schnappsEnv)
  assign(".SCHNAPPs_DEBUG", DEBUG, envir = .schnappsEnv)
  assign(".SCHNAPPs_DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
  assign("DEBUG", DEBUG, envir = .schnappsEnv)
  assign("DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
  assign("historyPath", historyPath, envir = .schnappsEnv)
  assign("defaultValues", defaultValues, envir = .schnappsEnv)
  # assign("historyFile", historyFile, envir = .schnappsEnv)
  
 
  # will be set during sourcing, but we need to define them, otherwise there will be a warning
  devscShinyApp <<- FALSE
  devscShinyApp <- FALSE


defaultValues = list()
defaultValues[["coEtgMinExpr"]] = 100
packagePath <<- packagePath
localContributionDir="."
assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = .schnappsEnv)

base::cat(file = stderr(), paste("\n\n\n", packagePath,"\n\n\n"))
source(paste0(packagePath,  "/ui.R"),local = T)
cat(file = stderr(), "here as\n")

source(paste0(packagePath,  "/server.R"),local = T)
# source("R/DotPlotwithModuleScore.R")
cat(file = stderr(), "here 1\n")

options("future.globals.maxSize")
options(future.globals.maxSize= 2024^3)
# options(shinyjqui.debug = TRUE)
cat(file = stderr(), "here 2\n")
shiny::addResourcePath(
  prefix = "www",
  directoryPath = paste0(packagePath, "/../www/")
)
cat(file = stderr(), "here 3\n")

shinyApp(ui = scShinyUI, server = scShinyServer, enableBookmarking = "server")

cat(file = stderr(), "here 4\n")

# runApp(app)
cat(file = stderr(), "here 5\n")

#schnapps(defaultValues = defaultValues, DEBUG = T, historyPath = "history", port = 3838)

