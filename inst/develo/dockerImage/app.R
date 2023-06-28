options(shiny.sanitize.errors = FALSE)
library(SCHNAPPs)
#
#schnappsLite(data = "epdc.rn-sham2.v2.lite.RData", DEBUG = T, historyPath = "history")
#schnappsLite(data = "rnAllEPDC.lite.RData", DEBUG = T, historyPath = "history")
defaultValueMultiGenes = "wt1,pdgfrb,agtr1a,col3a1,col1a1,postn,tbx18,scx,npr1,vegfa,npr2,notch2,tagln,acta2,ctgf, Rack1, Bmp4,  Eef2, col3a1,col1a1,tgfb3,postn"
defaultValueSingleGene = "wt1"


defaultValues = list()
defaultValues[["coEtgMinExpr"]] = 100
packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/")

base::cat(file = stderr(), paste("\n\n\n", packagePath,"\n\n\n"))
source(paste0(packagePath,  "/ui.R"))
source(paste0(packagePath,  "/server.R"))
source("R/DotPlotwithModuleScore.R")


options("future.globals.maxSize")
options(future.globals.maxSize= 2024^3)
options(shinyjqui.debug = TRUE)
shiny::addResourcePath(
  prefix = "www",
  directoryPath = "./inst/www/"
)

app <- shinyApp(ui = scShinyUI, server = scShinyServer, enableBookmarking = "server")


runApp(app, port=3838)

#schnapps(defaultValues = defaultValues, DEBUG = T, historyPath = "history", port = 3838)

