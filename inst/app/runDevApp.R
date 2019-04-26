#' this is used to run the app without installing it.
#'


devscShinyApp = TRUE
packagePath = "inst/app"
source(paste0(packagePath,  "/server.R"))
source(paste0(packagePath,  "/ui.R"))
shinyApp(ui = scShinyUI, server = scShinyServer)
