#' .schnappsEnv
#' Environment for package
.schnappsEnv <- new.env(parent = emptyenv())



#' shiny server relevant functions
#'
#' @details Shiny app for the analysis of single cell data
#'
#' @param data RdataFile data file with scEx, scEx_log, projections
#' @param localContributionDir path to the directory(ies) that contain additional functionality
#' @param defaultValueSingleGene single gene name to used as a default value.
#' @param defaultValueMultiGenes comma separated list of gene names to be used as a default value.
#' @param defaultValueRegExGene regular Expression used for gene selection
#' @param DEBUG TRUE/FALSE whether to show debugging information on the console
#' @param DEBUGSAVE TRUE/FALSE where or not save internal data (very time consuming)
#' @param historyPath location (directory) where history directories and data will be stored.
#' historyPath should be used to generate Report
#'
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar sidebarMenu dashboardBody
#' tabItem menuSubItem menuItem
#' @importFrom shinyTree renderTree shinyTree
#' @importFrom DT DTOutput dataTableProxy
#' @importFrom dplyr '%>%'
#' @importFrom shiny shinyApp runApp
#'
#' @export schnappsLite
#'
#' @examples
#' packPath = "schnappsPackage"
#' schnappsLite(data = paste0(packPath, "/data/scExLite.RData"))
schnappsLite <- function(data=RdataFile,
                         localContributionDir = "~/Rstudio/shHubgit/Dummy/",
                     defaultValueSingleGene = "CD52",
                     defaultValueMultiGenes = "CD52, S100A4, S100A9, S100A8",
                     defaultValueRegExGene = "", # tip: '^CD7$|^KIT$; genes with min expression
                     DEBUG = FALSE,
                     DEBUGSAVE = FALSE,
                     historyPath = NULL
                     # ,
                     # historyFile = NULL
                     
) {
  # on.exit({
  #   rm(list = c(".SCHNAPPs_locContributionDir",
  #        ".SCHNAPPs_defaultValueSingleGene",
  #        ".SCHNAPPs_defaultValueMultiGenes",
  #        ".SCHNAPPs_defaultValueRegExGene",
  #        ".SCHNAPPs_DEBUG",
  #        ".SCHNAPPs_DEBUGSAVE"),
  #      envir = .schnappsEnv)
  # })
  assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = .schnappsEnv)
  assign(".SCHNAPPs_defaultValueSingleGene", defaultValueSingleGene, envir = .schnappsEnv)
  assign(".SCHNAPPs_defaultValueMultiGenes", defaultValueMultiGenes, envir = .schnappsEnv)
  assign(".SCHNAPPs_defaultValueRegExGene", defaultValueRegExGene, envir = .schnappsEnv)
  assign(".SCHNAPPs_DEBUG", DEBUG, envir = .schnappsEnv)
  assign(".SCHNAPPs_DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
  assign("DEBUG", DEBUG, envir = .schnappsEnv)
  assign("DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
  assign("historyPath", historyPath, envir = .schnappsEnv)
  # assign("historyFile", historyFile, envir = .schnappsEnv)
  
  # will be set during sourcing, but we need to define them, otherwise there will be a warning
  scShinyUI <- NULL
  scShinyServer <- NULL
  packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/")
  source(paste0(packagePath, "/server.R"), local = TRUE)
  source(paste0(packagePath, "/ui.R"), local = TRUE)
  source(paste0(packagePath, "/server-lite.R"), local = TRUE)
  source(paste0(packagePath, "/ui-lite.R"), local = TRUE)
  
  # load data
  assign(".SCHNAPPs_LiteData", loadLiteData(file = data), envir = .schnappsEnv)
  
  
  if (is.null(.schnappsEnv$".SCHNAPPs_LiteData")) {
    # error loading
    exit()
  }
  
  app <- shinyApp(ui = scShinyUI, server = scShinyServer)
  runApp(app)
}

# library(SCHNAPPs)
#
# schnapps(localContributionDir = "~/Rstudio/__shHubgit/bjContributions",  defaultValueSingleGene = "cd52", defaultValueMultiGenes = "S100A4, CD52, S100A9, S100A8,")



#' Example data for schnapps app
#'
#' A dataset containing the prices and other attributes of almost 54,000
#' diamonds. The variables are as follows:
#'
#' * `scExLite`: singlecellExperiment object
#'
#' @format A data frame with 53940 rows and 10 variables
#' @source <http://www.diamondse.info/>
"scExLite"

#' * `ccol` : color definitions
#' @format a vector with colours used for clusters
"ccol"

#' 
#' * `scol` : color definitions for samples
#' @format a vector of colour values for samples
"scol"

#' 
#' * `pca` : pca results, Eigenvalues and vectors.
#' @format List with "x", "var_pcs", "rotation"
"pca"
