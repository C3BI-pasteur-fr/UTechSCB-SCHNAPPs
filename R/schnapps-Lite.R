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
#' @param maxCells maximal number of cells to use. Used to limit memory usage for server
#' @param historyPath location (directory) where history directories and data will be stored.
#' @param defaultValues list of default values to use for inputs
#' @param AllowClustering whether to include functionality to cluster cells.
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
#' # TODO fix:
#' # schnappsLite(data = paste0(packPath, "/data/scExLite.RData"))
schnappsLite <- function(data="RdataFile",
                         localContributionDir = "~/Rstudio/shHubgit/Dummy/",
                         defaultValueSingleGene = "CD52",
                         defaultValueMultiGenes = "CD52, S100A4, S100A9, S100A8",
                         defaultValueRegExGene = "", # tip: '^CD7$|^KIT$; genes with min expression
                         DEBUG = FALSE,
                         DEBUGSAVE = FALSE,
                         historyPath = NULL,
                         maxCells = 3000,
                         defaultValues = list(),
                         AllowClustering = FALSE
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
  assign("defaultValues", defaultValues, envir = .schnappsEnv)
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
  
  
  nCells = length(colnames(.schnappsEnv$.SCHNAPPs_LiteData$scEx))
  if (nCells > maxCells){
    cellIdx = unique(sort(sample(nCells, maxCells)))
    cells2keep = colnames(.schnappsEnv$.SCHNAPPs_LiteData$scEx)[cellIdx]
    for (na in names(.schnappsEnv$.SCHNAPPs_LiteData)){
      if (na == "pca") {
        .schnappsEnv$.SCHNAPPs_LiteData[[na]]$x = .schnappsEnv$.SCHNAPPs_LiteData[[na]]$x[cellIdx,]
        next()
      }
      if (na %in%  c("sampleCol", "clusterCol")) next()
      if (na == "tsne") {
        .schnappsEnv$.SCHNAPPs_LiteData[[na]] = .schnappsEnv$.SCHNAPPs_LiteData[[na]][cellIdx,]
        next()
      }
      if ( all(cells2keep %in% colnames(.schnappsEnv$.SCHNAPPs_LiteData[[na]]))){
        .schnappsEnv$.SCHNAPPs_LiteData[[na]] = .schnappsEnv$.SCHNAPPs_LiteData[[na]][,cells2keep]
      } else {
        if(all(cells2keep %in% rownames(.schnappsEnv$.SCHNAPPs_LiteData[[na]]))){
          .schnappsEnv$.SCHNAPPs_LiteData[[na]] = .schnappsEnv$.SCHNAPPs_LiteData[[na]][cells2keep,]
        } else {
          if(is.null(rownames(.schnappsEnv$.SCHNAPPs_LiteData[[na]]))) {
            .schnappsEnv$.SCHNAPPs_LiteData[[na]] = .schnappsEnv$.SCHNAPPs_LiteData[[na]][cellIdx]
          }else{
            cat(file = stderr(), paste("couldn't ", na, "\n"))
          }
        }
      }
    }
    
  }
  
  if (is.null(.schnappsEnv$".SCHNAPPs_LiteData")) {
    # error loading
    return(NULL)
  }
  # save(file = "global.RData", list = c(".schnappsEnv"))
  app <- shinyApp(ui = scShinyUI, server = scShinyServer)
  # runApp(app)
  
  
}

# library(SCHNAPPs)
#
# schnapps(localContributionDir = "~/Rstudio/__shHubgit/bjContributions",  defaultValueSingleGene = "cd52", defaultValueMultiGenes = "S100A4, CD52, S100A9, S100A8,")


