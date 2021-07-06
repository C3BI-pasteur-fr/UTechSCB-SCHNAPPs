#' .schnappsEnv
#' Environment for package
.schnappsEnv <- new.env(parent = emptyenv())



#' shiny server relevant functions
#'
#' @details Shiny app for the analysis of single cell data
#'
#'
#' @param localContributionDir path to the directory(ies) that contain additional functionality
#' @param defaultValueSingleGene single gene name to used as a default value.
#' @param defaultValueMultiGenes comma separated list of gene names to be used as a default value.
#' @param defaultValueRegExGene regular Expression used for gene selection
#' @param DEBUG TRUE/FALSE whether to show debugging information on the console
#' @param DEBUGSAVE TRUE/FALSE where or not save internal data (very time consuming)
#' @param historyPath location (directory) where history directories and data will be stored.
#' @param defaultValues list of default values to use for inputs
#'
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar sidebarMenu dashboardBody
#' tabItem menuSubItem menuItem
#' @importFrom shinyTree renderTree shinyTree
#' @importFrom DT DTOutput dataTableProxy
#' @importFrom dplyr '%>%'
#' @importFrom shiny shinyApp runApp
#'
#' @export schnapps
#'
#' @examples
#' # create example data
#' data("scEx", package = "SCHNAPPs")
#' save(file = "scEx.Rdata", list = "scEx")
#' # use "scEx.Rdata" with load data functionality within the shiny app
schnapps <- function(localContributionDir = "~/Rstudio/shHubgit/Dummy/",
                     defaultValueSingleGene = "CD52",
                     defaultValueMultiGenes = "CD52, S100A4, S100A9, S100A8",
                     defaultValueRegExGene = "", # tip: '^CD7$|^KIT$; genes with min expression
                     DEBUG = FALSE,
                     DEBUGSAVE = FALSE,
                     historyPath = NULL,
                     defaultValues = list()
                     # ,
                     # historyFile = NULL

) {
 # make sure libraries are loaded in right order
  library(tidySingleCellExperiment)
  library(Rtsne)
  library(pryr)
  library(dendsort)
  library(MASS)
  library(hms)
  library(uwot)
  library(evaluate)
  library(glue)
  library(profvis)
  library(dplyr)
  library(plyr)
  library(stringr)
  library(crayon)
  library(irlba)
  library(SeuratObject)
  library(Seurat)
  library(shinyjqui)
  library(shinycssloaders)
  library(shinyjs)
  library(edgeR)
  library(limma)
  library(shinydashboardPlus)
  library(shinydashboard)
  library(shinyBS)
  library(InteractiveComplexHeatmap)
  library(ComplexHeatmap)
  library(tidyr)
  library(psychTools)
  library(digest)
  library(magrittr)
  library(future)
  library(reactlog)
  library(kableExtra)
  library(gtools)
  library(debugme)
  library(BiocSingular)
  library(ggalluvial)
  library(scran)
  library(colourpicker)
  library(Matrix)
  library(kohonen)
  library(scater)
  library(scuttle)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Biobase)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(S4Vectors)
  library(BiocGenerics)
  library(MatrixGenerics)
  library(matrixStats)
  library(shinyWidgets)
  library(knitr)
  library(reshape2)
  library(RColorBrewer)
  library(pheatmap)
  library(DT)
  library(plotly)
  library(ggplot2)
  library(tibble)
  library(shinyTree)
  library(shiny)
  library(SCHNAPPs)
  
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
  devscShinyApp <<- FALSE
  devscShinyApp <- FALSE
  source(paste0(packagePath, "/server.R"), local = TRUE)
  source(paste0(packagePath, "/ui.R"), local = TRUE)
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
#' * `scEx`: singlecellExperiment object
#'
#' @format A data frame with 53940 rows and 10 variables
#' @source <http://www.diamondse.info/>
"scEx"
