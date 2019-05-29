#' shiny server relevant functions
#'
#' @details Shiny app for the analysis of single cell data
#'
#' @export schnapps
#' 
#' @param localContributionDir path to the directory(ies) that contain additional functionality
#' @param defaultValueSingleGene single gene name to used as a default value.
#' @param defaultValueMultiGenes comma separated list of gene names to be used as a default value.
#' @param defaultValueRegExGene regular Expression used for gene selection 
#' @param DEBUG TRUE/FALSE whether to show debugging information on the console
#' @param DEBUGSAVE TRUE/FALSE where or not save internal data (very time consuming)
#' 
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar sidebarMenu dashboardBody
#' tabItem menuSubItem menuItem
#' @importFrom  shiny actionButton setBookmarkExclude observeEvent updateRadioButtons shinyUI
#' checkboxInput htmlOutput downloadButton  tags withProgress
#' renderPrint renderUI exportTestValues renderText isolate downloadHandler fluidRow fileInput textInput
#' column numericInput textOutput icon uiOutput
#' radioButtons verbatimTextOutput wellPanel
#' @importFrom htmltools HTML br
#' @importFrom shinyTree renderTree shinyTree
#' @importFrom shinyMCE  tinyMCE
#' @importFrom magrittr %>%
#' @importFrom DT DTOutput dataTableProxy
#' @importFrom shinycssloaders withSpinner
#' @importFrom plotly event_data plot_ly
#' @importFrom shinyBS tipify bsAlert
#'
#'
#'


schnapps <- function(localContributionDir = "~/Rstudio/shHubgit/Dummy/",
                       defaultValueSingleGene = "CD52",
                       defaultValueMultiGenes = "CD52, S100A4, S100A9, S100A8",
                       defaultValueRegExGene = "", # tip: '^CD7$|^KIT$; genes with min expression
                       DEBUG = FALSE,
                       DEBUGSAVE = FALSE

                       ) {
  on.exit({
    rm(list = c(".SCHNAPPs_locContributionDir",
         ".SCHNAPPs_defaultValueSingleGene",
         ".SCHNAPPs_defaultValueMultiGenes",
         ".SCHNAPPs_defaultValueRegExGene",
         ".SCHNAPPs_DEBUG",
         ".SCHNAPPs_DEBUGSAVE"),
       envir = globalenv())
  })
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
#> [1] "scEx"
