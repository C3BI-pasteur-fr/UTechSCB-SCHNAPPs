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
#' @param workers number of CPUs to be used for normlization/clustering
#' @param launch.browser a function to call with the application's URL.
#'
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar sidebarMenu dashboardBody
#' tabItem menuSubItem menuItem
#' @importFrom shinyTree renderTree shinyTree
#' @importFrom DT DTOutput dataTableProxy
#' @importFrom dplyr '%>%'
#' @importFrom shiny shinyApp runApp
#'
#' @export setup_schnapps_env
#'
#' @examples
#' # create example data
#' data("scEx", package = "SCHNAPPs")
#' save(file = "scEx.Rdata", list = "scEx")
#' # use "scEx.Rdata" with load data functionality within the shiny app
# Create a setup function that extracts the initialization part
setup_schnapps_env <- function(localContributionDir = "~/Rstudio/shHubgit/Dummy/",
                               defaultValueSingleGene = "CD3g",
                               defaultValueMultiGenes = "cd3g, cd4, cd8b, ms4a1, TCF4, LILRA2, LYZ, cd79a, bcl11b, IL32, hbb, nkg7,MNDA",
                               defaultValueRegExGene = "",
                               DEBUG = FALSE,
                               DEBUGSAVE = FALSE,
                               historyPath = NULL,
                               defaultValues = list(),
                               workers = 2,
                               packagePath = NULL) {  # Add packagePath parameter
  
  # Initialize environment variables
  .schnappsEnv <- new.env(parent = emptyenv())
  
  # Assign variables to environment
  assign(".schnappsEnv", .schnappsEnv, envir = .GlobalEnv)  # Make it globally available
  
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
  
  # Install Wind if needed
  if (!"Wind" %in% installed.packages()) {
    devtools::install_github("haowulab/Wind", build_opts = c("--no-resave-data"))
  }
  
  # Set up parallel processing
  options(future.globals.maxSize = 2024^3)
  future::plan("multisession", workers = workers)
  
  # Handle packagePath
  if (is.null(packagePath)) {
    tryCatch({
      packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/")
    }, error = function(e) {
      # If SCHNAPPs package not found, ask for packagePath
      stop("SCHNAPPs package not found. Please provide packagePath parameter.")
    })
  }
  
  # Assign packagePath to environment
  assign("packagePath", packagePath, envir = .schnappsEnv)
  
  # Source required files
  message("Sourcing serverFunctions.R")
  source(paste0(packagePath, "/defaultValues.R"))
  source(paste0(packagePath, "/toolTips.R"))
  source(paste0(packagePath, "/shortCuts.def.R"))
  source(paste0(packagePath, "/serverFunctions.R"))
  
  
  # Set up BiocParallel
  BiocParallel::register(safeBPParam(workers))
  
  # Source server file
  message("Sourcing server.R")
  devscShinyApp <- FALSE
  source(paste0(packagePath, "/server.R"))
  
  # Create a new environment for functions
  fun_env <- new.env(parent = .GlobalEnv)
  
  # Source files into the function environment
  r_files = dir(path = normalizePath(paste0(packagePath, "/")),recursive = T,pattern = "reactives.R", full.names = T)
  for (fp in r_files){
    try(
      sys.source(fp, envir = fun_env)
    )
  }
  # Source files into the function environment
  r_files = dir(path = normalizePath(paste0(packagePath, "/")),recursive = T,pattern = "parameters.R", full.names = T)
  for (fp in r_files){
    try(
      sys.source(fp, envir = fun_env)
    )
  }
  
  # Attach the function environment or copy functions to .GlobalEnv
  attach(fun_env)  # Or use this to copy specific functions:
  # list2env(as.list(fun_env), envir = .GlobalEnv)
  
  return(list(
    env = .schnappsEnv,
    packagePath = packagePath,
    fun_env = fun_env
  ))
}

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
#' @param workers number of CPUs to be used for normlization/clustering
#' @param launch.browser a function to call with the application's URL.
#'
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar sidebarMenu dashboardBody
#' tabItem menuSubItem menuItem
#' @importFrom shinyTree renderTree shinyTree
#' @importFrom DT DTOutput dataTableProxy
#' @importFrom dplyr '%>%'
#' @importFrom shiny shinyApp runApp
#'
#' @export schnapps
#' @export setup_schnapps_env
#'
#' @examples
#' # create example data
#' data("scEx", package = "SCHNAPPs")
#' save(file = "scEx.Rdata", list = "scEx")
#' # use "scEx.Rdata" with load data functionality within the shiny app
# Create a setup function that extracts the initialization part
schnapps <- function(..., port = NULL, launch.browser = getOption("shiny.launch.browser", interactive())) {
  setup_result <- setup_schnapps_env(...)
  
  # Source UI and create app only if needed
  if (!is.null(launch.browser)) {
    source(paste0(setup_result$packagePath, "/ui.R"))
    app <- shinyApp(ui = scShinyUI, server = scShinyServer, enableBookmarking = "server")
    runApp(app, port = port, launch.browser = launch.browser)
  }
  
  invisible(setup_result)
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
"scExOrg"
