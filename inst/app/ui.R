
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#


require(shiny)
require(shinydashboard)
require(plotly)
require(shinythemes)
require(ggplot2)
require(DT)
require(edgeR)
require(pheatmap)
require(threejs)
require(shinyTree)
require(shinycssloaders)
require(shinyBS)

if (exists("devscShinyApp")) {
  if (devscShinyApp) {
    packagePath <- "inst/app"
  }
} else {
  packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE)
  packagePath <- paste0(packagePath, "/app/")
}
localContributionDir <- .SCHNAPPs_locContributionDir
defaultValueSingleGene <- .SCHNAPPs_defaultValueSingleGene
defaultValueMultiGenes <- .SCHNAPPs_defaultValueMultiGenes
defaultValueRegExGene <- .SCHNAPPs_defaultValueRegExGene
DEBUG <- .SCHNAPPs_DEBUG
DEBUGSAVE <- .SCHNAPPs_DEBUGSAVE

# source(paste0(packagePath,  "/ui.R"))

# this is where the general tabs are defined:
# if (file.exists(paste0(packagePath, "/defaultValues.R"))) {
#   source(paste0(packagePath, "/defaultValues.R"))
# }
# input, cell/gene selection tabs
# source('tabs.R',  local = TRUE)
source(paste0(packagePath, "/tabs.R"), local = FALSE)

# general tabs
allTabs <- list(
  inputTab,
  geneSelectionTab,
  cellSelectionTab,
  generalParametersTab
)
# parameters tab, includes basic normalization
source(paste0(packagePath, "/parameters.R"), local = TRUE)

# Basic menu Items
allMenus <- list(
  shinydashboard::menuItem("input",
    # id="inputID",
    tabName = "input", icon = icon("folder")
  ),
  shinydashboard::menuItem("Parameters",
    # id="parametersID",
    tabName = "parameters", icon = icon("gopuram"), parameterItems
  ),
  shinydashboard::menuItem(" Cell selection",
    # id="cellSelectionID",
    tabName = "cellSelection", icon = icon("ello")
  ),
  shinydashboard::menuItem("Gene selection",
    # id="geneSelectionID",
    tabName = "geneSelection", icon = icon("atom")
  )
)


# parse all ui.R files under contributions to include in application
uiFiles <- dir(path = c(paste0(packagePath, "/contributions"), localContributionDir), pattern = "ui.R", full.names = TRUE, recursive = TRUE)
for (fp in uiFiles) {
  menuList <- list()
  tabList <- list()
  source(fp, local = TRUE)

  for (li in menuList) {
    if (length(li) > 0) {
      # if(DEBUG)cat(file=stderr(), paste("menuList:", length(allMenus)," ", li$children, "\n"))
      allMenus[[length(allMenus) + 1 ]] <- li
    }
  }
  for (li in tabList) {
    if (length(li) > 0) {
      # if(DEBUG)cat(file=stderr(), paste(li$children[[1]], "\n"))
      allTabs[[length(allTabs) + 1]] <- li
    }
  }
}


mListNames <- c()
for (menuListItem in 1:length(allMenus)) {
  mListNames[menuListItem] = allMenus[[menuListItem]][3][[1]][[1]][3]$children[[2]]$children[[1]][1]
}
sollOrder <- c("input", "Parameters", "General QC", " Cell selection", "Gene selection", "Co-expression",
  "Data Exploration", "Subcluster analysis")
sollOrderIdx = c()
for (sIdx in 1:length(sollOrder)) {
  sollOrderIdx[sIdx] = which(sollOrder[sIdx] == mListNames)
}
sollOrderIdx = c(sollOrderIdx, which(!1:length(allMenus) %in%  sollOrderIdx))

allMenus <- allMenus[sollOrderIdx]

# todo
# parse all parameters.R files under contributions to include in application
# allTabs holds all tabs regardsless of their location in the GUI
parFiles <- dir(path = c(paste0(packagePath, "/contributions"), localContributionDir), pattern = "parameters.R", full.names = TRUE, recursive = TRUE)
for (fp in parFiles) {
  tabList <- list()
  source(fp, local = TRUE)

  for (li in tabList) {
    if (length(li) > 0) {
      # if(DEBUG)cat(file=stderr(), paste(li$children[[1]], "\n"))
      allTabs[[length(allTabs) + 1]] <- li
    }
  }
}

# search for parameter contribution submenu items (menuSubItem)
# parameterContributions = ""



scShinyUI <- shinyUI(
  shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = "SCHNAPPs"),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        id = "sideBarID",
        allMenus
      ),
      # shinyBS::tipify(
      #   checkboxInput("noStats", "don't display stats", FALSE),
      #   "check this if you are working on the cell/gene selection to avoid certain calculations"
      # ),

      shinyBS::tipify(
        htmlOutput("summaryStatsSideBar"),
        "<h3>Data summary</h3> <ul><li>medium UMI: shows how many genes are  expressed in log2 space of normalized data</li> </ul> ", "right"
      ),
      
      ### failed tests to change color of button text
      
      # <div id="scoped-content">
      #   <style type="text/css" scoped>
      #   h1 { color: red; } 
      # </style>
      #   
      #   <h1>Hello</h1>
      #   </div>
      #   
      # withTags({
      #   div(id="scoped-content",
      #       tags$style(type="text/css", scoped=NA, "downloadbutton a { color: red; }"),
      #       tags$a(id = "report", class = paste("btn btn-default shiny-download-link ", "downloadbutton"), href = "", target = "_blank", download = NA, 
      #              icon("download"), "Generate report"),
      #       downloadButton("report", "Generate report", class="downloadbutton")
      #   )
      # }),
      ## <div class="header" checked>
      ##   <p>Ready to take the Shiny tutorial? If so</p>
      ##   <a href="shiny.rstudio.com/tutorial">Click Here!</a>
      ## </div> 
      # tags$style(type="text/css", "downloadbutton a { color: #444; }"),
      downloadButton("report", "Generate report", class="butt"),
      tags$head(tags$style(".butt{color: black !important;}")), #  font color
      
      # commentetd out because currently no-one is using it
      # actionButton("goCalc", "Force Calculations"),
      
      # bookmarkButton(id = "bookmark1"),
      shinyBS::tipify(
        downloadButton("countscsv", "Download counts.csv", class="butt"),
        "<h3>download current normalized count data as CSV file</h3>"
      ),
      shinyBS::tipify(
        downloadButton("RDSsave", "Download Rds", class="butt"),
        "<h3>download current cell/gene configuration for reimport to this app</h3>"
      ),
      if (DEBUG) checkboxInput("DEBUGSAVE", "Save for DEBUG", FALSE),
      verbatimTextOutput("DEBUGSAVEstring")
    ), # dashboard side bar
    shinydashboard::dashboardBody(
      shinyBS::bsAlert("alert"),
      tags$div(
        allTabs,
        class = "tab-content"
      )
    ) # dashboard body
  ) # main dashboard
)
