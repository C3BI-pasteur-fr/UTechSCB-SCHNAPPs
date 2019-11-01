
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#


suppressMessages(require(shiny))
source(paste0(packagePath,  "/toolTips.R"), local = TRUE)
suppressMessages(require(shinydashboard))
suppressMessages(require(plotly))
suppressMessages(require(shinythemes))
suppressMessages(require(ggplot2))
suppressMessages(require(DT))
suppressMessages(require(edgeR))
suppressMessages(require(pheatmap))
suppressMessages(require(threejs))
suppressMessages(require(shinyTree))

if (exists("devscShinyApp")) {
  if (devscShinyApp) {
    packagePath <- "inst/app"
  }
} else {
  packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/")
}
localContributionDir <- get(".SCHNAPPs_locContributionDir", envir = .schnappsEnv)
defaultValueSingleGene <- get(".SCHNAPPs_defaultValueSingleGene", envir = .schnappsEnv)
defaultValueMultiGenes <- get(".SCHNAPPs_defaultValueMultiGenes", envir = .schnappsEnv)
defaultValueRegExGene <- get(".SCHNAPPs_defaultValueRegExGene", envir = .schnappsEnv)
DEBUG <- get(".SCHNAPPs_DEBUG", envir = .schnappsEnv)
DEBUGSAVE <- get(".SCHNAPPs_DEBUGSAVE", envir = .schnappsEnv)

# source(paste0(packagePath,  "/ui.R"))

# this is where the general tabs are defined:
# if (file.exists(paste0(packagePath, "/defaultValues.R"))) {
#   source(paste0(packagePath, "/defaultValues.R"))
# }
# input, cell/gene selection tabs
# source('tabs.R',  local = TRUE)
source(paste0(packagePath,  "/modulesUI.R"), local = FALSE)
source(paste0(packagePath, "/tabs.R"), local = TRUE)
# general tabs
allTabs <- list(
  inputTab,
  geneSelectionTab,
  cellSelectionTab,
  generalParametersTab,
  renameTab
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
  ),
  shinydashboard::menuItem("rename projections",
                           # id="geneSelectionID",
                           tabName = "renameProj", icon = icon("signature")
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
  mListNames[menuListItem] <- allMenus[[menuListItem]][3][[1]][[1]][3]$children[[2]]$children[[1]][1]
}
sollOrder <- c(
  "input", "Parameters", "General QC", " Cell selection", "Gene selection", "Co-expression",
  "Data Exploration", "Subcluster analysis"
)
sollOrderIdx <- c()
for (sIdx in 1:length(sollOrder)) {
  sollOrderIdx[sIdx] <- which(sollOrder[sIdx] == mListNames)
}
sollOrderIdx <- c(sollOrderIdx, which(!1:length(allMenus) %in% sollOrderIdx))

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
      htmlOutput("summaryStatsSideBar"),
      
      downloadButton("report", "Generate report", class = "butt"),
      tags$head(tags$style(".butt{color: black !important;}")), #  font color

      # bookmarkButton(id = "bookmark1"),
      downloadButton("countscsv", "Download (log) counts.csv", class = "butt"),
      downloadButton("RDSsave", "Download RData", class = "butt"),
      if (DEBUG) checkboxInput("DEBUGSAVE", "Save for DEBUG", FALSE),
      verbatimTextOutput("DEBUGSAVEstring"),
      if (exists("historyFile", envir = .schnappsEnv)){
        checkboxInput("save2History", "save to history file", FALSE)
      },
      verbatimTextOutput("save2Historystring")
    ), # dashboard side bar
    shinydashboard::dashboardBody(
      tags$div(
        allTabs,
        class = "tab-content"
      )
    ) # dashboard body
  ) # main dashboard
)
