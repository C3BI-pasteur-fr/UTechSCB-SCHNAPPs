

# Specific ui elements for the lite version of SCHNAPPs.
# we are overwriting essential reactives to remove any manipulation of the underlying data
# this allows for faster computations.
# 
# we require that ui.R has already been loaded.
#


### might remove the following if we include the source(ui.R)
suppressMessages(require(shiny))
source(paste0(packagePath, "/toolTips.R"), local = TRUE)
suppressMessages(require(shinydashboard))
suppressMessages(require(plotly))
suppressMessages(require(shinythemes))
suppressMessages(require(ggplot2))
suppressMessages(require(DT))
suppressMessages(require(edgeR))
suppressMessages(require(pheatmap))
suppressMessages(require(threejs))
suppressMessages(require(shinyTree))
suppressMessages(require(shinyjs))



introTab <- function(){
  shinydashboard::tabItem(
    "Intro",
    fluidPage(
      uiOutput('introRMD')
    )
  )
}


generalParametersTab <- function(){
  shinydashboard::tabItem(
    "generalParameters",
    fluidRow(div(h2("General parameters"), align = "center")),
    br(),
    box(
      title = "Colors", solidHeader = TRUE, width = 12, status = "primary",
      collapsible = F, collapsed = F,
      fluidRow(column(
        width = 12, offset = 1,
        actionButton("updateColors", "apply changes", width = "80%")
      )),
      br(),
      fluidRow(
        column(
          width = 6,
          uiOutput("sampleColorSelection")
        ),
        column(
          width = 6,
          uiOutput("clusterColorSelection")
        )
      )
    ),
    checkbsTT(item = "updateColors"),
    checkbsTT(item = "sampleColorSelection"),
    checkbsTT(item = "clusterColorSelection")
  )
}





scShinyUI <- function(request) {
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
  
  base::source(paste0(packagePath, "/serverFunctions.R"), local = TRUE)
  
  # source(paste0(packagePath,  "/ui.R"))
  
  # this is where the general tabs are defined:
  # if (file.exists(paste0(packagePath, "/defaultValues.R"))) {
  #   source(paste0(packagePath, "/defaultValues.R"))
  # }
  # input, cell/gene selection tabs
  # source('tabs.R',  local = TRUE)
  source(paste0(packagePath, "/modulesUI.R"), local = FALSE)
  # source(paste0(packagePath, "/tabs.R"), local = TRUE)
  
  # remove up to here ------------------
  
  # we still want to be able to change the colors
  introItem <- function(){
    shinydashboard::menuSubItem("Introduction", tabName = "Intro")
  }
  parameterItems <- list(
    
    shinydashboard::menuSubItem("General Parameters", tabName = "generalParameters"),
    shinydashboard::menuSubItem("Projections", tabName = "modifyProj")
  )
  

  
  # general tabs
  allTabs <- list(
    introTab(),
    generalParametersTab()
    # ,
    # modTab
  )
  # parameters tab, includes basic normalization
  source(paste0(packagePath, "/parameters.R"), local = TRUE)
  
  # Basic menu Items
  allMenus <- list(

    shinydashboard::menuItem("Introduction",
                             tabName = "Intro", icon = icon("dashboard")
    ),
    shinydashboard::menuItem("Parameters",
                             # id="parametersID",
                             tabName = "parameters", icon = icon("gopuram"), parameterItems
    )
    # ,
    # shinydashboard::menuItem("rename projections",
    #                          # id="geneSelectionID",
    #                          tabName = "modifyProj", icon = icon("signature")
    # )
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
    "Introduction", "Parameters", "General QC", "Co-expression",
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
  
  
  
  shinyUI(
    shinydashboard::dashboardPage(
      dheader(),
      # shinydashboard::dashboardHeader(title = paste("SCHNAPPs" , packageVersion("SCHNAPPs"))),
      shinydashboard::dashboardSidebar(
        shinydashboard::sidebarMenu(
          id = "sideBarID",
          allMenus
        ),
        htmlOutput("summaryStatsSideBar"),
        
        # downloadButton("report", "Generate report", class = "butt"),
        tags$head(tags$style(".butt{color: black !important;}")), #  font color; otherwise the text on these buttons is gray
        tags$head(tags$style(HTML("berndTest {background-color : rgb(255,34,22,0.1);}"))),
        
        # bookmarkButton(id = "bookmark1"),
        br(),
        downloadButton("countscsv", "Download (log) counts.csv", class = "butt"),
        br(),
        downloadButton("RDSsave", "Download RData", class = "butt"),
        br(),
        downloadButton("RmdSave", "Download History", class = "butt"),
        if (DEBUG) checkboxInput("DEBUGSAVE", "Save for DEBUG", FALSE),
        if (DEBUG) verbatimTextOutput("DEBUGSAVEstring"),
        if (exists("historyPath", envir = .schnappsEnv)){
          br()
        },
        if (exists("historyPath", envir = .schnappsEnv)){
          # checkboxInput("save2History", "save to history file", FALSE)
          actionButton("comment2History", "Add comment to history")
          
        },
        if (DEBUG) {
          br()
        },
        if (DEBUG) {
          actionButton("openBrowser", "open Browser")
        }
        # ,
        # verbatimTextOutput("save2Historystring")
        # ,verbatimTextOutput("currentTabInfo")
      ), # dashboard side bar
      shinydashboard::dashboardBody(
        shinyjs::useShinyjs(debug = TRUE),
        rintrojs::introjsUI(),
        inlineCSS(list(.red = "background-color: DarkSalmon; hover: red")),
        inlineCSS(list(.green = "background-color: lightgreen")),
        tags$div(
          allTabs,
          class = "tab-content"
        )
      ) # dashboard body
    ) # main dashboard
  )
}

if (DEBUG) cat(file = stderr(), "end: ui-lite\n")