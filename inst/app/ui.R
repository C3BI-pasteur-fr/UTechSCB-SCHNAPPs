# inst/app/ui.R
# 
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#


suppressMessages(require(shinyjqui))
suppressMessages(require(shiny))
source(normalizePath(paste0(packagePath, "/toolTips.R")), local = TRUE)
suppressMessages(require(shinydashboard))
suppressMessages(require(shinydashboardPlus))
suppressMessages(require(plotly))
# suppressMessages(require(shinythemes))
suppressMessages(require(ggplot2))
suppressMessages(require(DT))
suppressMessages(require(edgeR))
suppressMessages(require(pheatmap))
# suppressMessages(require(threejs))
suppressMessages(require(shinyTree))
# suppressMessages(require(shinyjs))

if (exists("devscShinyApp")) {
  if (devscShinyApp) {
    packagePath <- paste0("inst", .Platform$file.sep, "app")
    # setwd("~/Rstudio/UTechSCB-SCHNAPPs/")
  } else {
    packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/") %>% normalizePath()
  }
} else {
  packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/") %>% normalizePath()
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

# "introjsUI"
if ("rintrojs" %in% rownames(installed.packages())) {
  suppressMessages(require(rintrojs))
} else {
  introjsUI = function(...) {}
  cat(file = stderr(), "Please install introjsUI: install.packages('rintrojs')")
}


scShinyUI <- function(request) {
  library(shinyjqui)
  # jsCode <- 'shinyjs.hidemenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = "none"; x.classList.remove("menu-open");};shinyjs.showmenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = "block"; x.classList.add("menu-open");};'
  
  # browser()
  # load from history directory the old input variable that use defaultValues function
  dvFile = normalizePath(paste0(.schnappsEnv$historyPath, "/defaultValues.RData"))
  if(file.exists(dvFile)){
    cp = load(file=dvFile)
    if("defaultValues" %in% cp){
      # .schnappsEnv$defaultValues = defaultValues
      assign("defaultValues", defaultValues, envir = .schnappsEnv)
    }else{
      warning("defaultValues file exist but no defaultValues\n\n")
    }
  }
  source(normalizePath(paste0(packagePath, "/modulesUI.R")), local = FALSE)
  source(normalizePath(paste0(packagePath, "/tabs.R")), local = TRUE)
  # general tabs
  allTabs <- list(
    inputTab(),
    shortCutsTab(),
    geneSelectionTab(),
    cellSelectionTab(),
    clusterParametersTab() %>% checkAllowed(env = .schnappsEnv)
    # ,
    # renameTab()
  )
  
  # parameters tab, includes basic normalization
  source(normalizePath(paste0(packagePath, "/parameters.R")), local = TRUE)
  base::source(normalizePath(paste0(packagePath, "/serverFunctions.R")), local = TRUE)
  
  # Basic menu Items
  allMenus <- list(
    shinydashboard::menuItem("input",
                             # id="inputID",
                             tabName = "input", icon = icon("folder")
    ),
    shinydashboard::menuItem("short cuts",
                             tabName = "shortCuts", icon = icon("gopuram")
    ),
    shinydashboard::menuItem("Parameters",
                             # id="parametersID",
                             tabName = "parameters", icon = icon("gopuram"),
                             parameterItems()
    ),
    shinydashboard::menuItem(" Cell selection",
                             # id="cellSelectionID",
                             tabName = "cellSelection", icon = icon("ello")
    ),
    shinydashboard::menuItem("Gene selection",
                             # id="geneSelectionID",
                             tabName = "geneSelection", icon = icon("atom")
    )
    # ,
    # shinydashboard::menuItem("rename projections",
    #   # id="geneSelectionID",
    #   tabName = "modifyProj", icon = icon("signature")
    # )
  )
  
  
  # parse all ui.R files under contributions to include in application
  uiFiles <- dir(path = c(normalizePath(paste0(packagePath, "/contributions")), localContributionDir), pattern = "ui.R", full.names = TRUE, recursive = TRUE)
  for (fp in uiFiles) {
    # fp = uiFiles[5]
    menuList <- list()
    tabList <- list()
    source(fp, local = TRUE)
    
    for (li in menuList) {
      if (length(li) > 0) {
        # if(DEBUG)cat(file=stderr(), paste("menuList:", length(allMenus)," ", li$children, "\n"))
        allMenus[[length(allMenus) + 1]] <- li
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
    "input", "short cuts", "Parameters", "General QC", " Cell selection", "Gene selection", "Co-expression",
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
  parFiles <- dir(path = c(normalizePath(paste0(packagePath, "/contributions")), localContributionDir), 
                  pattern = "parameters.R",
                  full.names = TRUE, recursive = TRUE)
  for (fp in parFiles) {
    tabList <- list()
    source(fp, local = TRUE)
    
    for (li in tabList) {
      if (length(li) > 0) {
        # if(DEBUG)cat(file=stderr(), paste(li$children[[1]], "\n"))
        # allTabs[[length(allTabs) + 1]] <- li
      }
    }
  }
  
  mulist <- list(
    inputTab(),
    geneSelectionTab()
  )
  getallTabs <- function() {
    # tags$div(list(inputTab(),
    #               geneSelectionTab()),
    #          # inputTab(),
    #          class = "tab-content"
    # )  ,
    tags$div(mulist,
             # inputTab(),
             class = "tab-content"
    )
    tags$div(
      allTabs,
      # inputTab(),
      class = "tab-content"
    )
  }
  # search for parameter contribution submenu items (menuSubItem)
  # parameterContributions = ""
  
  getallMenus <- function() {
    allMenus
  }
  controlbarContext =NULL
  if(file.exists(paste0(packagePath, .Platform$file.sep, "controlbarContext.R"))) 
    source(file = paste0(packagePath, .Platform$file.sep,  "controlbarContext.R"), local = TRUE)
  
  # jsCode <- "shinyjs.hidemenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = 'none'; x.classList.remove('menu-open');}; shinyjs.showmenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = 'block'; x.classList.add('menu-open');};"
  # jsCode2 <- "shinyjs.pageCol = function(params){console.log(params);$('body').css('color', params);}"
  dep=htmltools::htmlDependency("jqueryui", "3.6.0", c(href="shared/jqueryui"), script = "jquery-ui.min.js")
  shinyUI(
    shinydashboardPlus::dashboardPage(
      dheader(),
      
      shinydashboardPlus::dashboardSidebar(
        shinydashboard::sidebarMenu(
          style = "height: 90vh; overflow-y: auto;", 
          id = "sideBarID",
          getallMenus(),
          htmlOutput("summaryStatsSideBar"),
          
          # downloadButton("report", "Generate report", class = "butt"),
          tags$head(tags$style(".butt{color: black !important;}")), #  font color; otherwise the text on these buttons is gray
          tags$head(tags$style(HTML("berndTest{background-color:rgba(255,34,22,0.1);}"))), # supposed to change the transparency of introjs area that is highlighted.
          
          # bookmarkButton(id = "bookmark1"),
          br(),
          downloadButton("countscsv", "Download (log) counts.csv", class = "butt"),
          br(),
          downloadButton("RDSsave", "Download RData", class = "butt"),
          br(),
          downloadButton("RmdSave", "Download History", class = "butt"),
          if (DEBUG) sc_checkboxInput("DEBUGSAVE", "Save for DEBUG", FALSE),
          if (DEBUG) verbatimTextOutput("DEBUGSAVEstring"),
          if (is.environment(.schnappsEnv)) {
            if (exists("historyPath", envir = .schnappsEnv)) {
              # sc_checkboxInput("save2History", "save to history file", FALSE)
              actionButton("comment2History", "Add comment to history")
            }
          },
          if (DEBUG) {
            actionButton("openBrowser", "open Browser")
          },
          # bookmarkButton(),
          actionButton("Quit", "quit")
        )
        # ,
        # verbatimTextOutput("save2Historystring")
        # ,verbatimTextOutput("currentTabInfo")
      ), # dashboard side bar
      shinydashboard::dashboardBody(
        tags$script(HTML("$('body').addClass('fixed');")),
        # shinyjs is not working jqui_ and since we need resize etc we cannot use shinyjs
        # understanding the java script - shiny interaction takes a bit too much time
        # https://shiny.rstudio.com/articles/communicating-with-js.html
        shinyjs::useShinyjs(debug = TRUE),
        # shinyjs::extendShinyjs(text = jsCode, functions = c("hidemenuItem", "showmenuItem")),
        # extendShinyjs(text = jsCode2, functions = c("pageCol")),
        tags$script("showmenuItem = function(targetid) {var x = document.getElementById(targetid); x.style.display = 'block'; x.classList.add('menu-open');};"),
        # alert('finished loading');
        tags$script("$(document).on('shiny:connected', function(){
  var startingUp = true;
  Shiny.setInputValue('startingUp', '1');
  function onTimeout() {
    if(startingUp){
      startingUp = 0;
      Shiny.setInputValue('startingUp', '0');

    };
  };
  function startIdleTimer() {
    if (idleTimer) clearTimeout(idleTimer);
    idleTimer = setTimeout(onTimeout, timeoutWarningMsecs);
  };
  $(document).on('shiny:idle', onTimeout);
});"),
        introjsUI(),
        shinyjs::inlineCSS(list(.red = "background-color: DarkSalmon; hover: red")),
        shinyjs::inlineCSS(list(.green = "background-color: lightgreen")),
        getallTabs(),
        
        ### !!!! https://github.com/Yang-Tang/shinyjqui/issues/87
        ### not working resize shinyjqui
        # if(DEBUG){
        #   library(profvis)
        #   library(shiny)
        #   profvis_ui("profiler")
        # },
        
        tags$head(tags$style(HTML("div.box-header {display: block;}"))),
        tags$head(tags$style(HTML("h3.box-title {display: block;}")))
        # tags$head(
        #   tags$script(version    = "1.12.1",
        #               src        = "www/shared/jqueryui",
        #               script     = "jquery-ui.min.js"
        #               ))
        # tags$div(list(inputTab(),
        #               geneSelectionTab()),
        #          # inputTab(),
        #          class = "tab-content"
        # )  ,
        # h4("Sum of all previous slider values:", textOutput("sum"))
      ), # dashboard body
      
      options = list(sidebarExpandOnHover = TRUE),
      controlbar = controlbarContext
    ) # main dashboard
  )
}

