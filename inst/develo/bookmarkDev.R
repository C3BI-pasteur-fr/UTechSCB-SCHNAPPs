library(shiny)
library(shinydashboard)

# base::source(file = "inst/develo/bookmarkUI.R", local = TRUE)
inputTab <- function(){
  shinydashboard::tabItem(
    tabName = "input",
    radioButtons("whichscLog",
                 label = "Compute normalizations?",
                 choices = c(
                   "disable log" = "disablescEx_log",
                   "use scEx from loaded data" = "useLog",
                   "calculate normalization here" = "calcLog"
                 ),
                 selected = "disablescEx_log"
    )
  )
}
# mulist <- list(inputTab())

# mulist2 <- function(mulist){mulist}
# mulist3 <- function() {list(inputTab())}
# mulist4 <- function() {mulist}

ui <- function(request) {
  mulist <- list(inputTab())
  
  shinyUI(
    shinydashboard::dashboardPage(
      shinydashboard::dashboardHeader(title = "SCHNAPPs"),
      shinydashboard::dashboardSidebar(
        # sideBar()
        shinydashboard::sidebarMenu(
          id = "sideBarID",
          shinydashboard::menuItem("input",
                                   # id="inputID",
                                   tabName = "input", icon = icon("folder")
          ),
        bookmarkButton()
        )), # dashboard side bar
      shinydashboard::dashboardBody(
        tags$div(
          class = "tab-content",
          mulist
          # list(inputTab())
        )
      )
   )
  )
  }

server <- function(input, output, session) {
}

sa = shinyApp(ui, server, enableBookmarking = "url")
runApp(sa)
