###
### Liana contributions UI elements
###

# Check if needed packages are installed, if not, don't show ui elements.
# This should prevent the reactives etc from ever being executed.

if(!"liana" %in% installed.packages()){
  cat(file = stderr(), paste("please install liana:
     remotes::install_github('saezlab/OmnipathR', dependencies = T)
     remotes::install_github('saezlab/liana', dependencies = T)
  "))
  return(NULL)
}
require(liana)



menuList <- list(
  shinydashboard::menuItem("Liana",
                           icon = icon("dashboard"),
                           tabName = "Liana", startExpanded = FALSE,
                           shinydashboard::menuSubItem("Liana basic", tabName = "LianaBasic")
  )
)

# Liana ----
tabList = list(
  shinydashboard::tabItem(
    "LianaBasic",
    box(
      title = "Ligand - receptor analysis using Liana", solidHeader = TRUE, width = 12, status = "primary",
      footer = "LIANA enables the use of any combination of ligand-receptor methods and resources, and their consensus.",
      fluidRow(
        column(
          width = 12, offset = 1,
          actionButton("updateLianaParameters", "apply changes", width = "80%")
        )
      ),
      br(),
      fluidRow(
        cellSelectionUI("Liana_dataInput"),
        box(
          fluidRow(
            column(width = 3, 
                   sc_selectizeInput(inputId = "Liana_idents_col", label = "the cell labels to be used.", 
                                     choices = c("dbCluster"), selected = defaultValue("Liana_idents_col", "dbCluster"))),
            column(width = 3, 
                   sc_selectizeInput(inputId = "Liana_resource", label = "resource(s) to be used by the methods", 
                                     choices = liana::show_resources(), selected = defaultValue("Liana_resource", "OmniPath"))),
            column(width = 3, 
                   sc_selectizeInput(inputId = "Liana_method", label = "method(s) to be run via liana",
                                     multiple = T,
                                     choices = liana::show_methods() , selected = defaultValue("Liana_method", "natmi"))),
            column(width = 3,
                   sc_numericInput("Liana_min_cells", "minimum cell per cell identity to be considered for analysis",
                                   defaultValue("Liana_min_cells", 5),
                                   min = 2, max = 10000
                   )
            )
            # parallelize
            # whether to parallelize cellphonedb-like
            # 
            # workers
            # number of workers to be called
          ),
          checkbsTT("Liana_idents_col"),
          checkbsTT("Liana_resource"),
          checkbsTT("Liana_method"),
          checkbsTT("Liana_min_cells")
        ),
        
      ),
      br(),
      fluidRow(column(
        12,
        jqui_resizable(plotlyOutput("Liana_dotPlot")) # %>% withSpinner()
      )),
      br(),
      actionButton("save2Hist_Liana_dotPlot", "save to history"),
      br(),
      fluidRow(column(
        12,
        jqui_resizable(plotOutput("Liana_Heatmap")) # %>% withSpinner()
      )),
      br(),
      actionButton("save2Hist_Liana_Heatmap", "save to history"),
      br(),
      fluidRow(
        column(
          width = 10,
          tableSelectionUi("Liana_raw_TableMod")
        )
      ),
      br(),
      fluidRow(column(width = 6, 
             sc_selectizeInput(inputId = "Liana_method_show", label = "results to show",
                               multiple = F,
                               choices = c("none") , selected = defaultValue("none", "none")))
             ),
      fluidRow(
        column(
          width = 10,
          tableSelectionUi("Liana_all_TableMod")
        )
      )
    )
    
  )
)