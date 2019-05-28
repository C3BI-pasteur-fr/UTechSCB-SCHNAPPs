library(magrittr)
menuList <- list(
  shinydashboard::menuItem("Subcluster analysis", 
                           icon = icon("bar-chart-o"),
           # id = 'subclusterID',
           tabName = "subcluster", startExpanded = FALSE,
    shinydashboard::menuSubItem("DGE analysis", tabName = "sCA_dge")
  )
)

tabList <- list(
  dgeTab = shinydashboard::tabItem(
    "sCA_dge",
    tags$ul(
      tags$li(
        strong("Subclustering"),
        paste(
          ":Select a group of cells in plot1 and a different group of cells in plot2",
          "for identifying differential features between these subclusters"
        )
      ),
      tags$li(
        strong("colors:"),
        paste("colored by cluster identity")
      ),
      tags$li(
        strong("selection hint:"),
        paste("you can also slect by groups you have defined in other plots.")
      ),
      tags$li(
        strong("selection hint:"),
        paste('also check out "Gene.count" to verify that number genes per cell.')
      )
    ),

    fluidRow(
      column(
        2,
        uiOutput("sCA_dgeClustersSelection")
      ),
      column(
        2,
        selectInput(
          "sCA_subscluster_x1",
          label = "X",
          choice = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne1"
        )
      ),
      column(
        2,
        selectInput(
          "sCA_subscluster_y1",
          label = "Y",
          choice = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne2"
        )
      )
    ),

    fluidRow(
      column(
        6,
        plotOutput("sCA_dge_plot1", brush = brushOpts(
          id =
            "db1"
        )) %>% shinycssloaders::withSpinner()
      ),
      column(
        6,
        plotOutput("sCA_dge_plot2", brush = brushOpts(
          id =
            "db2"
        )) %>% shinycssloaders::withSpinner()
      )
    ),
    shinydashboard::tabItem(
      "diffExpMethod",
      list(
        tags$h3("Method to use for differential gene expression analysis"),
        fluidRow(column(
          10,
          radioButtons(
            inputId = "sCA_dgeRadioButton",
            label = "Method to use",
            choices = "dgeChoices",
            selected = "DE_logNormalization",
            width = "100%"
          )
        ))
        # ,
        # fluidRow(column(10, verbatimTextOutput("sCA_dgeRadioButtonValue"))),
        # wellPanel(
        #   # This outputs the dynamic UI component
        #   uiOutput("dgeParametersDynamic")
        # )
      )
    ),
    
    # fluidRow(column(11,
    #   offset = 1,
    #   h4("Selected genes"), br(),
    #   textOutput("crSelectedGenes", inline = FALSE)
    # )), br(),
    fluidRow(column(11,
      offset = 0,
      h4("Differentially Expressed Genes", offset = 1),
      tableSelectionUi("sCA_dgeTable")
    ))
  )
)
