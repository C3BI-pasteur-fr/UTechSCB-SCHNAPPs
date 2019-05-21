require(shinyjqui)
library(magrittr)

# clusterUI -----------
# to select clusters from the list of available knn clusters
# Flexible 2 plot
# TODO select colour of points
clusterUI <- function(id) {
  if (DEBUG) {
    cat(file = stderr(), paste("clusterUI: ", NS(id)("clusters"), "\n"))
  }
  ns <- NS(id)
  tagList(
    fluidRow(
      column(6,
             offset = 0,
             textInput(ns("geneIds"), "comma separated list of genes for UmiCountPerGenes", value = "")),
      column(6,
             offset = 0,
             textInput(ns("geneIds2"), "comma separated list of genes for UmiCountPerGenes2", value = "")
      )
    ),
    fluidRow(
      # column(
      #   3,
      #   uiOutput(ns("clusters"))
      # ),
      column(
        4,
        selectInput(
          ns("dimension_x"),
          label = "X",
          choices = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne1"
        )
      ),
      column(
        4,
        selectInput(
          ns("dimension_y"),
          label = "Y",
          choices = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne2"
        )
      ),
      column(
        4,
        selectInput(
          ns("dimension_col"),
          label = "color",
          choices = c("Gene.count"),
          selected = "Gene.count"
        )
      )
    ),
    fluidRow(column(
      11,
      # jqui_resizable(plotOutput(ns("clusterPlot"), brush = brushOpts(id = ns("b1"))) )
      jqui_resizable(plotly::plotlyOutput(ns("clusterPlot")))
    )),
    fluidRow(
      checkboxInput(ns("moreOptions"), "show more options", FALSE),
      uiOutput(ns("additionalOptions"))
      # checkboxInput(ns("showCells"), "show cell names", FALSE),
      #
      # verbatimTextOutput(ns('cellSelection'))
    )
  )
}

# tableSelectionUi --------------
tableSelectionUi <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(div(
      h5("Selected itmes to be copied"),
      align = "left"
    )),
    fluidRow(
      verbatimTextOutput(ns("cellSelection"))
    ),
    fluidRow(
      downloadButton(ns("download_cellNameTable"), "Download table")
    ),
    fluidRow(
      h4("Cells", offset = 1),
      column(3,
             checkboxInput(ns("selectAll"), "Select all rows", FALSE)), 
      column(3,
             checkboxInput(ns("reorderCells"), "reorder cells by sum of selected genes", FALSE)), 
      br(),
      column(
        width = 12,
        DT::DTOutput(ns("cellNameTable")) %>% shinycssloaders::withSpinner(),
        style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
      )
    ),
    fluidRow()
  )
}

# pHeatMapUI --------------
pHeatMapUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(tags$h3("Heatmap plot")),
    fluidRow(column(11,
                    jqui_resizable(plotOutput(ns("pHeatMapPlot"),
                                              # height = "auto",
                                              brush = brushOpts(id = "crh1")
                    ),options = list( width="95%"))
    )
    ),
    fluidRow(
      checkboxInput(ns("moreOptions"), "show more options", FALSE),
      uiOutput(ns("additionalOptions")),
      downloadButton(ns("download_pHeatMapUI"), "Download PlotData")
      # checkboxInput(ns("showCells"), "show cell names", FALSE),
      #
      # verbatimTextOutput(ns('cellSelection'))
    )
  )
}
