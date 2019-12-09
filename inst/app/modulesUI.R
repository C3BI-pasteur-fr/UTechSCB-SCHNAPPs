suppressMessages(library(magrittr))

if ("shinycssloaders" %in% rownames(installed.packages())) {
  suppressMessages(library(shinycssloaders))
} else {
  withSpinner <-function(x) {
    x
  }
}
if ("shinyjqui" %in% rownames(installed.packages())) {
  suppressMessages(require(shinyjqui))
} else {
  jqui_resizable <- function(x, ...) {
    x
  }
}
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
      column(width = 6,
             textInput(ns("geneIds"), "comma separated list of genes for UmiCountPerGenes", value = "")),
      column(width = 6,
             textInput(ns("geneIds2"), "comma separated list of genes for UmiCountPerGenes2", value = "")
      )
    ),
    fluidRow(
      column(width = 4,
             selectInput(ns("dimension_x"),
                         label = "X",
                         choices = c("tsne1", "tsne2", "tsne3"),
                         selected = "tsne1"
             )
      ),
      column(width = 4,
             selectInput(ns("dimension_y"),
                         label = "Y",
                         choices = c("tsne1", "tsne2", "tsne3"),
                         selected = "tsne2"
             )
      ),
      column(width = 4,
             selectInput(ns("dimension_col"),
                         label = "color",
                         choices = c("Gene.count"),
                         selected = "Gene.count"
             )
      )
    ),
    box(
      title = "plot", solidHeader = TRUE, width = 12, status = 'primary', 
      fluidRow(
        column(width = 12,
               jqui_resizable(plotly::plotlyOutput(ns("clusterPlot")))
        )
      ),
      box(
        title = "additional options", solidHeader = TRUE, width = 12, status = 'primary', 
        collapsible = TRUE, collapsed = TRUE,
        fluidRow(
          column(width = 3,
                 checkboxInput(ns("logX"), "log transform X", value = FALSE)
          ),
          column(width = 3,
                 checkboxInput(ns("logY"), "log transform Y", value = FALSE)
          ),
          column(width = 3,
                 selectInput(ns("divideXBy"),
                             label = "Divide X by",
                             # choices = c("None", "Gene.count", "UMI.count"),
                             choices = c("None", "UmiCountPerGenes", "UmiCountPerGenes2"),
                             selected = "None"
                 )
          ),
          column(width = 3,
                 selectInput(ns("divideYBy"),
                             label = "Divide Y by",
                             # choices = c("None", "Gene.count", "UMI.count"),
                             choices = c("None", "UmiCountPerGenes", "UmiCountPerGenes2"),
                             selected = "None"
                 )
          )
        ),
        fluidRow(
          column(width = 12,
                 checkboxInput(ns("addToGroup"), "Add to group/otherwise overwrite", TRUE),
                 textInput(ns(id = "groupName"), label = "name group, also used in Plot to color selected cells red.", value = "cellGroupName"),
                 selectInput(ns("groupNames"),
                             label = "group names, !When modifying a group this list of cells is used as a reference!",
                             choices = c("plot"),
                             selected = "plot"
                 ),
                 verbatimTextOutput(ns("nCellsVisibleSelected")),
                 actionButton(ns("changeGroups"), "change current selection"),
                 checkboxInput(ns("showCells"), "show cell names", FALSE),
                 verbatimTextOutput(ns("cellSelection")),
                 actionButton(ns("save2Hist"), "save to history"),
                 uiOutput(ns("additionalOptions")) # TODO:is this still needed???
          )
        )
      )
    )
    # ,
    # fluidRow(column(width = 6,
    #                 checkboxInput(ns("moreOptions"), "show more options", FALSE),
    # )
    # )
  )
}

# tableSelectionUi --------------
tableSelectionUi <- function(id) {
  ns <- NS(id)
  tagList(
    box(width = 12,
        fluidRow(
          column(width = 12,
                 div(
                   h5("Selected itmes to be copied"),
                   align = "left"
                 ),
                 verbatimTextOutput(ns("rowSelection"))
          )
        ),
        fluidRow(
          column(width = 3,
          downloadButton(ns("download_cellNameTable"), "Download table")
          ),
          column(width = 3,
                 actionButton(ns("save2HistTabUi"), "Save to history"),
                 )
        ),
        fluidRow(
          h4("Cells", offset = 1),
          column(width = 3,
                 checkboxInput(ns("selectAll"), "Select all rows", FALSE)), 
          column(width = 3,
                 checkboxInput(ns("reorderCells"), "reorder cells by sum of selected genes", FALSE)), 
          br(),
          column(width = 12,
                 DT::DTOutput(ns("cellNameTable")) %>% withSpinner(),
                 style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
          )
        )
    )
  )
}

# pHeatMapUI --------------
pHeatMapUI <- function(id) {
  ns <- NS(id)
  tagList(
    box( width = 12,
         fluidRow(
           column(width = 12,
                  jqui_resizable(plotOutput(ns("pHeatMapPlot"),
                                            # height = "auto",
                                            brush = brushOpts(id = "crh1")
                  ),options = list( width="99%"))
           )
         ),
         box(
           title = "additional options", solidHeader = TRUE, width = 12, status = 'primary', 
           collapsible = TRUE, collapsed = TRUE,
           fluidRow(
             column(width = 12,
                    # checkboxInput(ns("moreOptions"), "show more options", FALSE),
                    checkboxInput(ns("showColTree"), label = "Show tree for cells", value = FALSE),
             )
           ),
           fluidRow(
             column(width = 6,
                    selectInput(ns("normRow"),
                                label = "scale by row (for color)",
                                choices = c("row", "column", "none"),
                                selected = "none"
                    ),
                    selectInput(
                      ns("ColNames"),
                      label = "group names",
                      choices = c(),
                      selected = "sampleNames",
                      multiple = TRUE
                    ),
                    
                    selectInput(
                      ns("orderNames"),
                      label = "order of columns",
                      choices = c(),
                      selected = "",
                      multiple = TRUE
                    )
             ),
             column(width = 6,
                    numericInput(
                      ns("heatmapWidth"),
                      label = "width of image in pixel",
                      min = 100, max = 20000, step = 10,
                      value = 800
                    ),
                    numericInput(
                      ns("heatmapHeight"),
                      label = "height of image in pixel",
                      min = 200, max = 20000, step = 10,
                      value = 300
                    ))
           ),
           fluidRow(
             column(width = 12,
                    # uiOutput(ns("additionalOptions")),
                    downloadButton(ns("download_pHeatMapUI"), "Download PlotData"),
                    actionButton(ns("save2HistHM"), "save to history")
                    
             )
           )
           
           # checkboxInput(ns("showCells"), "show cell names", FALSE),
           #
           # verbatimTextOutput(ns('cellSelection'))
           
         ) # box
    )
  )
}
