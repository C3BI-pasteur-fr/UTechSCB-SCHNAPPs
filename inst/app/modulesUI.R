suppressMessages(library(magrittr))

if ("shinyjqui" %in% rownames(installed.packages())) {
  suppressMessages(require(shinyjqui))
} else {
  jqui_resizable <- function(x, ...) {
    x
  }
}
if ("shinycssloaders" %in% rownames(installed.packages())) {
  suppressMessages(library(shinycssloaders))
} else {
  withSpinner <- function(x) {
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
      column(
        width = 6,
        sc_textInput(ns("geneIds"), "comma separated list of genes for UmiCountPerGenes", value = defaultValue(ns("geneIds"), ""))
      ),
      column(
        width = 6,
        sc_textInput(ns("geneIds2"), "comma separated list of genes for UmiCountPerGenes2", value = defaultValue(ns("geneIds2"), ""))
      )
    ),
    fluidRow(
      column(
        width = 4,
        {sc_selectInput(ns("dimension_x"),
                    label = "X",
                    choices = c(defaultValue(ns("dimension_x"), "tsne1"), "tsne2", "tsne3"),
                    selected = defaultValue(ns("dimension_x"), "tsne1")
        )}
      ),
      column(
        width = 4,
        sc_selectInput(ns("dimension_y"),
                    label = "Y",
                    choices = c("tsne1", defaultValue(ns("dimension_y"), "tsne2"), "tsne3"),
                    selected = defaultValue(ns("dimension_y"), "tsne2")
        )
      ),
      column(
        width = 4,
        sc_selectInput(ns("dimension_col"),
                    label = "color",
                    choices = c(defaultValue(ns("dimension_col"), "sampleNames")),
                    selected = defaultValue(ns("dimension_col"), "sampleNames")
        )
      )
    ),
    shinydashboardPlus::box(
      title = "plot", solidHeader = TRUE, width = 12, status = "primary",
      fluidRow(
        column(
          width = 12,
          
            plotly::plotlyOutput(ns("clusterPlot")) %>% jqui_resizable()
        )
      ),
      shinydashboardPlus::box(
        title = "additional options", solidHeader = TRUE, width = 12, status = "primary",
        # helpID = "twoDselectedAddOpt",
        id = ns("clusterAddOpt"),
        dropdown_icon = NULL,
        closable = FALSE,
        enable_dropdown = T,
        dropdown_menu = actionButton(inputId = ns("twoDselectedAddOpt"), label = "", icon = icon("fas fa-question")),
        collapsible = TRUE, collapsed = TRUE,
        fluidRow(
          column(
            width = 3,
            sc_checkboxInput(ns("logX"), "log transform X", value = defaultValue(ns("logX"), FALSE))
          ),
          column(
            width = 3,
            sc_checkboxInput(ns("logY"), "log transform Y", value = defaultValue(ns("logY"), FALSE))
          ),
          column(
            width = 3,
            sc_selectInput(ns("divideXBy"),
                        label = "Divide X by",
                        # choices = c("None", "Gene.count", "UMI.count"),
                        choices = c(defaultValue(ns("divideXBy"), "None"), "UmiCountPerGenes", "UmiCountPerGenes2"),
                        selected = defaultValue(ns("divideXBy"), "None")
            )
          ),
          column(
            width = 3,
            sc_selectInput(ns("divideYBy"),
                        label = "Divide Y by",
                        # choices = c("None", "Gene.count", "UMI.count"),
                        choices = c("None", "UmiCountPerGenes", "UmiCountPerGenes2", "normByCol"),
                        selected = defaultValue(ns("divideYBy"), "None")
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            sc_checkboxInput(ns("addToGroup"), "Add to group/otherwise overwrite", defaultValue(ns("groupNames"), TRUE)),
            list(sc_textInput(ns(id = "groupName"), label = "name group, also used in Plot to color selected cells red.", value = ""),
            sc_selectInput(ns("groupNames"),
                        label = "group names, !When modifying a group this list of cells is used as a reference!",
                        choices = c(defaultValue(ns("groupNames"), "plot"), "all", "none"),
                        selected = defaultValue(ns("groupNames"), "plot")
            ),
            verbatimTextOutput(ns("nCellsVisibleSelected")),
            actionButton(ns("changeGroups"), "change current selection")) %>% setId("groupTutorial"),
            sc_checkboxInput(ns("showCells"), "show cell names", defaultValue(ns("showCells"), FALSE)),
            verbatimTextOutput(ns("cellSelection")),
            actionButton(ns("save2Hist"), "save to history"),
            uiOutput(ns("additionalOptions")) # TODO:is this still needed???
          )
        )
      )
    )
    # ,
    # fluidRow(column(width = 6,
    #                 sc_checkboxInput(ns("moreOptions"), "show more options", FALSE),
    # )
    # )
  )
}

# tableSelectionUi --------------
tableSelectionUi <- function(id) {
  ns <- NS(id)
  tagList(
   # jqui_resizable(
     shinydashboardPlus::box(
      width = 12, height = 700,
      fluidRow(
        column(
          width = 12,
          div(
            h5("Selected itmes to be copied"),
            align = "left"
          ),
          verbatimTextOutput(ns("rowSelection"))
        )
      ),
      fluidRow(
        column(
          width = 6,
          downloadButton(ns("download_cellNameTable"), "Download table")
        ),
        column(
          width = 6,
          actionButton(ns("save2HistTabUi"), "Save to history"),
        )
      ),
      fluidRow(h4("Cells", offset = 1)),
      fluidRow(
        column(
          width = 3,
          sc_checkboxInput(ns("selectAll"), "Select all rows", FALSE)
        ),
        column(
          width = 3,
          sc_checkboxInput(ns("reorderCells"), "reorder cells by sum of selected genes", FALSE)
        ),
        column(
          width = 3,
          sc_checkboxInput(ns("showAllCells"), "show all columns", defaultValue(ns("showAllCells"), FALSE))
        ),
        column(
          width = 3,
          sc_checkboxInput(ns("showRowNames"),"show row names", defaultValue(ns("showRowNames"), FALSE))
        )
      ),
      br(),
      fluidRow(column(
        width = 12,
        DT::DTOutput(ns("cellNameTable")) ,
        style = "overflow-y: scroll;overflow-x: scroll;"
      )
      )
    )
    # )
     )
}

# pHeatMapUI --------------
pHeatMapUI <- function(id) {
  ns <- NS(id)
  tagList(
    shinydashboardPlus::box(
      width = 12,
      fluidRow(
        column(
          width = 12,
          # jqui_resizable(plotOutput(ns("pHeatMapPlot"),
          #                           # height = "auto",
          #                           brush = brushOpts(id = "crh1")
          # ), options = list(width = "99%"))
          # originalHeatmapOutput(heatmap_id = ns("pHeatMapPlot"), title = NULL) %>% jqui_resizable()
          plotOutput(ns("pHeatMapPlot"), click = ns("heatmap_click"), 
                     brush = ns("heatmap_brush")) %>% jqui_resizable()
        )),
      fluidRow(
        column(width = 12,
               verbatimTextOutput(ns("pHeatMapPlotSelection"))
               )
      ),
      shinydashboardPlus::box(
        title = "additional options", solidHeader = TRUE, width = 12, status = "primary",
        collapsible = TRUE, collapsed = TRUE,
        id=ns("heatmapAddOpt"),
        fluidRow(
          column(
            width = 12,
            # sc_checkboxInput(ns("moreOptions"), "show more options", FALSE),
            # sc_checkboxInput(ns("showColTree"), label = "Show tree for cells", value = defaultValue(ns("showColTree"),FALSE)),
            sc_selectInput(
              ns("sortingCols"),
              label = "order columns by",
              choices = c("dendrogram", "list", "gene (click)"),
              selected = defaultValue(ns("sortingCols"), "list"),
              multiple = FALSE
            ),
            sc_selectInput(
              ns("sortingRows"),
              label = "order rows by",
              choices = c("dendrogram", "list"),
              selected = defaultValue(ns("sortingRows"), "dendrogram"),
              multiple = FALSE
            )
          )
        ),
        fluidRow(
          column(
            width = 6,
            sc_selectInput(ns("normRow"),
                        label = "scale by row (for color)",
                        choices = c("row", "column", "row_order", "col_order", defaultValue(ns("normRow"), "none")),
                        selected = defaultValue(ns("normRow"), "none")
            ),
            sc_selectInput(
              ns("ColNames"),
              label = "group names",
              choices = defaultValue(ns("ColNames"), "sampleNames"),
              selected = defaultValue(ns("ColNames"), c("dbCluster", "sampleNames")),
              multiple = TRUE
            ),
            
            sc_selectInput(
              ns("orderNames"),
              label = "order of columns",
              choices = defaultValue(ns("orderNames"), "sampleNames"),
              selected = defaultValue(ns("orderNames"), c("dbCluster", "sampleNames")),
              multiple = TRUE
            ),
            sliderInput(
              ns("heatmapMinMaxValue"),
              label = "min/max value for heatmap",
              min = -10000000,
              max = 100000000,
              value = defaultValue(ns("heatmapMinMaxValue"), c(-10000000,100000000))
            )
          ),
          column(
            width = 6,
            # sc_numericInput(
            #   ns("heatmapWidth"),
            #   label = "width of image in pixel",
            #   min = 100, max = 20000, step = 10,
            #   value = defaultValue(ns("heatmapWidth"),800)
            # ),
            # sc_numericInput(
            #   ns("heatmapHeight"),
            #   label = "height of image in pixel",
            #   min = 200, max = 20000, step = 10,
            #   value = defaultValue(ns("heatmapHeight"), 300)
            # ),
            sc_textInput(ns("heatMapGrpName"),"name selected cells"),
            actionButton(ns("heatMapGrpNameButton"),"create projection"),
            br(),
            sc_selectInput(
              ns("colPal"),
              label = "color palette to choose from",
              choices = c("none", "Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", 
                          "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", 
                          "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
              selected = defaultValue(ns("colPal"), "none"),
              multiple = FALSE
            ),
            # uiOutput(ns("addColours")),
            # actionButton(ns("heatmapTest"),"test"),
            # verbatimTextOutput(ns("heatmapTestout")),
            sliderInput(
              ns("heatmapCellGrp"),
              label = "number of cell groups",
              min = 2,
              max = 1000,
              value = defaultValue(ns("heatmapCellGrp"), 2)
            )
          )
        ),
        
        fluidRow(
          column(
            width = 12,
            div(
              h5("Selected itmes to be copied"),
              align = "left"
            ),
            verbatimTextOutput(ns("heatmapSelectedGenes")))),
        fluidRow(
          column(
            width = 12,
            # uiOutput(ns("additionalOptions")),
            downloadButton(ns("download_pHeatMapUI"), "Download PlotData"),
            actionButton(ns("save2HistHM"), "save to history")
          )
        )
        
        # sc_checkboxInput(ns("showCells"), "show cell names", FALSE),
        #
        # verbatimTextOutput(ns('cellSelection'))
      ) # box
    )
  )
}


# cellSelectionUI --------------
cellSelectionUI <- function(id) {
  ns <- NS(id)
  shinydashboardPlus::box(
    width = 6, title = "Selection of input cells",
    fluidRow(
      column(
        width = 6,
        # uiOutput("DE_clusterSelectionPanelPlot")
        sc_selectInput(
          inputId = ns("Mod_clusterPP"), label = "Clusters/Factor to use",
          choices = c(defaultValue(ns("Mod_clusterPP"), "dbCluster"), "sampleNames"),
          selected = defaultValue(ns("Mod_clusterPP"), "dbCluster")
        )
      ),
      column(
        width = 6,
        sc_selectInput(
          inputId = ns("Mod_PPGrp"), label = "Values to use",
          choices = defaultValue(ns("Mod_PPGrp"), "1"), selected = defaultValue(ns("Mod_PPGrp"), "1"), multiple = TRUE
        )
      )
    )
  )
}
