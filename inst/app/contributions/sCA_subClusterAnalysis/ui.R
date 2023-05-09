suppressMessages(library(magrittr))
source(paste0(packagePath, "/modulesUI.R"), local = TRUE)

menuList <- list(
  shinydashboard::menuItem("Subcluster analysis",
                           icon = icon("biohazard"),
                           # id = 'subclusterID',
                           tabName = "subcluster", startExpanded = FALSE,
                           shinydashboard::menuSubItem("DGE analysis", tabName = "sCA_dge")
  )
)
# sCA_dge ----
tabList <- list(
  dgeTab = shinydashboard::tabItem(
    "sCA_dge",
    tabBox(
      title = "Differential expression", width = 12, 
      tabPanel(title = "Select data", 
               width = 12, id = "modProjDGE",
               footer = {
                 tags$ul(
                   tags$li(
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
                   
                 )},
               cellSelectionUI("sCA_dataInput"),
               shinydashboardPlus::box(width = 6,
                                   title = "Dimensions for plot",
                                   fluidRow(
                                     column(width = 6,
                                            sc_selectInput(
                                              "sCA_subscluster_x1",
                                              label = "X",
                                              choices = c(defaultValue("sCA_subscluster_x1", "tsne1"), "tsne2", "tsne3"),
                                              selected = defaultValue("sCA_subscluster_x1", "tsne1")
                                            )
                                     ),
                                     column(width = 6,
                                            sc_selectInput(
                                              "sCA_subscluster_y1",
                                              label = "Y",
                                              choices = c("tsne1", defaultValue("sCA_subscluster_y1", "tsne2"), "tsne3"),
                                              selected = defaultValue("sCA_subscluster_y1", "tsne2")
                                            )
                                     )
                                   )
               ),
               fluidRow(
                 column(
                   width = 10, offset = 1,
                   sc_selectizeInput(
                     inputId = "sCA_dgeRadioButton",
                     label = "Method to use",
                     choices = defaultValue("sCA_dgeRadioButton", "none"),
                     selected = defaultValue("sCA_dgeRadioButton", "none"),
                     width = "100%"
                   )
                 )
               ),
               fluidRow(
                 column(width = 12, offset = 1,
                        actionButton("updateDGEParameters", "apply changes", width = "80%")
                 )
               ),
               
               
               br(),
               shinydashboardPlus::box(width = 12,
                                   fluidRow(
                                     
                                     column(width = 6,
                                            # plotly::plotlyOutput("sCA_dge_plot1")
                                            plotOutput("sCA_dge_plot1", brush = brushOpts(
                                              id = "db1"
                                            )) %>% withSpinner()
                                     )
                                     ,
                                     column(width = 6,
                                            # plotly::plotlyOutput("sCA_dge_plot2")
                                            plotOutput("sCA_dge_plot2", brush = brushOpts(
                                              id = "db2"
                                            )) %>% withSpinner()
                                     )
                                   )
               ),
               br(),
               shinydashboardPlus::box(
                 title = "additional options", solidHeader = TRUE, width = 12, status = "primary",
                 dropdown_icon = NULL,
                 closable = FALSE,
                 enable_dropdown = T,
                 collapsible = TRUE, collapsed = TRUE,
                 id = "dgeAddOptions",
                 shinydashboardPlus::box(title = "scDEA parameters", width = 12,solidHeader = TRUE, 
                     fluidRow(
                       column(width = 3,
                              sc_checkboxInput("scDEA_parallel", "use parallel implementations", value = defaultValue("scDEA_parallel", FALSE)),
                       )
                     ),
                     fluidRow(
                       column(width = 3,
                              sc_checkboxInput("scDEA_BPSC", "BPSC", value = defaultValue("scDEA_BPSC", TRUE)),
                              sc_checkboxInput("scDEA_DEsingle", "DEsingle", value = defaultValue("scDEA_DEsingle", TRUE)),
                              sc_checkboxInput("scDEA_DESeq2", "DESeq2", value = defaultValue("scDEA_DESeq2", TRUE)),
                              sc_checkboxInput("scDEA_edgeR", "edgeR", value = defaultValue("scDEA_edgeR", TRUE))),
                       column(width = 3,
                              sc_checkboxInput("scDEA_MAST", "MAST", value = defaultValue("scDEA_MAST", TRUE)),
                              sc_checkboxInput("scDEA_monocle", "monocle", value = defaultValue("scDEA_monocle", TRUE)),
                              sc_checkboxInput("scDEA_scDD", "scDD", value = defaultValue("scDEA_scDD", TRUE)),
                              sc_checkboxInput("scDEA_Ttest", "Ttest", value = defaultValue("scDEA_Ttest", TRUE))),
                       column(width=3,
                              sc_checkboxInput("scDEA_Wilcoxon", "Wilcoxon", value = defaultValue("scDEA_Wilcoxon", TRUE)),
                              sc_checkboxInput("scDEA_limma", "limma", value = defaultValue("scDEA_limma", TRUE)),
                              sc_checkboxInput("scDEA_Seurat", "Seurat", value = defaultValue("scDEA_Seurat", TRUE)),
                              sc_checkboxInput("scDEA_zingeR.edgeR", "zingeR.edgeR", value = defaultValue("scDEA_zingeR.edgeR", TRUE))
                       )
                     )
                     
                 )
               )
      ),

      if ("manhattanly" %in% rownames(installed.packages()))
        tabPanel(
          title = "Volcano plot",  width = 12, 
          collapsible = FALSE, collapsed = FALSE,
          fluidRow(
            column(
              width = 6,
              sc_numericInput(
                inputId = "sCA_volc_effectLimit", label = "x-axis threshold", 
                value = defaultValue("sCA_volc_effectLimit", 1), min = 0.0, max = 10000,
                step = 0.1
              )
            ),
            column(
              width = 6,
              sc_numericInput(
                inputId = "sCA_volc_pval", label = "y-axis threshold", 
                value = defaultValue("sCA_volc_pval", 5), min = 0.0, max = 110000,
                step = 0.1
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              verbatimTextOutput("sCA_volc_selected")
            )
          ),
          fluidRow(
            column(
              width = 12,
              # https://stackoverflow.com/questions/44412382/clicking-same-plotly-marker-twice-does-not-trigger-events-twice
              # useShinyjs(),
              # code to reset plotlys event_data() to NULL -> executed upon action button click
              # note that "A" needs to be replaced with plotly source string if used
              # extendShinyjs(text = "shinyjs.sCA_volcanoPlot_resetClick = function() { Shiny.onInputChange('plotly_selected-A', 'null'); }", functions = "sCA_volcanoPlot_resetClick"),
              
              if ("manhattanly" %in% rownames(installed.packages()))
                plotly::plotlyOutput("sCA_volcanoPlot") %>% jqui_resizable()
            )
          ),
          br(),
          actionButton("save2HistVolc", "save to history")
        ),
      tabPanel(
        title = "Differentially Expressed Genes",  width = 12, 
        collapsible = FALSE, collapsed = TRUE,
        fluidRow(
          column(
            12,
            tableSelectionUi("sCA_dgeTable")
          )
        )
      )
    )
  )
)
