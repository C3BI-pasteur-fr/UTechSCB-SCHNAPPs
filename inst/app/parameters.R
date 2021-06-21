
if (DEBUG) cat(file = stderr(), paste("parameters:", length(allTabs), " ", "\n"))


# here we add content to the page on the rigth (main visualization window)
allTabs[[length(allTabs) + 1]] <- list(
  shinydashboard::tabItem(
    "genParams",
    
    # shinydashboard::box(
    #   title = "Sample of normalized values", solidHeader = TRUE, width = 12, status = "primary",
    #   collapsible = FALSE, collapsed = TRUE,
    #   
    # )
    shinydashboard::box(
      title = "Colors", solidHeader = TRUE, width = 12, status = "primary",
      collapsible = TRUE, collapsed = TRUE,
      fluidRow(column(
        width = 12, offset = 1,
        actionButton("updateColors", "apply changes", width = "80%")
      )),
      br(),
      # tabBox(title = "modify colors", width = 12, id = "modCols",
      uiOutput("ColorSelection")
      # tabPanel(
      #   title = "Sample", solidHeader = TRUE, width = 12, value = "SampleColorPanel",
      #   fluidRow(
      #     column(
      #       width = 6,
      #       uiOutput("sampleColorSelection")
      #     ))),
      # tabPanel(
      #   title = "Cluster", solidHeader = TRUE, width = 12, value = "ClusterColorPanel",
      #   fluidRow(
      #     column(
      #       width = 6,
      #       uiOutput("clusterColorSelection")
      #     )))
      # )
    ),
    checkbsTT(item = "updateColors"),
    checkbsTT(item = "sampleColorSelection"),
    checkbsTT(item = "clusterColorSelection")
  )
)
if (DEBUG) {
  cat(file = stderr(), paste("end: parameters.R\n"))
}
