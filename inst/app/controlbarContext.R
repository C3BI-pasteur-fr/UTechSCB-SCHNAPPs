controlbarContext <- shinydashboardPlus::dashboardControlbar(
  width = 500, collapsed = F,
  id = "controlbar",
  controlbarMenu(
    id = "menu",
    controlbarItem(
      "Workflow 1",
      shinydashboard::box(
        title = "Loading data",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 div(h4(a(id = "wkfl1.LoadData.click", "Load data"))),
                 div(h5("Load the data count matrix, annotations, and GMT information without activating the Log-counts. ")),
                 div(h5("The count matrix can be in text form (csv/tsv), an RData file with a singlecellExpression object.")),
                 div(h5("Since the first steps of QC will use only the raw counts the potentially expensive normalization procedure is not needed.")),
                 div(h5("Optionally, the data can be subsampled.")),
                 div(h5("The \"regular expression to count genes/cell\" defines genes for which a collective count should be retained, 
                 even if those genes are removed during the filtering process. This way, one can keep track of ribosomal counts. 
                 The regular expression '^MT-' searches of all genes whose gene name starts (^) with 'MT-'")),
                 div(h5("")),
                 div(h5("")),
                 div(h5("")),
                 div(h5(""))
          )),
        br()),
      shinydashboard::box(
        title = "QC",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 div(h4("Counts"), align = "center"),
                 div(a(id="wkfl1.Gene.selection.click", "check feature counts"))
          )),
        br()),
      shinydashboard::box(
        title = "QC",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=10, offset = 1,
                 div(h4("Counts"), align = "center"),
                 div(a(id="wkfl1.Gene.selection.click", "check feature counts"))
          )),
        br())
    )
    
  )
)
