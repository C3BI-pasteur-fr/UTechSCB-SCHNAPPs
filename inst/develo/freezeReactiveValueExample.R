
ui <- fluidPage(
  selectInput("data", "Data Set", c("mtcars", "pressure")),
  checkboxGroupInput("cols", "Columns (select 2)", character(0)),
  plotOutput("plot")
)

server <- function(input, output, session) {
  observe({
    data <- get(input$data)
    # Sets a flag on input$cols to essentially do req(FALSE) if input$cols
    # is accessed. Without this, an error will momentarily show whenever a
    # new data set is selected.
    freezeReactiveValue(input, "cols")
    updateCheckboxGroupInput(session, "cols", choices = names(data))
  })
  
  output$plot <- renderPlot({
    # When a new data set is selected, input$cols will have been invalidated
    # above, and this will essentially do the same as req(FALSE), causing
    # this observer to stop and raise a silent exception.
    cols <- input$cols
    data <- get(input$data)
    browser()
    if (length(cols) == 2) {
      plot(data[[ cols[1] ]], data[[ cols[2] ]])
    }
  })
}

shinyApp(ui, server)
