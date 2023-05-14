## Only run examples in interactive R sessions
if (interactive()) {
  
  ui <- fluidPage(
    actionButton("minus", "-1"),
    actionButton("plus", "+1"),
    actionButton("nothing", "+1"),
    br(),
    textOutput("value")
  )
  
  # The comments below show the equivalent logic using reactiveValues()
  server <- function(input, output, session) {
    value <- reactiveVal(0)       # rv <- reactiveValues(value = 0)
    
    observeEvent(input$minus, {
      cat(file = stderr(), "min\n")
      newValue <- value() - 1     # newValue <- rv$value - 1
      value(newValue)             # rv$value <- newValue
    })
    
    observeEvent(input$plus, {
      cat(file = stderr(), "max\n")
      newValue <- value() + 1     # newValue <- rv$value + 1
      value(newValue)             # rv$value <- newValue
    })
    observeEvent(input$nothing, {
      cat(file = stderr(), "nothing\n")
      newValue <- value()    # newValue <- rv$value + 1
      value(newValue)             # rv$value <- newValue
    })
    
    output$value <- renderText({
      cat(file = stderr(), "render\n")
      value()                     # rv$value
    })
  }
  
  shinyApp(ui, server)
  
}

