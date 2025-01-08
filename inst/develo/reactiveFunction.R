library(shiny)

# Wrapper function for reactive with debugging
reactive_debug <- function(expr) {
  reactive({
    start_time <- Sys.time()
    
    # Call the original reactive expression
    output <- expr()
    
    # Calculate and print execution time
    end_time <- Sys.time()
    exec_time <- end_time - start_time
    cat("Execution time:", exec_time, "seconds\n")
    
    # Return the result of the reactive expression
    return(output)
  })
}

# Define UI
ui <- fluidPage(
  titlePanel("Reactive Debug Example"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("num", "Choose a number:", 
                  min = 1, max = 1000, value = 500)
    ),
    mainPanel(
      textOutput("result"),
      verbatimTextOutput("debug_info")
    )
  )
)

# Define server logic
server <- function(input, output) {
  # Use the reactive_debug function
  exFun = function() {
    Sys.sleep(1)  # Simulate a time-consuming computation
    input$num^2
  }
  number_squared <- reactive_debug(exFun)

  output$result <- renderText({
    num_sq <- number_squared()
    paste("The square of", input$num, "is", num_sq)
  })
  
  output$debug_info <- renderPrint({
    num_sq <- number_squared()
    cat("Debug info: Input number squared =", num_sq, "\n")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)