library(shiny)
library(parallel)

#
# reactive variables
# 
detachedProc <- reactiveValues()
detachedProc$process <- NULL
detachedProc$msg <- NULL
detachedProc$obs <- NULL
# counter <- 0
# results <- list()
dfEmpty <- data.frame(results = numeric(0))


#
# Long computation
#
analyze <- function() {
    out <- lapply(1:50, function(x) {
        Sys.sleep(1)
        rnorm(1)
    })
    data.frame(results = unlist(out))
}

#
# Shiny app
#
shinyApp(
    ui = fluidPage(
        column(6,
               wellPanel(
                   tags$label("Press start and wait 5 seconds for the process to finish"),
                   actionButton("start", "Start", class = "btn-primary"),
                   actionButton("stop", "Stop", class = "btn-danger"),
                   textOutput('msg'),
                   tableOutput('result')
               )
        ),
        column(6,
               wellPanel(
                   sliderInput(
                       "inputTest",
                       "Shiny is responsive during computation",
                       min = 10,
                       max = 100,
                       value = 40
                   ),
                   plotOutput("testPlot")
               ))),
    server = function(input, output, session)
    {
        #
        # Add something to play with during waiting
        #
        output$testPlot <- renderPlot({
            plot(rnorm(input$inputTest))
        })
        
        #
        # Render messages
        #
        output$msg <- renderText({
            detachedProc$msg
        })
        
        #
        # Render results
        #
        output$result <- renderTable({
            print(detachedProc$result)
            detachedProc$result
        })
        
        #
        # Start the process
        #
        observeEvent(input$start, {
            if (!is.null(detachedProc$process))
                return()
            detachedProc$result <- dfEmpty
            detachedProc$process <- mcparallel({
                analyze()
            })
            
            detachedProc$msg <- sprintf("%1$s started", detachedProc$process$pid)
            
        })
        
        
        #
        # Stop the process
        #
        observeEvent(input$stop, {
            detachedProc$result <- dfEmpty
            if (!is.null(detachedProc$process)) {
                tools::pskill(detachedProc$process$pid)
                detachedProc$msg <- sprintf("%1$s killed", detachedProc$process$pid)
                detachedProc$process <- NULL
                
                if (!is.null(detachedProc$obs)) {
                    detachedProc$obs$destroy()
                }
            }
        })
        
        #
        # Handle process event
        #
        observeEvent(detachedProc$process, {
            detachedProc$obs <- observe({
                invalidateLater(500, session)
                isolate({
                    result <- mccollect(detachedProc$process, wait = FALSE)
                    if (!is.null(result)) {
                        detachedProc$result <- result
                        detachedProc$obs$destroy()
                        detachedProc$process <- NULL
                    }
                })
            })
        })
    }
)
