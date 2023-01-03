library(shiny)
library(parallel)
library(future)
library("BiocParallel")

plan(multisession, workers = 5)
# f <- future({
#   cat("Hello world!\n")
#   3.14
#    })
# v <- value(f)
# v
# availableCores()
# v %<-% {
#   cat("Hello world!\n")
#   3.14
# }
# v
# v
# 


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
  # this is not showing because it is run in a different process
  cat(file = stderr(), "analyze function\n")
  out <- lapply(1:10, function(x) {
    Sys.sleep(1)
    rnorm(1)
  })
  data.frame(results = unlist(out))
}

# f = future({
#   # detachedProc$process$pid = Sys.getpid()
#   analyze()
# },seed=NULL)
# 
# value(f)


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
    
    actrivateObserver <- reactiveVal(0)
    #
    # Add something to play with during waiting
    #
    output$testPlot <- renderPlot({
      cat(file = stderr(), "output$testPlot\n")
      plot(rnorm(input$inputTest))
    })
    
    #
    # Render messages
    #
    output$msg <- renderText({
      cat(file = stderr(), "output$msg\n")
      detachedProc$msg
    })
    
    #
    # Render results
    #
    output$result <- renderTable({
      cat(file = stderr(), "output$result\n")
      
      print(detachedProc$result)
      detachedProc$result
    })
    
    #
    # Start the process
    #
    observeEvent(input$start, {
      cat(file = stderr(), "input$start\n")
      
      if (!is.null(detachedProc$process))
        return()
      detachedProc$result <- dfEmpty
      #span the process/function call
      detachedProc$process <- future({
        # detachedProc$process$pid = Sys.getpid()
        analyze()
      },seed=NULL)
      actrivateObserver(1)
      detachedProc$msg <- sprintf("%1$s started", detachedProc$process$pid)
      
    })
    
    
    #
    # Stop the process
    #
    observeEvent(input$stop, {
      detachedProc$result <- dfEmpty
      if (!is.null(detachedProc$process)) {
        # need to get process id and kill using #For windows
        # system(sprintf("taskkill /F /PID %s", v[[i]]))
        #For Linux
        # system(sprintf("kill -9 %s", v[[i]]))
        
        # tools::pskill(detachedProc$process$pid)
        
        detachedProc$msg <- sprintf("%1$s killed", detachedProc$process$pid)
        detachedProc$process <- NULL
        
        # if (!is.null(detachedProc$obs)) {
        # detachedProc$obs$destroy()
        # }
      }
    })
    
    #
    # Handle process event
    #
    # observeEvent(detachedProc$process, {
      # cat(file = stderr(), "detachedProc$process\n")
      detachedProc$obs <- observe({
        cat(file = stderr(), "detachedProc$obs\n")
        # this will re-execute the collection process of mcparallel
        # if(!is.null(detachedProc$process))
        if(actrivateObserver()>0)
          invalidateLater(500, session)
        isolate({
          if(resolved(detachedProc$process))
            if(!is.null(detachedProc$process)){
              detachedProc$result <- value(detachedProc$process)
              detachedProc$process <- NULL
              actrivateObserver(0)
            }
          
        })
      })
    # })
  }
)
