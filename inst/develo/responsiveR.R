library(shiny)
library(parallel)
library(future)
library("BiocParallel")
library("future.callr")

plan(callr, workers = 5)
# plan(multisession, workers = 5)
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
detachedProc$PID <- NULL
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
  return(Sys.getpid())
}

getWorkerPIDs <- function(){
  v <- listenv::listenv()  # requires listenv package
  for (ii in 1:works) {
    v[[ii]] %<-% {
      Sys.getpid()
    }
  }
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
    
    activateObserver <- reactiveVal(0)
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
      pl=plan()
      
      #span the process/function call
      detachedProc$process <- future({
        # detachedProc$process$pid = Sys.getpid()
        analyze()
      },seed=NULL,
      packages = NULL,
      globals = list(analyze=analyze), # we specify all variables with the function call
      lazy = FALSE, #start immediatly
      stdout = structure(TRUE, drop = TRUE))
      activateObserver(1)
      cat(file = stderr(), paste("input$start",detachedProc$process$process$get_pid(),"me:",Sys.getpid(),"\n"))
      if("callr" %in% class(pl)){
        detachedProc$PID = detachedProc$process$process$get_pid()
      }else{
        # if("multisession" %in% class(pl)){
        #   currentWorkerPIDs = getWorkerPIDs()
        # }
        #
        # please use callr otherwise we cannot kill process (for now)
        # 
      }
      
      # browser()
      detachedProc$msg <- sprintf("%1$s started", detachedProc$process$pid)
      
    })
    
    
    #
    # Stop the process
    #
    observeEvent(input$stop, {
      # we can only kill when using future.callr
      cat(file = stderr(), "input$stop\n")
      if (!is.null(detachedProc$PID)) {
        if("running" == detachedProc$process$state){
          #For windows
          # system(sprintf("taskkill /F /PID %s", v[[i]]))
          
          #For Linux
          system(sprintf("kill -9 %s", detachedProc$PID))
          activateObserver(0)
          detachedProc$PID = NULL
          detachedProc$process = NULL
        }
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
        if(activateObserver()>0)
          invalidateLater(500, session)
        isolate({
          if(resolved(detachedProc$process))
            if(!is.null(detachedProc$process)){
              detachedProc$result <- value(detachedProc$process)
              detachedProc$process <- NULL
              activateObserver(0)
            }
          
        })
      })
    # })
  }
)
