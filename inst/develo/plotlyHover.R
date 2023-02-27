

require(shiny)
require(plotly)
require(shinyjqui)
library(ComplexHeatmap)
library("InteractiveComplexHeatmap")
#---- plotly variables ----

my_bar_color       <- '#60809f'
my_font_axis_color <- '#8c8c73'

my_light_line_color <- "#c6c6b9"
my_highlight_line_color <- "#60809f"

my_line_skinny <- .75
my_line_reg    <- 1
my_line_thick  <- 2

my_ticks <- list(
  size = 10,
  color = my_font_axis_color
)

my_title_font <- list(
  size = 12,
  color = my_font_axis_color
)

#---- ui ----

ui <- fluidPage(
  br(),
  fluidRow(
    column(
      width = 6,
      h4("Hover over a bar and it's data will be highlighted in the line graph on the right! But how does it work?"),
      tags$head(tags$style(".butt{color: black !important;}")) #  font color; otherwise the text on these buttons is gray
      
    ), actionButton(inputId = "inputHelp", label = "help", icon = icon("fas fa-question"))
  ),
  br(),
  fluidRow(
    column(
      width = 6,
      plotlyOutput("graph1")
    ),
    column(
      width = 6,
      plotlyOutput("graph2") %>% shinyjqui::jqui_resizable()
    )
  ),
  br(),
  br(),
  fluidRow(
    column(
      width = 6,
      h6("Just the ouput of the plotly hover event in case you need to see it"),
      verbatimTextOutput("hover_stuff"),
      InteractiveComplexHeatmapOutput("ht2")
    )
  )
)

#---- server ----
#
m = matrix(rnorm(100*100), nrow = 100)
ht = ComplexHeatmap::pheatmap(m)
# ht = draw(ht)

server <- function(input, output, session) {
  
  output$graph1 <- renderPlotly({
    plot_ly(
      y = c("search", "forms", "admin"),
      x = c(20, 14, 23),
      marker = list(color = my_bar_color),
      name = "Transactions",
      type = "bar",
      source = "bar_plot"
    )
    
  })
  
  
  
  
  search <- rnorm(25, mean = 1)
  forms <- rnorm(25, mean = 1)
  admin <- rnorm(25, mean = 1)
  x <- c(1:25)
  data <- data.frame(x, search, forms, admin)
  
  my_traces <- c("search", "forms", "admin")
  
  output$graph2 <- renderPlotly({
    
    plot_ly(data, x = ~x,
            y = ~search,
            name = 'search',
            type = 'scatter',
            mode = 'lines',
            line = list(color = my_light_line_color,
                        width = my_line_skinny)) %>%
      add_trace(y = ~forms,
                name = 'forms',
                mode = 'lines',
                line = list(color = my_light_line_color,
                            width = my_line_skinny)) %>%
      add_trace(y = ~admin,
                name = 'admin',
                mode = 'lines',
                line = list(color = my_light_line_color,
                            width = my_line_skinny)) %>% 
      event_register('plotly_unhover')
    
  })
  
  output$hover_stuff <- renderPrint({
    
    event_data("plotly_hover",
               source = "bar_plot")
    
    
  })
  
  
  observeEvent(event_data("plotly_hover",
                          source = "bar_plot"), {
                            
                            plotlyProxy("graph2", session) %>%
                              
                              plotlyProxyInvoke(
                                method = "restyle",
                                list(                           # this format updates every trace, why?
                                  line = list(
                                    color = my_light_line_color,
                                    width = my_line_skinny
                                  )
                                )
                              ) %>% 
                              
                              plotlyProxyInvoke(
                                method = "restyle",
                                "line",                        # this format targets one trace, why?
                                list(
                                  color = my_highlight_line_color,
                                  width = my_line_thick
                                ),
                                as.integer(match(event_data("plotly_hover", source = "bar_plot")[["y"]],
                                                 my_traces)-1)
                              )
                          })
  
  
}

#---- run the app ----

shinyApp(ui, server)
