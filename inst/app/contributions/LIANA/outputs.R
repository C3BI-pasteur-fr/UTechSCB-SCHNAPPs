require(liana)


Liana_dataInput <- callModule(
  cellSelectionModule,
  "Liana_dataInput"
)


# observer of button Color SOM ----
observe(label = "ob_LianaParameter", 
        {
          if (DEBUG) cat(file = stderr(), "ob_LianaParameter\n")
          # browser()
          input$updateLianaParameters
          setRedGreenButtonCurrent(
            vars = list(
              c("Liana_resource", input$Liana_resource),
              c("Liana_idents_col", input$Liana_idents_col),
              c("Liana_method", input$Liana_method),
              c("Liana_min_cells", input$Liana_min_cells),
              c("Liana_dataInput-Mod_PPGrp", input$'Liana_dataInput-Mod_PPGrp'),
              c("Liana_dataInput-Mod_clusterPP", input$'Liana_dataInput-Mod_clusterPP')
            )
          )
          updateButtonColor(buttonName = "updateLianaParameters", parameters = c(
            "coE_geneSOM", "Liana_idents_col", "Liana_min_cells", "Liana_method",
            "Liana_dataInput-Mod_PPGrp", "Liana_dataInput-Mod_clusterPP"
          ))
          
        })

observe(label = "somxy",{
  .schnappsEnv$defaultValues[["Liana_resource"]] = input$Liana_resource
  .schnappsEnv$defaultValues[["Liana_idents_col"]] = input$Liana_idents_col
  .schnappsEnv$defaultValues[["Liana_method"]] = input$Liana_method
  .schnappsEnv$defaultValues[["Liana_min_cells"]] = input$Liana_min_cells
  .schnappsEnv$defaultValues[["Liana_dataInput-Mod_PPGrp"]] = input$'Liana_dataInput-Mod_PPGrp'
  .schnappsEnv$defaultValues[["Liana_dataInput-Mod_clusterPP"]] = input$'Liana_dataInput-Mod_clusterPP'
})

output$Liana_dotPlot <- renderPlotly({
  liana_scEx = liana_scExReact()
  if (is.null(liana_scEx)) return(NULL)
  
  liana_scEx %>%
    liana::liana_dotplot(source_groups = unique(liana_scEx$source),
                  target_groups = unique(liana_scEx$target),
                  ntop = 20) %>% ggplotly()
})

