require(liana)


Liana_dataInput <- callModule(
  cellSelectionModule,
  "Liana_dataInput"
)

callModule(
  tableSelectionServer,
  "Liana_raw_TableMod"
)
callModule(
  tableSelectionServer,
  "Liana_raw_TableMod",
  Liana_aggr_Table, caption = "Aggregate Table"
)
callModule(
  tableSelectionServer,
  "Liana_all_TableMod",
  Liana_all_Table, caption = "Aggregate Table"
)

observeEvent(
  label = "li1",
  eventExpr = liana_scExReact(),
  handlerExpr = {
    liana_scEx = liana_scExReact()
    req(Liana_all_Table)
    # browser()
    updateSelectInput(session, "Liana_method_show",
                      choices = names(liana_scEx)
    )
  })

Liana_all_Table <- reactive({
  if (DEBUG) cat(file = stderr(), "Liana_all_Table started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "Liana_all_Table")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "Liana_all_Table")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("Liana_all_Table", id = "Liana_all_Table", duration = NULL)
  }
  liana_scEx = liana_scExReact()
  mShow = input$Liana_method_show
  req(mShow)
  if (is.null(liana_scEx)) {
    if (DEBUG) if (is.null(liana_scEx)) cat(file = stderr(), "Liana_all_Table liana_scEx null.\n")
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/Liana_all_Table.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/Liana_all_Table.RData")
  # browser()
  if(! mShow %in% names(liana_scEx)) return(NULL)
  liana_scEx[[mShow]]
})


Liana_aggr_Table <- reactive({
  if (DEBUG) cat(file = stderr(), "Liana_aggr_Table started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "Liana_aggr_Table")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "Liana_aggr_Table")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("Liana_aggr_Table", id = "Liana_aggr_Table", duration = NULL)
  }
  liana_scEx = liana_aggr()
  if (is.null(liana_scEx)) {
    if (DEBUG) if (is.null(liana_scEx)) cat(file = stderr(), "Liana_aggr_Table liana_scEx null.\n")
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/Liana_aggr_Table.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/Liana_aggr_Table.RData")
  liana_scEx
})



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

output$Liana_dotPlot <- plotly::renderPlotly({
  if (DEBUG) cat(file = stderr(), "Liana_dotPlot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "Liana_dotPlot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "Liana_dotPlot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("Liana_dotPlot", id = "Liana_dotPlot", duration = NULL)
  }
  
  liana_scEx = liana_aggr()
  if (is.null(liana_scEx)) return(NULL)
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/Liana_dotPlot.Rdata", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/Liana_dotPlot.Rdata")
  if(!is(liana_scEx,"tbl_df")){
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("liana_scEx is not a tbl_df", id = "Liana_dotPlotError", duration = NULL, type = "error")
    }
    cat(file = stderr(), "ERROR: liana_scEx is not a tbl_df.\n")
    browser()
    return(NULL)
  }
  retVal <- liana_scEx %>%
    liana::liana_dotplot(source_groups = unique(liana_scEx$source),
                         target_groups = unique(liana_scEx$target),
                         ntop = 20) 
  if(is.null(retVal)){
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("Liana_dotPlotError", id = "Liana_dotPlot", duration = NULL,type = "error")
    }
    cat(file = stderr(), "ERROR in Liana_dotPlot \n")
    browser()
    return(NULL)
  }
  
  af = liana::liana_dotplot
  # remove env because it is too big
  specEnv = emptyenv()
  environment(af) = new.env(parent = specEnv)
  .schnappsEnv[["coE_dotPlot_GeneSets"]] <- list(plotFunc = af,
                                                 source_groups = unique(liana_scEx$source),
                                                 target_groups = unique(liana_scEx$target),
                                                 ntop = 20)
  printTimeEnd(start.time, "Liana_dotPlot")
  return(retVal %>% ggplotly())
  
  
})


output$Liana_Heatmap <- renderPlot({
  if (DEBUG) cat(file = stderr(), "Liana_Heatmap started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "Liana_Heatmap")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "Liana_Heatmap")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("Liana_Heatmap", id = "Liana_Heatmap", duration = NULL)
  }
  
  liana_scEx = liana_aggr()
  if (is.null(liana_scEx)) return(NULL)
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/Liana_Heatmap.Rdata", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/Liana_Heatmap.Rdata")
  if(!is(liana_scEx,"tbl_df")){
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("liana_scEx is not a tbl_df", id = "Liana_HeatmapError", duration = NULL, type = "error")
    }
    cat(file = stderr(), "ERROR: liana_scEx is not a tbl_df.\n")
    
  }
  
  liana_truncscEx <- liana_scEx %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01) # note that these pvals are already corrected
  
  # how to get this to work???
  heat_freq(liana_truncscEx) 
  
})
