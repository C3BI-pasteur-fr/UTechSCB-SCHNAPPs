

liana_scExReact <- reactive({
  if (DEBUG) cat(file = stderr(), "liana_scExReact started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "liana_scExReact")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "liana_scExReact")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("liana_scExReact", id = "liana_scExReact", duration = NULL)
  }
  if(input$updateLianaParameters < 1){return(NULL)}
  scEx = isolate(scEx())
  scEx_log = isolate(scEx_log())
  proj = isolate(projections())
  resource = isolate(input$Liana_resource)
  method = isolate(input$Liana_method)
  idents_col = isolate(input$Liana_idents_col)
  min_cells = isolate(input$Liana_min_cells)
  selectedCells = isolate(Liana_dataInput())
  cellNs <- isolate(selectedCells$cellNames())
  
  if(is.null(scEx)) return(NULL)
  if(is.null(scEx_log)) return(NULL)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/liana_scExReact.Rdata", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/liana_scExReact.Rdata")
  
  # browser()
  proj = proj[cellNs,]
  proj[,idents_col] = droplevels(proj[,idents_col])
  liana_scEx <- liana_scExFunc(scEx=scEx[, cellNs], scEx_log = scEx_log[, cellNs], proj, resource, idents_col, method, min_cells)
})

liana_aggr <- reactive({
  liana_scEx = liana_scExReact()
  method = isolate(input$Liana_method)
  req(method, liana_scEx)
  
  if(length(method) > 1){
    liana_scEx <- liana_scEx %>% liana_aggregate()
  }
  return(liana_scEx)
})

liana_scExFunc <- function(scEx, scEx_log, proj, resource, idents_col, method, min_cells){
  require(liana)
  # browser()
  # scEx %>% dplyr::glimpse()
  rownames(scEx) = toupper(rownames(scEx))
  rownames(scEx_log) = toupper(rownames(scEx_log))
  colData(scEx) = as(proj, "DFrame")
  # assays(scEx)
  assays(scEx)[["logcounts"]] = as(assays(scEx_log)[["logcounts"]],"CsparseMatrix")
  
  # what if only one method is used? do we still need to aggregate?
  liana_scEx <- liana_wrap(scEx, idents_col = idents_col, assay="logcounts",
                           base = 2 , # log expression base
                           method = method,
                           resource = resource)
  
  return(liana_scEx)
  
}


