

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
  input$updateLianaParameters
  scEx = scEx()
  scEx_log = scEx_log()
  proj = projections()
  if(is.null(scEx)) return(NULL)
  if(is.null(scEx_log)) return(NULL)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/liana_scExReact.Rdata", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/liana_scExReact.Rdata")
  
  browser()
  liana_scEx <- liana_scExFunc(scEx, scEx_log, proj)
 })

liana_scExFunc <- function(scEx, scEx_log, proj){
  # scEx %>% dplyr::glimpse()
  rownames(scEx) = toupper(rownames(scEx))
  rownames(scEx_log) = toupper(rownames(scEx_log))
  colData(scEx) = as(proj, "DFrame")
  # assays(scEx)
  assays(scEx)[["logcounts"]] = as(assays(scEx_log)[["logcounts"]],"CsparseMatrix")
  
  liana_scEx <- liana_wrap(scEx, idents_col = "dbCluster", assay="logcounts",
                           base = 2 , # log expression base
                           resource = c( "OmniPath"  ))
  return(liana_scEx)
  
}


