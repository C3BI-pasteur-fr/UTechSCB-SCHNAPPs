# local helper function
add2workflowObsList <- function(li, wkfl, workflowObsList){
  if(!wkfl %in% names(workflowObsList)) workflowObsList[[wkfl]] = list()
  workflowObsList[[wkfl]][[li$id]] =   list(id = li$id, fun=li$Func)
  return(workflowObsList)
}

# shortCutsTab ----
workflowObsList = list()

# wkfl1 ----
wkfl = "wkfl1"
workflowObsList[[wkfl]] = list()

## LoadData ----
workflowObsList = list(
  id = "LoadData",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "input"
    )
    updateRadioButtons(session = session, "whichscLog", selected = "disablescEx_log")
    updateBox("addOptInput", action = "toggle")
    
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

# QC
## gQC_sampleHist ----
workflowObsList = list(
  id = "gQC_sampleHist",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "gQC_sampleHist"
    )
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## gQC_umiHist ----
workflowObsList = list(
  id = "gQC_umiHist",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "gQC_umiHist"
    )
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## nFeatureViolin ----
workflowObsList = list(
  id = "nFeatureViolin",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpression"
    )
    # browser()
    # js$showmenuItem("coexpressionSelected")
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "Feature.count"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-divideYBy",
                      selected = "None"
    )
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## nFeatureSelection ----
workflowObsList = list(
  id = "nFeatureSelection",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "barcode"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "Feature.count"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-divideYBy",
                      selected = "None"
    )
    updateBox("coE_selected-clusterAddOpt", action = "toggle")
    updateCheckboxInput(session = session, "coE_selected-showCells",value = TRUE)
    
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## nCountViolin ----
workflowObsList = list(
  id = "nCountViolin",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "UMI.count"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-divideYBy",
                      selected = "None"
    )
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## nCountSelection ----
workflowObsList = list(
  id = "nCountSelection",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "barcode"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "UMI.count"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-divideYBy",
                      selected = "None"
    )
    updateCheckboxInput(session = session, "coE_selected-showCells",value = TRUE)
    
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## npMTViolin ----
workflowObsList = list(
  id = "npMTViolin",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "percent.mt"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-divideYBy",
                      selected = "None"
    )
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## npMTSelection ----
workflowObsList = list(
  id = "npMTSelection",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "barcode"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "percent.mt"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-divideYBy",
                      selected = "None"
    )
    updateCheckboxInput(session = session, "coE_selected-showCells",value = TRUE)
    
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## countFeature ----
workflowObsList = list(
  id = "countFeature",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "nCount_RNA"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "Feature.count"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-divideYBy",
                      selected = "None"
    )
    updateCheckboxInput(session = session, "coE_selected-showCells",value = TRUE)
    
  }
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## wkfl1.countMt.click ----
workflowObsList = list(
  id = "countMt",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "nCount_RNA"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "percent.mt"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-divideYBy",
                      selected = "None"
    )
    updateCheckboxInput(session = session, "coE_selected-showCells",value = TRUE)
    
  }
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## go2CellSelection ----
workflowObsList = list(
  id = "go2CellSelection",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "cellSelection"
    )
    updateBox("cellSelectionParameters", action = "toggle")
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## combineVars1 ----
workflowObsList = list(
  id = "combineVars1",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "modifyProj"
    )
    updateTabsetPanel(
      session = session,
      "modProj",
      selected = "gQC_combProj"
    )
    # "combine.Proj.Tab"
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## renameLevels1 ----
workflowObsList = list(
  id = "renameLevels1",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "modifyProj"
    )
    updateTabsetPanel(
      session = session,
      "modProj",
      selected = "gQC_renameLev"
    )
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)


# Normalization
## selectNormParameters ----
workflowObsList = list(
  id = "selectNormParameters",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "clusterParameters"
    )
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)



## selectPCAParameters ----
workflowObsList = list(
  id = "setPCAparameters",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "clusterParameters"
    )
    updateNumericInput(
      session = session,
      inputId = "pcaRank",
      value = 50
    )
    updateNumericInput(
      session = session,
      inputId = "pcaN",
      value = 2000
    )
    updateSelectizeInput(
      session = session,
      inputId = "hvgSelection", 
      selected = "getTopHVGs"
    )
    updateCheckboxInput(
      session = session,
      inputId = "pcaCenter",
      value = TRUE
    )
    updateCheckboxInput(
      session = session,
      inputId = "pcaScale",
      value = FALSE
    )
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)


## selectPCAParameters ----
workflowObsList = list(
  id = "setClusterParameters",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "clusterParameters"
    )
    updateNumericInput(
      session = session,
      inputId = "seurClustDims", 
      value = "50"
    )
    updateNumericInput(
      session = session,
      inputId = "seurClustk.param", 
      value = "20"
    )
  }
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## geneSelection ----
workflowObsList = list(
  id = "geneSelection",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "geneSelection"
    )
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## input2 ----
workflowObsList = list(
  id = "input2",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "input"
    )
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## umap ----
workflowObsList = list(
  id = "umap",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "gQC_umapPlot"
    )
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## setUMAPparameters ----

workflowObsList = list(
  id = "setUMAPparameters",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "gQC_umapPlot"
    )
    updateBox("addUMAPoptions", action = "toggle")
    updateNumericInput(
      session = session,
      inputId = "gQC_um_n_neighbors", 
      value = "20"
    )
    updateSelectizeInput(
      session = session,
      inputId = "gQC_um_spread", 
      selected = "7"
    )
    updateNumericInput(
      session = session,
      inputId = "gQC_um_metric", 
      value = "cosine"
    )
  }
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## viewUMAP2D ----
workflowObsList = list(
  id = "viewUMAP2D",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "gQC_umapPlot"
    )
    updateSelectInput(session, "gQC_umap_main-dimension_x",
                      selected = "UMAP1"
    )
    updateSelectInput(session, "gQC_umap_main-dimension_y",
                      selected = "UMAP2"
    )
    updateSelectInput(session, "gQC_umap_main-dimension_col",
                      selected = "dbCluster"
    )
  }
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

## viewUMAP3D ---- 
workflowObsList = list(
  id = "viewUMAP3D",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "gQC_tsnePlot"
    )
    updateSelectizeInput(session, "gQC_dim3D_x",
                      selected = "UMAP1"
    )
    updateSelectizeInput(session, "gQC_dim3D_y",
                      selected = "UMAP2"
    )
    updateSelectizeInput(session, "gQC_dim3D_z",
                      selected = "UMAP3"
    )
    updateSelectizeInput(session, "gQC_col3D",
                      selected = "dbCluster"
    )
  }
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)



## renameCluster ----
workflowObsList = list(
  id = "renameCluster",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "Projections"
    )
    updateTabsetPanel(
      session = session,
      "modProj",
      selected = "renameProj"
    )
    updateSelectizeInput(
      session = session,
      inputId = "oldPrj", 
      selected = "dbCluster"
    )
    updateTextInput(
      session = session,
      inputId = "newPrj",
      value = "dbClusterLog"
    )
  }
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)



# tSNE
# gQC_tsnePerplexity = 30
# gQC_tsneSeed = 2104
# 
#  resolution=0.5
# pcaRank=30
# seurClustDims = 30
# gQC_um_n_neighbors = 30


# div(h5("Then combine these two new projections using a temp variable.", a(id="combineVars1","combine projections."))),
# div(h5("Then rename the levels of this new temp variable.", a(id="renameLevels1","rename levels."))),




# UMI.count

# gQC_variancePC
# updateTabItems(
#   session = session,
#   "sideBarID",
#   selected = "coexpressionSelected"
# )
# updateSelectInput(session, "coE_selected-dimension_x",
#                   selected = "tsne1"
# )
# updateSelectInput(session, "coE_selected-dimension_y",
#                   selected = "tsne2"
# )
# updateSelectInput(session, "coE_selected-dimension_col",
#                   selected = "cellDensity"
# )
# updateSelectInput(session, "coE_selected-divideYBy",
#                   selected = "None"
# )

##### generate Observers =====

# for loops are not working with observes/onclick/layzy stuff
lapply(names(workflowObsList),FUN = function(nam){
  lapply(names(workflowObsList[[nam]]),FUN = function(namItem){
    fun = workflowObsList[[nam]][[namItem]][["fun"]]
    divFunc = function(x, session, ...){
      x(session)
    }
    cat(file = stderr(), paste0("=+++++++======== ", nam,".", make.names(namItem),".click\n"))
    shinyjs::onclick(paste0(nam,".", make.names(namItem),".click"),divFunc(x=fun, session=session))
  })
})

