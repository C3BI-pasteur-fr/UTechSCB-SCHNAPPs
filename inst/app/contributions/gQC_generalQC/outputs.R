# source("moduleServer.R", local = TRUE)
# source("reactives.R", local = TRUE)
# TODO: verify that this anything and then integrate in DUMMY
myZippedReportFiles <- c("gqcProjections.csv")

require(stringr)

.schnappsEnv$gQC_X1 <- "tsne1"
.schnappsEnv$gQC_X2 <- "tsne2"
.schnappsEnv$gQC_X3 <- "tsne3"
.schnappsEnv$gQC_col <- "sampleNames"
observe(label = "ob2", {
  if (DEBUG) cat(file = stderr(), "observe: gQC_dim3D_x\n")
  .schnappsEnv$gQC_X1 <- input$gQC_dim3D_x
})
observe(label = "ob3", {
  if (DEBUG) cat(file = stderr(), "observe: gQC_dim3D_y\n")
  .schnappsEnv$gQC_X2 <- input$gQC_dim3D_y
})
observe(label = "ob4", {
  if (DEBUG) cat(file = stderr(), "observe: gQC_dim3D_z\n")
  .schnappsEnv$gQC_X3 <- input$gQC_dim3D_z
})
observe(label = "ob5", {
  if (DEBUG) cat(file = stderr(), "observe: gQC_col3D\n")
  .schnappsEnv$gQC_col <- input$gQC_col3D
})

observe(label = "obs_gQC_tsnePerplexity", x = {
  .schnappsEnv$defaultValues[["gQC_tsnePerplexity"]] = input$gQC_tsnePerplexity
})
observe(label = "obs_gQC_tsneSeed", x = {
  .schnappsEnv$defaultValues[["gQC_tsneSeed"]] = input$gQC_tsneSeed
})
observe(label = "obs_gQC_tsneTheta", x = {
  .schnappsEnv$defaultValues[["gQC_tsneTheta"]] = input$gQC_tsneTheta
})
observe(label = "obs_gQC_tsneDim", x = {
  .schnappsEnv$defaultValues[["gQC_tsneDim"]] = input$gQC_tsneDim
})

observe(label = "obs_gQC_dim3D_x", x = {
  .schnappsEnv$defaultValues[["gQC_dim3D_x"]] = input$gQC_dim3D_x
})
observe(label = "obs_gQC_dim3D_y", x = {
  .schnappsEnv$defaultValues[["gQC_dim3D_y"]] = input$gQC_dim3D_y
})
observe(label = "obs_gQC_dim3D_z", x = {
  .schnappsEnv$defaultValues[["gQC_dim3D_z"]] = input$gQC_dim3D_z
})
observe(label = "obs_gQC_col3D", x = {
  .schnappsEnv$defaultValues[["gQC_col3D"]] = input$gQC_col3D
})
observe(label = "obs_gQC_um_spread", x = {
  .schnappsEnv$defaultValues[["gQC_um_spread"]] = input$gQC_um_spread
})
observe(label = "obs_gQC_um_n_components", x = {
  .schnappsEnv$defaultValues[["gQC_um_n_components"]] = input$gQC_um_n_components
})
observe(label = "obs_gQC_um_n_neighbors", x = {
  .schnappsEnv$defaultValues[["gQC_um_n_neighbors"]] = input$gQC_um_n_neighbors
})
observe(label = "obs_gQC_um_init", x = {
  .schnappsEnv$defaultValues[["gQC_um_init"]] = input$gQC_um_init
})

observe(label = "obs_gQC_um_negative_sample_rate", x = {
  .schnappsEnv$defaultValues[["gQC_um_negative_sample_rate"]] = input$gQC_um_negative_sample_rate
})

observe(label = "obs_gQC_um_randSeed", x = {
  .schnappsEnv$defaultValues[["gQC_um_randSeed"]] = input$gQC_um_randSeed
})

observe(label = "obs_gQC_um_local_connectivity", x = {
  .schnappsEnv$defaultValues[["gQC_um_local_connectivity"]] = input$gQC_um_local_connectivity
})
observe(label ="obs_gQC_um_bandwidth", x = {
  .schnappsEnv$defaultValues[["gQC_um_bandwidth"]] = input$gQC_um_bandwidth
})

observe(label ="obs_gQC_um_n_epochs", x = {
  .schnappsEnv$defaultValues[["gQC_um_n_epochs"]] = input$gQC_um_n_epochs
})

observe(label ="obs_gQC_um_set_op_mix_ratio", x = {
  .schnappsEnv$defaultValues[["gQC_um_set_op_mix_ratio"]] = input$gQC_um_set_op_mix_ratio
})

observe(label ="obs_gQC_um_metric", x = {
  .schnappsEnv$defaultValues[["gQC_um_metric"]] = input$gQC_um_metric
})

observe(label ="obs_gQC_um_min_dist", x = {
  .schnappsEnv$defaultValues[["gQC_um_min_dist"]] = input$gQC_um_min_dist
})
observe(label ="obs_gQC_binSize", x = {
  .schnappsEnv$defaultValues[["gQC_binSize"]] = input$gQC_binSize
})


# gQC_update3DInput ----
#' gQC_update3DInput
#' update axes for tsne display
gQC_update3DInput <- reactive({
  if (DEBUG) cat(file = stderr(), "gQC_update3DInput started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_update3DInput")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_update3DInput")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_update3DInput", id = "gQC_update3DInput", duration = NULL)
  }
  
  projections <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  # choices = colnames(projections)[unlist(lapply(colnames(projections), function(x) !is.factor(projections[,x])))]
  choices <- colnames(projections)
  # Can also set the label and select items
  updateSelectInput(session, "gQC_dim3D_x",
                    choices = choices,
                    selected = .schnappsEnv$gQC_X1
  )
  
  updateSelectInput(session, "gQC_dim3D_y",
                    choices = choices,
                    selected = .schnappsEnv$gQC_X2
  )
  updateSelectInput(session, "gQC_dim3D_z",
                    choices = choices,
                    selected = .schnappsEnv$gQC_X3
  )
  updateSelectInput(session, "gQC_col3D",
                    choices = colnames(projections),
                    selected = .schnappsEnv$gQC_col
  )
})

# observer of UMAP button ----
observe(label = "ob_UMAPParams", {
  # save(file = "~/SCHNAPPsDebug/ob_UMAPParams.RData", list = c(ls(), ".schnappsEnv"))
  # load("~/SCHNAPPsDebug/updateButtonColor.RData")
  deepDebug()
  if (DEBUG) cat(file = stderr(), "observe umapVars\n")
  
  input$activateUMAP
  setRedGreenButtonCurrent(
    vars = list(
      c("gQC_um_randSeed", input$gQC_um_randSeed),
      c("gQC_um_n_neighbors", input$gQC_um_n_neighbors),
      c("gQC_um_n_components", input$gQC_um_n_components),
      c("gQC_um_n_epochs", input$gQC_um_n_epochs),
      # c("um_alpha", input$um_alpha),
      c("gQC_um_init", input$gQC_um_init),
      c("gQC_um_min_dist", input$gQC_um_min_dist),
      c("gQC_um_set_op_mix_ratio", input$gQC_um_set_op_mix_ratio),
      c("gQC_um_local_connectivity", input$gQC_um_local_connectivity),
      c("gQC_um_bandwidth", input$gQC_um_bandwidth),
      c("um_gamma", input$um_gamma),
      c("gQC_um_negative_sample_rate", input$gQC_um_negative_sample_rate),
      c("gQC_um_metric", input$gQC_um_metric),
      c("gQC_um_spread", input$gQC_um_spread)
    )
  )
  
  updateButtonColor(buttonName = "activateUMAP", parameters = c(
    "gQC_um_randSeed", "gQC_um_n_neighbors", "gQC_um_n_components", "gQC_um_n_epochs", 
    "gQC_um_init", "gQC_um_min_dist", "gQC_um_set_op_mix_ratio", 
    "gQC_um_local_connectivity", "gQC_um_bandwidth", "um_gamma", 
    "gQC_um_negative_sample_rate", "gQC_um_metric", "gQC_um_spread"
  ))
})



# observe: cellNameTable_rows_selected ----
observe(label = "ob_tsneParams", {
  deepDebug()
  if (DEBUG) cat(file = stderr(), "observe tsneVars\n")
  out <- tsne()
  if (is.null(out)) {
    .schnappsEnv$calculated_gQC_tsneDim <- "NA"
  }
  input$updatetsneParameters
  
  setRedGreenButtonCurrent(
    vars = list(
      c("gQC_tsneDim", input$gQC_tsneDim),
      c("gQC_tsnePerplexity", input$gQC_tsnePerplexity),
      c("gQC_tsneTheta", input$gQC_tsneTheta),
      c("gQC_tsneSeed", input$gQC_tsneSeed)
    )
  )
  updateButtonColor(buttonName = "updatetsneParameters", parameters = c(
    "gQC_tsneDim", "gQC_tsnePerplexity",
    "gQC_tsneTheta", "gQC_tsneSeed"
  ))
})

# gQC_tsne_main ----
output$gQC_tsne_main <- plotly::renderPlotly({
  if (DEBUG) cat(file = stderr(), "gQC_tsne_main started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_tsne_main")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_tsne_main")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_tsne_main", id = "gQC_tsne_main", duration = NULL)
  }
  
  upI <- gQC_update3DInput()
  projections <- projections()
  dimX <- input$gQC_dim3D_x
  dimY <- input$gQC_dim3D_y
  dimZ <- input$gQC_dim3D_z
  dimCol <- input$gQC_col3D
  scols <- sampleCols$colPal
  ccols <- clusterCols$colPal
  
  if (is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$gQC_tsne_main:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/gQC_tsne_main.RData", list = c(ls()))
  }
  # cp =load(file="~/SCHNAPPsDebug/gQC_tsne_main.RData")
  
  retVal <- tsnePlot(projections, dimX, dimY, dimZ, dimCol, scols, ccols)
  
  exportTestValues(tsnePlot = {
    str(retVal)
  })
  layout(retVal)
})

# gQC_umap_main 2D plot ----
callModule(
  clusterServer,
  "gQC_umap_main",
  projections
)

# gQC_projectionTableMod ----
callModule(
  tableSelectionServer,
  "gQC_projectionTableMod",
  projectionTable, caption = "Table with projections"
)

# gQC_projectionCombTableMod ----
callModule(
  tableSelectionServer,
  "gQC_projCombTableMod",
  projectionTable, caption = "Tables with projections"
)

# gQC_plotUmiHist ----
output$gQC_plotUmiHist <- plotly::renderPlotly({
  if (DEBUG) cat(file = stderr(), "gQC_plotUmiHist started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_plotUmiHist")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_plotUmiHist")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_plotUmiHist", id = "gQC_plotUmiHist", duration = NULL)
  }
  
  scEx <- scEx()
  scols <- sampleCols$colPal
  binSize <- input$gQC_binSize
  
  if (is.null(scEx)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/gQC_plotUmiHist.RData", list = c(ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/gQC_plotUmiHist.RData")
  
  dat <- data.frame(counts = Matrix::colSums(assays(scEx)[["counts"]]))
  dat$sample <- colData(scEx)$sampleNames
  
  
  fig <- plotly::plot_ly(alpha = 1,
                         nbinsx = binSize)
  # dat[dat$sample == levels(colData(scEx)$sampleNames)[[1]],],
  # x = ~counts, 
  # # y = ~counts,
  # type="histogram")
  lev = levels(colData(scEx)$sampleNames)
  for (idx in seq_along(lev)) {
    fig <- fig %>% add_trace(
      type = 'histogram', color = I(scols[idx]), name = lev[idx],
      x = dat[dat$sample == levels(colData(scEx)$sampleNames)[[idx]],"counts"]
    )
  }
  fig <- fig %>% layout(
    barmode="stack",
    bargap=0.1,
    title = "Histogram of UMIs",
    yaxis = list(title = "Number of cells"),
    xaxis = list(title = "UMI count"))
  
  fig
  # marker = list(color = scols))
  
  
  
  # retVal <- ggplot(data = dat, aes(counts, fill = sample)) +
  #   geom_histogram(bins = 50) +
  #   labs(title = "Histogram for raw counts", x = "count", y = "Frequency") +
  #   scale_fill_manual(values = scols, aesthetics = "fill")
  # 
  .schnappsEnv[["gQC_plotUmiHist"]] <- fig
  return(fig)
})

# gQC_plotSampleHist -----
output$gQC_plotSampleHist <- plotly::renderPlotly({
  if (DEBUG) cat(file = stderr(), "gQC_plotSampleHist started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_plotSampleHist")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_plotSampleHist")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_plotSampleHist", id = "gQC_plotSampleHist", duration = NULL)
  }
  
  sampleInf <- sampleInfo()
  scols <- sampleCols$colPal
  
  if (is.null(sampleInf)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/sampleHist.RData", list = c(ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/sampleHist.RData")
  retVal <- gQC_sampleHistFunc(sampleInf, scols)
  .schnappsEnv[["gQC_plotSampleHist"]] <- retVal
  return(retVal)
})

output$gQC_variancePCA <- renderPlot({
  if (DEBUG) cat(file = stderr(), "gQC_variancePCA started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_variancePCA")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_variancePCA")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_variancePCA", id = "gQC_variancePCA", duration = NULL)
  }
  pca <- pcaReact()
  if (is.null(pca)) {
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/gQC_variancePCA.RData", list = c(ls()))
  }
  # load(file = "~/SCHNAPPsDebug/gQC_variancePCA.RData")
  
  # h2("Variances of PCs")
  
  
  
  df <- data.frame(var = pca$var_pcs, pc = 1:length(pca$var_pcs))
  retVal <- plotHistVarPC(df, pc, var) 
  
  af = plotHistVarPC
  # remove env because it is too big
  environment(af) = new.env(parent = emptyenv())
  .schnappsEnv[["gQC_variancePCA"]] <- list(plotFunc = af,
                                            df, pc, var)
  
  # .schnappsEnv[["gQC_variancePCA"]] <- retVal
  return(retVal)
  # barplot(pca$var_pcs, main = "Variance captured by first PCs")
})

# rename projections observers ----

observeEvent(
  label = "ob30",
  eventExpr = input$updatePrjsButton,
  handlerExpr = {
    deepDebug()
    if (DEBUG) cat(file = stderr(), "updatePrjsButton\n")
    oldPrj <- input$oldPrj
    newPrj <- input$newPrj
    projections <- projections()
    acn = allCellNames()
    newPrjs <- projectionsTable$newProjections
    
    if (is.null(projections)) {
      return(NULL)
    }
    
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/updatePrjsButton.RData",
           list = c("normaliztionParameters", ls())
      )
    }
    # cp = load(file="~/SCHNAPPsDebug/updatePrjsButton.RData")
    # deepDebug()
    if (newPrj %in% colnames(projections)) {
      showNotification(
        "New column name already used",
        type = "error",
        duration = NULL
      )
      return(NULL)
    }
    if (ncol(newPrjs) == 0) {
      newPrjs = data.frame(row.names = acn)
      newPrjs[,newPrj] = NA
      newPrjs[rownames(projections),newPrj] <- projections[, oldPrj, drop = FALSE]
    } else {
      # newPrjs <- cbind(newPrjs[rownames(projections), , drop = FALSE], projections[, oldPrj, drop = FALSE])
      # deepDebug()
      newPrjs <- dplyr::full_join(
        tibble::rownames_to_column(newPrjs), 
        tibble::rownames_to_column(projections[, oldPrj, drop = FALSE]), 
        by='rowname')
      rownames(newPrjs) = newPrjs[,1]
      newPrjs = newPrjs[,-1]
    }
    colnames(newPrjs)[ncol(newPrjs)] <- newPrj
    projectionsTable$newProjections <- newPrjs
  }
)

# # rename projections
# observe(label = "ob27", {
#   projections <- projections()
#   
# })

observe(label = "ob28", {
  deepDebug()
  input$newPrj
  updateTextInput(session, "newPrj", value = make.names(input$newPrj, unique = TRUE))
})

observeEvent(
  label = "ob29",
  eventExpr = input$delPrjsButton,
  handlerExpr = {
    deepDebug()
    if (DEBUG) cat(file = stderr(), "updatePrjsButton\n")
    newPrjs <- projectionsTable$newProjections
    delPrj <- input$delPrj
    if (is.null(projections)) {
      return(NULL)
    }
    if (!delPrj %in% colnames(newPrjs)) {
      return(NULL)
    }
    # deepDebug()
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/delPrjsButton.RData",
           list = c("normaliztionParameters", ls())
      )
    }
    # load(file="~/SCHNAPPsDebug/delPrjsButton.RData")
    
    projectionsTable$newProjections <- newPrjs[, -which(colnames(newPrjs) == delPrj), drop = FALSE]
  }
)

# combine projections observers ----

observeEvent(
  label = "gQC_updateCombPrjsButton",
  eventExpr = input$gQC_updateCombPrjsButton,
  handlerExpr = {
    if (DEBUG) cat(file = stderr(), "gQC_updateCombPrjsButton\n")
    prj1 <- input$gQC_combPrj1
    prj2 <- input$gQC_combPrj2
    newPrj <- make.names(input$gQC_newCombPrj)
    projections <- projections()
    newPrjs <- projectionsTable$newProjections
    acn = allCellNames()
    
    if (is.null(projections)) {
      return(NULL)
    }
    if (!all(c(prj1, prj2) %in% colnames(projections))) {
      return(NULL)
    }
    
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/gQC_updateCombPrjsButton.RData",
           list = c("normaliztionParameters", ls())
      )
    }
    # cp=  load(file="~/SCHNAPPsDebug/gQC_updateCombPrjsButton.RData")
    if (newPrj %in% colnames(projections)) {
      showNotification(
        "New column name already used",
        type = "error",
        duration = NULL
      )
      return(NULL)
    }
    # deepDebug()
    combProjections = data.frame(row.names = rownames(projections), 
                                 newPrj=  paste(projections[,prj1], projections[,prj2], sep = " - ") %>% as.factor())
    if (length(levels(combProjections[,1])) > 100) {
      out = showModal(verifyLevelModal(NLevel = length(levels(combProjections[,1]))))
      if (DEBUG) cat(file = stderr(), paste("gQC_updateCombPrjsButton modal out:", out, "\n"))
      # deepDebug()
    }
    if (ncol(newPrjs) == 0) {
      newPrjs <- data.frame(row.names = acn)
      newPrjs[rownames(combProjections),newPrj] = combProjections[,1]
      # rownames(newPrjs) = rownames(projections)
    } else {
      # newPrjs <- cbind(newPrjs[rownames(projections), , drop = FALSE], combProjections)
      # deepDebug()
      newPrjs <- dplyr::left_join(
        tibble::rownames_to_column(newPrjs), 
        tibble::rownames_to_column(combProjections), 
        by='rowname')
      rownames(newPrjs) = newPrjs[,1]
      newPrjs = newPrjs[,-1]
    }
    colnames(newPrjs)[ncol(newPrjs)] <- newPrj
    projectionsTable$newProjections  <- newPrjs
  }
)

# output$gQC_orgLevels ----
output$gQC_orgLevels = renderText({
  rnProj = input$gQC_rnProj
  projections = projections()
  shiny::req(rnProj)
  shiny::req(projections)
  if (! rnProj %in% colnames(projections)) return(NULL)
  # deepDebug()
  paste(levels(factor(projections[,rnProj])), collapse = ", ")
})

# rename projection levels ----


verifyLevelModal <- function(NLevel, failed = FALSE) {
  modalDialog(
    span(paste(
      "There are ", NLevel, "new levels, are you sure you want to do this?\n")
    )
    ,
    footer = tagList(
      modalButton("Cancel"),
      actionButton("commentok", "OK")
    )
  )
}

# gQC_rearrange levels ----
observeEvent(eventExpr = input$gQC_raProj,
             label = "raLevBtn",
             handlerExpr = {
               deepDebug()
               projections <- projections()
               projFactors <- projFactors()
               if(is.null(projections)) return()
               if(is.null(projFactors)) return()
               if(input$gQC_raProj %in% projFactors){
                 
                 projLevels = levels(projections[,input$gQC_raProj])
                 updateOrderInput(
                   session,
                   'gQC_newRaLev',
                   items = projLevels
                 )
               }else{
                 updateOrderInput(
                   session,
                   'gQC_newRaLev',
                   items = "not yet"
                 )
               }
               
             }
)
observeEvent(eventExpr = input$gQC_rearrangeLevButton,
             label = "raLevBtn",
             handlerExpr = {
               deepDebug()
               newProjName = make.names(input$gQC_newRaPrj)
               newLevelOrder = input$gQC_newRaLev
               projections = projections()
               raProj = input$gQC_raProj
               acn = allCellNames()
               newPrjs <- projectionsTable$newProjections
               if (is.null(projections)) {
                 return(NULL)
               }
               if (.schnappsEnv$DEBUGSAVE) {
                 save(file = "~/SCHNAPPsDebug/gQC_rearrangeButton.RData",
                      list = c("normaliztionParameters", ls())
                 )
               }
               # cp=  load(file="~/SCHNAPPsDebug/gQC_rearrangeButton.RData")
               orgLevelNames = levels(factor(projections[,raProj]))
               newLbVec = stringr::str_trim(str_split(newLevelOrder, ","))
               names(newLbVec) = orgLevelNames
               # deepDebug()
               # sampe projections as displayed, i.e. only those available for the cells
               # otherwise the diplay (output$gQC_orgLevels) has to be changed as well
               projections[,raProj] =  factor(projections[,raProj])
               if(is.null(
                 tryCatch({
                   
                   if (ncol(newPrjs) == 0) {
                     newPrjs = data.frame(row.names = acn)
                     newPrjs[,newProjName] = "NA"
                     # drop = TRUE: we re interested in the vector not the data frame
                     newPrjs[rownames(projections),newProjName] <- as.character(projections[, raProj, drop = TRUE])
                   } else {
                     # deepDebug()
                     newPrjs <- dplyr::full_join(
                       tibble::rownames_to_column(newPrjs), 
                       tibble::rownames_to_column(projections[, raProj, drop = FALSE]), 
                       by='rowname')
                     rownames(newPrjs) = newPrjs[,1]
                     newPrjs = newPrjs[,-1]
                     # newPrjs <- cbind(newPrjs[rownames(projections), , drop = FALSE], projections[,raProj])
                   }
                   
                   newPrjs[,ncol(newPrjs)] = as.factor(newPrjs[,ncol(newPrjs)])
                   
                   # in case there was NA introduced by hidden cells
                   if ("NA" %in% levels(newPrjs[,ncol(newPrjs)]) & !"NA" %in% newLevelOrder ){
                     newLevelOrder = c(newLevelOrder, "NA")
                   }
                   if (!length(levels(newPrjs[,ncol(newPrjs)])) == length(newLevelOrder) ){
                     cat(file = stderr(), paste("number of levels not correct\n\nold levels:\n"))
                     cat(file = stderr(), levels(newPrjs[,ncol(newPrjs)]))
                     cat(file = stderr(), paste("\n\n\nnew levels:\n"))
                     cat(file = stderr(), newLevelOrder)
                     cat(file = stderr(), paste("\n"))
                     showNotification("number of levels not correct. See console", id = "renameProbl", duration = NULL, type = "error")
                     return(NULL)
                   }
                   newPrjs[,ncol(newPrjs)] = factor(newPrjs[,ncol(newPrjs)], levels = newLevelOrder)
                 }, error=function(w){
                   # deepDebug()
                   cat(file = stderr(), paste("something went wrong during releveling", w,"\n"))
                   showNotification("problem with names", id = "renameProbl", duration = NULL, type = "error")
                   return(NULL)
                 }))) return(NULL)
               # newProjName = make.unique(c(colnames(projections),newProjName))[length(c(colnames(projections),newProjName))]
               # updateTextInput(session, "gQC_newRnPrj", value = newProjName)
               colnames(newPrjs)[ncol(newPrjs)] <- newProjName
               projectionsTable$newProjections  <- newPrjs
               
             })


# gQC_renameLevButton ----
observeEvent(eventExpr = input$gQC_renameLevButton,
             label = "rnLevBtn",
             handlerExpr = {
               deepDebug()
               newLables = input$gQC_renameLev
               rnProj = input$gQC_rnProj
               newProjName = make.names(input$gQC_newRnPrj)
               projections = projections()
               acn = allCellNames()
               newPrjs <- projectionsTable$newProjections
               if (is.null(projections)) {
                 return(NULL)
               }
               orgLevelNames = levels(factor(projections[,rnProj]))
               newLbVec = stringr::str_trim(str_split(newLables, ",")[[1]])
               names(newLbVec) = orgLevelNames
               
               if (.schnappsEnv$DEBUGSAVE) {
                 save(file = "~/SCHNAPPsDebug/gQC_renameLevButton.RData",
                      list = c("normaliztionParameters", ls())
                 )
               }
               # cp=  load(file="~/SCHNAPPsDebug/gQC_renameLevButton.RData")
               # deepDebug()
               # sampe projections as displayed, i.e. only those available for the cells
               # otherwise the diplay (output$gQC_orgLevels) has to be changed as well
               proj2Add =  projections[,rnProj,drop=FALSE]
               proj2Add[,rnProj] = factor(proj2Add[,rnProj])
               
               if(is.null(
                 tryCatch({
                   
                   if (ncol(newPrjs) == 0) {
                     newPrjs = data.frame(row.names = acn)
                     newPrjs[,newProjName] = "NA"
                     # drop = TRUE: we re interested in the vector not the data frame
                     newPrjs[rownames(projections),newProjName] <- as.character(projections[, rnProj, drop = TRUE])
                   } else {
                     # deepDebug()
                     newPrjs <- dplyr::full_join(
                       tibble::rownames_to_column(newPrjs), 
                       tibble::rownames_to_column(projections[, rnProj, drop = FALSE]), 
                       by='rowname')
                     rownames(newPrjs) = newPrjs[,1]
                     newPrjs = newPrjs[,-1]
                     # newPrjs <- cbind(newPrjs[rownames(projections), , drop = FALSE], projections[,rnProj])
                   }
                   
                   newPrjs[,ncol(newPrjs)] = as.factor(newPrjs[,ncol(newPrjs)])
                   
                   
                   if ("NA" %in% levels(newPrjs[,ncol(newPrjs)]) & !"NA" %in% stringr::str_trim(newLbVec) ){
                     newLevelNames =  levels(newPrjs[,ncol(newPrjs)])
                     naPos = which ("NA" == newLevelNames)
                     newLbVec = newLbVec[newLevelNames]
                     newLbVec[which(is.na(newLbVec))] = "NA"
                   }
                   if (!length(levels(newPrjs[,ncol(newPrjs)])) == length(stringr::str_trim(newLbVec)) ){
                     cat(file = stderr(), paste("number of levels not correct\n\nold levels:\n"))
                     cat(file = stderr(), levels(newPrjs[,ncol(newPrjs)]))
                     cat(file = stderr(), paste("\n\n\nnew levels:\n"))
                     cat(file = stderr(), stringr::str_trim(newLbVec))
                     cat(file = stderr(), paste("\n"))
                     showNotification("number of levels not correct. See console", id = "renameProbl", duration = NULL, type = "error")
                     return(NULL)
                   }
                   levels(newPrjs[,ncol(newPrjs)]) = stringr::str_trim(newLbVec)
                 }, error=function(w){
                   # deepDebug()
                   cat(file = stderr(), paste("something went wrong during releveling", w,"\n"))
                   showNotification("problem with names", id = "renameProbl", duration = NULL, type = "error")
                   return(NULL)
                 }))) return(NULL)
               newProjName = make.unique(c(colnames(projections),newProjName))[length(c(colnames(projections),newProjName))]
               updateTextInput(session, "gQC_newRnPrj", value = newProjName)
               colnames(newPrjs)[ncol(newPrjs)] <- newProjName
               projectionsTable$newProjections  <- newPrjs
               
             })

observeEvent(eventExpr = input$gQC_rnProj,
             label = "gqc1",
             handlerExpr = {
               deepDebug()
               rnProj = input$gQC_rnProj
               projections = projections()
               shiny::req(rnProj)
               shiny::req(projections)
               if (! rnProj %in% colnames(projections)) return(NULL)
               
               # deepDebug()
               updateTextAreaInput(session, inputId = "gQC_renameLev", value = paste(as.character(levels(factor(projections[,rnProj]))), 
                                                                                     collapse = ", "))
             })


# rename projections
.schnappsEnv$gQC_combPrj1 <- "tsne1"
.schnappsEnv$gQC_combPrj2 <- "tsne1"
.schnappsEnv$gQC_rnProj <- "tsne1"

observe(label = "ob27b", {
  deepDebug()
  projections <- projections()
  projFactors <- projFactors()
  
  # only factorials?
  updateSelectInput(session, "gQC_raProj",
                    choices = projFactors,
                    selected = .schnappsEnv$gQC_raProj
  )
  updateSelectInput(session, "gQC_combPrj1",
                    choices = colnames(projections),
                    selected = .schnappsEnv$gQC_combPrj1
  )
  updateSelectInput(session, "gQC_combPrj2",
                    choices = colnames(projections),
                    selected = .schnappsEnv$gQC_combPrj2
  )
  updateSelectInput(session, "gQC_rnProj",
                    choices = colnames(projections),
                    selected = .schnappsEnv$gQC_rnProj
  )
  updateSelectInput(session, "oldPrj",
                    choices = c(colnames(projections)),
                    selected = .schnappsEnv$oldPrj
  )
  updateSelectInput(session, "delPrj",
                    choices = c(colnames(projectionsTable$newProjections)),
                    selected = .schnappsEnv$delPrj
  )
  updateSelectInput(session, "gQC_windProj",
                    choices = projFactors,
                    selected = .schnappsEnv$gQC_windProj
  )
  
  
})
observe(label = "ob27c", {
  if (DEBUG) cat(file = stderr(), "observe: gQC_combPrj1\n")
  .schnappsEnv$gQC_combPrj1 <- input$gQC_combPrj1
})
observe(label = "ob27d", {
  if (DEBUG) cat(file = stderr(), "observe: gQC_combPrj2\n")
  .schnappsEnv$gQC_combPrj2 <- input$gQC_combPrj2
})
observe(label = "ob27e", {
  if (DEBUG) cat(file = stderr(), "observe: gQC_rnProj\n")
  .schnappsEnv$gQC_rnProj <- input$gQC_rnProj
})
observe(label = "ob27f", {
  if (DEBUG) cat(file = stderr(), "observe: oldPrj\n")
  .schnappsEnv$oldPrj <- input$oldPrj
})
observe(label = "ob27g", {
  if (DEBUG) cat(file = stderr(), "observe: delPrj\n")
  .schnappsEnv$delPrj <- input$delPrj
})
observe(label = "ob27h", {
  if (DEBUG) cat(file = stderr(), "observe: gQC_windProj\n")
  .schnappsEnv$gQC_windProj <- input$gQC_windProj
  .schnappsEnv$defaultValues[["gQC_windProj"]] <- input$gQC_windProj
})
observe(label = "ob27i", {
  if (DEBUG) cat(file = stderr(), "observe: gQC_raProj\n")
  .schnappsEnv$gQC_raProj <- input$gQC_raProj
})


# rename levels

output$gQC_renameLev <- renderText({"text"})


# WIND ----

output$gQC_windHC <- renderPlot({
  require(Wind)
  # remotes::install_github("renozao/xbioc")
  # library(xbioc)
  if ("xbioc" %in% rownames(installed.packages())){
    require(xbioc)
  }  else {
    is_logscale <- function(x) {return(T)
      cat(file = stderr(), "Please install xbioc: remotes::install_github('renozao/xbioc')")
    }
  }
  if (DEBUG) cat(file = stderr(), "gQC_windHC started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "gQC_windHC")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gQC_windHC")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_windHC", id = "gQC_windHC", duration = NULL)
  }
  
  scEx <- scEx()
  projections <- projections()
  pca = pcaReact()
  gQC_windProj <- input$gQC_windProj
  # browser()
  if (is.null(projections) | is.null(scEx) | !gQC_windProj %in% colnames(projections)) {
    return(NULL)
  }
  if (length(levels(projections[,gQC_windProj]))<3) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("Projections have less than 3 levels", id = "gQC_windHCPR", duration = 20, type = "warning")
    }
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/gQC_windHC.RData", list = c(ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/gQC_windHC.RData")
  Y <- as.matrix(assays(scEx)[[1]])
  # if(is_logscale(Y)) {
  #   Y = exp(Y)
  # }
  
  trueclass <- projections[,gQC_windProj]
  ctStruct = tryCatch({
    createRef(Y, classes = trueclass)
    },error = function(e) {
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("Problem with WIND", type = "warning", duration = NULL)
      }
      cat(file = stderr(), paste("\n+++++ Error in WIND\n\t", e, "\n"))
      return(NULL)
    }
  )
  if(is.null(ctStruct)) return(NULL)
  plot(ctStruct$hc, xlab="", axes=FALSE, ylab="", ann=FALSE)
})
