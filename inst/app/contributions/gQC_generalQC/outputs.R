# source("moduleServer.R", local = TRUE)
# source("reactives.R", local = TRUE)

# TODO: verify that this anything and then integrate in DUMMY
myZippedReportFiles <- c("gqcProjections.csv")



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
  # browser()
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
    yaxis = list(title = "Number of samples"),
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
      save(
        file = "~/SCHNAPPsDebug/updatePrjsButton.RData",
        list = c("normaliztionParameters", ls())
      )
    }
    # cp = load(file="~/SCHNAPPsDebug/updatePrjsButton.RData")
    # browser()
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
      # browser()
      newPrjs <- dplyr::left_join(
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
  input$newPrj
  updateTextInput(session, "newPrj", value = make.names(input$newPrj, unique = TRUE))
})

observeEvent(
  label = "ob29",
  eventExpr = input$delPrjsButton,
  handlerExpr = {
    if (DEBUG) cat(file = stderr(), "updatePrjsButton\n")
    newPrjs <- projectionsTable$newProjections
    delPrj <- input$delPrj
    if (is.null(projections)) {
      return(NULL)
    }
    if (!delPrj %in% colnames(newPrjs)) {
      return(NULL)
    }
    # browser()
    if (.schnappsEnv$DEBUGSAVE) {
      save(
        file = "~/SCHNAPPsDebug/delPrjsButton.RData",
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
      save(
        file = "~/SCHNAPPsDebug/gQC_updateCombPrjsButton.RData",
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
    # browser()
    combProjections = data.frame(row.names = rownames(projections), 
                                 paste(projections[,prj1], projections[,prj2], sep = " - ") %>% as.factor())
    if (length(levels(combProjections)) > 100) {
      out = showModal(verifyLevelModal(NLevel = length(levels(combProjections))))
      # browser()
    }
    if (ncol(newPrjs) == 0) {
      newPrjs <- data.frame(row.names = acn)
      newPrjs[rownames(combProjections),newPrj] = combProjections
      # rownames(newPrjs) = rownames(projections)
    } else {
      # newPrjs <- cbind(newPrjs[rownames(projections), , drop = FALSE], combProjections)
      # browser()
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
  # browser()
  paste(levels(projections[,rnProj]), collapse = ", ")
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

# gQC_renameLevButton ----
observeEvent(eventExpr = input$gQC_renameLevButton,
             label = "rnBtn",
             handlerExpr = {
               newLables = input$gQC_renameLev
               rnProj = input$gQC_rnProj
               newProjName = make.names(input$gQC_newRnPrj)
               projections = projections()
               acn = allCellNames()
               newPrjs <- projectionsTable$newProjections
               if (is.null(projections)) {
                 return(NULL)
               }
               
               if (.schnappsEnv$DEBUGSAVE) {
                 save(
                   file = "~/SCHNAPPsDebug/gQC_renameLevButton.RData",
                   list = c("normaliztionParameters", ls())
                 )
               }
               # cp=  load(file="~/SCHNAPPsDebug/gQC_renameLevButton.RData")
               # browser()
               
               if(is.null(
                 tryCatch({
                   newLbVec = str_split(newLables, ",")[[1]]
                   if (ncol(newPrjs) == 0) {
                     newPrjs = data.frame(row.names = acn)
                     newPrjs[,newPrj] = NA
                     newPrjs[rownames(projections),newPrj] <- projections[, rnProj, drop = FALSE]
                   } else {
                     # browser()
                     newPrjs <- dplyr::left_join(
                       tibble::rownames_to_column(newPrjs), 
                       tibble::rownames_to_column(projections[, rnProj, drop = FALSE]), 
                       by='rowname')
                     rownames(newPrjs) = newPrjs[,1]
                     newPrjs = newPrjs[,-1]
                     # newPrjs <- cbind(newPrjs[rownames(projections), , drop = FALSE], projections[,rnProj])
                   }
                   newPrjs[,ncol(newPrjs)] = as.factor(newPrjs[,ncol(newPrjs)])
                   levels(newPrjs[,ncol(newPrjs)]) = stringr::str_trim(newLbVec)
                 }, error=function(w){
                   # browser()
                   cat(file = stderr(), paste("something went wrong during releveling", w,"\n"))
                   showNotification("problem with names", id = "renameProbl", duration = NULL)
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
               rnProj = input$gQC_rnProj
               projections = projections()
               shiny::req(rnProj)
               shiny::req(projections)
               if (! rnProj %in% colnames(projections)) return(NULL)
               # browser()
               updateTextAreaInput(session, inputId = "gQC_renameLev", value = paste(levels(projections[,rnProj]), collapse = ", "))
             })


# rename projections
.schnappsEnv$gQC_combPrj1 <- "tsne1"
.schnappsEnv$gQC_combPrj2 <- "tsne1"
.schnappsEnv$gQC_rnProj <- "tsne1"

observe(label = "ob27b", {
  projections <- projections()
  projFactors <- projFactors()
  
  # only factorials?
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
})


# rename levels

output$gQC_renameLev <- renderText({"text"})


# WIND ----

output$gQC_windHC <- renderPlot({
  require(Wind)
  # remotes::install_github("renozao/xbioc")
  # library(xbioc)
  if ("xbioc" %in% rownames(installed.packages())){
    require(xbioc)}  else {
      is_logscale <- function(x) {return(T)}
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
  
  scEx_log <- scEx_log()
  projections <- projections()
  pca = pca()
  gQC_windProj <- input$gQC_windProj
  
  if (is.null(projections) | is.null(scEx_log)) {
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
  Y <- as.matrix(assays(scEx_log)[[1]])
  if(is_logscale(Y)) {
    Y = exp(Y)
  }
  
  trueclass <- projections[,gQC_windProj]
  ctStruct = createRef(Y, trueclass)
  plot(ctStruct$hc, xlab="", axes=FALSE, ylab="", ann=FALSE)
})
