suppressMessages(require(Rsomoclu))

coE_SOM_dataInput <- callModule(
  cellSelectionModule,
  "coE_SOM_dataInput"
)

if(!is.null(.schnappsEnv$historyPath)){
  mC_metaCell_dataInput <- callModule(
    cellSelectionModule,
    "mC_metaCell_dataInput"
  )
}


# SOM heatmap module -----
callModule(
  pHeatMapModule,
  "coE_heatmapSOM",
  coE_heatmapSOMReactive
)

# observer of button Color SOM ----
observe(label = "ob_somParameter", 
        {
          if (DEBUG) cat(file = stderr(), "ob_somParameter\n")
          # browser()
          input$updateSOMParameters
          setRedGreenButtonCurrent(
            vars = list(
              c("coE_geneSOM", input$coE_geneSOM),
              c("coE_dimSOM", input$coE_dimSOM),
              c("coE_SOM_dataInput-Mod_PPGrp", input$'coE_SOM_dataInput-Mod_PPGrp'),
              c("coE_SOM_dataInput-Mod_clusterPP", input$'coE_SOM_dataInput-Mod_clusterPP')
            )
          )
          updateButtonColor(buttonName = "updateSOMParameters", parameters = c(
            "coE_geneSOM", "coE_dimSOM",
            "coE_SOM_dataInput-Mod_PPGrp", "coE_SOM_dataInput-Mod_clusterPP"
          ))
          
        })

observe(label = "somxy",{
  .schnappsEnv$defaultValues[["coE_dimSOMX"]] = input$coE_dimSOMX
  .schnappsEnv$defaultValues[["coE_dimSOMY"]] = input$coE_dimSOMY
})

output$coE_SOMcodebook <- renderPlot({
  sommap = coE_somMapReact()
  if (is.null(sommap)) return(NULL)
  plot(sommap, type="codes", main = "Codes")
})


output$coE_SOMcomponents <- renderPlot({
  sommap = coE_somMapReact()
  if (is.null(sommap)) return(NULL)
  
  plot(sommap, type = "property", property = sommap$codes[[1]][,1],
       main = colnames(sommap$codes)[1])
})
output$coE_SOMuMat <- renderPlot({
  sommap = coE_somMapReact()
  if (is.null(sommap)) return(NULL)
  plot(sommap, type="dist.neighbours")
})

output$coE_somInfo <- renderText({
  
  idxx = input$coE_dimSOMX - 1
  idxy = input$coE_dimSOMY - 1
  res2 = coE_somTrainReact()
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_somInfo.RData", list = c(ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/coE_somInfo.RData")
  
  if (is.null(res2)) {
    return("Genes in neuron\n")
  }
  
  paste("Genes in neuron (", length(which(res2$globalBmus[,1]==idxx & res2$globalBmus[,2]==idxy)), ")\n",
        paste(names(which(res2$globalBmus[,1]==idxx & res2$globalBmus[,2]==idxy)), collapse = ", ", sep = ",")
  )
  
})


output$coE_somInfoSymbol <- renderText({
  
  idxx = input$coE_dimSOMX - 1
  idxy = input$coE_dimSOMY - 1
  res2 = coE_somTrainReact()
  scEx_log = scEx_log()
  
  featureData <- rowData(scEx_log)
  geneName = geneName2Index(genesin, featureData)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_somInfo.RData", list = c(ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/coE_somInfo.RData")
  
  if (is.null(res2)) {
    return("Genes in neuron\n")
  }
  
  paste("Genes in neuron (", length(which(res2$globalBmus[,1]==idxx & res2$globalBmus[,2]==idxy)), ")\n",
        paste(featureData[names(which(res2$globalBmus[,1]==idxx & res2$globalBmus[,2]==idxy)), "symbol"], collapse = ", ", sep = ",")
  )
  
})


# observer of button Color metaCell ----
observe(label = "ob_metaCellParameter", 
        {
          if (DEBUG) cat(file = stderr(), "ob_metaCellParameter\n")
          # browser()
          input$updateMetaCellParameters
          setRedGreenButtonCurrent(
            vars = list(
              c("coE_geneSOM", input$coE_geneSOM),
              c("coE_dimSOM", input$coE_dimSOM),
              c("coE_SOM_dataInput-Mod_PPGrp", input$'coE_SOM_dataInput-Mod_PPGrp'),
              c("coE_SOM_dataInput-Mod_clusterPP", input$'coE_SOM_dataInput-Mod_clusterPP')
            )
          )
          updateButtonColor(buttonName = "updateMetaCellParameters", parameters = c(
            "coE_geneSOM", "coE_dimSOM",
            "coE_SOM_dataInput-Mod_PPGrp", "coE_SOM_dataInput-Mod_clusterPP"
          ))
          
        })

observe(label = "somxy",{
  .schnappsEnv$defaultValues[["coE_dimSOMX"]] = input$coE_dimSOMX
  .schnappsEnv$defaultValues[["coE_dimSOMY"]] = input$coE_dimSOMY
})

### fixed function mcell_mc_hierarchy ----
mcell_mc_hierarchy_bj <- function (mc_id, mc_hc, T_gap) {
  mc = scdb_mc(mc_id)
  if (is.null(mc)) {
    stop("undefined meta cell object ", mc_id)
  }
  mc_ord = rep(0, length(mc_hc$order))
  mc_ord[mc_hc$order] = 1:length(mc_hc$order)
  n_mc = ncol(mc@mc_fp)
  lfp = log2(mc@mc_fp)
  parent = rep(-1, n_mc - 1)
  cells = rep(0, times = n_mc - 1)
  mcs = list()
  gaps = rep(0, times = n_mc - 1)
  for (i in 1:nrow(mc_hc$merge)) {
    left = mc_hc$merge[i, 1]
    right = mc_hc$merge[i, 2]
    if (left < 0) {
      cs = names(which(mc@mc == -left))
      cell_left = length(cs)
      mc_left = -left
    }
    else {
      cell_left = cells[left]
      mc_left = mcs[left]
      parent[left] = i
      gaps[left] = mc_hc$height[i] - mc_hc$height[left]
    }
    if (right < 0) {
      cs = names(which(mc@mc == -right))
      cell_right = length(cs)
      mc_right = -right
    }
    else {
      cell_right = cells[right]
      mc_right = mcs[right]
      parent[right] = i
      gaps[right] = mc_hc$height[i] - mc_hc$height[right]
    }
    mcs[[i]] = c(unlist(mc_right), unlist(mc_left))
    cells[i] = cell_right + cell_left
  }
  n_min_outcells = 300
  hits = list()
  hit_i = 1
  sup_x = c()
  for (i in which(gaps > T_gap)) {
    j = parent[i]
    mincells = cells[i] + n_min_outcells
    while (j != -1 && cells[j] < mincells) {
      j = parent[j]
    }
    if (j != -1) {
      N = cells[j]
      n = cells[i]
      mcs_in = mcs[[i]]
      mcs_out = setdiff(mcs[[j]], mcs[[i]])
      lfp_avg_in = apply(lfp[, mcs_in], 1, mean)
      lfp_min_in = apply(lfp[, mcs_in], 1, min)
      lfp_max_in = apply(lfp[, mcs_in], 1, max)
      if (length(mcs_out) > 1) {
        lfp_max_out = apply(lfp[, mcs_out], 1, max)
        lfp_min_out = apply(lfp[, mcs_out], 1, min)
        lfp_avg_out = apply(lfp[, mcs_out], 1, mean)
      }
      else {
        if (length(mcs_out) == 0) {
          message("zero length mcs out??")
          stop("boom")
        }
        lfp_max_out = lfp[, mcs_out]
        lfp_min_out = lfp[, mcs_out]
        lfp_avg_out = lfp[, mcs_out]
      }
      e_marks = tail(sort(lfp_avg_in), 20)
      sep_marks = tail(sort(lfp_min_in), 20)
      marks_gap = tail(sort(lfp_avg_in - lfp_avg_out), 
                       20)
      marks_gap_anti = head(sort(lfp_avg_in - lfp_avg_out), 
                            20)
      x_ord = mean(mc_ord[mcs[[i]]])
      sup_x = c(sup_x, x_ord)
      hits[[hit_i]] = list(marks = e_marks, min_marks = sep_marks, 
                           marks_gap = marks_gap, marks_gap_anti = marks_gap_anti, 
                           mcs = mcs[[i]], x_ord = x_ord, sup_mcs = mcs[[j]])
      hit_i = hit_i + 1
    }
  }
  hits = hits[order(sup_x)]
  return(hits)
}


observeEvent(input$updateMetaCellParameters,{
  if(input$updateMetaCellParameters <1 ) return(NULL)
  if (DEBUG) cat(file = stderr(), "metaCellSetup started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "metaCellSetup")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "metaCellSetup")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("metaCellSetup", id = "metaCellSetup", duration = NULL)
  }
  
  # if (!is.null(getDefaultReactiveDomain())) {
  #   removeNotification(id = "metaCellSetup")
  # }
  
  scEx <- isolate(scEx())
  projections <- isolate(projections())
  # input$updateSOMParameters
  sampCol <- isolate(sampleCols$colPal)
  ccols <- isolate(clusterCols$colPal)
  # coE_updateInputSOMt()
  
  if(is.null(scEx))return(NULL)
  
  # genesin <- isolate(input$coE_geneSOM)
  selectedCells <- isolate(mC_metaCell_dataInput())
  
  mC_recalc <- isolate(input$mC_recalc)
  mC_bin_for_cutoff = isolate(input$mC_bin_for_cutoff)
  mC_T_vm <- isolate(input$mC_T_vm)
  mC_T_tot <- isolate(input$mC_T_tot)
  mC_T_top3 <- isolate(input$mC_T_top3)
  mC_K <- isolate(input$mC_K)
  mC_dsamp <- isolate(input$mC_dsamp)
  mC_min_mc_size <- isolate(input$mC_min_mc_size)
  mC_resamp_n <- isolate(input$mC_resamp_n)
  mc_min_mc_size2 <- isolate(input$mc_min_mc_size2)
  mc_alpha <- isolate(input$mc_alpha)
  mc_K2 <- isolate(input$mc_K2)
  T_lfc <- isolate(input$mc_T_lfc)
  mc_T_gap <- isolate(input$mc_T_gap)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/metaCellSetup.RData", list = c(ls(), ".schnappsEnv$historyPath", ".schnappsEnv$DEBUG"))
  }
  # cp = load(file = "~/SCHNAPPsDebug/metaCellSetup.RData")
  
  library("metacell")
  metacellDir = file.path(.schnappsEnv$historyPath,"metaDataDbscEx")
  metacellFigsDir = file.path(metacellDir,"metaDataScExFigs")
  # browser()
  
  # metaDataScExFigs
  if(!dir.exists(metacellDir)) {
    dir.create(metacellDir,recursive = T)
    if(exists(".scdb")) {rm(".scdb", envir = globalenv())} # this can happen if meta cell was already calculated
    scdb_init(metacellDir, force_reinit=F)
  }else{
    if(mC_recalc){
      unlink(metacellDir, recursive=TRUE)
      dir.create(metacellDir,recursive = T)
      # if(exists(".scdb")) {rm(".scdb", envir = globalenv())}
      scdb_init(metacellDir, force_reinit=T)
    } else {
      # take what was already done
    }
  }
  if (DEBUG) cat(file = stderr(), "scdb_init\n")
  
  if (DEBUG) cat(file = stderr(), "scm_import_sce_to_mat\n")
  scm = scm_import_sce_to_mat(scEx, counts_slot = "counts")
  if (DEBUG) cat(file = stderr(), "scdb_add_mat\n")
  scdb_add_mat("metaDataDb", scm)
  #> remote mode
  #> summing up total of 0 paralog genes into 0 unique genes
  #> [1] TRUE
  if (DEBUG) cat(file = stderr(), "scdb_mat\n")
  mat = scdb_mat("metaDataDb")
  if (.schnappsEnv$DEBUG) cat(file = stderr(), paste("metaCell Dim: ",dim(mat@mat),"\n"))
  
  if(!dir.exists(metacellFigsDir)) dir.create(metacellFigsDir,recursive = T)
  scfigs_init(metacellFigsDir)
  
  
  
  if (DEBUG) cat(file = stderr(), "mcell_add_gene_stat\n")
  mcell_add_gene_stat(gstat_id="metaDataDb", mat_id="metaDataDb", force=T)
  if (DEBUG) cat(file = stderr(), "mcell_gset_filter_varmean\n")
  mcell_gset_filter_varmean(gset_id="metaDataDb_feats", gstat_id="metaDataDb", 
                            T_vm=mC_T_vm, force_new=T)
  if (DEBUG) cat(file = stderr(), "mcell_gset_filter_cov\n")
  mcell_gset_filter_cov(gset_id = "metaDataDb_feats", gstat_id="metaDataDb", 
                        T_tot=mC_T_tot, T_top3=mC_T_top3)
  
  if (DEBUG) cat(file = stderr(), "mcell_plot_umis_per_cell\n")
  mcell_plot_umis_per_cell(mat_id = "metaDataDb", bin_for_cutoff = mC_bin_for_cutoff)
  
  
  # This adds to the database a new cgraph object named test_graph. The K=100 parameter is important, 
  # as it affects the size distribution of the derived metacells. Note that constructing the graph can 
  # become computationally intensive if going beyond 20-30,000 cells. The system is currently limited by 
  # memory, and we have generated a graph on 160,000 cells on machines with 0.5TB RAM. For more modest data 
  # sets (e.g. few 10x lanes or MARS-seq experiments), things will run very quickly.
  
  if (DEBUG) cat(file = stderr(), "mcell_add_cgraph_from_mat_bknn\n")
  mcell_add_cgraph_from_mat_bknn(mat_id="metaDataDb", 
                                 gset_id = "metaDataDb_feats", 
                                 graph_id="metaDataDb_graph",
                                 K=mC_K,
                                 dsamp=mC_dsamp)
  # The next step will use the cgraph to sample five hundred metacell partitions, each covering 75% of the cells and organizing them in dense subgraphs:
  
  if (DEBUG) cat(file = stderr(), "mcell_coclust_from_graph_resamp\n")
  mcell_coclust_from_graph_resamp(
    coc_id="metaDataDb_coc500", 
    graph_id="metaDataDb_graph",
    min_mc_size=mC_min_mc_size, 
    p_resamp=0.75, n_resamp=mC_resamp_n)
  
  if (DEBUG) cat(file = stderr(), "mcell_mc_from_coclust_balanced\n")
  mcell_mc_from_coclust_balanced(
    coc_id="metaDataDb_coc500", 
    mat_id= "metaDataDb",
    mc_id= "metaDataDb_mc", 
    K=mc_K2, min_mc_size=mc_min_mc_size2, alpha=mc_alpha)
  
  if (DEBUG) cat(file = stderr(), "mcell_plot_outlier_heatmap\n")
  mcell_plot_outlier_heatmap(mc_id="metaDataDb_mc", mat_id = "metaDataDb", 
                             T_lfc=T_lfc,
                             max_genes_to_plot = 500)
  if (DEBUG) cat(file = stderr(), "mcell_mc_split_filt\n")
  mcell_mc_split_filt(new_mc_id="metaDataDb_mc_f", 
                      mc_id="metaDataDb_mc", 
                      mat_id="metaDataDb",
                      T_lfc=T_lfc, plot_mats=T)
  if (DEBUG) cat(file = stderr(), "mcell_gset_from_mc_markers\n")
  mcell_gset_from_mc_markers(gset_id="metaDataDb_markers", mc_id="metaDataDb_mc_f")
  
  if (DEBUG) cat(file = stderr(), "mcell_mc2d_force_knn\n")
  mcell_mc2d_force_knn(mc2d_id="metaDataDb_2dproj",
                       mc_id="metaDataDb_mc", 
                       graph_id="metaDataDb_graph")
  tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
  tgconfig::set_param("mcell_mc2d_width",2000, "metacell")
  if (DEBUG) cat(file = stderr(), "mcell_mc2d_plot\n")
  mcell_mc2d_plot(mc2d_id="metaDataDb_2dproj",
                  cell_outline = T)
  if (DEBUG) cat(file = stderr(), "mcell_mc_hclust_confu\n")
  mc_hc = mcell_mc_hclust_confu(mc_id="metaDataDb_mc_f", 
                                graph_id="metaDataDb_graph")
  if (DEBUG) cat(file = stderr(), "mcell_mc_hierarchy\n")
  mc_sup = mcell_mc_hierarchy_bj(mc_id="metaDataDb_mc_f",
                                 mc_hc=mc_hc, T_gap=mc_T_gap)
  if (DEBUG) cat(file = stderr(), "mcell_mc_plot_hierarchy\n")
  mcell_mc_plot_hierarchy(mc_id="metaDataDb_mc_f", 
                          graph_id="metaDataDb_graph", 
                          mc_order=mc_hc$order, 
                          sup_mc = mc_sup, 
                          width=2000, heigh=2000, min_nmc=2)
  mcUpdated(isolate(mcUpdated()) + 1)
  
})

# output$mc_metaDataDb_2dproj.2d_graph_proj <- renderUI({
#   mcUpdated()
#   div(id = "shortCut.geneSelection.click",
#       tags$img(src = "metacell/metaDataDb_2dproj.2d_graph_proj.png", 
#                width = 200, height = 200)
#   )
# })

output$mc_metaDataDb_2dproj.2d_graph_proj <- renderImage(deleteFile = T,{
  # mcell_plot_umis_per_cell(mat_id = "metaDataDb", bin_for_cutoff = 5)
  
  # react on reactive value
  mcUpdated()
  
  metacellDir = file.path(.schnappsEnv$historyPath,"metaDataDbscEx")
  metacellFigsDir = file.path(metacellDir,"metaDataScExFigs")
  pngfile <- file.path(metacellFigsDir,"metaDataDb_2dproj.2d_graph_proj.png")
  outfile <- file.path(tempdir(),"metaDataDb_2dproj.2d_graph_proj.png")
  file.copy(from = pngfile, to = outfile,overwrite = T)
  retVal <- list(
    src = normalizePath(outfile, mustWork = FALSE),
    contentType = "image/png",
    width = 100,
    height = 100,
    alt = "metaDataDb_2dproj.2d_graph_proj.png should be here"
  )
})

output$mc_outliers <- renderUI({
  # react on reactive value
  mcu = mcUpdated()
  if (mcu == 0) return(NULL)
  metacellDir = file.path(.schnappsEnv$historyPath,"metaDataDbscEx")
  metacellFigsDir = file.path(metacellDir,"metaDataScExFigs","metaDataDb_mc.outliers")
  pngFiles = dir(path = metacellFigsDir,pattern = "mc.*_out.png")
  if (length(pngFiles)<1) return(NULL)
  
  image_output_list <- 
    lapply(1:length(pngFiles),
           function(i)
           {
             imagename = paste0("mc.",i, "_out_image")
             imageOutput(imagename)             
             
           })
  do.call(tagList, image_output_list)
  
})

observeEvent(input$updateMetaCellOutliers, {
  if (DEBUG) cat(file = stderr(), "updateMetaCellOutliers\n")
  
  metacellDir = file.path(.schnappsEnv$historyPath,"metaDataDbscEx")
  metacellFigsDir = file.path(metacellDir,"metaDataScExFigs","metaDataDb_mc.outliers")
  pngFiles = dir(path = metacellFigsDir,pattern = "mc.*_out.png")
  if (length(pngFiles)<1) return(NULL)
  tdir = tempdir()
  outList = list()
  for(idx in seq(length(pngFiles))){
    pngFile = file.path(metacellFigsDir, paste0("mc",idx,"_out.png"))
    outfile = file.path(tdir, paste0("mc",idx,"_out.png"))
    file.copy(from = pngFile, to = outfile,overwrite = T)
    local({
      output[[paste0("mc.",idx, "_out_image")]] = 
        renderImage(deleteFile = F,{
          list(src = normalizePath(outfile, mustWork = FALSE),
               contentType = "image/png",
               width = 1000,
               height = 1000,
               alt = paste("mc",idx,"_out.png should be here")
          )
        }) 
    })
  }
})


# ## Show 'codebook'
# plot(sommap, type="codes", main = "Codes")
# ## Show 'component planes'
# plot(sommap, type = "property", property = sommap$codes[[1]][,1],
#      main = colnames(sommap$codes)[1])
# ## Show 'U-Matrix'
# plot(sommap, type="dist.neighbours")
