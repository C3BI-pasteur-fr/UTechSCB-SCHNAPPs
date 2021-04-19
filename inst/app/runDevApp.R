  #' this is used to run the app without installing it.
  #'
  #'
  library(reactlog)
  # if (!exists(".schnappsEnv")) {
  .schnappsEnv <- new.env(parent=emptyenv())
  # }
  # 
  library(future)
  plan("multiprocess", workers = 6)
  
  library("BiocParallel")
  register(MulticoreParam(6))
  
  localContributionDir = "~/Rstudio/SCHNAPPsContributions/"
  # localContributionDir = ""
  defaultValueSingleGene = "CDH2" # CD52
  defaultValueMultiGenes = "CAV1, MIR205HG, KCNE1B, ANKRD66, SCGB3A2, SCGB3A1, CDH2" # itgae, cd69, itga1" # CD52, S100A9, S100A4
  # defaultValueMultiGenes = "prf1, Gzmb, IFNG, PDCD1, HAVCR2, LAG3, TSC22D3,ZFP36L2"
  defaultValueRegExGene = "" # tip: '^CD7$|^KIT$; genes with min expression
  DEBUG = T
  DEBUGSAVE = F
  historyPath = "/Volumes/Oct2020/RStudio_history/"
  historyPath = NULL
  
  assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = .schnappsEnv)
  assign(".SCHNAPPs_defaultValueSingleGene", defaultValueSingleGene, envir = .schnappsEnv)
  assign(".SCHNAPPs_defaultValueMultiGenes", defaultValueMultiGenes, envir = .schnappsEnv)
  assign(".SCHNAPPs_defaultValueRegExGene", defaultValueRegExGene, envir = .schnappsEnv)
  assign(".SCHNAPPs_DEBUG", DEBUG, envir = .schnappsEnv)
  assign(".SCHNAPPs_DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
  assign("localContributionDir", localContributionDir, envir = .schnappsEnv)
  assign("defaultValueSingleGene", defaultValueSingleGene, envir = .schnappsEnv)
  assign("defaultValueMultiGenes", defaultValueMultiGenes, envir = .schnappsEnv)
  assign("defaultValueRegExGene", defaultValueRegExGene, envir = .schnappsEnv)
  assign("DEBUG", DEBUG, envir = .schnappsEnv)
  assign("DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
  assign("historyPath", historyPath, envir = .schnappsEnv)
  ls(.schnappsEnv)
  
  
  # # Scran parameters
  defaultValues = list()
  
  
  # Seurat parameters
  # defaultValues = list()
  
  # defaultValues[["selectIds"]] = ""
  defaultValues[["pcaN"]] = 500
  defaultValues[["pcaScale"]] = TRUE
  defaultValues[["sampleInput"]] =FALSE
  defaultValues[["hvgSelection"]] = "vst"
  # defaultValues[["alluiv1"]] = "seurartCluster"
  defaultValues[["alluiv2"]] = "dbCluster"
  # defaultValues[["tabsetCluster"]] = "seurat_Clustering"
  defaultValues[["minGenesGS"]] = 1
  defaultValues[["minGenes"]] = 1
  defaultValues[["maxGenes"]] = 50000
  defaultValues[["seurClustDims"]] = 15
  defaultValues[["seurClustk.param"]] = 15
  defaultValues[["cellPatternRM"]] = "-s1|-s2"
  defaultValues[["gQC_binSize"]] = 200
  defaultValues[["selectIds"]] = "^MT-|^RP|^MRP|MALAT1|B2M|EEF1A1"
  defaultValues[["selectIds"]] = ""
  defaultValues[["whichscLog"]] = "calcLog"
  defaultValues[["whichscLog"]] = "disablescEx_log"
  
  defaultValues[["gQC_um_n_neighbors"]] = 20 
  defaultValues[["gQC_um_spread"]] = 6 
  defaultValues[["gQC_um_local_connectivity"]] = 2
  # defaultValues[["useSeuratPCA"]] = TRUE
  
  # defaultValues[["DE_panelplotids"]] = c("CD8A", "CD4", "CD8B", "FCER1G", "CCR7", "GZMK", "FoxP3", "GZMK", "GZMB", "CCR7", "LEF1", "CCL5", "VIM", "CCL5", "TCF7", "NKG7", "LGALs1", "NKG7", "SELL", "CST7", "ANXA2", "CST7", "IL7R", "HLA-DRB1", "KLRG1", "CD27", "CTLA4A")
  assign("defaultValues", defaultValues, envir = .schnappsEnv)
  
  devscShinyApp = TRUE
  packagePath <<- "inst/app"
  source(paste0(packagePath,  "/ui.R"))
  source(paste0(packagePath,  "/server.R"))
  
  app <- shinyApp(ui = scShinyUI, server = scShinyServer, enableBookmarking = "server")
  options(shiny.reactlog=TRUE)
  
  # 
  # 
  options(keep.source=TRUE)
  p <- profvis::profvis({
    runApp(app)
  })
  htmlwidgets::saveWidget(p, '~/profvis1.html')
  
  # schnapps(Ï
  # defaultValueMultiGenes = "IL7R, CCR7,CD14, LYZ ,IL7R, S100A4,MS4A1 ,CD8A,FCGR3A, MS4A7 ,GNLY, NKG7,FCER1A, CST3,PPBP",
  # defaultValueSingleGene = "MS4A1", DEBUG=TRUE
  # )
  # 
  # 
  # sctkEx = SCtkExperiment(assays=list(counts=as.matrix(assays(scEx)[['counts']]),
  #                                     logcounts = as.matrix(assays(scEx)[['logcounts']])),
  #                         colData = colData(scEx),
  #                         rowData = rowData(scEx))
  # singleCellTK(inSCE = sctkEx)
  # 
  # 
  # 
  # library(SCHNAPPs)
  # schnapps(DEBUG = T, historyPath = "/Volumes/Oct2020/RStudio/history/celia/")
  # 