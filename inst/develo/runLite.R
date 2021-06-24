.schnappsEnv <- new.env(parent=emptyenv())
localContributionDir = "~/Rstudio/SCHNAPPsContributions/"
localContributionDir = ""
defaultValueSingleGene = "itgae" # CD52
defaultValueMultiGenes = "wt1,pdgfrb,agtr1a,col3a1,col1a1,postn,tbx18,scx,npr1,vegfa,npr2,notch2,tagln,acta2,ctgf, Rack1, Bmp4,  Eef2, col3a1,col1a1,tgfb3,postn"
# defaultValueMultiGenes = "CD52, S100A9, S100A4" # itgae, cd69, itga1" # CD52, S100A9, S100A4
# defaultValueMultiGenes = "prf1, Gzmb, IFNG, PDCD1, HAVCR2, LAG3, TSC22D3,ZFP36L2"
defaultValueRegExGene = "" # tip: '^CD7$|^KIT$; genes with min expression
DEBUG = T
DEBUGSAVE = F
historyPath = "~/Rstudio/Schnapps/history"
# historyPath = NULL
AllowClustering = F

.schnappsEnv 

assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueSingleGene", defaultValueSingleGene, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueMultiGenes", defaultValueMultiGenes, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueRegExGene", defaultValueRegExGene, envir = .schnappsEnv)
assign(".SCHNAPPs_DEBUG", DEBUG, envir = .schnappsEnv)
assign(".SCHNAPPs_DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
assign("DEBUG", DEBUG, envir = .schnappsEnv)
assign("DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
assign("historyPath", historyPath, envir = .schnappsEnv)
# assign("historyFile", historyFile, envir = .schnappsEnv)

defaultValueMultiGenes = "wt1,pdgfrb,agtr1a,col3a1,col1a1,postn,tbx18,scx,npr1,vegfa,npr2,notch2,tagln,acta2,ctgf, Rack1, Bmp4,  Eef2, col3a1,col1a1,tgfb3,postn"
defaultValueSingleGene = "wt1"

defaultValues = list()
defaultValues[["coEtgMinExpr"]] = 100
defaultValues[["coE_selected-groupNames"]] = "plot"
defaultValues[["coE_heatmap_geneids"]] = defaultValueMultiGenes
defaultValues[["coExpHeatmapModule-heatmapMinMaxValue"]] = c(0,25)
defaultValues[["coE_selected-dimension_x"]] = "UMAP1"
defaultValues[["coE_selected-dimension_y"]] = "UMAP2"
defaultValues[["coE_selected-dimension_col"]] = "sampleNames"
defaultValues[["coE_geneGrpVioIds2"]] = defaultValueMultiGenes
defaultValues[["coE_dimension_xVioiGrp2"]] = c("seuratCluster", "sampleNames")
defaultValues[["coE_geneGrpVioIds"]]  = defaultValueMultiGenes
defaultValues[["coE_dimension_xVioiGrp"]] = "seuratCluster"
defaultValues[["alluiv1"]] = "sampleNames"
defaultValues[["alluiv2"]] = "seuratCluster"
defaultValues[["DE_Exp_dataInput-Mod_clusterPP"]] = "seuratCluster"
defaultValues[["DE_Exp_dataInput-Mod_PPGrp"]] = as.character(c(0,1,2,3,4,5,6))
defaultValues[["DE_gene_id"]] = defaultValueSingleGene
defaultValues[["DE_expclusters-dimension_x"]] = "UMAP1"
defaultValues[["DE_expclusters-dimension_y"]] = "UMAP2" 
defaultValues[["DE_expclusters-dimension_col"]] = "sampleNames" 
defaultValues[["DE_PanelPlotCellSelection-Mod_clusterPP"]] = "seuratCluster"
defaultValues[["DE_PanelPlotCellSelection-Mod_PPGrp"]] = as.character(c(0,1,2,3,4,5,6))
defaultValues[["DE_dim_x"]] = "UMAP1"
defaultValues[["DE_dim_y"]] = "UMAP2"
defaultValues[["DE_nCol"]] = 3
defaultValues[["sCA_dataInput-Mod_clusterPP"]] = "seuratCluster"
defaultValues[["sCA_dataInput-Mod_PPGrp"]] = as.character(c(0,1,2,3,4,5,6))
defaultValues[["sCA_subscluster_x1"]] = "UMAP1"
defaultValues[["sCA_subscluster_y1"]] = "UMAP2"
defaultValues[["DEBUGSAVE"]] = DEBUGSAVE


assign("defaultValues", defaultValues, envir = .schnappsEnv)


# will be set during sourcing, but we need to define them, otherwise there will be a warning
scShinyUI <- NULL
scShinyServer <- NULL
# packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/")

devscShinyApp = T   # TRUE = look for local sources
packagePath <- "inst/app"

source(paste0(packagePath, "/serverFunctions.R"), local = TRUE)
source(paste0(packagePath, "/server.R"), local = TRUE)
source(paste0(packagePath, "/ui.R"), local = T)
source(paste0(packagePath, "/server-lite.R"), local = TRUE)
source(paste0(packagePath, "/ui-lite.R"), local = T)

# load data
maxCells = 300000
# file = "inst/develo/testApp/HPVC.lite.RData"
file = "../scShinyHubData/mca_Seurat_afterClust_CtrMem.schnapps.RData"
# data = "~/Rstudio/UTechSCB-SCHNAPPs/data/scExLite.RData"
assign(".SCHNAPPs_LiteData", loadLiteData(fileName = file), envir = .schnappsEnv)


nCells = length(colnames(.schnappsEnv$.SCHNAPPs_LiteData$scEx))
if (nCells > maxCells){
  cellIdx = unique(sort(base::sample(nCells, maxCells)))
  cells2keep = colnames(.schnappsEnv$.SCHNAPPs_LiteData$scEx)[cellIdx]
  for (na in names(.schnappsEnv$.SCHNAPPs_LiteData)){
    if (na == "pca") {
      .schnappsEnv$.SCHNAPPs_LiteData[[na]]$x = .schnappsEnv$.SCHNAPPs_LiteData[[na]]$x[cellIdx,]
      next()
    }
    if (na %in%  c("sampleCol", "clusterCol")) next()
    if (na == "tsne") {
      .schnappsEnv$.SCHNAPPs_LiteData[[na]] = .schnappsEnv$.SCHNAPPs_LiteData[[na]][cellIdx,]
      next()
    }
    if ( all(cells2keep %in% colnames(.schnappsEnv$.SCHNAPPs_LiteData[[na]]))){
      .schnappsEnv$.SCHNAPPs_LiteData[[na]] = .schnappsEnv$.SCHNAPPs_LiteData[[na]][,cells2keep]
    } else {
      if(all(cells2keep %in% rownames(.schnappsEnv$.SCHNAPPs_LiteData[[na]]))){
        .schnappsEnv$.SCHNAPPs_LiteData[[na]] = .schnappsEnv$.SCHNAPPs_LiteData[[na]][cells2keep,]
      } else {
        if(is.null(rownames(.schnappsEnv$.SCHNAPPs_LiteData[[na]]))) {
          .schnappsEnv$.SCHNAPPs_LiteData[[na]] = .schnappsEnv$.SCHNAPPs_LiteData[[na]][cellIdx]
        }else{
          cat(file = stderr(), paste("couldn't ", na, "\n"))
        }
      }
    }
  }
  
}

if (is.null(.schnappsEnv$".SCHNAPPs_LiteData")) {
  # error loading
  stop(".schnappsEnv$.SCHNAPPs_LiteData not given\n")
}

app <- shinyApp(ui = scShinyUI, server = scShinyServer)
# options(shiny.reactlog=TRUE) 
runApp(app)
# 

