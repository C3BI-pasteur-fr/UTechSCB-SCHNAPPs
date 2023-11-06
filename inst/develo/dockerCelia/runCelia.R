library(dplyr)

.schnappsEnv <- new.env(parent=emptyenv())
localContributionDir = NULL
defaultValueSingleGene = "NANOG" # CD52
defaultValueMultiGenes = "RTN4, NEUROD1, ROBO1, NRG1, DLL1, SLIT2, NRP2, B2M, FLNA, PTN, HES1, FGF13, SOX1, CDK5RAP2, GPI, PAX6, NRCAM, DLG4, TJP1, NANOG"
defaultValueRegExGene = "" # tip: '^CD7$|^KIT$; genes with min expression
DEBUG = T
DEBUGSAVE = F
historyPath = NULL

AllowClustering = F

assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueSingleGene", defaultValueSingleGene, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueMultiGenes", defaultValueMultiGenes, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueRegExGene", defaultValueRegExGene, envir = .schnappsEnv)
assign(".SCHNAPPs_DEBUG", DEBUG, envir = .schnappsEnv)
assign(".SCHNAPPs_DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
assign("DEBUG", DEBUG, envir = .schnappsEnv)
assign("DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
assign("historyPath", historyPath, envir = .schnappsEnv)

defaultValues = list()
defaultValues[["coEtgMinExpr"]] = 100
defaultValues[["coE_selected-groupNames"]] = "plot"
defaultValues[["coE_heatmap_geneids"]] = defaultValueMultiGenes
defaultValues[["coExpHeatmapModule-heatmapMinMaxValue"]] = c(0,25)
defaultValues[["coE_selected-dimension_x"]] = "tsne1"
defaultValues[["coE_selected-dimension_y"]] = "tsne2"
defaultValues[["coE_selected-dimension_col"]] = "types.cell"
defaultValues[["coE_geneGrpVioIds2"]] = defaultValueMultiGenes
defaultValues[["coE_dimension_xVioiGrp2"]] = c("dbCluster", "types.cell")
defaultValues[["coE_geneGrpVioIds"]]  = defaultValueMultiGenes
defaultValues[["coE_dimension_xVioiGrp"]] = "types.cell"
defaultValues[["alluiv1"]] = "types.cell"
defaultValues[["alluiv2"]] = "dbCluster"
defaultValues[["DE_Exp_dataInput-Mod_clusterPP"]] = "dbCluster"
defaultValues[["DE_Exp_dataInput-Mod_PPGrp"]] = as.character(c(0,1,2,3,4,5,6))
defaultValues[["DE_gene_id"]] = defaultValueSingleGene
defaultValues[["DE_expclusters-dimension_x"]] = "tsne1"
defaultValues[["DE_expclusters-dimension_y"]] = "tsne2" 
defaultValues[["DE_expclusters-dimension_col"]] = "types.cell" 
defaultValues[["DE_PanelPlotCellSelection-Mod_clusterPP"]] = "sampleNames"
defaultValues[["DE_PanelPlotCellSelection-Mod_PPGrp"]] = "All"
defaultValues[["DE_dim_x"]] = "types.cell"
defaultValues[["DE_dim_y"]] = "UMI.count"
defaultValues[["DE_nCol"]] = 3
defaultValues[["sCA_dataInput-Mod_clusterPP"]] = "types.cell"
defaultValues[["sCA_dataInput-Mod_PPGrp"]] = c("stem.like", "mature", "NPC")
defaultValues[["sCA_subscluster_x1"]] = "types.cell"
defaultValues[["sCA_subscluster_y1"]] = "UMI.count"
defaultValues[["DEBUGSAVE"]] = FALSE
defaultValues[["gQC_windProj"]] = "dbCluster"

defaultValues[["sCA_dgeRadioButton"]] = "seurat:t-test"

assign("defaultValues", defaultValues, envir = .schnappsEnv)


# will be set during sourcing, but we need to define them, otherwise there will be a warning
scShinyUI <- NULL
scShinyServer <- NULL

packagePath <- find.package("SCHNAPPs", lib.loc = NULL, quiet = TRUE) %>% paste0("/app/")

cat(file = stderr(), packagePath,"\n")

devscShinyApp = F   # TRUE = look for local sources
# packagePath <- "inst/app"

source(paste0(packagePath, "/serverFunctions.R"), local = TRUE)
source(paste0(packagePath, "/server.R"), local = TRUE)
source(paste0(packagePath, "/ui.R"), local = T)
source(paste0(packagePath, "/server-lite.R"), local = TRUE)
source(paste0(packagePath, "/ui-lite.R"), local = T)

# load data
maxCells = 300000
# file = "inst/develo/testApp/HPVC.lite.RData"
# file = "/Users/bernd/Rstudio/UTechSCB-SCHNAPPs/inst/develo/dockerCelia/celia.schnapps.RData"
file = "/srv/shiny-server/schnapps/celia.schnapps.RData"
# data = "~/Rstudio/UTechSCB-SCHNAPPs/data/scExLite.RData"
assign(".SCHNAPPs_LiteData", loadLiteData(fileName = file), envir = .schnappsEnv)

# GMT DATA !!!!
# if(!isEmpty(gmtFiles)){
#   fileInfo = gmtFiles %>% file.info()
#   latestFile = fileInfo %>% pull("ctime") %>% order()  %>% last()
#   tempEnv =  new.env(parent=emptyenv())
#   cp = load(rownames(fileInfo)[latestFile], envir = tempEnv)
#   if("gmtData" %in% cp) {
#     gmtData(tempEnv$gmtData$gmtData) 
#   }
# }
gmtData = list()
gmtData[["Mature Neurons"]]$name = 'Mature Neurons'
gmtData[["Mature Neurons"]]$desc = "Mature Neurons"
gmtData[["Mature Neurons"]]$genes = c("ELAVL4", "CUX2", "TUBB3", "MAP2", "RTN4", "NEUROD1", "ROBO1", "NRG1", "DLL1")

gmtData[["neuronal precursor"]]$name = 'neuronal precursor'
gmtData[["neuronal precursor"]]$desc = "neuronal precursor"
gmtData[["neuronal precursor"]]$genes = c("NES", "FABP7", "SLIT2", "NRP2", "B2M", "FLNA", "PTN", "HES1", "FGF13", "SOX1", "CDK5RAP2", "PAX6", "SOX2")

gmtData[["precusor cells"]]$name = 'precusor cells'
gmtData[["precusor cells"]]$desc = "precusor cells"
gmtData[["precusor cells"]]$genes = c("TJP1", "NANOG", "POU5F1")

gmtData[["neuro epitel"]]$name = 'neuro epitel'
gmtData[["neuro epitel"]]$desc = "NEP"
gmtData[["neuro epitel"]]$genes = c("PAX5", "EMX1", "DLX2", "PAX2")

assign(".SCHNAPPs_GMTData", gmtData, envir = .schnappsEnv)


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

shiny::addResourcePath(
  prefix = "www",
  directoryPath = "~/www/"
)


app <- shinyApp(ui = scShinyUI, server = scShinyServer)
# options(shiny.reactlog=TRUE) 
runApp(app, port = 3838,host = "0.0.0.0")
# 

