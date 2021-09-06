options(shiny.sanitize.errors = FALSE)
library(SCHNAPPs)
#
#schnappsLite(data = "epdc.rn-sham2.v2.lite.RData", DEBUG = T, historyPath = "history")
#schnappsLite(data = "rnAllEPDC.lite.RData", DEBUG = T, historyPath = "history")
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
defaultValues[["coE_dimension_xVioiGrp2"]] = list("seuratCluster", "sampleNames")
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


schnapps(defaultValues = defaultValues, 
         defaultValueMultiGenes = "wt1,pdgfrb,agtr1a,col3a1,col1a1,postn,tbx18,scx,npr1,vegfa,npr2,notch2,tagln,acta2,ctgf, Rack1, Bmp4,  Eef2, col3a1,col1a1,tgfb3,postn" , defaultValueSingleGene = "wt1", 
         localContributionDir = "/SCHNAPPsContributions/", DEBUG = T, historyPath = "history", port = 3839)

