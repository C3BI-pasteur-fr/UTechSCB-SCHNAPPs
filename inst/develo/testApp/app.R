options(shiny.sanitize.errors = FALSE)
library(SCHNAPPs)
.schnappsEnv <- new.env(parent=emptyenv())
#
#schnappsLite(data = "epdc.rn-sham2.v2.lite.RData", DEBUG = T, historyPath = "history")
#schnappsLite(data = "rnAllEPDC.lite.RData", DEBUG = T, historyPath = "history")
defaultValueMultiGenes = "wt1,pdgfrb,agtr1a,col3a1,col1a1,postn,tbx18,scx,npr1,vegfa,npr2,notch2,tagln,acta2,ctgf, Rack1, Bmp4,  Eef2, col3a1,col1a1,tgfb3,postn"
defaultValueSingleGene = "wt1"
localContributionDir = "./tmp"
assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = .schnappsEnv)


defaultValues = list()
defaultValues[["coEtgMinExpr"]] = 100
defaultValues[["coE_selected-groupNames"]] = "plot"
defaultValues[["coE_heatmap_geneids"]] = defaultValueMultiGenes
defaultValues[["coExpHeatmapModule-heatmapMinMaxValue"]] = c(0,25)
defaultValues[["coE_selected-dimension_x"]] = "UMAP1"
defaultValues[["coE_selected-dimension_y"]] = "UMAP2"
defaultValues[["coE_selected-dimension_col"]] = "sampleNames"
defaultValues[["coE_geneGrpVioIds2"]] = defaultValueMultiGenes
defaultValues[["coE_dimension_xVioiGrp2"]] = c("Seurat", "sampleNames")
defaultValues[["coE_geneGrpVioIds"]]  = defaultValueMultiGenes
defaultValues[["coE_dimension_xVioiGrp"]] = "Seurat"
defaultValues[["alluiv1"]] = "sampleNames"
defaultValues[["alluiv2"]] = "Seurat"
defaultValues[["DE_Exp_dataInput-Mod_clusterPP"]] = "Seurat"
defaultValues[["DE_Exp_dataInput-Mod_PPGrp"]] = as.character(c(0,1,2,3,4,5,6,7,8))
defaultValues[["DE_gene_id"]] = defaultValueSingleGene
defaultValues[["DE_expclusters-dimension_x"]] = "UMAP1"
defaultValues[["DE_expclusters-dimension_y"]] = "UMAP2" 
defaultValues[["DE_expclusters-dimension_col"]] = "sampleNames" 
defaultValues[["DE_PanelPlotCellSelection-Mod_clusterPP"]] = "Seurat"
defaultValues[["DE_PanelPlotCellSelection-Mod_PPGrp"]] = as.character(c(0,1,2,3,4,5,6,7,8))
defaultValues[["DE_dim_x"]] = "UMAP1"
defaultValues[["DE_dim_y"]] = "UMAP2"
defaultValues[["sCA_dataInput-Mod_clusterPP"]] = "Seurat"
defaultValues[["sCA_dataInput-Mod_PPGrp"]] = as.character(c(0,1,2,3,4,5,6,7,8))
defaultValues[["sCA_subscluster_x1"]] = "UMAP1"
defaultValues[["sCA_subscluster_y1"]] = "UMAP2"

devscShinyApp=FALSE
rm("packagePath")
app = schnappsLite(defaultValues = defaultValues, 
                   # data = "inst/develo/testApp/HPVC.lite.RData", 
                   data = "/Users/bernd/Rstudio/UTechSCB-SCHNAPPs/inst/develo/testApp/HPVC.lite.RData", 
                   defaultValueMultiGenes = "wt1,pdgfrb,agtr1a,col3a1,col1a1,postn,tbx18,scx,npr1,vegfa,npr2,notch2,tagln,acta2,ctgf, Rack1, Bmp4,  Eef2, col3a1,col1a1,tgfb3,postn" , 
             defaultValueSingleGene = "wt1", 
             DEBUG = T)
# runApp(app)
