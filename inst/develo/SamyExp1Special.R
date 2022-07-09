# require(purrr)
# require(plotly)
# require(Seurat)
# require(shiny)
# 
# ui <- fluidPage(
#   fluidRow(
#     column(width = 2,
#   numericInput(inputId = "lowFeature", label = "lowFeature", value = 2000)),
#   column(width = 2,
#   numericInput(inputId = "pctMT", label = "pctMT", value  = 5)),
#   column(width = 2,
#   numericInput(inputId = "highFeature",  label = "highFeature", value = 7000)),
#   column(width = 2,
#   numericInput(inputId = "nMT",  label = "nMT", value = 5000)),
#   column(width = 2,
#   numericInput(inputId = "umiMTRatio",  label = "umiMTRatio" , value= 100))
#   ),
#   
#   actionButton(inputId = "button", label = "apply"),
#   plotlyOutput("plotHist"),
#   fluidRow(column(width = 6,
#                   plotlyOutput("plotUmap")),
#                   column(width = 6,
#                          plotlyOutput("plotUmap2")) ),
#   fluidRow(
#     column(width = 12,
#            textOutput("numbers"))
#   )
# )
# 
# 
# server <- function(input, output, session) {
#   load("~/Rstudio/scShinyHubData/samy.1.RData")
#   # load("~/Rstudio/UTechSCB-SCHNAPPs/example.SCEX.RData")
#   # sc = scEx[,1:200]
#   # save(file = "example.SCEX.RData", list = c("sc"))
#   
#   genesOfInterest = list(AT1 = c("AGER", "EMP2", "CAV1", "CLDN18", "CAV2"), 
#                          AT2 = c("ALOX15B", "LRRK2", "ROS1", "SFTPA1", "CSF3R"), 
#                          BASAL = c("MIR205HG", "KRT5", "EYA2", "CYP24A1", "KRT17"),
#                          Ciliated = c( "KCNE1B", "TCTEX1D4", "ANKRD66", "DYDC2", "C22orf15"), 
#                          Club = c("SCGB3A2", "CYP2B7P", "ITGA9", "MGP", "ATL2"), 
#                          Goblet = c("BPIFB1", "MUC5B", "LTF", "SCGB3A1", "SCGB1A1"), 
#                          PNEC = c("GRP", "CHGA", "SCG2", "ASCL1", "CALCA"), 
#                          Ionocyte = c("ASCL3", "FOXI1", "ATP6V1G3", "BSND", "HEPACAM2"), 
#                          Mesothelial = c("CALB2", "VCAM1", "MEDAG", "GAS1", "HAS1"), 
#                          AberrantBasaloid = c("PRSS2", "CDH2", "CAMK2N1", "MMP7", "EPHB2"))
#   
#   
#   seuratdataOrg <- UpdateSeuratObject(as.Seurat(scEx, assay = "RNA", data=NULL))
#   
#   # seuratdata <- CreateSeuratObject(
#   #   counts = counts, min.cells = min.cells, min.features = min.genes,
#   #   project = "test"
#   # )
#   
#   seuratdataOrg[["nCount_RNA"]] =  colSums(seuratdataOrg, slot= 'counts')
#   seuratdataOrg[["nFeature_RNA"]] =  colSums(GetAssayData(seuratdataOrg, slot= 'counts') >0 )
#   
#   
#   mito.genes <- grep(pattern = "^MT-", ignore.case = TRUE, x = rownames(x = seuratdataOrg), value = TRUE)
#   RPS.genes <- grep(pattern = "^RPS", ignore.case = TRUE, x = rownames(x = seuratdataOrg), value = TRUE)
#   
#   seuratdataOrg[["nMito_RNA"]] =  colSums(seuratdataOrg[mito.genes,], slot= 'counts')
#   seuratdataOrg[["nRps_RNA"]] =  colSums(seuratdataOrg[RPS.genes,], slot= 'counts')
#   seuratdataOrg[["umiMT.ratio"]] = seuratdataOrg[["nCount_RNA"]]/seuratdataOrg[["nMito_RNA"]]
#   
#   
#   potentialFactorials = c()
#   potentialNumericals = c()
#   for (cn in colnames(seuratdataOrg@meta.data)) {
#     nUniq = length(unique(seuratdataOrg@meta.data[,cn]))
#     if (nUniq <40 & nUniq >1) {
#       potentialFactorials = c(potentialFactorials, cn)
#     }
#     if (is(seuratdataOrg@meta.data[,cn], "numeric")){
#       potentialNumericals = c(potentialNumericals, cn)
#     }
#   }
#   
#   if (is(genesOfInterest, "list")) {
#     for (idx in 1:length(genesOfInterest)) {
#       gio = genesOfInterest[[idx]]
#       gio = gio[gio %in% rownames(scEx)]
#       genesOfInterest[[idx]] = gio
#     }
#   } else {
#     genesOfInterest = genesOfInterest[genesOfInterest %in% rownames(scEx)]
#   }
#   seuratdataOrg = PercentageFeatureSet(object = seuratdataOrg, pattern = "^mt-|^MT-", col.name = "percent.mt")
#   
#   
#   #############################################
#   
#   
#   
#   seuratdata =  reactiveVal(seuratdataOrg)
#   
#   observeEvent(input$button, {
#     
#     sc <- subset(seuratdataOrg, subset = 
#                            nFeature_RNA > input$lowFeature & 
#                            nFeature_RNA < input$highFeature & 
#                            percent.mt < input$pctMT & 
#                            nMito_RNA < input$nMT &
#                            umiMT.ratio < input$umiMTRatio 
#                          # &
#                          #   phases == "G1"
#     )
#     all.genes <- rownames(sc)
#     sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
#     sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
#     sc <- ScaleData(sc, features = all.genes)
#     
#     sc <- RunPCA(sc, ndims.print = 1, nfeatures.print = 1)
#     sc <- FindNeighbors(sc, dims = 1:30)
#     sc <- FindClusters(sc, resolution = 0.5)
#     
#     sc <- RunUMAP(sc, dims = 1:30)
#     
#     seuratdata(sc)
#     })
#   
#   output$plotHist <- renderPlotly({
#     binSize = 20
#     scols = c("#00AA00", "#AA0000")
#     
#     sc = seuratdata()
#     tabl = table(sc[["cellTypes"]][,1], sc[["sampleNames"]][,1])
#     sc
#     tabl[,"SAC"] =  tabl[,"SAC"] / sum( tabl[,"SAC"]) * 100
#     tabl[,"WNA"] =  tabl[,"WNA"] / sum( tabl[,"WNA"]) * 100
#     fig <- plotly::plot_ly( alpha = 1)
#     # browser()
#     lev = unique(sc[["sampleNames"]][,1])
#     for (idx in seq_along(lev)) {
#       fig <- fig %>% add_trace(
#         type = 'bar', color = I(scols[idx]), 
#         x = rownames(tabl),
#         y = tabl[,lev[idx]],
#         name = lev[idx]
#       )
#     }
#     fig <- fig %>% layout(
#       barmode="group",
#       title = paste("percent sample"),
#       bargap=0.1)
#     cat(file = stderr(), "fig1 done\n")
#     fig
#     
#   })
#   
#   output$plotUmap = renderPlotly({
#     seuratdata = seuratdata()
#     cn = "cellTypes"
#     DimPlot(seuratdata, reduction = "umap", group.by = cn, label = T) + ggtitle(cn)
#   }
#   )
#   output$plotUmap2 = renderPlotly({
#     seuratdata = seuratdata()
#     cn = "sampleNames"
#     DimPlot(seuratdata, reduction = "umap", group.by = cn, label = T) + ggtitle(cn)
#   }
#   )
#   output$numbers <- renderPrint({ 
#     sc = seuratdata()
#     str(table(sc[["sampleNames"]][,1]))
#   })
# }
# 
# 
# 
# 
# shinyApp(ui, server)
# 
# 
# 
# 
# # marker = list(color = scols))
# 
# 
# 
# 
# 
