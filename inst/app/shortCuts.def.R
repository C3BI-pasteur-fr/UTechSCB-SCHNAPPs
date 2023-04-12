# local helper function
add2ShortCutList <- function(li, name, group, shortCutsList){
  assign(paste0(name,"Head"), li$Head, envir = parent.env(environment()))
  assign(paste0(name,"Descr"), li$Descr, envir = parent.env(environment()))
  assign(paste0(name,"Explain"), li$Explain, envir = parent.env(environment()))
  assign(paste0(name,"URL"), li$URL, envir = parent.env(environment()))
  assign(paste0(name,"Func"), li$Func, envir = parent.env(environment()))
  shortCutsList[[group]][[li$Head]] <- list(descr = li$Descr, 
                                            explain = li$Explain, 
                                            url=li$URL, 
                                            fun=li$Func)
  return(shortCutsList)
}

# shortCutsList = list(
#   Head = ,
#   Descr = ,
#   Explain = ,
#   URL = "www/images/",
#   Func = 
# ) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)

# shortCutsTab ----
shortCutsList = list()

# should be moved in a language dependent file to be sourced
.generalDescr = "Plots shown are representative and not related to the data used.
<b>Most of the plots need transformed data. Please make sure under input you selected 'calculate logcounts using SCHNAPPs'</b>"
.generalExplain = "Click on a figure to get to the that representation in the app"

# gene/cell manipulationa
.geneCellHeader = "Gene/cell manipulations"
group = .geneCellHeader
shortCutsList[[group]] = list()



shortCutsList = list(
  Head = "Gene selection",
  Descr = "gene selection panel",
  Explain = "go directly to the gene selection panel where you can modify the regular expression for filtering out genes, set the minimum number of UMIs and the genes you want to keep.",
  URL = "www/images/geneSelection.png",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "geneSelection"
    )}
  
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)

shortCutsList = list(
  Head = "Cell selection",
  Descr = "cell selection panel",
  Explain = "go directly to the cell selection panel where you can modify the regular expression for filtering out genes, set the minimum number of UMIs and the genes you want to keep.",
  URL = "www/images/cellSelection.png",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "cellSelection"
    )}
  
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)



# 2D plots
.twoDHeader = "2D plots"
group = .twoDHeader
shortCutsList[[group]] = list()

shortCutsList = list(
  Head = "Clustered violin plot",
  Descr = "Violin plot",
  Explain = "Violin plot with two categorical values: X axis has a factorial; Y axis is numerical; if color is factorial then one gets a split violin plot",
  URL = "www/images/2dViolin.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "dbCluster"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "UMI.count"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "sampleNames"
    )
    
  }
  
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)



shortCutsList = list(
  Head = "sorted values",
  Descr = "sort values",
  Explain = "sort values: useful for selecting cells based on a threshold. You can zoom-in to define an accurate threshold.",
  URL = "www/images/barCodeVsGene.png",
  Func = function(session){
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
    
  }
  
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


shortCutsList = list(
  Head = "Heatmap",
  Descr = "heatmap ",
  Explain = "heatmap with many ways to customize and select. Here, gene names is set to being empty, which triggers a one vs rest comparison based on cluster assignments. The parameters 'number of marker genes', 'min log fold change', and direction govern what is shown in the heatmap.",
  URL = "www/images/heatmap.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionAll"
    )
    # open additional options
    # set id="coExpHeatmapModule-ColNames"
    updateSelectizeInput(session, "coExpHeatmapModule-ColNames",selected = c( "dbCluster", "sampleNames"))
    updateSelectizeInput(session, "coExpHeatmapModule-orderNames",selected = c( "dbCluster", "sampleNames"))
    updateTextInput(session, "coE_heatmap_geneids",
                    value = "")
    updateBox(id="coExpHeatmapModule-heatmapAddOpt", session = session,action = "toggle")
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)

shortCutsList = list(
  Head = "2d tSNE plot",
  Descr = "2d tSNE",
  Explain = "heatmap with many ways to customize and select: Set x/y axis to tsne1/tsne2; color by anything",
  URL = "www/images/tsne2D.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "tsne1"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "tsne2"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "dbCluster"
    )
    
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


shortCutsList = list(
  Head = "dotPlotSelected",
  Descr = "dot plot",
  Explain = "Dot plot: when x/y are factorial. The size of the dots represents the number of cells.",
  URL = "www/images/dotPlotSelected.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "dbCluster"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "dbCluster"
    )
    
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


shortCutsList = list(
  Head = "histogramSelected",
  Descr = "histogram",
  Explain = "histogram of factorial values. If the color is also factorial, each bar is divided by the respective amount.",
  URL = "www/images/HistogramSelected.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "histogram"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "dbCluster"
    )
    
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


shortCutsList = list(
  Head = "histogramNormBy",
  Descr = "histogram -2",
  Explain = "histogram of factorial values. If the color is also factorial, each bar is divided by the respective amount. Now the sum of heights for each color is 1. This is achieved by setting divide Y by to 'normByCol'.",
  URL = "www/images/HistogramNormby.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "sampleNames"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "histogram"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "dbCluster"
    )
    updateSelectInput(session, "coE_selected-divideYBy",
                      selected = "normByCol"
    )
    
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)

shortCutsList = list(
  Head = "cellDensity",
  Descr = "cell density",
  Explain = "color by cell density. Especially useful when working with many cells. Only works if x/y are both continous variables.",
  URL = "www/images/cellDensity.png",
  Func =function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "coexpressionSelected"
    )
    updateSelectInput(session, "coE_selected-dimension_x",
                      selected = "tsne1"
    )
    updateSelectInput(session, "coE_selected-dimension_y",
                      selected = "tsne2"
    )
    updateSelectInput(session, "coE_selected-dimension_col",
                      selected = "cellDensity"
    )
    updateSelectInput(session, "coE_selected-divideYBy",
                      selected = "None"
    )
    
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


shortCutsList = list(
  Head = "dotPlotGeneSetsMscore",
  Descr = "Gene list dot plot with M-score",
  Explain = "Group cells by genes sets, calculate the Seurat M-score and plot in a dot plot. A GMT has to be loaded or gene sets have to be manually defined beforehand.",
  URL = "www/images/dotPlotGeneSetsMscore.png",
  Func =function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "geneSets"
    )
    updateTabsetPanel(
      session = session,
      "geneSetPlotsTabBox",
      selected = "dotPlotMscore"
    )
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)

shortCutsList = list(
  Head = "dotPlotGeneSet",
  Descr = "Gene list dot plot",
  Explain = "Group cells by genes sets, calculate the Seurat M-score and plot in a dot plot. A GMT has to be loaded or gene sets have to be manually defined beforehand.",
  URL = "www/images/dotPlotGeneSetsMscore.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "geneSets"
    )
    updateTabsetPanel(
      session = session,
      "geneSetPlotsTabBox",
      selected = "dotPlotGeneSet"
    )
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


shortCutsList = list(
  Head = "UMI histogram",
  Descr = "UMI histogram.",
  Explain = "UMI histogram.",
  URL = "www/images/UMIhistogram.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "gQC_umiHist"
    )
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)

shortCutsList = list(
  Head = "scaterPlot",
  Descr = "scater plot, summarizing the 50 most highly variable genes.",
  Explain = "Click on 'apply changes' and potentially increase the memory if you run into problems (see console).",
  URL = "www/images/scaterPlot.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "DE_scaterQC"
    )
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


shortCutsList = list(
  Head = "WIND",
  Descr = "Dendrogram showing the relation of leaves for a factorial.",
  Explain = "Dendrogram showing the relation of leaves for a factorial.",
  URL = "www/images/WINDdendrogram.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "modifyProj"
    )
    updateTabsetPanel(
      session = session,
      "modProj",
      selected = "gQC_wind"
    )
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


# 
shortCutsList = list(
  Head = "UMAP",
  Descr = "UMAP",
  Explain = "To activate UMAP projection and optimize parameters click on apply changes",
  URL = "www/images/UMAPwParameters.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "gQC_umapPlot"
    )
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)

# VariancePCs
shortCutsList = list(
  Head = "PCEigenvalues",
  Descr = "PCA Eigenvalues",
  Explain = "Show the Eigenvalues of a PCA to estimate the number of PCs needed for downstream calculations.",
  URL = "www/images/VariancePCs.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "gQC_variancePC"
    )
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


shortCutsList = list(
  Head = "ViolinCobinations",
  Descr = "Violin plot of combinations",
  Explain = "Show violin plot of combination of genes. Set max 5 genes to investigate.",
  URL = "www/images/ViolinCobinations.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "CoExpressionViolin"
    )
    updateTabsetPanel(
      session = session,
      "violinPlots",
      selected = "permViol"
    )
    updateCheckboxInput(session = session, "coE_showExpression",value = FALSE)
    updateCheckboxInput(session = session, "coE_showPermutations",value = TRUE)
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)

# PanelPlotUMI
shortCutsList = list(
  Head = "PanelPlotUMI",
  Descr = "UMI counts as panel plot",
  Explain = "Given a factorial (x-axis) show a box plot of a set of genes individually.",
  URL = "www/images/PanelPlotUMI.png",
  Func =function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "DE_panelPlot"
    )
    updateSelectInput(session, "DE_dim_x",
                      selected = "dbCluster"
    )
    updateSelectInput(session, "DE_dim_y",
                      selected = "UMI.count"
    )
    
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)

# PanelPlottSNE
shortCutsList = list(
  Head = "PanelPlotUMI",
  Descr = "tSNE projection as panel plot",
  Explain = "Plot individual genes on a tSNE projection in parallel",
  URL = "www/images/PanelPlottSNE.png",
  Func =function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "DE_panelPlot"
    )
    updateSelectInput(session, "DE_dim_x",
                      selected = "tsne1"
    )
    updateSelectInput(session, "DE_dim_y",
                      selected = "tsne2"
    )
    
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


# ViolinPlotGeneCount
shortCutsList = list(
  Head = "ViolinPlotGeneCount",
  Descr = "Violin plot of number of given expressed.",
  Explain = "For each cell the number of genes that are expressed within a given range is counted and plotted in a violin plot.",
  URL = "www/images/ViolinPlotGeneCount.png",
  Func =function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "DE_panelPlot"
    )
    updateTabsetPanel(
      session = session,
      "violinPlots",
      selected = "permViol"
    )
    updateCheckboxInput(session = session, "coE_showExpression",value = FALSE)
    updateCheckboxInput(session = session, "coE_showPermutations",value = FALSE)
    
    }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)

###################################

# VolcanoPlot
shortCutsList = list(
  Head = "VolcanoPlot",
  Descr = "Volcano plot of DGE",
  Explain = "Volcano plot of differentially expressed gene. To show this a differential gene expression analysis has to be performed. The link goes to this page. Cells selected on the left panel will be compared to cells selected in the right panel. If no cells are seleted in the right panel all other cells will be used. The display of the volcano plot needs a few seconds to build.",
  URL = "www/images/VolcanoPlot.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "sCA_dge"
    )
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


# 3DExpressionGene
shortCutsList = list(
  Head = "3DExpressionGene",
  Descr = "3D expression of gene",
  Explain = "Show the expression intensity of a gene or the sum of a list of genes in a 3D representation.",
  URL = "www/images/3DExpressionGene.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "DE_expression"
    )
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


###################################


# 3DUMAP
shortCutsList = list(
  Head = "3DUMAP",
  Descr = "3D UMAP",
  Explain = "given that the UMAP has been calulated, show it in 3D",
  URL = "www/images/3DUMAP.png",
  Func = function(session){
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
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)

# 3DClusterTsne
shortCutsList = list(
  Head = "3DClusterTsne",
  Descr = "tSNE in 3D",
  Explain = "show the tSNE projection in 3D",
  URL = "www/images/3DClusterTsne.png",
  Func = function(session){
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "gQC_tsnePlot"
    )
    updateSelectizeInput(session, "gQC_dim3D_x",
                         selected = "tsne1"
    )
    updateSelectizeInput(session, "gQC_dim3D_y",
                         selected = "tsne2"
    )
    updateSelectizeInput(session, "gQC_dim3D_z",
                         selected = "tsne3"
    )
    updateSelectizeInput(session, "gQC_col3D",
                         selected = "dbCluster"
    )
  }
) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)





# shortCutsList = list(
#   Head = ,
#   Descr = ,
#   Explain = ,
#   URL = "www/images/",
#   Func = 
# ) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


