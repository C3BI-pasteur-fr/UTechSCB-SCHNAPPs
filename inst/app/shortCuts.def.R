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


# shortCutsList = list(
#   Head = ,
#   Descr = ,
#   Explain = ,
#   URL = "www/images/",
#   Func = 
# ) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


# shortCutsList = list(
#   Head = ,
#   Descr = ,
#   Explain = ,
#   URL = "www/images/",
#   Func = 
# ) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


# shortCutsList = list(
#   Head = ,
#   Descr = ,
#   Explain = ,
#   URL = "www/images/",
#   Func = 
# ) %>% add2ShortCutList(name = ".shortCut", group = group, shortCutsList)


