controlbarContext <- shinydashboardPlus::dashboardControlbar(
  width = 500, collapsed = T,
  id = "controlbar",
  controlbarMenu(
    id = "menu",
    # Workflow 1 ----
    controlbarItem(
      "counts 2 transformed data",
      ## Loading data ----
      shinydashboardPlus::box(
        title = "Loading data",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 div(h4(a(id = "wkfl1.LoadData.click", "Load data"))),
                 div(h6(HTML("[<em>Input</em>]"))),
                 div(h5("Load the data count matrix, annotations, and GMT information without activating the Log-counts. ")),
                 div(h5("The count matrix can be in text form (csv/tsv), an RData file with a singlecellExpression object.")),
                 div(h5("Since the first steps of QC will use only the raw counts the potentially expensive normalization 
                        procedure is not needed.")),
                 div(h5("For this workflow it is assumed that an RData file with a singleCellExperiment object is loaded that 
                        has been been processed with Seurat's CellCycleScoring, singleR from celldex/singleR has been applied 
                        and descriptions from ENSEMBL have been added.")),
                 div(h5("Optionally, the data can be subsampled.")),
                 div(h5("The \"regular expression to count genes/cell\" defines genes for which a collective count should be retained, 
                 even if those genes are removed during the filtering process. This way, one can keep track of ribosomal counts. 
                 The regular expression '^MT-' searches of all genes whose gene name starts (^) with 'MT-'")),
                 div(h5("")),
                 div(h5("")),
                 div(h5("")),
                 div(h5(""))
          )),
        br()),
      ## QC cells ----
      shinydashboardPlus::box(
        title = "QC - cells",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 div(h4("Counts"), align = "center"),
                 br(),
                 {
                   withTags(
                     nav(div(),ol(
                       li(
                         div(a(id="wkfl1.gQC_sampleHist.click", "check sample histogram")),
                         div(h6(HTML("[<em>General QC - Sample histogram</em>]"))),
                         div(h5("check that the samples have a similar distribution. An imbalance in samples can cause artefacts
                            in later analyses. In case of imbalances, either work on individual samples, or subsample to roughly 
                            the lowest sample.")),
                         div(h5("This graph can be saved to the history RMD file by clicking on the 'save to history' button."))
                       ),
                       li(
                         div(a(id="wkfl1.gQC_umiHist.click", "check UMI histogram")),
                         div(h6(HTML("[<em>General QC - UMI histogram</em>]"))),
                         div(h5("check that the there is only one obvious distribution. Multiple peaks may indicate contaminations")),
                         div(h5("take note of where to potentially set a threshold.")),
                         div(h5("This graph can be saved to the history RMD file by clicking on the 'save to history' button."))
                       ),
                       li(
                         div(a(id="wkfl1.nFeatureViolin.click", "visualize the number of features per cell in a violin plot per sample.")),
                         div(h6(HTML("[<em>Co-expression - 2D plot</em>] x = sampleNames; y=Feature.Count; color = sampleNames"))),
                       ),
                       li(
                         div(a(id="wkfl1.nFeatureSelection.click", "select cells to be removed based on the number of features per cell")),
                         div(h6(HTML("[<em>Co-expression - 2D plot</em>] x = barcode; y=Feature.Count; color = sampleNames"))),
                         div(h5("set 'name group, also used in Plot to color selected cells red.' to 'rmCells'.")),
                         div(h5("Select the cells to be removed and click 'change current selection' button, while 'Add to group/otherwise overwrite' is checked.")),
                         div(h5("repeat until all cells are selected. Using the plotly functionality it might be necessary to zoom into the plot and select within the zoomed region."))
                       ),
                       li(
                         div(a(id="wkfl1.nCountViolin.click", "view number of UMIs per cell")),
                         div(h6(HTML("[<em>Co-expression - 2D plot</em>] x = sampleNames; y=UMI.Count; color = sampleNames"))),
                         div(h5("Proceed similar to 3."))
                       ),
                       li(
                         div(a(id="wkfl1.nCountSelection.click", "select cells to be removed based on the number of UMIs per cell")),
                         div(h6(HTML("[<em>Co-expression - 2D plot</em>] x = barcode; y=UMI.Count; color = sampleNames"))),
                         div(h5("Proceed similar to 4."))
                       ),
                       li(
                         div(a(id="wkfl1.npMTViolin.click", "view  number of mitochondrial genes per cell")),
                         div(h6(HTML("[<em>Co-expression - 2D plot</em>] x = sampleNames; y=percent.mt; color = sampleNames"))),
                         div(h5("Proceed similar to 3."))
                       ),
                       li(
                         div(a(id="wkfl1.npMTSelection.click", "select cells to be removed based on the number of mitochondrial genes per cell")),
                         div(h6(HTML("[<em>Co-expression - 2D plot</em>] x = barcode; y=percent.mt; color = sampleNames"))),
                         div(h5("Proceed similar to 4."))
                       ),
                       li(
                         div(a(id="wkfl1.countFeature.click", "view the selected cells in the count vs Feature projection")),
                         div(h6(HTML("[<em>Co-expression - 2D plot</em>] x = nCount_RNA; y=Feature.count; color = sampleNames"))),
                         div(h5("test the effect of log-transforming the data.")),
                         div(h5("Proceed similar to 4."))
                       ),
                       li(
                         div(a(id="wkfl1.countMt.click", "view the selected cells in the count vs mitochondrial percentage projection")),
                         div(h6(HTML("[<em>Co-expression - 2D plot</em>] x = nCount_RNA; y=percent.mt; color = sampleNames"))),
                         div(h5("test the effect of log-transforming the data.")),
                         div(h5("Proceed similar to 4."))
                       ),
                       li(
                         div(h5("copy cell names")),
                         div(h6(HTML("[<em>Co-expression - 2D plot</em>] show cell names = checked; tripple-click"))),
                         div(h5("Selected cell names are shown at the bottom of the page (tripple click). Those can be copied into memory and then pasted 
                              in the ", a("'Cell selection'",id="wkfl1.go2CellSelection.click")," tab. ")),
                         div(h5("paste the cell names under 'Cells to be removed' and click 'apply changes' to update the underlying count 
                            matrix.")),
                         div(h6(HTML("[<em>Cell selection</em>] Cells to be removed")))
                         
                       ),
                       li(
                         div(h5("Validate using the previous figures with updated data.")),
                         div(h5("Optionally, define cells as high/normal/low with respect to mitochondrial content.")),
                         div(h5("To do so, define two new groups similar to 4. using the mitochondrial vs. nCount projection.")),
                         div(h5("Then combine these two new projections using a temp variable.", a(id="wkfl1.combineVars1.click","combine projections."))),
                         div(h6(HTML("[<em>Parameters - Projections - combine projections</em>]"))),
                         div(h5("Then rename the levels of this new temp variable.", a(id="wkfl1.renameLevels1.click","rename levels."))),
                         div(h6(HTML("[<em>Parameters - Projections - rename levels</em>]"))),
                         div(h5("Visualize the result in a plot similar to 10 with the color set to the new projection name.")),
                       )
                     )))
                 }
          ),
          br()
        )),
      ## Normalization ----
      shinydashboardPlus::box(
        title = "Normalization",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 div(h4("Processed counts"), align = "center"),
                 div(h5("Before starting the computations set the parameters for the time consuming processes."), align = "left"),
                 withTags(
                   nav(div(),ol(
                     li(
                       div(h5(a(id="wkfl1.selectNormParameters.click","select normalization parameters")), align = "left"),
                       div(h6(HTML("[<em>Parameters - Cluster parameters</em>]"))),
                       div(h5("Set the normalization method to 'scEx_log'"), align = "left"),
                       div(h5("scEX_log transforms using log2(+1) and devides by the total 
                              number of cells, which leaves values in a range between 0 and 1.
                              'scale by' defines a value to scale up to pseudo counts that are easier to handle.
                              When set to '0' the lowest value above zero will be chosen.")),
                       div(h5("Click on 'apply changes' button"), align = "left")
                     ),
                     li(
                       div(h5("select PCA parameters"), align = "left"),
                       div(h6(HTML("[<em>Parameters - Cluster parameters</em>]"))),
                       div(h5("Below the normalization method the parameters for the PCA can be set."), align = "left"),
                       div(h5("Set the 'Number of components', the number of variable features to use, 
                              whether to center and scale the data, which method to use and if to apply the Seurat 
                              implementation of not."), align = "left"),
                       div(h5(a(id="wkfl1.setPCAparameters.click", "Do this for me"), align = "left")),
                       div(h5("Click on 'apply changes' button"), align = "left")
                     ),
                     li(
                       div(h5("select Clustering parameters"), align = "left"),
                       div(h6(HTML("[<em>Parameters - Cluster parameters</em>]"))),
                       div(h5("Below the PCA parameters for the cluster parameters can be set."), align = "left"),
                       div(h5("There are several clustering methods available. SIMLR takes a very long time 
                              to run and should not be used as a first option."), align = "left"),div(h5(a(id="wkfl1.setClusterParameters.click", "Do this for me"), align = "left")),
                       div(h5("The Seurat implementation for clustering allows setting the number of PCs to use that 
                              were previously defined. K for the k-nearest neighbor algorithm can be set as well."), align = "left"),
                       div(h5(a(id="wkfl1.setClusterParameters.click", "Do this for me"), align = "left"))
                     ),
                     li(
                       div(h5(a(id="wkfl1.geneSelection.click", "Gene selection")), align = "left"),
                       div(h6(HTML("[<em>Gene selection</em>]"))),
                       div(h5("The min expression threshold is already applied durint the QC phase.")),
                       div(h5("Verify that the 'regular expression for selection of genes to be removed' is correct.")),
                       div(h5("Click on 'apply changes' button"), align = "left"),
                       div(h5("Validate that the interesting genes are still in the 'Gene kept' table and which 
                              genes are in the 'genes removed' table.")),
                     ),
                     li(
                       div(h5(a(id = "wkfl1.umap.click", "Activate the UMAP")), align = "left"),
                       div(h6(HTML("[<em>Parameters - Umap</em>]"))),
                       div(h5("Set 'N Neighbors' to 20; metric to cosine; spread to 7.")),
                       div(h5(a(id="wkfl1.setUMAPparameters.click", "Do this for me"), align = "left")),
                       div(h5("Click on 'apply changes' button"), align = "left")
                     ),
                     li(
                       div(h5(a(id = "wkfl1.input2.click", "Activate the calculations")), align = "left"),
                       div(h6(HTML("[<em>input</em>]"))),
                       div(h5("On the input select 'calculate logcounts using SCHNAPPs'"), align = "left"),
                       div(h5("This starts the calculations and can take a few minutes.")),
                     ),
                     li(
                       div(h5(a(id = "wkfl1.viewUMAP2D.click", "Visualize clusters in UMAP projection")), align = "left"),
                       div(h6(HTML("[<em>parameters - UMAP</em>]"))),
                       div(h5("This is only one of the possibilities to view the 2D projections")),
                     ),
                     li(
                       div(h5(a(id = "wkfl1.viewUMAP3D.click", "Visualize clusters in 3D UMAP projection")), align = "left"),
                       div(h6(HTML("[<em>parameters - tSNE</em>]"))),
                       div(h5("This is only possibility to view the 3D projections")),
                     ),
                     li(
                       div(h5(a(id = "wkfl1.renamedbCluster.click", "rename some Variables")), align = "left"),
                       div(h6(HTML("[<em>Parameters - Projections - Rename projections</em>]"))),
                       div(h5("SCHNAPPs uses the variables UMAPx, tsnex, dbCluster internally and overwrites them
                              whenever the calculations need to be updated. To be able to compare different settings,
                              for example when calculating cluster assignments, those projections can be renamed."), align = "left"),
                       div(h5("Select 'dbCluster' in 'projections to copy + rename' and set new name of Projection to 
                              'dbClusterLog'.")),
                       div(h5(a(id="wkfl1.renameCluster.click", "Do this for me"), align = "left")),
                       div(h5("Click on the 'rename' button.")),
                       div(h5("repeat this procedure for UMAP1,2,3 and PC_1"), align = "left"),
                     ),
                     li(
                       div(h5("save RData file"), align = "left"),
                       div(h6(HTML("[Button Download RData</em>]"))),
                       div(h5("click on the button and wait. It takes some time to compile all the data.
                              Save to a location of your liking."))
                     )
                   )
                   )
                 )
          )),
        br()),
      
    ),
    
    
    # wkfl2: specific questions ----
    controlbarItem(
      "specifiic analyses",
      ## Loading data ----
      shinydashboardPlus::box(
        title = "Loading data",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 div(h4(a(id = "wkfl2.LoadData2.click", "Load data"))),
                 div(h6(HTML("[<em>Input</em>]"))),
                 div(h5("Restart the app (this is not necessary, but it assures a 'cleaner' environment)."), align = "left"),
                 div(h5("Load data and GMT file"), align = "left")
          )
        )
        ,
        br()),
      shinydashboardPlus::box(
        title = "Heatmap",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 withTags(
                   nav(div(),ol(
                     li(
                       div(h5("show standard heatmap"), align = "left"),
                       div(h6(HTML("[<em>Co-expression -  heatmap</em>]"))),
                       div(h5(""), align = "left")
                     ),
                     li(
                       div(h5("change order"), align = "left")),
                     li(
                       div(h5("select cells / genes"), align = "left")),
                     li(
                       div(h5("up regulated genes between clusters"), align = "left")),
                     li(
                       div(h5("order by expression of genes"), align = "left"))
                   )))
          )),
        br()),
      ## DoubletFinder ----
      
      if("DoubletFinder" %in% installed.packages()){
        shinydashboardPlus::box(
        title = "QC - DoubletFinder",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 div(h5("DoubletFinder"), align = "center"),
                 div(h6(HTML("[<em>General QC - DoubletFinder</em>]"))),
          ))
        
      )},
      ## Investigate PCA
      
      shinydashboardPlus::box(
        title = "Investigate PCA",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 withTags(
                   nav(div(),ol(
                     li(
                       div(h4("PCA"), align = "center"),
                       div(h6(HTML("[<em>Parameters - Cluster Parameters - Loadings</em>]"))),
                       div(h5("check 'show row names'"), align = "left")),
                     div(h5("select top and bottom 10 genes of PC1"), align = "left")),
                     li(
                       div(h5("define gene set"), align = "left"),
                       div(h6(HTML("[<em>Parameters - Gene sets - edit gene set</em>]"))),
                     )
                     
                   ))))),
      br()),
    shinydashboardPlus::box(
      title = "compare two clusterings",background = "navy",
      width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
      fluidRow(
        column(width=11, offset = 0,
               withTags(
                 nav(div(),ol(
                   li(
                     div(h5("load data with annotation from Seb"), align = "left")),
                   li(
                     div(h5("DGE of outlier cells"), align = "left"))
                 )))
        )),
      br()),
    shinydashboardPlus::box(
      title = "Gene sets",background = "navy",
      width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
      fluidRow(
        column(width=11, offset = 0,
               withTags(
                 nav(div(),ol(
                   li(
                     div(h5("work with c7 gene set"), align = "left")),
                   li(
                     div(h5("load gene sets"), align = "left")),
                   li(
                     div(h5("visualize in Dot plot"), align = "left")),
                   li(
                     div(h5("dot plot with summary"), align = "left")),
                   li(
                     div(h5(""), align = "left"))
                 )))
        )),
      br()
    ),
    controlbarItem(
      "advanced features",
      ## Loading data ----
      shinydashboardPlus::box(
        title = "Trajectories",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 withTags(
                   nav(div(),ol(
                     li(
                       div(h5("trajectory: scorpius"), align = "left")),
                     li(
                       div(h5("trajectory: tempora"), align = "left")),
                     li(
                       div(h5("trajectory: elpigraph"), align = "left"))
                   )))
          )
        )
        ,
        br()
      ),
      shinydashboardPlus::box(
        title = "common features of cells",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 withTags(
                   nav(div(),ol(
                     li(
                       div(h5("heatmap"), align = "left")),
                     li(
                       div(h5("coefficient of variance"), align = "left")),
                     li(
                       div(h5("correlation coefficients"), align = "left"))
                   )))
          )
        ),
        br()
      ),
      shinydashboardPlus::box(
        title = "common features of cells",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 withTags(
                   nav(div(),ol(
                     li(
                       div(h5("heatmap"), align = "left")),
                     li(
                       div(h5("coefficient of variance"), align = "left")),
                     li(
                       div(h5("correlation coefficients"), align = "left"))
                   )))
          )
        )
        ,
        br())
    )
  )
)
