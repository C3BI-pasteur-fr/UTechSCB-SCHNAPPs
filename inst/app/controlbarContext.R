controlbarContext <- shinydashboardPlus::dashboardControlbar(
  width = 500, collapsed = F,
  id = "controlbar",
  controlbarMenu(
    id = "menu",
    controlbarItem(
      "Workflow 1",
      shinydashboardPlus::box(
        title = "Loading data",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=11, offset = 0,
                 div(h4(a(id = "wkfl1.LoadData.click", "Load data"))),
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
                         div(h5("check that the samples have a similar distribution. An imbalance in samples can cause artefacts
                            in later analyses. In case of imbalances, either work on individual samples, or subsample to roughly 
                            the lowest sample.")),
                         div(h5("This graph can be saved to the history RMD file by clicking on the 'save to history' button."))
                       ),
                       li(
                         div(a(id="wkfl1.gQC_umiHist.click", "check UMI histogram")),
                         div(h5("check that the there is only one obvious distribution. Multiple peaks may indicate contaminations")),
                         div(h5("take note of where to potentially set a threshold.")),
                         div(h5("This graph can be saved to the history RMD file by clicking on the 'save to history' button."))
                       ),
                       li(
                         div(a(id="wkfl1.nFeatureViolin.click", "visualize the number of features per cell in a violin plot per sample."))
                       ),
                       li(
                         div(a(id="wkfl1.nFeatureSelection.click", "select cells to be removed based on the number of features per cell")),
                         div(h5("set 'name group, also used in Plot to color selected cells red.' to 'rmCells'.")),
                         div(h5("Select the cells to be removed and click 'change current selection' button, while 'Add to group/otherwise overwrite' is checked.")),
                         div(h5("repeat until all cells are selected. Using the plotly functionality it might be necessary to zoom into the plot and select within the zoomed region."))
                       ),
                       li(
                         div(a(id="wkfl1.nCountViolin.click", "view number of UMIs per cell")),
                         div(h5("Proceed similar to 3."))
                       ),
                       li(
                         div(a(id="wkfl1.nCountSelection.click", "select cells to be removed based on the number of UMIs per cell")),
                         div(h5("Proceed similar to 4."))
                       ),
                       li(
                         div(a(id="wkfl1.npMTViolin.click", "view  number of mitochondrial genes per cell")),
                         div(h5("Proceed similar to 3."))
                       ),
                       li(
                         div(a(id="wkfl1.npMTSelection.click", "select cells to be removed based on the number of mitochondrial genes per cell")),
                         div(h5("Proceed similar to 4."))
                       ),
                       li(
                         div(a(id="wkfl1.countFeature.click", "view the selected cells in the count vs Feature projection")),
                         div(h5("test the effect of log-transforming the data.")),
                         div(h5("Proceed similar to 4."))
                       ),
                       li(
                         div(a(id="wkfl1.countMt.click", "view the selected cells in the count vs mitochondrial percentage projection")),
                         div(h5("test the effect of log-transforming the data.")),
                         div(h5("Proceed similar to 4."))
                       ),
                       li(
                         div(h5("Selected cell names are shown at the bottom of the page (tripple click). Those can be copied into memory and then pasted 
                              in the ", a("'Cell selection'",id="wkfl1.go2CellSelection.click")," tab. ")),
                         div(h5("paste the cell names under 'Cells to be removed' and click 'apply changes' to update the underlying count 
                            matrix."))
                       ),
                       li(
                         div(h5("Validate using the previous figures with updated data.")),
                         div(h5("Optionally, define cells as high/normal/low with respect to mitochondrial content.")),
                         div(h5("To do so, define two new groups similar to 4. using the mitochondrial vs. nCount projection.")),
                         div(h5("Then combine these two new projections using a temp variable.", a(id="combineVars1","combine projections."))),
                         div(h5("Then rename the levels of this new temp variable.", a(id="renameLevels1","rename levels."))),
                       )
                     )))
                 }
          )),
        br()),
      shinydashboardPlus::box(
        title = "Normalization",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=10, offset = 1,
                 div(h4("Counts"), align = "center"),
                 div(h5("select normalization parameters"), align = "left"),
                 div(h5("select PCA parameters"), align = "left"),
                 div(h5("select Clustering parameters"), align = "left"),
                 div(h5("select normalization parameters"), align = "left"),
                 div(h5("save RData file"), align = "left"),
                 
          )),
        br()),
      shinydashboardPlus::box(
        title = "QC - genes",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=10, offset = 1,
                 div(h4("PCA"), align = "center"),
                 div(h5("select top and bottom 10 genes of PC1"), align = "left"),
                 div(h5("select PCA parameters"), align = "left"),
                 
          )),
        br()),
      shinydashboardPlus::box(
        title = "QC - DoubletFinder",background = "navy",
        width = 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
        fluidRow(
          column(width=10, offset = 1,
                 div(h4("Counts"), align = "center"),
                 
          )),
        br())
    )
    
  )
)
