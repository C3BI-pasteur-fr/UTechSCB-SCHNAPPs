# shortCutsTab ----

# should be moved in a language dependent file to be sourced
.geneCellHeader = "Gene/cell manipulations"
.twoDHeader = "2D plots"
.generalDescr = "Plots shown are representative and not related to the data used."
.generalExplain = "Click on a figure to get to the that representation in the app"

# specific text
.geneSelectionHead = "Gene selection"
.geneSelectionDescr = "description"



# shortCutsTab is used in tabs.R 
shortCutsTab <- function() {
  shinydashboard::tabItem(
    tabName = "shortCuts",
    fluidRow(div(h3(.generalExplain), 
                 align = "center")),
    div(h6(.generalDescr)),
    
    # Gene/cell manipulations ----
    fluidRow(div(h4(.geneCellHeader), 
                 align = "left")),
    fluidRow(
      column(width = 4,
             div(h4(.geneSelectionHead)),
             # shortCut.geneSelection is the output variable described in shortCuts.R
             # click relates to the onclick event in surtCuts.R
             uiOutput("shortCut.geneSelection", click = "shortCut.geneSelection.click"),
             div(h6(.geneSelectionDescr))
      )
    ),
    
    # 2D plots ----
    fluidRow(div(h4(.twoDHeader), 
                 align = "left")),
    
    fluidRow(
    )
  )
}
