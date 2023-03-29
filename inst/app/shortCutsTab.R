source(paste0(packagePath, "/shortCuts.def.R"), local=T)
# shortCutsList[[.twoDHeader]]  = list()

shortCutItem <- function(scLset, namItem){
  # cat(file = stderr(),paste0("++++++++sCOut.", make.names(namItem),"\n"))
  userBox(
    title = userDescription(
      title = scLset[[namItem]][["descr"]],
      subtitle = NULL,
      type = 2,
      image = scLset[[namItem]][["url"]],
    ),
    collapsed = F,
    collapsible = TRUE,
    status='navy',
    background = "blue",
    width=12,
    boxToolSize = 'lg',
    headerBorder=F,
    # imageElevation = 2,
    id = paste0("sCOut.", make.names(namItem),".plus"),
    uiOutput(paste0("sCOut.", make.names(namItem)), 
             click = paste0("sCOut.", make.names(namItem),".click")),
    div(h6(scLset[[namItem]][["explain"]]))
  )
}

# shortCutsTab is used in tabs.R 
shortCutsTabFunc <- function(shortCutsList) {
  # outTag = 
  rows = list()
  
  for (nam in names(shortCutsList)){
    namItemList = list()
    rows = append(rows, fluidRow(div(h4(nam), 
                                     align = "left")) %>% list())
    idx = 1
    for (namItem in names(shortCutsList[[nam]])){
      namItemList = append(namItemList, 
                           column(width = 4, shortCutItem(shortCutsList[[nam]], namItem)) %>% list()
      )
      if(idx%%3 == 0){
        rows = append(rows, fluidRow(
          namItemList
        ) %>% list()
        )
        namItemList = list()
      }
      idx = idx + 1
    }
    if(length(namItemList)>0) {
      rows = append(rows, fluidRow(
        namItemList
      ) %>% list()
      )
    }
  }
  return(rows)
}

shortCutsTab <- function(){
  shinydashboard::tabItem(
    tabName = "shortCuts",
    fluidRow(div(h3(.generalExplain), 
                 align = "center")),
    div(h6(.generalDescr)),
    actionButton("expandAllShortCuts","toggle all short cuts (not working)"),
    shortCutsTabFunc(shortCutsList)
  )
}
# shortCutsTab()
