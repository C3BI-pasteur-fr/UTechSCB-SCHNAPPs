# Definition of the short cuts tab images and on-click events.
source(paste0(packagePath, "/shortCuts.def.R"), local=T)

# output images ----
# output$shortCut.geneSelection <- renderUI({
#   div(id = "shortCut.geneSelection.click",
#       tags$img(src = "www/images/schnappsLogo.png", 
#                width = 200, height = 200)
#       #   contentType = "image/png",
#       #   width = 500,
#       #   height = 500,
#       #   alt = "Scater plot will be here when 'apply changes' is checked"
#   )
# })

# for loops are not working with observes/onclick/layzy stuff
lapply(names(shortCutsList),FUN = function(nam){
  lapply(names(shortCutsList[[nam]]),FUN = function(namItem){
    id = paste0("sCOut.", make.names(namItem),".click")
    url = shortCutsList[[nam]][[namItem]][["url"]]
    fun = shortCutsList[[nam]][[namItem]][["fun"]]
    output[[paste0("sCOut.", make.names(namItem))]] <-  renderUI({
      div(id = id,
          tags$img(src = url, 
                   width = 400, height = 400))})
    divFunc = function(x, session, ...){
      x(session)
    }
    cat(file = stderr(), paste0("========= sCOut.", make.names(namItem),".click\n"))
    shinyjs::onclick(paste0("sCOut.", make.names(namItem),".click"),divFunc( x=fun, session=session))
  })
})




observeEvent(eventExpr = input$expandAllShortCuts, label = "expandAllShortCuts", {
  clicked  = input$expandAllShortCuts
  cat(file = stderr(),  "=============toggle all\n")
  lapply(names(shortCutsList),FUN = function(nam){
    lapply(names(shortCutsList[[nam]]),FUN = function(namItem){
      cat(file = stderr(),  paste("====", paste0("sCOut.", make.names(namItem),".plus"),"\n"))
      shinydashboardPlus::updateBox(paste0("sCOut.", make.names(namItem),".plus"), action = "toggle")
    })
  })
  # updateBox("sCOut.Gene.selection.plus", action = "remove", session = session)
  # updateBox("sCOut.Gene.selection.plus", action = "toggle", session = session)
})


# for (nam in names(shortCutsList)){
#   for (namItem in names(shortCutsList[[nam]])){
#     # since this is layzy evaluation we need to pass the variables in a new env
#     env=new.env()
#     env$id = paste0("sCOut.", make.names(namItem),".click")
#     env$url = shortCutsList[[nam]][[namItem]][["url"]]
#     env$fun = shortCutsList[[nam]][[namItem]][["fun"]]
#     output[[paste0("sCOut.", make.names(namItem))]] <-  renderUI(env = env, {
#       div(id = id,
#           tags$img(src = url, 
#                    width = 200, height = 200))})
#     divFunc = function(x, session, ...){
#       browser()
#       # cat(file=stderr(), paste("=========", id, url, fun, "\n"))
#       x(session)
#       }
#     environment(divFunc) = env
#     cat(file=stderr(), paste("========= create onclick", env$id, env$url, "\n"))
#     onclick(paste0("sCOut.", make.names(namItem),".click"),divFunc(event, x=env$fun, session=session)
#     #         {
#     #   cat(file = stderr(), paste("========onclick"), namItem, "\n")
#     # }
#     )
#   }
# }



# on click events ----
# onclick(
#   "shortCut.geneSelection.click", 
#   { 
#     updateTabItems(
#       session = session,
#       "sideBarID",
#       selected = "geneSelection"
#     )
#   }
# )