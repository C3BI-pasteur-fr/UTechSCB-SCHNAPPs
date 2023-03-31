# local helper function
add2workflowObsList <- function(li, wkfl, workflowObsList){
  if(!wkfl %in% names(workflowObsList)) workflowObsList[[wkfl]] = list()
  workflowObsList[[wkfl]][[li$id]] =   list(id = li$id, fun=li$Func)
  return(workflowObsList)
}

# shortCutsTab ----
workflowObsList = list()


wkfl = "wkfl1"
workflowObsList[[wkfl]] = list()

# load data link
workflowObsList = list(
  id = "LoadData",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "input"
    )
    updateRadioButtons(session = session, "whichscLog", selected = "disablescEx_log")
    updateBox("addOptInput", action = "toggle")

  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)

# show same distribution
workflowObsList = list(
  id = "showSamples",
  Func = function(session) {
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "input"
    )
    updateRadioButtons(session = session, "whichscLog", selected = "disablescEx_log")
    updateBox("addOptInput", action = "toggle")
    
  }
  
) %>% add2workflowObsList(wkfl = wkfl, workflowObsList)



# for loops are not working with observes/onclick/layzy stuff
lapply(names(workflowObsList),FUN = function(nam){
  lapply(names(workflowObsList[[nam]]),FUN = function(namItem){
    fun = workflowObsList[[nam]][[namItem]][["fun"]]
    divFunc = function(x, session, ...){
      x(session)
    }
    cat(file = stderr(), paste0("=+++++++======== ", nam,".", make.names(namItem),".click\n"))
    onclick(paste0(nam,".", make.names(namItem),".click"),divFunc(x=fun, session=session))
  })
})

