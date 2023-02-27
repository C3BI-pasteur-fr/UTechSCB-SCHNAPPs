# Definition of the short cuts tab images and on-click events.

# output images ----
output$shortCut.geneSelection <- renderUI({
  div(id = "shortCut.geneSelection.click",
      tags$img(src = "www/images/schnappsLogo.png", 
               width = 200, height = 200)
      #   contentType = "image/png",
      #   width = 500,
      #   height = 500,
      #   alt = "Scater plot will be here when 'apply changes' is checked"
  )
})

# on click events ----
onclick(
  "shortCut.geneSelection.click", 
  { 
    updateTabItems(
      session = session,
      "sideBarID",
      selected = "geneSelection"
    )
  }
)