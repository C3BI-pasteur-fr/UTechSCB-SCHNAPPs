oldname = "renderPlotly"
newName = "plotly::renderPlotly"
system(paste0("find . -name \"*.R\"  -exec sed -i '' -e 's/", oldname, "/", newName,"/g' {} +"))


system(paste0("find . -name \"*.Rmd\"  -exec sed -i '' -e 's/", oldname, "/", newName,"/g' {} +"))
