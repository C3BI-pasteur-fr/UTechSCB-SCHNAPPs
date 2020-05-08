oldname = ", ls(envir = globalenv())"
newName = ""
system(paste0("find . -name \"*.R\"  -exec sed -i '' -e 's/", oldname, .Platform$file.sep, newName,"/g' {} +"))
system(paste0("find . -name \"*.Rmd\"  -exec sed -i '' -e 's/", oldname, .Platform$file.sep, newName,"/g' {} +"))
