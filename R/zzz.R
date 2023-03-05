.onLoad <- function(libname, pkgname) {
  cat(file = stderr(),"\n\non Load\n\n")
  shiny::addResourcePath(
    prefix = "www",
    directoryPath = system.file(
      "www",
      package = "SCHNAPPs"
    )
  )
}

.onUnload <- function(libname, pkgname) {
  shiny::removeResourcePath("www")
}
