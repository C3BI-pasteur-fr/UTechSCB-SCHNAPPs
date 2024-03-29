.onLoad <- function(libname, pkgname) {
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
