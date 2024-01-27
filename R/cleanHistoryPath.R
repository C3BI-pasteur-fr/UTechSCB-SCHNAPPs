#' clean SCHNAPPs history path of RData files that are not referenced in any Rmd file in that directory.
#'
#' @param path History path that holds min. one Rmd file
#' @param verbose Print what is/would be executed
#' @param dryrun don't execute the removal commands.
#'
#' @return nothing
#' @export cleanHistoryPath
#'
#' @examples
cleanHistoryPath <- function(path, verbose=TRUE, dryrun=TRUE){
  if(dryrun) verbose = TRUE
  if(!dir.exists(path)){
    stop("path is not a directory")
  }
  rmdFiles = dir(path = path, pattern = "*.Rmd", full.names = TRUE)
  if(length(rmdFiles)<1){
    stop("no Rmd files found.")
  }
  if(verbose){
    message(paste("rmd files found:", paste(rmdFiles,collapse = ", ")))
  }
  # load text from Rmd Files
  rmdText = c()
  for (fp in rmdFiles){
    rmdText = c(rmdText, read.delim2(fp, header = F, sep = "\n", quote = NULL))
  }
  # check if rdata file exists in rmd files
  rdataFiles = dir(path = path, pattern = ".*\\..*\\.RData", full.names = TRUE)
  for(fp in rdataFiles){
    if(!any(grepl(basename(fp), rmdText))){
      if(verbose)
        message(paste("not found: ", basename(fp)))
      if(!dryrun)
        file.remove(fp)
    } else {
      if (verbose)
      message(paste("found: ", basename(fp)))
    }
  }
  
}

# path = "/Volumes/LaCie2022/RStudio_history/celia/hist_2023-May-26.15.18/"
# 
# cleanHistoryPath(path, verbose=TRUE, dryrun=F)
  