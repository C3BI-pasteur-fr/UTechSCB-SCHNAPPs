fps <- list.files(pattern = ".RData$|.RDs$",
                  path = "F:/Rstudio/data/",
                  recursive = TRUE,
                  full.names = T,
                  ignore.case = T)
for (fp in fps) {
  print(fp)
  cont = T
  scEx = NULL
  cp = tryCatch(load(file = fp),
           error = function(x){
             # print(x)
             cont = F
           }
  ) 
  if (cont) {
    if (is.null(scEx)){
      print(paste(fp, "no scEX"))
    } else {
      if (startsWith(rownames(scEx)[1], "ENS")) {
        if("symbol" %in% colnames(rowData(scEx))) {
          print(head(rownames(scEx)))
          rownames(scEx) = make.unique(rowData(scEx)$symbol)
          save(file = paste0(fp,".noENSG.RData"), list = c("scEx"))
        }
      }
    }
    
  }
}

rownames(scEx) = make.unique(rowData(scEx)$symbol)