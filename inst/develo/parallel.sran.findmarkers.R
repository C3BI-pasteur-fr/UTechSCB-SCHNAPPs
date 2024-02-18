

library(furrr)
library(pryr)
# library(lineprof)
# source(file = "~/Rstudio/scran/R/findMarkers.R")
cp = load("/Users/bernd/SCHNAPPsDebug/clusterServerreturnValues.RData")
direction="up"
lfc = 1
object_size(scEx_log, 
            projections$dbCluster,
            direction = direction,
            lfc = lfc)
# 4.18 GB
register(MulticoreParam(6))
{
  start.time <- base::Sys.time()
  wmarkers <- tryCatch({
    prof <- utils::Rprof(filename = "Rprof.6.out",memory.profiling=T, filter.callframes = T)
    scran::findMarkers(scEx_log, 
                       projections$sampleNames,
                       direction = direction,
                       lfc = lfc,
                       BPPARAM = BiocParallel::bpparam())
    
    
    
  }, error = function(e) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("coE_scranFindMarkerTableReact", id = "coE_scranFindMarkerTableReactError", duration = NULL, type = "error")
    }
    browser()
    cat(file = stderr(), unlist(e))
    return(NULL)
  })
  cat(file = stderr(), base::Sys.time() - start.time)
}
# 47 with 2
# ~30 min with 4 cpus
# 25.77998 with 6
# 34 min with 8
# kills machine with 16
# 
summaryRprof(filename = "Rprof.8.out",
             memory = c( "both")) # both", "tseries", "stats"
# tseries
# "Rprof.2.out"
# vsize.small vsize.large     nodes duplications                         stack:2
# 0.02   107826376 10671913536  756221144          739 "tryCatch":"scran::findMarkers"
# 0.04        4096     2693960     117544           44 "tryCatch":"scran::findMarkers"
# "Rprof.4.out"
# vsize.small vsize.large     nodes duplications                         stack:2
# 0.02   116030640 12796485800  987884184          288 "tryCatch":"scran::findMarkers"
# 0.04       29544        4568    1101016          319 "tryCatch":"scran::findMarkers"
# 0.06       19128     6735648     717528          176 "tryCatch":"scran::findMarkers"
# "Rprof.8.out"
# vsize.small vsize.large     nodes duplications                         stack:2
# 0.02   107873280 10670562584  757334536          727 "tryCatch":"scran::findMarkers"
# 0.04        5208     4040936     148008           56 "tryCatch":"scran::findMarkers"
# "Rprof.16.out" # killed
# vsize.small vsize.large      nodes duplications                         stack:2
# 0.02   118854744 10720463736 1019764144          726 "tryCatch":"scran::findMarkers"
# 0.04        5376     4040936     148344           57 "tryCatch":"scran::findMarkers"

