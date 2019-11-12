load("../scShinyHubData/scCourse2019.RData")


params <- list(
  min.mean = 1,
  min.size = 20,
  method = "igraph",
  use.ranks = TRUE,
  d = 10,
  BPPARAM = SerialParam()
)
params$assay.type <- "counts"
params$x <- scEx

register(MulticoreParam(
  workers = ifelse(detectCores() > 1, detectCores() - 1, 1)
),
default = TRUE
)

register(MulticoreParam(
  workers = 13
),
default = TRUE
)
system.time(retVal <- tryCatch(
  {
    suppressMessages(do.call("quickCluster", params))
  },
  error = function(e) {
    cat(file = stderr(), paste("\nProblem with clustering", e, "\n\n"))
    return(NULL)
  },
  warning = function(e) {
    if (DEBUG) cat(file = stderr(), paste("\nclustering produced Warning:\n", e, "\n"))
    # return(suppressMessages(do.call("quickCluster", params)))
  }
))


system.time(t <- bplapply(1:40000, FUN))
system.time(t <- bplapply(1:4000000, FUN))
system.time(t <- bplapply(1:40000, FUN))
system.time(t <- bplapply(1:400000, FUN))
system.time(t <- bplapply(1:400000, FUN), BPPARAM = SerialParam())
system.time(t <- bplapply(1:400000, FUN, BPPARAM = SerialParam()))
system.time(t <- bplapply(1:400000, FUN, BPPARAM = MultiCoreParam()))
system.time(t <- bplapply(1:400000, FUN, BPPARAM = MulticoreParam()))
system.time(t <- bplapply(1:400000, FUN, BPPARAM = SnowParam()))
system.time(t <- bplapply(1:400000, FUN))
