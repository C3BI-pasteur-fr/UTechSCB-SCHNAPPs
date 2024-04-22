# update to the bleeding edge modifications 

if (!"xbioc" %in% rownames(installed.packages())) {
  remotes::install_github("renozao/xbioc")
}

devtools::install_github('C3BI-pasteur-fr/UTechSCB-SCHNAPPs', dependencies = TRUE)

BiocManager::valid()

library(SCHNAPPs)

library(future)
# plan("multiprocess", workers = 6)
plan(sequential)

library("BiocParallel")
register(safeBPParam(6))

schnapps(DEBUG=T)

##### 
# Different ways of starting SCHNAPPs
##### 

# reproducible work, allows creating and modifying a history of important event during the analysis.
# this take time and disc space, so it is not used by default.
schnapps(historyPath = "./history")


# adding functionality, e.g. trajectory inference
# see https://github.com/baj12/SCHNAPPsContributions for further information
schnapps(localContributionDir = "someplace on the disc with more tools")


