remotes::install_github('saezlab/liana')
library(liana)
library(tidyverse)
library(magrittr)
library(SingleCellExperiment)

if(!require("circlize")){
  install.packages("circlize", quiet = TRUE,
                   repos = "http://cran.us.r-project.org")
}
require("circlize")

show_methods() # multiple
show_resources() # only one

liana_path <- system.file(package = "liana")
testdata <-
  readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

testdata %>% dplyr::glimpse()
liana_test <- liana_wrap(testdata)
liana_test %>% dplyr::glimpse()


cp = load("inst/develo/puceal.1-3.project.2023-11-21.RData")
cp

scEx %>% dplyr::glimpse()
rownames(scEx) = toupper(rownames(scEx))
assays(scEx)
assays(scEx)[["logcounts"]] = as(assays(scEx)[["logcounts"]],"CsparseMatrix")

liana_scEx <- liana_wrap(scEx, idents_col = "dbCluster", assay="logcounts",
                         base = 2 , # log expression base
           resource = c( "MouseConsensus"  ))

liana_scEx %>% dplyr::glimpse()

liana_scEx <- liana_scEx %>%
  liana_aggregate()
dplyr::glimpse(liana_scEx)

liana_scEx %>%
  liana_dotplot(source_groups = unique(liana_scEx$source)[1:3],
                target_groups = unique(liana_scEx$target)[1:3],
                ntop = 20)

liana_truncscEx <- liana_scEx %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

# how to get this to work???
heat_freq(liana_truncscEx)
liana_truncscEx %>% dplyr::glimpse()
liana_truncscEx$source %>% table()
liana_truncscEx$target %>% table()
freqs <- liana_truncscEx %>% liana:::.get_freq()
liana_heatmap(mat = freqs)

p <- chord_freq(liana_truncscEx,
                source_groups = unique(liana_scEx$source)[1:3],
                target_groups = unique(liana_scEx$target)[1:3])

