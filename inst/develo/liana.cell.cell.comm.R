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

liana_test <- liana_wrap(scEx, idents_col = "dbCluster", assay="logcounts",
                         base = 2 , # log expression base
           resource = c( "MouseConsensus"  ))

liana_test %>% dplyr::glimpse()

liana_test <- liana_test %>%
  liana_aggregate()
dplyr::glimpse(liana_test)

liana_test %>%
  liana_dotplot(source_groups = unique(liana_test$source)[1:3],
                target_groups = unique(liana_test$target)[1:3],
                ntop = 20)

liana_trunc <- liana_test %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

heat_freq(liana_trunc)

p <- chord_freq(liana_trunc,
                source_groups = unique(liana_test$source)[1:3],
                target_groups = unique(liana_test$target)[1:3])
