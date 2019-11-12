library(SingleCellExperiment)
cp <- load("~/Downloads/mv1-2-e16-avc_for_trajectory.Rds")
outfile <- "~/Downloads/mv1-2-e16-avc_for_trajectory.v2.Rds"
featureData(gbm)
featuredata$symbol <- featuredata$Associated.Gene.Name
featuredata$id <- rownames(featuredata)
scEx <- SingleCellExperiment(
  assay = list(counts = gbm@assayData$exprs),
  colData = phenoData(gbm)@data,
  rowData = featuredata
)

save(file = outfile, list = c("scEx"))
