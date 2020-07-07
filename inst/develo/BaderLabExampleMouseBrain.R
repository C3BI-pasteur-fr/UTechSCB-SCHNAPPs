install.packages("scClustViz")
devtools::install_github("BaderLab/scClustViz")
BiocManager::install("org.Mm.eg.db")
devtools::install_github('immunogenomics/presto')


library("scClustViz")
library('presto')
cp = load("~/Downloads/AMB/AllCellTypes.RData")

renameCols = function(scEx, oldName, newName) {
  colnames(colData(scEx)) = stringr::str_replace(colnames(colData(scEx)),oldName,newName)
  scEx
}
scEx = renameCols(scEx, "sample_order","sampleNames")
colData(scEx)$barcode = rownames(colData(scEx))

rowData(scEx)$Description = ""
rowData(scEx)$id = rownames(rowData(scEx))
rowData(scEx)$symbol = rownames(rowData(scEx))

#ident is cell type
## colData(scEx)[names(Clusters(AMB_sCV[[1]])),"cellType"] = as.character(Clusters(AMB_sCV[[1]]))
## colData(scEx)$cellType = as.factor(colData(scEx)$cellType)


save(file = "~/Downloads/AMB/AllCellTypes.schnapps.RData", list = c("scEx"))
