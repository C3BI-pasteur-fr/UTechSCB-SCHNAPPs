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

load(file = "~/Downloads/AMB/AllCellTypes.schnapps.RData")
library(Seurat)
library(ggplot2)
seuratdata <- UpdateSeuratObject(as.Seurat(scEx, assay = "RNA", data=NULL))
seuratdata[["nCount_RNA"]] =  colSums(seuratdata, slot= 'counts')
seuratdata[["nFeature_RNA"]] =  colSums(GetAssayData(seuratdata, slot= 'counts') >0 )
seuratdata = PercentageFeatureSet(object = seuratdata, pattern = "^mt-|^MT-", col.name = "percent.mt")


mito.genes <- grep(pattern = "^MT-", ignore.case = TRUE, x = rownames(x = seuratdata), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", ignore.case = TRUE, x = rownames(x = seuratdata), value = TRUE)

seuratdata[["nMito_RNA"]] =  colSums(seuratdata[mito.genes,], slot= 'counts')
seuratdata[["nRps_RNA"]] =  colSums(seuratdata[RPS.genes,], slot= 'counts')
seuratdata[["umiMT.ratio"]] = seuratdata[["nCount_RNA"]]/seuratdata[["nMito_RNA"]]
FeatureScatter(object = seuratdata, feature1 = "nCount_RNA", feature2 = "percent.mt",  pt.size = 0.01)+ geom_density_2d()
FeatureScatter(object = seuratdata, feature1 = "percent.mt", feature2 = "nMito_RNA",  pt.size = 0.01)+ geom_density_2d()
FeatureScatter(object = seuratdata, feature1 = "umiMT.ratio", feature2 = "nMito_RNA",  pt.size = 0.01)+ geom_density_2d()
FeatureScatter(object = seuratdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",  pt.size = 0.01)+ geom_density_2d()

