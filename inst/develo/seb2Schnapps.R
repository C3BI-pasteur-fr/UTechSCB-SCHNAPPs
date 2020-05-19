# install scater https://bioconductor.org/packages/release/bioc/html/scater.html
library(scater)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
# library(loomR)
library(Seurat)
library(patchwork)


# file from Seb's pipeline (changed the extension to RData)
cp = load("/Users/bernd/Downloads/FBJ_SC01_GFPc_so.RData")
# creates "FBJ_SC01_GFPc_so" a Seurat object

scEx = as.SingleCellExperiment(FBJ_SC01_GFPc_so)

cp = load("~/Downloads/rnEPDC.lite.RData")
# cp
# "scEx"            "pca"             "scol"            "ccol"            "geneCount"       "umiCount"        "beforeFilterPrj" "tsne"            "dbCluster"      
# "umapReact"


# problems:
# the data is in two files from Seb. 
# pca etc are not calculated
# but the data is normalized.

# Solution:
# 1. save as singleCellExperiment object
# 2. load both files with standard SCHNAPPs and use normalization


cp = load("/Users/bernd/Downloads/FBJ_SC01_GFPc_so.RData")
scEx = as.SingleCellExperiment(FBJ_SC01_GFPc_so)
colnames(colData(scEx))[1] = "sampleNames"
colData(scEx)$barcode = rownames(colData(scEx))
# rownames(rowData(scEx))
rowData(scEx)$symbol = rownames(rowData(scEx))
rowData(scEx)$id = rownames(rowData(scEx))
rowData(scEx)$Description = rownames(rowData(scEx))
save(file = "FBJ_SC01_GFPc_so.scEx.RData", list = c("scEx"))

cp = load("/Users/bernd/Downloads/FBJ_SC01_BufferC_so.RData")
scEx = as.SingleCellExperiment(get(cp))
colnames(colData(scEx))[1] = "sampleNames"
colData(scEx)$barcode = rownames(colData(scEx))
# rownames(rowData(scEx))
rowData(scEx)$symbol = rownames(rowData(scEx))
rowData(scEx)$id = rownames(rowData(scEx))
rowData(scEx)$Description = rownames(rowData(scEx))
save(file = "FBJ_SC01_BufferC_so.scEx.RData", list = c("scEx"))

