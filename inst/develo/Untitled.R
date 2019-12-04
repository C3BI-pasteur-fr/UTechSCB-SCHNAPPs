library(rhdf5)
h5ls("~/Downloads/GSE115189/suppl/GSM3169075_filtered_gene_bc_matrices_h5.h5")
mydata <- h5read("~/Downloads/GSE115189/suppl/GSM3169075_filtered_gene_bc_matrices_h5.h5", "/Homo_sapiens.GRCh38.v90.cellranger/data")


load(file = "GSE3611.RData")

as(gse[[2]], "SingleCellExperiment")
assayData(gse[[1]])$exprs

colnames(gse[[1]])
pData(gse[[1]])
gse[[1]]
SingleCellExperiment(assays = assayData(gse[[1]])$exprs)

# need to set rownames before

scEx <- SingleCellExperiment(
  counts = assayData(gse[[1]])$exprs,
  colData = pData(phenoData(gse[[1]]))
)

example_sce <- SingleCellExperiment(
  assays = list(counts = sc_example_counts),
  colData = sc_example_cell_info
)


# gsm holds data
# gse just describes the objects
# still don't know if individual cells are also individual gsm and how to link annoation
# also unclear on how to tranform to scEx


gpl <- getGEO("GPL17021")

Meta(gse[[1]])
