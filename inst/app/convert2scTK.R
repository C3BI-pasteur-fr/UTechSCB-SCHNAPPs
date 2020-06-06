library(singleCellTK)
library(SCHNAPPs)

cd = load("~/Downloads/rnEPDC.lite.RData")
cd
counts_mat <- assay(scEx, "counts")
sample_annot <- colData(scEx)
row_annot <- rowData(scEx)
newSCE <- createSCE(assayFile = counts_mat, annotFile = sample_annot, 
                    featureFile = row_annot, assayName = "counts",
                    inputDataFrames = TRUE, createLogCounts = TRUE)


singleCellTK(newSCE)

