# sc RNAseq PhD course 2021 
# data preparation for SCHNAPPs
# 

library(hdf5r)
library(SingleR)
library(Seurat)
library(SCHNAPPs)
library(Matrix)
library(SingleCellExperiment)
library(celldex)
library(AnnotationDbi)
library(org.Mm.eg.db)

filename = "../mouseHeartDatasetAnalysis/data/heart_10k_v3_filtered_feature_bc_matrix.h5"

# wrapper script to apply singleR to data
# 
applySingleR <- function(scEx, hpca.se, colName, label = "label.main") {
  if (!"logcounts" %in% assayNames(scEx)){
    scEx <- scater::logNormCounts(scEx)
  }
  
  common <- intersect(rownames(scEx), rownames(hpca.se))
  hpca.se <- hpca.se[common,]
  scExC <- scEx[common,]
  pred.hpca <- SingleR(test = scExC, ref = hpca.se, labels = colData(hpca.se)[,label])
  to.remove <- pruneScores(pred.hpca)
  new.pruned <- pred.hpca$labels
  new.pruned[pruneScores(pred.hpca, nmads=5)] <- NA
  cdata = colData(scEx)
  cdata[, colName] <- NA
  cdata[colnames(scExC), colName]  <- new.pruned
  colData(scEx) <- cdata
  return(scEx)
}
# similar to Seurat function, but appends sample name to cell identifiers
# 
schnappsRead10X_h5 <- function(filename, use.names = TRUE, unique.features = TRUE, sampleName = "1") {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = "r")
  genomes <- names(x = infile)
  output <- list()
  if (!infile$attr_exists("PYTABLES_FORMAT_VERSION")) {
    if (use.names) {
      feature_slot <- "features/name"
    }
    else {
      feature_slot <- "features/id"
    }
  }  else {
    if (use.names) {
      feature_slot <- "gene_names"
    }
    else {
      feature_slot <- "genes"
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, "/data")]]
    indices <- infile[[paste0(genome, "/indices")]]
    indptr <- infile[[paste0(genome, "/indptr")]]
    shp <- infile[[paste0(genome, "/shape")]]
    features <- infile[[paste0(genome, "/", feature_slot)]][]
    barcodes <- infile[[paste0(genome, "/barcodes")]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1, p = indptr[],
      x = as.numeric(x = counts[]), dims = shp[], giveCsparse = FALSE
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = "CsparseMatrix")
    if (infile$exists(name = paste0(genome, "/features"))) {
      types <- infile[[paste0(genome, "/features/feature_type")]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(X = types.unique, FUN = function(x) {
          return(sparse.mat[which(x = types == x), ])
        }, simplify = FALSE, USE.NAMES = TRUE)
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  
  if (length(x = output) == 1) {
    if (class(output[[genome]]) == "list") {
      for (le in 1:length(output[[genome]])) {
        colnames(output[[genome]][[le]]) <- paste0(sub("(.*)-(.*)", "\\1", colnames(output[[genome]][[le]])), "-", sampleName)
      }
    }else{
      colnames(output[[genome]]) <- paste0(sub("(.*)-(.*)", "\\1", colnames(output[[genome]])), "-", sampleName)
    }
    return(output[[genome]])
  }
  else {
    # colnames(output) <- paste0(sub("(.*)-(.*)", "\\1", colnames(output)), "-", sampleName)
    return(output)
  }
}

seuratObj = schnappsRead10X_h5(filename = filename)

dim(seuratObj)

# remove genes with no expression
seuratObj = seuratObj[rowSums(as.matrix(seuratObj))>0,]

dim(seuratObj)

# cells/genes with equal expression pose a problem
# these are genes with min. expression.
which(duplicated(as.matrix(seuratObj)))
rowSums(as.matrix(seuratObj)[duplicated(as.matrix(seuratObj)),])

seuratObj = unique(as.matrix(seuratObj))

dim(seuratObj)


## cell cycle prediction
## 
seuratObj = CreateSeuratObject(counts = seuratObj)
seuratObj <- NormalizeData(seuratObj)

# genes defined in Seurat package
# Mouse genes are first character upper case rest lower case
# HS is all upper case
s.genes <- stringr::str_to_title(cc.genes$s.genes[cc.genes$s.genes %in% toupper(rownames(seuratObj))])
g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes[cc.genes$g2m.genes %in% toupper(rownames(seuratObj))])
seuratObj = CellCycleScoring(seuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^mt-")

## cell type prediction
## 
pd <- data.frame(
  barcode = sub("(.*)-(.*)", "\\1", colnames(seuratObj)),
  sampleNames = sub(".*-(.*)", "\\1", colnames(seuratObj))
)
pd$barcode <- as.character(pd$barcode)
rownames(pd) <- colnames(seuratObj)


scEx <- SingleCellExperiment(
  assay = list(counts = GetAssayData(object = seuratObj, slot = "counts")),
  colData = pd,
  rowData = NULL
)

scEx = scater::logNormCounts(scEx)


system.time(
  scExTemp <- applySingleR(scEx = scEx, hpca.se = celldex::MouseRNAseqData(), colName = "cellTypes", label = "label.fine")
)

table(colData(scExTemp)$cellTypes)

if (all(colnames(scEx) == rownames(seuratObj[[]]))){
  scEx[["cellTypes"]] = as.factor(colData(scExTemp)$cellTypes)
  scEx[["Phase"]] = as.factor(seuratObj[[]]$Phase)
  scEx[["S.Score"]] = seuratObj[[]]$S.Score
  scEx[["G2M.Score"]] = seuratObj[[]]$G2M.Score
  scEx[["percent.mt"]] = seuratObj[[]]$"percent.mt"
}else{
  print ("not true")
}

rownames(rowData(scEx)) 
gaAnno = mapIds(org.Mm.eg.db, rownames(rowData(scEx)), keytype = "SYMBOL", "GOALL")
description = mapIds(org.Mm.eg.db, rownames(rowData(scEx)), keytype = "SYMBOL", "GENENAME")
ontology = mapIds(org.Mm.eg.db, rownames(rowData(scEx)), keytype = "SYMBOL", "ONTOLOGY")


rowData(scEx)[["id"]] = rownames(scEx) # needed
rowData(scEx)[["symbol"]] = rownames(scEx) # needed 
rowData(scEx)[["GOALL"]] = gaAnno[rownames(scEx)]
rowData(scEx)[["Description"]] = description[rownames(scEx)] # needed
rowData(scEx)[["Ontology"]] = ontology[rownames(scEx)]

fname = "data/scCoursePhD.2021"
outfile <- paste0(fname, ".RData")
save(file = outfile, list = c("scEx"))





