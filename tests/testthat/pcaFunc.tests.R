library(testthat)
library(shiny)
# Load your function (assuming it's in a file named 'your_functions.R')
source("inst/app/reactives.R")

# Create mock data for testing
library(SingleCellExperiment)
# scEx <- SingleCellExperiment(
#   assays = list(logcounts = matrix(rnorm(100), nrow = 10)),
#   colData = DataFrame(sampleNames = factor(rep(c("A", "B"), each = 5))),
#   rowData = DataFrame(geneNames = paste("Gene", 1:10))
# )

counts <- matrix(rpois(10000, 5), nrow = 100, ncol = 100)
colnames(counts) =  paste0("Cell", 1:100)
rownames(counts) = paste0("gene",1:100)
scEx <- SingleCellExperiment(assays = list(counts = counts))
rowData(scEx)$sampleNames <- paste0("Sample", 1:100)

colData(scEx)$sampleNames <- c(rep("sample1",30), rep("sample2",30), rep("sample3",40))


scEx_log <- scEx
assayNames(scEx_log) <- "logcounts"
rank <- 3
center <- FALSE
scale <- FALSE
useSeuratPCA <- TRUE
pcaGenes <- colnames(assays(scEx)[["logcounts"]])
rmGenes <- c()
featureData <- rowData(scEx)
pcaN <- 10
maxGenes <- 1000
hvgSelection <- "vst"
inputNormalization <- "DE_something"

DEBUG <<- "TRUE"
if(is.null(.schnappsEnv)){
  .schnappsEnv =  new.env(parent=emptyenv())
  .schnappsEnv$DEBUGSAVE = FALSE
}


test_that("Test for pcaFunc with valid input", {
  result <- pcaFunc(
    scEx,
    scEx_log,
    rank,
    center,
    scale,
    useSeuratPCA,
    pcaGenes,
    rmGenes,
    featureData,
    pcaN,
    maxGenes,
    hvgSelection,
    inputNormalization
  )
  
  # Perform assertions on the result
  expect_true(!is.null(result), "Result should not be NULL")
  expect_type(result, "list")
  expect_true("x" %in% names(result), "Result should contain 'x'")
  expect_true("var_pcs" %in% names(result), "Result should contain 'var_pcs'")
  expect_true("rotation" %in% names(result), "Result should contain 'rotation'")
  # Add more specific assertions based on the expected behavior of pcaFunc
})

center <- T
scale <- T
useSeuratPCA <- F
inputNormalization <- "DE_seuratSCTnorm"


test_that("Test for pcaFunc with valid input", {
  result <- pcaFunc(
    scEx,
    scEx_log,
    rank,
    center,
    scale,
    useSeuratPCA,
    pcaGenes,
    rmGenes,
    featureData,
    pcaN,
    maxGenes,
    hvgSelection,
    inputNormalization
  )
  
  # Perform assertions on the result
  expect_true(!is.null(result), "Result should not be NULL")
  expect_type(result, "list")
  expect_true("x" %in% names(result), "Result should contain 'x'")
  expect_true("var_pcs" %in% names(result), "Result should contain 'var_pcs'")
  expect_true("rotation" %in% names(result), "Result should contain 'rotation'")
  # Add more specific assertions based on the expected behavior of pcaFunc
})



test_that("Test for pcaFunc with invalid input", {
  # Test with invalid input that should return NULL
  invalid_result <- pcaFunc(
    scEx,
    scEx_log,
    rank = 0,  # Invalid rank
    center,
    scale,
    useSeuratPCA,
    pcaGenes,
    rmGenes,
    featureData,
    pcaN,
    maxGenes,
    hvgSelection,
    inputNormalization
  )
  
  # Perform assertions on the invalid result
  expect_true(is.null(invalid_result), "Result should be NULL for invalid input")
})
