source("../../inst/app/serverFunctions.R")
library(SingleCellExperiment)
library(shiny)
library(testthat)
library(SingleCellExperiment)
library(plotly)

# # Define a test case for the pltHighExp function
# test_that("pltHighExp plots the highest expressed genes", {
#   # Create a mock SingleCellExperiment object
#   # For simplicity, we will use a matrix for the expression data
#   scaterReads <- SingleCellExperiment(assays = list(counts = matrix(1:100, nrow = 10, ncol = 10)))
#   
#   # Specify the number of highest expressed genes and colors for the plot
#   n <- 5
#   scols <- c("red", "blue", "green")
#   
#   # Call the pltHighExp function
#   result <- pltHighExp(scaterReads, n, scols)
#   
#   # Check if the result is a ggplot2 plot object
#   expect_is(result, "ggplot")
#   
#   # You can also add further expectations to validate the plot content
#   # For example, you can test if the plot is colored correctly
#   
#   # Example expectation: the plot should have a custom color scale
#   expect_identical(result$scales$scales[[1]]$palette, scols)
# })

# Define a test case for the pltHighExp function
test_that("pltHighExp plots the highest expressed genes", {
 
  
  # Generate a mock SingleCellExperiment object with sampleNames
  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  colnames(counts) =  paste0("Cell", 1:10)
  scaterReads <- SingleCellExperiment(assays = list(counts = counts))
  rowData(scaterReads)$sampleNames <- paste0("Sample", 1:10)
  colData(scaterReads)$sampleNames <- c(rep("sample1",3), rep("sample2",3), rep("sample3",4))
  # Specify the number of highest expressed genes and colors for the plot
  n <- 5
  scols <- c("red", "blue", "green")
  
  # Call the pltHighExp function
  result <- pltHighExp(scaterReads, n, scols)
  
  # Check if the result is a ggplot2 plot object
  expect_s3_class(result, "ggplot")
  
  # You can also add further expectations to validate the plot content
  # For example, you can test if the plot is colored correctly
  
  # Extract the colors used in the plot
  colors_used <- unique(ggplot2::ggplot_build(result)$data[[1]]$colour)
  
  # Example expectation: the plot should have the specified colors
  expect_identical(colors_used, scols)
})

# Mock featureData for testing
features <- data.frame(symbol = c("GeneA", "GeneB", "GeneC"), row.names = c("Row1", "Row2", "Row3"))

# Test for known gene names
test_that("Known gene names are correctly indexed", {
  expected_indices <- c("Row1", "Row3")
  test_indices <- geneName2Index(c("GeneA, GeneC"), features)
  expect_equal(test_indices, expected_indices)
})

# Test for case-insensitivity and space removal
test_that("Function handles case-insensitivity and spaces", {
  expected_indices <- c("Row1", "Row2")
  test_indices <- geneName2Index(c(" genea ,GENE B"), features)
  expect_equal(test_indices, expected_indices)
})

# Test for gene names not in featureData
test_that("Warning for gene names not in featureData", {
  expect_warning(geneName2Index(c("GeneD, GeneE"), features))
})

# Test for null input
test_that("Handle null input correctly", {
  expect_null(geneName2Index(NULL, features))
})

# Run the tests

# Mock featureData for testing
features <- data.frame(symbol = c("GeneA", "GeneB", "GeneC"), row.names = c("Row1", "Row2", "Row3"))

# Test for known gene names
test_that("Known gene names are correctly indexed", {
  expected_indices <- c("Row1", "Row3")
  test_indices <- geneName2Index(c("GeneA, GeneC"), features)
  expect_equal(test_indices, expected_indices)
})

# Test for case-insensitivity and space removal
test_that("Function handles case-insensitivity and spaces", {
  expected_indices <- c("Row1", "Row2")
  test_indices <- geneName2Index(g_id = c(" genea, GENE B"), featureData = features)
  expect_equal(test_indices, expected_indices)
})

DEBUG <<-T
# Test for gene names not in featureData
test_that("Warning for gene names not in featureData", {
  expect_warning(geneName2Index(c("GeneD, GeneE"), features))
})

# Test for null input
test_that("Handle null input correctly", {
  expect_null(geneName2Index(NULL, features))
})

# Run the tests
# test_dir("path/to/your/tests")
# test_file("tests/testthat/serverFunctionstests.R")

#### update umi


test_that("geneName2Index() returns correct gene names", {
  # Define sample inputs and expected output
  geneNames <- c("gene1, gene2")
  geneNames2 <- c("gene3, gene4")
  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  colnames(counts) = paste0("Cell", 1:10)

  scEx <- SingleCellExperiment(assays = list(counts = counts))
  rownames(scEx) <- paste0("gene", 1:10)
  rowData(scEx)$symbol = paste0("gene", 1:10)
  colData(scEx)$sampleNames <- c(rep("sample1",3), rep("sample2",3), rep("sample3",4))
  projections <- list()
  
  # Call the function
  result <- updateProjectionsWithUmiCount(geneNames, geneNames2, scEx, projections)
  
  # Check the output
  expect_equal(result$UmiCountPerGenes, colSums(counts[1:2, ]))
  expect_equal(result$UmiCountPerGenes2, colSums(counts[3:4, ]))
})

test_that("geneName2Index() returns 0 when geneNames and geneNames2 are empty", {
  # Define sample inputs and expected output
  g_id <- character(0)
  scEx <- SingleCellExperiment(assays = list(counts = counts))
  rownames(scEx) <- paste0("gene", 1:10)
  rowData(scEx)$symbol = paste0("gene", 1:10)
  colData(scEx)$sampleNames <- c(rep("sample1", 3), rep("sample2", 3), rep("sample3", 4))
  projections <- rowData(scEx)
  
  # Call the function
  result <- geneName2Index(g_id, projections)
  
  # Check the output
  expect_equal(result, NULL)
})



# Sample data for testing
set.seed(123)
scEx_log <- SingleCellExperiment(assays = list(counts = matrix(rnorm(1000), ncol = 10)))
projections <- data.frame(PC1 = rnorm(10), PC2 = rnorm(10), sampleNames = sample(letters[1:3], 10, replace = TRUE))
featureData <- data.frame(symbol = rownames(scEx_log))
geneNames <- c("GeneA", "GeneB")
geneNames2 <- c("GeneC", "GeneD")
dimX <- "PC1"
dimY <- "PC2"
clId <- "Cluster1"
grpN <- NULL
legend.position <- "topright"
grpNs <- NULL
logx <- FALSE
logy <- FALSE
divXBy <- "None"
divYBy <- "None"
dimCol <- "sampleNames"
colors <- c("a" = "blue", "b" = "red", "c" = "green")

# Test cases
test_that("plot2DprojectionNew returns a plotly object", {
  p <- plot2DprojectionNew(
    scEx_log, projections, "Sample Group 1", featureData,
    geneNames, geneNames2, dimX, dimY, clId, grpN, legend.position, grpNs,
    logx, logy, divXBy, divYBy, dimCol, colors
  )
  expect_s3_class(p, "plotly")
})

test_that("plot2DprojectionNew returns NULL for invalid projections", {
  invalid_projections <- data.frame(PC1 = rnorm(10), PC2 = rnorm(10))
  p <- plot2DprojectionNew(
    scEx_log, invalid_projections, "Sample Group 1", featureData,
    geneNames, geneNames2, dimX, dimY, clId, grpN, legend.position, grpNs,
    logx, logy, divXBy, divYBy, dimCol, colors
  )
  expect_null(p)
})

test_that("plot2DprojectionNew handles empty projections", {
  empty_projections <- data.frame()
  p <- plot2DprojectionNew(
    scEx_log, empty_projections, "Sample Group 1", featureData,
    geneNames, geneNames2, dimX, dimY, clId, grpN, legend.position, grpNs,
    logx, logy, divXBy, divYBy, dimCol, colors
  )
  expect_null(p)
})

test_that("plot2DprojectionNew handles missing columns in projections", {
  missing_col_projections <- data.frame(PC1 = rnorm(10))
  p <- plot2DprojectionNew(
    scEx_log, missing_col_projections, "Sample Group 1", featureData,
    geneNames, geneNames2, dimX, dimY, clId, grpN, legend.position, grpNs,
    logx, logy, divXBy, divYBy, dimCol, colors
  )
  expect_null(p)
})

test_that("plot2DprojectionNew handles NULL projections", {
  p <- plot2DprojectionNew(
    scEx_log, NULL, "Sample Group 1", featureData,
    geneNames, geneNames2, dimX, dimY, clId, grpN, legend.position, grpNs,
    logx, logy, divXBy, divYBy, dimCol, colors
  )
  expect_null(p)
})

