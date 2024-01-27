source("../../inst/app/serverFunctions.R")
library(SingleCellExperiment)
library(shiny)
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
library(testthat)

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
library(testthat)

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

context("Example Code")

test_that("geneName2Index() returns correct gene names", {
  # Define sample inputs and expected output
  dimX <- 1
  dimY <- 2
  geneNames <- c("gene1", "gene2")
  geneNames2 <- c("gene3", "gene4")
  scEx <- SingleCellExperiment(assays = list(counts = counts))
  rownames(rowData(scEx)) <- paste0("gene", 1:10)
  colData(scEx)$sampleNames <- c(rep("sample1",3), rep("sample2",3), rep("sample3",4))
  projections <- list()
  
  # Call the function
  result <- updateProjectionsWithUmiCount(dimX, dimY, geneNames, geneNames2, scEx, projections)
  
  # Check the output
  expect_equal(result$UmiCountPerGenes, c(1, 1))
  expect_equal(result$UmiCountPerGenes2, 0)
})

test_that("geneName2Index() returns 0 when geneNames and geneNames2 are empty", {
  # Define sample inputs and expected output
  dimX <- 1
  dimY <- 1
  geneNames <- character(0)
  geneNames2 <- character(0)
  scEx <- matrix(nrow = 1, ncol = 1)
  projections <- list(UmiCountPerGenes = 1, UmiCountPerGenes2 = 1)
  
  # Call the function
  result <- example_function(dimX, dimY, geneNames, geneNames2, scEx, projections)
  
  # Check the output
  expect_equal(result$UmiCountPerGenes, 0)
  expect_equal(result$UmiCountPerGenes2, 0)
})



