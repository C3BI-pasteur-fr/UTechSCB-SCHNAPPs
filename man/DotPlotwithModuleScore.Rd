% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DotPlotwithModuleScore.R
\name{DotPlotwithModuleScore}
\alias{DotPlotwithModuleScore}
\title{DotPlotwithModuleScore}
\usage{
DotPlotwithModuleScore(
  object,
  assay = NULL,
  features,
  featureDat = featureDat,
  cols = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  clusters = "sampleNames",
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
)
}
\arguments{
\item{object}{A Seurat object.}

\item{assay}{The assay to be used for plotting. Defaults to NULL, in which case the default assay of the object is used.}

\item{features}{A list of character vectors of gene names in "symbol" annotation (SCHNAPPs specific).}

\item{featureDat}{A data frame containing information about features, particularly the "symbol" column.}

\item{cols}{A vector of colors for plotting.}

\item{col.min}{Minimum value for color scaling.}

\item{col.max}{Maximum value for color scaling.}

\item{dot.min}{Minimum value for dot scaling.}

\item{dot.scale}{Scaling factor for dot size.}

\item{idents}{Identifiers for cells.}

\item{group.by}{Grouping variable.}

\item{split.by}{Variable to split the plot by.}

\item{cluster.idents}{Boolean indicating whether to cluster identifiers.}

\item{clusters}{Column name for grouping identifiers.}

\item{scale}{Logical indicating whether to scale the data.}

\item{scale.by}{The type of scaling to be applied ('size' or 'radius').}

\item{scale.min}{Minimum value for scaling.}

\item{scale.max}{Maximum value for scaling.}
}
\value{
A ggplot2 object representing the DotPlot with module scores.
}
\description{
Adapted from Seurat::DotPlot but adds the module score per group.
}
