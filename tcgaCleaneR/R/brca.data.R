# TCGA Breast Cancer Data

#' @title TCGA Breast Cancer Data
#'
#' @description This dataset is a condensed version of the original TCGA Breast Cancer Data. The datset
#' is a \code{SummarizedExperiment} class data object which contains information about Genes, Samples,
#' list of row count and normalised row count at gene-level and gene IDs.
#'
#' @format The S4 data object is a matrix container with information on 100 genes. The rows represent
#' features of interest (e.g. genes) and columns represent samples. The objects contain one or more
#' assays, each represented by a matrix-like object of numeric or other mode. Information about these
#' features is stored in a DataFrame object, accessible using the function \code{rowData()}. The
#' sample information is accessible using the function \code{colData()}. The assays contain read counts
#' per gene with set of raw and normalised values.
"brca.data"
