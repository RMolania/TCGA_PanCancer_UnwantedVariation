# Purity filter function

#' @title Filter Samples based on Tumor Purity
#'
#' @description This function is a part of the data wrangling functionality of `tcgaCleaneR`.
#' It allows user to handle sample purity in \code{SummarizedExperiment} S4 class Cancer Dataset (e.g. TCGA dataset) by
#' filtering out the samples that are above a specific purity threshold.
#'
#'
#' @param data SummarizedExperiment S4 class Dataset. E.g. TCGA Dataset.
#' @param purity_cutoff numeric: Sample purity cutoff
#'
#' @return S4 data object
#' @export
#'
#' @examples
#'
#' filterSamplesByPurity(data= brca.data,purity_cutoff= 0.496)
#'
filterSamplesByPurity <- function(data,purity_cutoff){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  keep.samples <- sample.info$Purity_singscore > purity_cutoff
  brca.se.filtered <- data[ ,keep.samples]
  return(brca.se.filtered)
}

