# Removing samples based on library size

#' @title Filter Samples based on Library Size
#'
#' @description This function is a part of the data wrangling functionality of `tcgaCleaneR`.
#' It allows user to handle the bias in \code{SummarizedExperiment} S4 class Cancer Dataset (e.g. TCGA dataset) due to
#' library size by filtering out the samples with sample size greater than the threshold. Using \code{plotLibSize},
#' user can determine the threshold.
#'
#' @param data SummarizedExperiment S4 class Dataset. E.g. TCGA Dataset.
#' @param ls_cutoff numeric: library size threshold
#'
#' @return S4 data object
#' @export
#'
#' @examples
#'
#' filterSamplesByLibSize(data = brca.data, ls_cutoff = 17.5)
#'
filterSamplesByLibSize <- function(data,ls_cutoff){
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  library_size <- log2(colSums(raw.count))
  keep.samples <- library_size > ls_cutoff
  brca.se.filtered <- data[ , keep.samples]
  return(brca.se.filtered)
}


