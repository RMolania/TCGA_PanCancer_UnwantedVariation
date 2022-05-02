# Removing lowly expressed genes

#' @title Low Count Genes Filter
#'
#' @description This function is a part of the data wrangling functionality of `tcgaCleaneR`.
#' It allows user to input the \code{SummarizedExperiment} S4 class Cancer Dataset (e.g. TCGA dataset) and the threshold
#' for the minimum gene count and sample count.
#'
#' @param data SummarizedExperiment S4 class Dataset. E.g. TCGA Dataset.
#' @param gene_count numeric: gene count threshold
#' @param sample_size numeric: sample size threshold
#'
#' @return S4 data object
#' @export
#'
#' @examples
#'
#' filterLowExprGenes(data=brca.data,gene_count = 20,sample_size = 200)
#'
filterLowExprGenes <- function(data,gene_count,sample_size){
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  keep.high <- apply(
    raw.count,
    1,
    function(x) length(x[x>gene_count])>=sample_size
  )
  data.filtered <- data[keep.high , ]
  return(data.filtered)
}
