# Filtering Genes Function

#' @title Gene Filter
#'
#' @description This function is a part of the data wrangling functionality of `tcgaCleaneR`.
#' It allows user to input the \code{SummarizedExperiment} S4 class Cancer Dataset (e.g. TCGA dataset) and the required
#' genes to filter data based on genes.
#'
#' @usage filterGenesByBiotypes(data,gene.type)
#'
#' @param data SummarizedExperiment S4 class Dataset. E.g. TCGA Dataset.
#' @param gene.type A character vector of items.
#'
#' @return S4 data object
#' @export
#'
#' @examples
#' filterGenesByBiotypes(data=brca.data,gene.type=c("protein.coding"))
#' \dontrun{
#' filterGenesByBiotypes(data=brca.data,gene.type=c("protein.coding","snRNA"))
#'}
filterGenesByBiotypes <- function(data,gene.type){
  gene.annot.rm <-  as.data.frame(SummarizedExperiment::rowData(data))
  keep.genes.rm <- gene.annot.rm$Gene_BioType %in% gene.type
  data.filtered <- data[keep.genes.rm , ]
  return(data.filtered)
}

