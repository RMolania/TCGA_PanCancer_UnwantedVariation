# Generate PCA function

#' @title Generate PCA
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`. It performs PCA using SVD
#' algorithm (\code{runSVD()}) on the \code{SummarizedExperiment} class TCGA Cancer Data across all \code{assays()} and
#' generate PCs.
#'
#' @param data S4 data object
#' @param nPcs numeric: Number of PCs that needs to be generated
#' @param is.log logical: Checks if the S4 data has log values. If 'False', it converts data to log scale.
#'
#' @return A List of S4 class containing n PCs. For all three \code{assays()}, two values are returned. \code{sing.val} contains singular values with \code{u} containing PCs. \code{variation} contains the variation of each PC.
#' @export
#'
#' @examples
#' \dontrun{
#' computePCA(data = brca.data, nPcs = 10, is.log = FALSE)
#' }
#'

computePCA <- function(data, nPcs, is.log){
  .pca <- function(data, nPcs, is.log) {
    if(is.log){
      data <- data
    }else{
      data <- log2(data + 1)
    }
    svd <- BiocSingular::runSVD(
      x = t(data),
      k = nPcs,
      BSPARAM = BiocSingular::bsparam(),
      center = TRUE,
      scale = FALSE
    )
    percent <- svd$d^2/sum(svd$d^2)*100
    percent <-
      sapply(
        seq_along(percent),
        function(i) {round(percent[i], 1)})
    return(list(
      sing.val = svd,
      variation = percent))
  }
  tcga.harmonized <- names(SummarizedExperiment::assays(data))
  pca.cancer.tcga  <- lapply(
    tcga.harmonized,
    function(x){
      .pca(
        data = as.matrix(SummarizedExperiment::assay(data, x)),
        nPcs = nPcs,
        is.log = is.log)
    })
  names(pca.cancer.tcga) <- tcga.harmonized
  return(pca.cancer.tcga)
}
