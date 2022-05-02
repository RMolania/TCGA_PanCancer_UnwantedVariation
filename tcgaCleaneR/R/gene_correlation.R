#  Gene Correlation analysis function

#' @title Gene Correlation analysis
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`.
#' It helps to run correlation analysis between gene expression level and unwanted variation effect because of Library
#' Size or Purity. This function can help to quantify the association between an individual gene’s expression level
#' and library size or tumor purity for a given Cancer type.
#'
#' @param data S4 data object
#' @param is.log logical: Checks if the S4 data has log values. It 'False', it converts data to log scale.
#' @param type character: Variation variable to perform correlation with. type included are 'librarysize', 'purity_HTseq_counts', 'purity_HTseq_FPKM' and 'purity_HTseq_FPKM.UQ'.
#' @param cor.method a character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated. Default is spearman.
#' @param n.cores The number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be at least one, and parallelization requires at least two cores.
#'
#' @return A S3 data frame. The output contains the correlation test output containing pvalue, adj p-value and
#' Spearman's rank correlation coefficient. Along with the data frame output the function also returns a histogram
#' for Spearman's rank correlation coefficient for easy analysis of the test results.
#' @export
#'
#' @examples
#' \dontrun{
#' df <- computeCorr(data = brca.data,is.log = FALSE,type = "purity",cor.method = 'spearman',n.cores = 1)
#' computeCorr(data = brca.data,is.log = FALSE,type = "librarysize",cor.method = 'pearson',n.cores = 1)
#' }

computeCorr <- function(data, is.log, type, cor.method, n.cores){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  .corr.gene.variable <- function(expr.data, is.log, variable, method, n.cores){
    if(is.log){
      expr.data <- expr.data
    }else{
      expr.data <- log2(expr.data + 1)
    }
    rho <- parallel::mclapply(
      1:nrow(expr.data),
      function(x){
        round(cor.test(x = expr.data[x, ], y = variable, method = method)[[4]], 6)},
      mc.cores = n.cores)

    pval <- parallel::mclapply(
      1:nrow(expr.data),
      function(x){
        cor.test(x = expr.data[x, ], y = variable, method = method)[[3]]},
      mc.cores = n.cores)

    results <- data.frame(
      genes = row.names(expr.data),
      rho = unlist(rho),
      pvalue = unlist(pval),
      adj.pvalue = p.adjust(unlist(pval), 'BH')
    )
  }
  sample.info$ls <- log2(colSums(raw.count))
  raw.count <- as.matrix(raw.count)
  if(type == "librarysize"){
    genes.df <- .corr.gene.variable(
      expr.data = raw.count,
      is.log = is.log,
      variable = sample.info$ls,
      method = cor.method,
      n.cores = 1
    )
    hist(genes.df$rho)
  } else
    if(type == "purity"){
      genes.df <- .corr.gene.variable(
        expr.data = raw.count,
        is.log = is.log,
        variable = sample.info$Purity_singscore,
        method = cor.method,
        n.cores = 1
      )
      hist(genes.df$rho)
    }
  return(genes.df)
}
