# Anova Function

#' @title Anova test for batch effect on individual gene expressions
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`. It use ANOVA F statistics to
#' summarize the effects of a qualitative source of unwanted variation (e.g. batches) on the expression levels of
#' individual genes. Unwanted Variations such as Plate effect and Time effect can be analysed by this test.
#'
#'
#' @param data S4 data object
#' @param variable character: The predictor variable to \code{lm} model. The variables included are 'Time' and 'Plate'
#' @param is.log logical: Checks if the S4 data has log values. It 'False', it converts data to log scale.
#' @param n.cores The number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be
#' at least one, and parallel computing requires at least two cores.
#'
#' @return A S3 data frame. The output contains the Anova test (F test) scores corresponding to all genes in S4 data object.
#' @export
#'
#' @examples
#' \dontrun{
#' computeANOVA(data = brca.data, variable = "Time", is.log = FALSE, n.cores = 1)
#' computeANOVA(data = brca.data, variable = "Plate", is.log = FALSE, n.cores = 1)
#' }

computeANOVA <- function(data, variable, is.log, n.cores){
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  sample.info$ls <- log2(colSums(raw.count))
  raw.count <- as.matrix(raw.count)
  average.exp <- log2(rowMeans(raw.count))
  if(is.log == TRUE){
    raw.count <- raw.count
  }else{
    raw.count <- log2(raw.count + 1)
  }
  if (variable == "Plate"){
    anova_test <- parallel::mclapply(
      1:nrow(raw.count),
      function(x) {
        MASS::dropterm(lm(raw.count[x , ] ~ sample.info$Plates), test = 'F')[c(5:6)]
      }
      , mc.cores = n.cores)
  } else
    if (variable == "Time"){
      anova_test <- parallel::mclapply(
        1:nrow(raw.count),
        function(x) {
          MASS::dropterm(lm(raw.count[x , ] ~ sample.info$Year), test = 'F')[c(5:6)]
        }
        , mc.cores = n.cores)
    }
  test.values <- data.frame(
    Genes = row.names(raw.count),
    FValue = round(unlist(lapply(anova_test, function(x) x$`F Value`[2])), digits = 4) ,
    PValue = unlist(lapply(anova_test, function(x) x$`Pr(F)`[2])),
    Adj.PValue = p.adjust(unlist(lapply(anova_test, function(x) x$`Pr(F)`[2])), method = 'BH'),
    Mean = round(average.exp, digits = 2)
  )
  return(test.values)
}
