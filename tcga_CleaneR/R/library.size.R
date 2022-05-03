# determining library size

#' @title Library Size
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`.
#' It allows user to analyze the Library Size bias, a technical bias in the TCGA RNA-seq.
#' The user can input the \code{SummarizedExperiment} S4 class Cancer Dataset (e.g. TCGA dataset) and the type of plot
#' to analyse variations in library size across years and samples.
#'
#' @param data SummarizedExperiment S4 class Dataset. E.g. TCGA Dataset.
#' @param plot_type character: Plot type
#'
#' @return Scatter Plot, Box plot
#' @export
#'
#' @examples
#' plotLibSize(data = brca.data, plot_type = "Scatterplot")
#' \dontrun{
#' plotLibSize(data = brca.data, plot_type = "Boxplot")
#' }
plotLibSize <- function(data,plot_type){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  library_size <- log2(colSums(raw.count))
  to.plot.ls <- data.frame(ls = library_size, samples = 1:length(library_size))
  if (plot_type == "Boxplot"){
    ggplot2::ggplot(to.plot.ls, aes(x=sample.info$Year, y=ls, fill = factor(sample.info$Year))) +
      geom_boxplot(alpha=0.3) +
      ylab('Log2 library size (total counts)') +
      xlab('Year') +
      labs(fill='Year') +
      ggtitle('Library size')
  } else
    if (plot_type == "Scatterplot"){
      ggplot2::ggplot(to.plot.ls, aes(x = samples, y = ls, color = factor(sample.info$Year))) +
        geom_point() +
        ylab('Log2 library size (total counts)') +
        xlab('Samples') +
        ggplot2::labs(color='Year') +
        ggtitle('Library size') +
        theme(
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black', size = 1),
          plot.title = element_text(size = 18),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          strip.text.x = element_text(size = 10),
          strip.text = element_text(size = 14))
    }
}
