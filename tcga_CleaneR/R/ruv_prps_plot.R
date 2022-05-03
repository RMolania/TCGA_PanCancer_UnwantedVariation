# RUV-III - PRPS((Pseudo replicate of pseudo sample)) map

#' @title RUV-III PRPS Map
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`. It helps to map the pseudo samples
#' for different biology in a cancer type replicated across batches. The replication of of pseudo samples across batches
#' are known as Pseudo Replicates.
#'
#' @param data S4 data object
#' @param n numeric: Minimum no. of samples needed to make pseudo sample. By default it consider 3 samples.
#'
#' @return PRPS Map which can be used to show how the biological variation are distributed across batches.
#' @export
#' @examples
#' \dontrun{
#' plotPRPS(data = brca.data)
#' }

plotPRPS <- function(data,n=3){
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  #sample.info$biology <- sample(letters[1:4], nrow(sample.info), replace = TRUE)
  sample.info$new.batch <- paste0(
    sample.info$Year, #sample.info$Year,
    '_',
    sample.info$Plates #sample.info$Plates
  )


  df_count <- sample.info %>%
    dplyr::count(new.batch, Subtypes)

  if (n!=3){
    df_count$use <- 'Un-selected'
    df_count$use[df_count$n > n-1] <- 'Selected'
  } else {
    df_count$use <- 'Un-selected'
    df_count$use[df_count$n > 2] <- 'Selected'
  }


  ggplot(df_count, aes(x = new.batch, y = Subtypes)) +
    geom_count(aes(color = use)) +
    geom_text(aes(
      label = n,
      hjust = 0.5,
      vjust = 0.5
    )) +
    xlab('Years-plates') +
    ylab('Biological groups') +
    theme_bw()+
    theme(
      axis.line = element_line(colour = 'black', size = .85),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(
        size = 12,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      strip.text.x = element_text(size = 10),
      legend.position = 'none'
    )
}
