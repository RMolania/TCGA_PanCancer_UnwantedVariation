#  PCs and library size/ purity/ time - Regression Analysis (regression + vector correlation for time)

#' @title Regression Analysis between PCs and unwanted variation
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`.
#' R2 values of fitted linear models are used to quantity the strength of the (linear) relationships between a single
#' quantitative source of unwanted variation such as sample (log) library size or tumor purity and global sample
#' summary statistics such as the first k PCs.
#' The function runs linear regression between unwanted variations in TCGA RNA-seq like
#' library size and purity with PCs from \code{computePCA}.
#' The output is a linear plot that compares the three \code{assays} in \code{SummarizedExperiment} TCGA Cancer data
#' across n PCs and R-sq.
#' For variable like Time which is not continuous, dummy variables are used to run
#' vector correlation \code{stats::cancor} between dummy time variable and n PCs with the same
#' linear output.
#'
#' @param pca.data list: PCA output from \code{computePCA}.
#' @param data S4 data object
#' @param type character: The response variable to \code{lm} model. groups included are 'librarysize', 'purity' and 'time'.
#' @param nPCs numeric: Number of PCs that needs to be used for regression
#'
#' @return Linear Plot the compares the correlation between library size (or Purity, time) and PCs across three datasets. When output is stored in a object the user can also access values used to plot the linear graphs.
#' @export
#'
#' @examples
#' \dontrun{
#' plotPCsVar(pca.data, data = brca.data, type = "purity", nPCs = 10)
#' df <- plotPCsVar(pca.data, data = brca.data, type = "time", nPCs = 8)
#' df
#' }

plotPCsVar <- function(pca.data, data, type, nPCs){
  raw.count <- as.data.frame(SummarizedExperiment::assay(data, 'HTseq_counts'))
  library.size <- log2(colSums(raw.count))
  data.set.names <- names(SummarizedExperiment::assays(data))
  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  if(type == "librarysize"){
    corr.cancer.tcga <- lapply(
      data.set.names,
      function(x){
        pcs <- pca.data[[x]]$sing.val$u
        tcga.ls.rSquared <- sapply(
          1:nPCs,
          function(y) {
            lm.ls <- summary(lm(library.size ~ pcs[, 1:y]))$r.squared
          })
      })
  } else
    if(type == "purity"){
      corr.cancer.tcga <- lapply(
        data.set.names,
        function(x){
          pcs <- pca.data[[x]]$sing.val$u
          tcga.ls.rSquared <- sapply(
            1:nPCs,
            function(y) {
              purity.ls <- summary(lm(sample.info$Purity_singscore ~ pcs[, 1:y]))$r.squared
            })
        })
    } else
      if(type == "time"){
        time.years <- fastDummies::dummy_cols(sample.info$Year)
        time.years <- time.years[, c(2:ncol(time.years))]
        corr.cancer.tcga <-
          lapply(
            data.set.names,
            function(x){
              sapply(
                1:nPCs,
                function(y) {
                  cca.pam50 <- stats::cancor(
                    x = pca.data[[x]]$sing.val$u[, 1:y, drop = FALSE],
                    y = time.years)
                  1 - prod(1 - cca.pam50$cor ^ 2)
                })

            })
      } else
        if(type == "biology"){
          biology.types <- fastDummies::dummy_cols(sample.info$Subtypes)
          biology.types <- biology.types[, c(2:ncol(biology.types))]
          corr.cancer.tcga <-
            lapply(
              data.set.names,
              function(x){
                sapply(
                  1:nPCs,
                  function(y) {
                    cca.pam50 <- stats::cancor(
                      x = pca.data[[x]]$sing.val$u[, 1:y, drop = FALSE],
                      y = biology.types)
                    1 - prod(1 - cca.pam50$cor ^ 2)
                  })

              })
        }
  names(corr.cancer.tcga) <- data.set.names
  if (length(SummarizedExperiment::assays(data)) == 4){
    # Visualize
    corr.normAssess <- data.frame(
      Raw.counts = corr.cancer.tcga$HTseq_counts,
      FPKM = corr.cancer.tcga$HTseq_FPKM,
      FPKM.UQ = corr.cancer.tcga$HTseq_FPKM.UQ,
      RUV_III = corr.cancer.tcga$RUV_III,
      pcs = c(1:nPCs)
    )

    # Visualize
    dataSets.colors <- wesanderson::wes_palette(
      n = 4,
      name = "GrandBudapest1")[c(1,2,4,3)]
    names(dataSets.colors) <- c(
      'Raw counts',
      'FPKM',
      'FPKM.UQ',
      'RUV_III'
    )

    corr.normAssess.p <- corr.normAssess %>%
      tidyr::pivot_longer(
        -pcs,
        names_to = 'Datasets',
        values_to = 'CorrValue') %>%
      dplyr::mutate(Datasets = replace(
        Datasets,
        Datasets == 'Raw.counts', 'Raw counts')) %>%
      dplyr::mutate(
        Datasets = factor(
          x = Datasets,
          levels = c('Raw counts', 'FPKM', 'FPKM.UQ', 'RUV_III'))) %>%
      data.frame(.)

    ggplot2::ggplot(corr.normAssess.p, ggplot2::aes(x = pcs, y = CorrValue, group = Datasets)) +
      ggplot2::geom_line(ggplot2::aes(color = Datasets), size = 1) +
      ggplot2::geom_point(ggplot2::aes(color = Datasets), size = 3) +
      ggplot2::xlab('PCs') +
      {if (type == 'purity')
      {
        ggplot2::ylab (expression("R"^"2"))
      } else if (type == 'librarysize'){
        ggplot2::ylab (expression("R"^"2"))
      } else if (type == 'time'){
        ggplot2::ylab ('Vector Correlation')
      } else if (type == 'biology'){
        ggplot2::ylab ('Vector Correlation')
      }} +
      ggplot2::scale_color_manual(
        values = c(dataSets.colors[1:4]),
        labels = c('Raw counts', 'FPKM','FPKM.UQ', 'RUV_III')) +
      ggplot2::scale_x_continuous(breaks = (1:nPCs), labels = c('PC1', paste0('PC1:', 2:nPCs)) ) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = 'black', size = 1),
        axis.title.x = ggplot2::element_text(size = 14),
        axis.title.y = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 12),
        legend.text = ggplot2::element_text(size = 10),
        legend.title = ggplot2::element_text(size = 14),
        strip.text.x = ggplot2::element_text(size = 10)
      ) + {if (type == 'purity')
      {
        ggplot2::ggtitle('Purity : PCs Regression Plot')
      } else if (type == 'librarysize'){
        ggplot2::ggtitle('Library Size : PCs Regression Plot')
      } else if (type == 'time'){
        ggplot2::ggtitle('Time : PCs Correlation Plot')
      } else if (type == 'biology'){
        ggplot2::ggtitle('Biology : PCs Correlation Plot')
      } }
  } else
    if(length(SummarizedExperiment::assays(data)) == 3){
  # Visualize
  corr.normAssess <- data.frame(
    Raw.counts = corr.cancer.tcga$HTseq_counts,
    FPKM = corr.cancer.tcga$HTseq_FPKM,
    FPKM.UQ = corr.cancer.tcga$HTseq_FPKM.UQ,
    pcs = c(1:nPCs)
  )

  # Visualize
  dataSets.colors <- wesanderson::wes_palette(
    n = 4,
    name = "GrandBudapest1")[c(1,2,4)]
  names(dataSets.colors) <- c(
    'Raw counts',
    'FPKM',
    'FPKM.UQ'
  )

  corr.normAssess.p <- corr.normAssess %>%
    tidyr::pivot_longer(
      -pcs,
      names_to = 'Datasets',
      values_to = 'CorrValue') %>%
    dplyr::mutate(Datasets = replace(
      Datasets,
      Datasets == 'Raw.counts', 'Raw counts')) %>%
    dplyr::mutate(
      Datasets = factor(
        x = Datasets,
        levels = c('Raw counts', 'FPKM', 'FPKM.UQ'))) %>%
    data.frame(.)

  ggplot2::ggplot(corr.normAssess.p, aes(x = pcs, y = CorrValue, group = Datasets)) +
    geom_line(aes(color = Datasets), size = 1) +
    geom_point(aes(color = Datasets), size = 3) +
    xlab('PCs') +
    {if (type == 'purity')
    {
      ggplot2::ylab (expression("R"^"2"))
    } else if (type == 'librarysize'){
      ggplot2::ylab (expression("R"^"2"))
    } else if (type == 'time'){
      ggplot2::ylab ('Vector Correlation')
    } else if (type == 'biology'){
      ggplot2::ylab ('Vector Correlation')
    }} +
    scale_color_manual(
      values = c(dataSets.colors[1:3]),
      labels = c('Raw counts', 'FPKM','FPKM.UQ')) +
    scale_x_continuous(breaks = (1:nPCs), labels = c('PC1', paste0('PC1:', 2:nPCs)) ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = 'black', size = 1),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 14),
      strip.text.x = element_text(size = 10)
    ) + {if (type == 'purity')
    {
      ggtitle('Purity : PCs Regression Plot')
    } else if (type == 'librarysize'){
      ggtitle('Library Size : PCs Regression Plot')
    } else if (type == 'time'){
      ggtitle('Time : PCs Correlation Plot')
    } else if (type == 'biology'){
      ggplot2::ggtitle('Biology : PCs Correlation Plot')
    }}
    }
  #return(corr.normAssess)
}
