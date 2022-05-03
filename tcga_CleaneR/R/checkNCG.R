# Check Negative Control Gene function

#' @title Check Unwanted variation in Negative Control Genes
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`. It helps to visualize the presence
#' of unwanted variation in the Negative Control Genes (NCG) for RUV-III analysis. Using PCA analysis any relation between
#' the negative control sets in the TCGA data and the unwanted variation type can be visualized. The plot between PCs and
#' unwanted variation can help verify the presence of unwanted variation in NCG sets.
#'
#' @param data S4 data object
#' @param ncg_set character: A character vector of NCG sets from TCGA data. At present TCGA data includes the following NCG variables i.e. "scRNAseq_HK","RNAseq_HK","Microrray_HK","NanostringPanCancer_HK" and "sinscorePanCancer_HK".
#' @param group character: Color code PCs based on group. groups included are 'Time', 'Tissue', 'Plate', 'TSS', 'Center'
#' @param plot_type character: Plot type. Includes 'DensityPlot' and 'BoxPlot'.
#' @param nPcs numeric: Number of PCs that needs to be computed
#' @param npcs numeric: Number of PCs that needs to be plotted
#' @param is.log logical: Checks if the S4 data has log values. If 'False', it converts data to log scale.
#'
#' @return Density Plot, Box Plot
#' @export
#'
#' @examples
#' \dontrun{
#' checkNegCtrlGenes(data =brca.data, ncg_set= c("Microrray_HK"), group='Time', plot_type="DensityPlot", nPcs=10, npcs = 3, is.log=FALSE)
#' checkNegCtrlGenes(data =brca.data, ncg_set= c("Microrray_HK","scRNAseq_HK"), group='Tissue', plot_type="BoxPlot", nPcs=7, npcs = 2, is.log=FALSE)
#' }
checkNegCtrlGenes <- function(data, ncg_set, group, plot_type, nPcs, npcs, is.log){
  gene.annot.rm <-  as.data.frame(SummarizedExperiment::rowData(data))
  hk_genes <- c("scRNAseq_HK","RNAseq_HK","Microrray_HK","NanostringPanCancer_HK","sinscorePanCancer_HK")
  nhk_genes <- setdiff(names(gene.annot.rm), hk_genes)
  names.use <- hk_genes[(hk_genes %in% ncg_set)]
  keep.ncg.genes <- c(nhk_genes,names.use)
  SummarizedExperiment::rowData(data) <- SummarizedExperiment::rowData(data)[names(SummarizedExperiment::rowData(data)) %in% keep.ncg.genes]
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
  tcga.harmonized <- tcga.harmonized[1]
  pca.cancer.tcga  <- lapply(
    tcga.harmonized,
    function(x){
      .pca(
        data = as.matrix(SummarizedExperiment::assay(data, x)),
        nPcs = nPcs,
        is.log = FALSE)
    })
  names(pca.cancer.tcga) <- tcga.harmonized

  .scatter.density.pc <- function(
    pcs,
    pc.var,
    group.name,
    group,
    color,
    strokeSize,
    pointSize,
    strokeColor,
    alpha){
    pair.pcs <- utils::combn(ncol(pcs), 2)
    pList <- list()
    for(i in 1:ncol(pair.pcs)){
      if(i == 1){
        x <- pair.pcs[1,i]
        y <- pair.pcs[2,i]
        p <- ggplot2::ggplot(mapping = ggplot2::aes(
          x = pcs[,x],
          y = pcs[,y],
          fill = group)) +
          ggplot2::xlab(paste0('PC', x, ' (', pc.var[x], '%)')) +
          ggplot2::ylab(paste0('PC', y, ' (', pc.var[y], '%)')) +
          ggplot2::geom_point(
            ggplot2::aes(fill = group),
            pch = 21,
            color = strokeColor,
            stroke = strokeSize,
            size = pointSize,
            alpha = alpha) +
          ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
          ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
          ggplot2::theme(
            legend.position = "right",
            panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black", size = 1.1),
            legend.background = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 12),
            legend.title = ggplot2::element_text(size = 14),
            legend.key = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 10),
            axis.text.y = ggplot2::element_text(size = 10),
            axis.title.x = ggplot2::element_text(size = 14),
            axis.title.y = ggplot2::element_text(size = 14)) +
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 4))) +
          ggplot2::scale_fill_manual(name = group.name, values = color)

        le <- ggpubr::get_legend(p)
      }else{
        x <- pair.pcs[1,i]
        y <- pair.pcs[2,i]
        p <- ggplot2::ggplot(mapping = ggplot2::aes(
          x = pcs[,x],
          y = pcs[,y],
          fill = group)) +
          ggplot2::xlab(paste0('PC', x, ' (',pc.var[x],  '%)')) +
          ggplot2::ylab(paste0('PC', y, ' (',pc.var[y], '%)')) +
          ggplot2::geom_point(
            ggplot2::aes(fill = group),
            pch = 21,
            color = strokeColor,
            stroke = strokeSize,
            size = pointSize,
            alpha = alpha) +
          ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
          ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
          ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black", size = 1.1),
            legend.position = "none",
            axis.text.x = ggplot2::element_text(size = 10),
            axis.text.y = ggplot2::element_text(size = 10),
            axis.title.x = ggplot2::element_text(size = 14),
            axis.title.y = ggplot2::element_text(size = 14)) +
          ggplot2::scale_fill_manual(values = color, name = group.name)
      }
      p <- p + ggplot2::theme(legend.position = "none")
      xdens <- cowplot::axis_canvas(p, axis = "x")+
        ggplot2::geom_density(
          mapping = ggplot2::aes(
            x = pcs[,x],
            fill = group),
          alpha = 0.7,
          size = 0.2
        ) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::scale_fill_manual(values = color)

      ydens <- cowplot::axis_canvas(
        p,
        axis = "y",
        coord_flip = TRUE) +
        ggplot2::geom_density(
          mapping = ggplot2::aes(
            x = pcs[,y],
            fill = group),
          alpha = 0.7,
          size = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::scale_fill_manual(name = group.name, values = color) +
        ggplot2::coord_flip()

      p1 <- insert_xaxis_grob(
        p,
        xdens,
        grid::unit(.2, "null"),
        position = "top"
      )
      p2 <- insert_yaxis_grob(
        p1,
        ydens,
        grid::unit(.2, "null"),
        position = "right"
      )
      pList[[i]] <- ggdraw(p2)
    }
    pList[[i+1]] <- le
    return(pList)
  }

  sample.info <-  as.data.frame(SummarizedExperiment::colData(data))
  data.set.names <- tcga.harmonized
  pca.data <- pca.cancer.tcga

  if (group == "Time"){
    if(plot_type == "DensityPlot"){
      currentCols <- RColorBrewer::brewer.pal(8, "Dark2")[-5]
      years.colors <- currentCols[1:length(unique(sample.info$Year))]
      pca.plots.time <- lapply(
        data.set.names,#tcga.harmonized,
        function(x){
          pcs <- pca.data[[x]]
          p <- .scatter.density.pc(
            pcs = pcs$sing.val$u[,1:npcs],
            pc.var = pcs$var,
            group.name = 'Time',
            group = sample.info$Year,
            color = years.colors,#unique(sample.info$Year),
            strokeSize = .2,
            pointSize = 3,
            strokeColor = 'gray30',
            alpha = .6)
          do.call(
            gridExtra::grid.arrange,
            c(p,
              ncol = 4, top = x)
          )
        })
    } else if (plot_type == "BoxPlot"){
      for (i in 1:npcs){
        boxplot(pca.data$HTseq_counts$sing.val$u[,i] ~ sample.info$Year)
      }
    }
  }
  else
    if (group == "Tissue"){
      if(plot_type == "DensityPlot"){
        pca.plots.time <- lapply(
          data.set.names,
          function(x){
            pcs <- pca.data[[x]]
            p <- .scatter.density.pc(
              pcs = pcs$sing.val$u[,1:npcs],
              pc.var = pcs$var,
              group.name = 'Tissue',
              group = sample.info$Tissues,
              color = c("#252525","#D9D9D9"),#c('red', 'blue'),
              strokeSize = .2,
              pointSize = 3,
              strokeColor = 'gray30',
              alpha = .6)
            do.call(
              gridExtra::grid.arrange,
              c(p,
                ncol = 4,
                top = x)
            )
          })
      } else if (plot_type == "BoxPlot"){
        for (i in 1:npcs){
          boxplot(pca.data$HTseq_counts$sing.val$u[,i] ~ sample.info$Tissues)
        }
      }
    } else
      if (group == "Plate"){
        if(plot_type == "DensityPlot"){
          pca.plots.time <- lapply(
            data.set.names,
            function(x){
              pcs <- pca.data[[x]]
              p <- .scatter.density.pc(
                pcs = pcs$sing.val$u[,1:npcs],
                pc.var = pcs$var,
                group.name = 'Plate',
                group = sample.info$Plates,
                color = unique(factor(sample.info$Plates)),
                strokeSize = .2,
                pointSize = 3,
                strokeColor = 'gray30',
                alpha = .6)
              do.call(
                gridExtra::grid.arrange,
                c(p,
                  ncol = 4,
                  top = x)
              )
            })
        } else if (plot_type == "BoxPlot"){
          for (i in 1:npcs){
            boxplot(pca.data$HTseq_counts$sing.val$u[,i] ~ sample.info$Plates)
          }
        }
      } else
        if (group == "TSS"){
          if(plot_type == "DensityPlot"){
            pca.plots.time <- lapply(
              data.set.names,
              function(x){
                pcs <- pca.data[[x]]
                p <- .scatter.density.pc(
                  pcs = pcs$sing.val$u[,1:npcs],
                  pc.var = pcs$var,
                  group.name = 'TSS',
                  group = sample.info$TSS,
                  color = unique(factor(sample.info$TSS)),
                  strokeSize = .2,
                  pointSize = 3,
                  strokeColor = 'gray30',
                  alpha = .6)
                do.call(
                  gridExtra::grid.arrange,
                  c(p,
                    ncol = 4,
                    top = x)
                )
              })
          } else if (plot_type == "BoxPlot"){
            for (i in 1:npcs){
              boxplot(pca.data$HTseq_counts$sing.val$u[,i] ~ sample.info$TSS)
            }
          }
        } else
          if (group == "Center"){
            if(plot_type == "DensityPlot"){
              pca.plots.time <- lapply(
                data.set.names,
                function(x){
                  pcs <- pca.data[[x]]
                  p <- .scatter.density.pc(
                    pcs = pcs$sing.val$u[,1:npcs],
                    pc.var = pcs$var,
                    group.name = 'Center',
                    group = sample.info$Center,
                    color = unique(sample.info$Center),
                    strokeSize = .2,
                    pointSize = 3,
                    strokeColor = 'gray30',
                    alpha = .6)
                  do.call(
                    gridExtra::grid.arrange,
                    c(p,
                      ncol = 4,
                      top = x)
                  )
                })
            } else if (plot_type == "BoxPlot"){
              for (i in 1:npcs){
                boxplot(pca.data$HTseq_counts$sing.val$u[,i] ~ sample.info$Center)
              }
            }
          }
}
