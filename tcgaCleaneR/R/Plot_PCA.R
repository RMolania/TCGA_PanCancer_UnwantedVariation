# Plot PCA function

#' @title PCA Visualization
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`. It helps to visualize the PCs from \code{get.pca}.
#'
#' @param pca.data list: PCA output from \code{computePCA}.
#' @param data S4 data object
#' @param group character: Color code PCs based on group. groups included are 'Time', 'Tissue', 'Plate', 'TSS', 'Center'
#' @param plot_type character: Plot type
#' @param pcs.no numeric vector: PCs that needs to be plotted
#'
#' @return Density Plot, Box Plot
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom cowplot axis_canvas
#' @importFrom cowplot insert_xaxis_grob
#' @importFrom cowplot insert_yaxis_grob
#'
#' @examples
#' \dontrun{
#' plotPC(pca.data, data = brca.data, group = "Time", plot_type = "DensityPlot", pcs.no = c(1,2,3))
#' plotPC(pca.data = df6, data = brca.data, group = "Plate", plot_type = "BoxPlot", pcs.no = c(1,2,4))
#' }
plotPC <- function(pca.data, data, group, plot_type, pcs.no){
  .scatter.density.pc <- function(
    pcs,
    pc.var,
    pcs.no,
    group.name,
    group,
    color,
    strokeSize,
    pointSize,
    strokeColor,
    alpha){
    colnames(pcs) <- paste0('pc', pcs.no)
    pair.pcs <- utils::combn(ncol(pcs), 2)
    pcs.var <- utils::combn(pcs.no, 2)
    pList <- list()
    for(i in 1:ncol(pair.pcs)){
      if(i == 1){
        x <- pair.pcs[1,i]
        y <- pair.pcs[2,i]
        a <- pcs.var[1,i]
        b <- pcs.var[2,i]
        p <- ggplot(mapping = aes(
          x = pcs[,x],
          y = pcs[,y],
          fill = group)) +
          xlab(paste0('PC', a, ' (', pc.var[a], '%)')) +
          ylab(paste0('PC', b, ' (', pc.var[b], '%)')) +
          geom_point(
            aes(fill = group),
            pch = 21,
            color = strokeColor,
            stroke = strokeSize,
            size = pointSize,
            alpha = alpha) +
          scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
          theme(
            legend.position = "right",
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size = 1.1),
            legend.background = element_blank(),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.key = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14)) +
          guides(fill = guide_legend(override.aes = list(size = 4))) +
          scale_fill_manual(name = group.name, values = color)

        le <- ggpubr::get_legend(p)
      }else{
        x <- pair.pcs[1,i]
        y <- pair.pcs[2,i]
        a <- pcs.var[1,i]
        b <- pcs.var[2,i]
        p <- ggplot(mapping = aes(
          x = pcs[,x],
          y = pcs[,y],
          fill = group)) +
          xlab(paste0('PC', a, ' (',pc.var[a],  '%)')) +
          ylab(paste0('PC', b, ' (',pc.var[b], '%)')) +
          geom_point(
            aes(fill = group),
            pch = 21,
            color = strokeColor,
            stroke = strokeSize,
            size = pointSize,
            alpha = alpha) +
          scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
          theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size = 1.1),
            legend.position = "none",
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14)) +
          scale_fill_manual(values = color, name = group.name)
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
  data.set.names <- names(SummarizedExperiment::assays(data))

  if (group == "Time"){
    if(plot_type == "DensityPlot"){
      currentCols <- RColorBrewer::brewer.pal(8, "Dark2")[-5]
      years.colors <- currentCols[1:length(unique(sample.info$Year))]
      pca.plots.time <- lapply(
        data.set.names,#tcga.harmonized,
        function(x){
          pcs <- pca.data[[x]]
          p <- .scatter.density.pc(
            pcs = pcs$sing.val$u[,pcs.no],
            pc.var = pcs$var,
            pcs.no = pcs.no,
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

      for (i in data.set.names){
        to.plot.pc <- pca.data[[i]]$sing.val$u
        to.plot.pc <- as.data.frame(to.plot.pc)
        to.plot.pc <- to.plot.pc[,pcs.no]
        colnames(to.plot.pc) <- paste0('pc', pcs.no)
        to.plot.pc$variable = sample.info$Year
        to.plot.pc <- to.plot.pc %>%
          tidyr::pivot_longer(-variable, names_to = 'pcs', values_to = 'var') %>%
          data.frame(.)
        p <- ggplot(to.plot.pc, aes(x = variable, y = var)) +
          geom_boxplot()+
          facet_wrap(~pcs) +
          ylab('PC') +
          ggtitle(i) +
          theme(
            panel.background = element_blank(),
            plot.title = element_text(size = 22),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            strip.text.x = element_text(size = 20),
            strip.text = element_text(size = 22))
        print(p)
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
              pcs = pcs$sing.val$u[,pcs.no],
              pc.var = pcs$var,
              pcs.no = pcs.no,
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

        for (i in data.set.names){
          to.plot.pc <- pca.data[[i]]$sing.val$u
          to.plot.pc <- as.data.frame(to.plot.pc)
          to.plot.pc <- to.plot.pc[,pcs.no]
          colnames(to.plot.pc) <- paste0('pc', pcs.no)
          to.plot.pc$variable = sample.info$Tissues
          to.plot.pc <- to.plot.pc %>%
            tidyr::pivot_longer(-variable, names_to = 'pcs', values_to = 'var') %>%
            data.frame(.)
          p <- ggplot(to.plot.pc, aes(x = variable, y = var)) +
            geom_boxplot()+
            facet_wrap(~pcs) +
            ylab('PC') +
            ggtitle(i) +
            theme(
              panel.background = element_blank(),
              plot.title = element_text(size = 22),
              axis.line = element_line(colour = 'black', size = 1),
              axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14),
              strip.text.x = element_text(size = 20),
              strip.text = element_text(size = 22))
          print(p)
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
                pcs = pcs$sing.val$u[,pcs.no],
                pc.var = pcs$var,
                pcs.no = pcs.no,
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

          for (i in data.set.names){
            to.plot.pc <- pca.data[[i]]$sing.val$u
            to.plot.pc <- as.data.frame(to.plot.pc)
            to.plot.pc <- to.plot.pc[,pcs.no]
            colnames(to.plot.pc) <- paste0('pc', pcs.no)
            to.plot.pc$variable = sample.info$Plates
            to.plot.pc <- to.plot.pc %>%
              tidyr::pivot_longer(-variable, names_to = 'pcs', values_to = 'var') %>%
              data.frame(.)
            p <- ggplot(to.plot.pc, aes(x = variable, y = var)) +
              geom_boxplot()+
              facet_wrap(~pcs) +
              ylab('PC') +
              ggtitle(i) +
              theme(
                panel.background = element_blank(),
                plot.title = element_text(size = 22),
                axis.line = element_line(colour = 'black', size = 1),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                strip.text.x = element_text(size = 20),
                strip.text = element_text(size = 22))
            print(p)
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
                  pcs = pcs$sing.val$u[,pcs.no],
                  pc.var = pcs$var,
                  pcs.no = pcs.no,
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

            for (i in data.set.names){
              to.plot.pc <- pca.data[[i]]$sing.val$u
              to.plot.pc <- as.data.frame(to.plot.pc)
              to.plot.pc <- to.plot.pc[,pcs.no]
              colnames(to.plot.pc) <- paste0('pc', pcs.no)
              to.plot.pc$variable = sample.info$TSS
              to.plot.pc <- to.plot.pc %>%
                tidyr::pivot_longer(-variable, names_to = 'pcs', values_to = 'var') %>%
                data.frame(.)
              p <- ggplot(to.plot.pc, aes(x = variable, y = var)) +
                geom_boxplot()+
                facet_wrap(~pcs) +
                ylab('PC') +
                ggtitle(i) +
                theme(
                  panel.background = element_blank(),
                  plot.title = element_text(size = 22),
                  axis.line = element_line(colour = 'black', size = 1),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16),
                  axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  legend.text = element_text(size = 12),
                  legend.title = element_text(size = 14),
                  strip.text.x = element_text(size = 20),
                  strip.text = element_text(size = 22))
              print(p)
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
                    pcs = pcs$sing.val$u[,pcs.no],
                    pc.var = pcs$var,
                    pcs.no = pcs.no,
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

              for (i in data.set.names){
                to.plot.pc <- pca.data[[i]]$sing.val$u
                to.plot.pc <- as.data.frame(to.plot.pc)
                to.plot.pc <- to.plot.pc[,pcs.no]
                colnames(to.plot.pc) <- paste0('pc', pcs.no)
                to.plot.pc$variable = sample.info$Center
                to.plot.pc <- to.plot.pc %>%
                  tidyr::pivot_longer(-variable, names_to = 'pcs', values_to = 'var') %>%
                  data.frame(.)
                p <- ggplot(to.plot.pc, aes(x = variable, y = var)) +
                  geom_boxplot()+
                  facet_wrap(~pcs) +
                  ylab('PC') +
                  ggtitle(i) +
                  theme(
                    panel.background = element_blank(),
                    plot.title = element_text(size = 22),
                    axis.line = element_line(colour = 'black', size = 1),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 12),
                    legend.title = element_text(size = 14),
                    strip.text.x = element_text(size = 20),
                    strip.text = element_text(size = 22))
                print(p)
              }
            }
          }
  #return(data.set.names)
}
