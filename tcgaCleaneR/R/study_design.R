# Study Design Function

#' @title Study Design Plot
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`. It helps to plot
#' a Combined Heat Map \code{ComplexHeatmap} for some known unwanted variations and present a dashboard level analysis
#' of the TCGA data before any analysis.
#'
#' @param data S4 data object: Input TCGA or similar \code{SummarizedExperiment} Dataset
#'
#' @return Heat Map Plot
#' @export
#'
#' @examples
#' \dontrun{
#' plotStudyOutline(data = brca.data)
#' }
#'
plotStudyOutline <- function(data){
  data$ls <- log2(colSums(SummarizedExperiment::assay(data, 'HTseq_counts')))
  sample.info <- as.data.frame(SummarizedExperiment::colData(data))

  currentCols <-  c(
    RColorBrewer::brewer.pal(8, "Dark2")[-5],
    RColorBrewer::brewer.pal(10, "Paired"),
    RColorBrewer::brewer.pal(12, "Set3"),
    RColorBrewer::brewer.pal(9, "Blues")[c(8, 3, 7, 4, 6, 9, 5)],
    RColorBrewer::brewer.pal(9, "Oranges")[c(8, 3, 7, 4, 6, 9, 5)],
    RColorBrewer::brewer.pal(9, "Greens")[c(8, 3, 7, 4, 6, 9, 5)],
    RColorBrewer::brewer.pal(9, "Purples")[c(8, 3, 7, 4, 6, 9, 5)],
    RColorBrewer::brewer.pal(9, "Reds")[c(8, 3, 7, 4, 6, 9, 5)],
    RColorBrewer::brewer.pal(9, "Greys")[c(8, 3, 7, 4, 6, 9, 5)],
    RColorBrewer::brewer.pal(9, "BuGn")[c(8, 3, 7, 4, 6, 9, 5)],
    RColorBrewer::brewer.pal(9, "PuRd")[c(8, 3, 7, 4, 6, 9, 5)],
    RColorBrewer::brewer.pal(9, "BuPu")[c(8, 3, 7, 4, 6, 9, 5)],
    RColorBrewer::brewer.pal(9, "YlGn")[c(8, 3, 7, 4, 6, 9, 5)],
    RColorBrewer::brewer.pal(10, "Paired")
  )


  cols <- c(
    'Year',
    'Plates',
    'TSS',
    'Tissues',
    'Center',
    'ls',
    'Purity_singscore',
    'Tumor.stage'
  )

  sample.info <- sample.info[ , cols]
  years.colors <- currentCols[1:length(unique(sample.info$Year))]
  ### Year
  H.time <- ComplexHeatmap::Heatmap(
    rev(sample.info$Year),
    cluster_columns  = FALSE,
    column_names_gp = grid::gpar(fontsize = 18),
    col =  years.colors,
    name = 'Time (years)',
    heatmap_legend_param = list(
      color_bar = "discrete" ,
      ncol = 2,
      title_gp = grid::gpar(fontsize = 18)
    )
  )

  ### Plates
  n.plates <- length(unique(sample.info$Plates)) # 38
  colfunc <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, 'PRGn')[-6])
  color.plates <- colfunc(n.plates)

  H.plate <- ComplexHeatmap::Heatmap(
    rev(sample.info$Plates),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_gp = grid::gpar(fontsize = 18),
    col = color.plates,
    name = 'Plates',
    heatmap_legend_param = list(
      color_bar = "discrete" ,
      ncol = 4,
      title_gp = grid::gpar(fontsize = 18)
    )
  )

  ### TSS
  n.tss <- length(unique(sample.info$TSS)) # 40
  colfunc <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, 'BrBG')[-6])
  color.tss <- colfunc(n.tss)
  H.tss <- ComplexHeatmap::Heatmap(
    rev(sample.info$TSS),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_gp = grid::gpar(fontsize = 18),
    col = color.tss,
    name = 'Tissue source sites',
    heatmap_legend_param = list(
      color_bar = "discrete" ,
      ncol = 4,
      title_gp = grid::gpar(fontsize = 18)
    )
  )

  ### Tissue
  H.tissue <- ComplexHeatmap::Heatmap(
    rev(sample.info$Tissues),
    cluster_rows = FALSE,
    column_names_gp = grid::gpar(fontsize = 18),
    col = c("#252525", "#D9D9D9"),
    name = 'Tissues',
    heatmap_legend_param = list(
      color_bar = "discrete" ,
      direction = "vertical",
      ncol = 1,
      title_gp = grid::gpar(fontsize = 18),
      labels = c(
        'Cancer',
        'Normal')
    )
  )
  ### Purity
  H.purity <- ComplexHeatmap::Heatmap(
    rev(sample.info$Purity_singscore),
    column_names_gp = grid::gpar(fontsize = 18),
    cluster_rows = FALSE,
    name = 'Tumor purity score',
    col = viridis::plasma(n = 10),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 18)
    )
  )

  ### library size
  H.ls <- ComplexHeatmap::Heatmap(
    rev(sample.info$ls),
    cluster_rows = FALSE,
    name = 'Library size',
    column_names_gp = grid::gpar(fontsize = 18),
    col = viridis::viridis(n = 10),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 18)
    )
  )

  ### Tumor Stage
  tumor.stage.colors <- currentCols[1:length(unique(sample.info$Tumor.stage))]
  H.tumor.stage <- ComplexHeatmap::Heatmap(
    rev(sample.info$Tumor.stage),
    cluster_columns  = FALSE,
    column_names_gp = grid::gpar(fontsize = 18),
    col =  tumor.stage.colors,
    name = 'Tumor Stage',
    heatmap_legend_param = list(
      color_bar = "discrete" ,
      ncol = 2,
      title_gp = grid::gpar(fontsize = 18)
    )
  )

  ComplexHeatmap::draw(
    H.time +
      H.plate +
      H.tss +
      H.tissue +
      H.purity +
      H.ls +
      H.tumor.stage,
    merge_legends = FALSE,
    heatmap_legend_side = 'left'
  )
}
