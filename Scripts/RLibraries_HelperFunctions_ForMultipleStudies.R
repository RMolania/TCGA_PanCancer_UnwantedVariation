

library(cowplot)
library(ggplot2)
library(mclust)
library(ruv)


.rle.comp <- function(expr.data, is.log) {
  if (is.log)
    expr.data <- expr.data
  else
    expr.data <- log2(expr.data + 1)
  rle.data <- expr.data - matrixStats::rowMedians(expr.data)
  rle.med <- matrixStats::colMedians(rle.data)
  rle.iqr <- matrixStats::colIQRs(rle.data)
  return(list(
    rle = rle.data,
    rle.med = rle.med,
    rle.iqr = rle.iqr
  ))
}


.CreatePseudoSamplesForLsPurityBatch <- function(
  expr.data,
  sample.info,
  librarySize,
  biology,
  batch,
  purity,
  include.ls = FALSE,
  include.purity = FALSE,
  minSamplesPerBatchPS = 3,
  minSamplesForPurityPerBiology = 12,
  minSamplesForPurityPS = 3,
  minSamplesForLibrarySizePerBatch = 10,
  minSamplesForLibrarySizePS = 3
){
  ### checl
  if(include.purity & minSamplesForPurityPS > minSamplesForPurityPerBiology){
    stop('error: minSamplesForPurityPS can not be smaller than minSamplesForPurityPerBiology')
  } else if(include.purity & minSamplesForPurityPerBiology < 2*minSamplesForPurityPS){
    stop('error: minSamplesForPurityPerBiology should be at least two times larger than minSamplesForPurityPS')
  } else if(include.ls & minSamplesForLibrarySizePS > minSamplesForLibrarySizePerBatch) {
    stop('error: minSamplesForLibrarySizePerBatch can not be smaller than minSamplesForLibrarySizePS')
  } else if(include.ls & minSamplesForLibrarySizePerBatch < 2*minSamplesForLibrarySizePS ){
    stop('error: minSamplesForLibrarySizePerBatch should be at least two times larger than minSamplesForLibrarySizePS')
  }
  ### Biology
  row.names(sample.info) <- colnames(expr.data)
  sample.info$biology <- apply(
    sample.info[ , biology, drop = FALSE],
    1,
    paste,
    collapse = "-"
  )
  ### Biology - Batch
  sample.info$biology.batch <- apply(
    sample.info[, c(biology, batch)],
    1,
    paste,
    collapse = "_"
  )
  ### removing batch effects
  # create PS per biology/batch
  selected.biology.ps.batch <- unlist(lapply(
    unique(sample.info$biology), 
    function(x){
      index <- sample.info$biology == x
      if(sum( table(sample.info$biology.batch[index] ) >= minSamplesPerBatchPS) > 1 ){
        x
      }
    }))
  if(length(selected.biology.ps.batch) > 0){
    message('PRPS are generated for batch effects')
  }else{
    message('error: there are not enough samples to create 
            pseudo-samples for batch effects removal, you may want to lower minSamplesPerBatchPS')
  }
  sample.info.ps.batch <- sample.info[sample.info$biology %in% selected.biology.ps.batch , ]
  expr.data.ps.batch <- expr.data[, row.names(sample.info.ps.batch)]
  ### sort samples
  selected.batches <- names(which(table(sample.info.ps.batch$biology.batch) >= minSamplesPerBatchPS))
  ps.batch <- sapply(
    selected.batches,
    function(x) {
      index <- sample.info.ps.batch$biology.batch == x
      Matrix::rowMeans(expr.data.ps.batch[, index])
    })
  
  if(include.ls){
    selected.batches.ls <- names(
      which(table(sample.info$biology.batch) >= minSamplesForLibrarySizePerBatch)
    )
    if(length(selected.batches.ls) > 0){
      message('PRPS are generated for library size effects')
      sample.info <- sample.info[
        with(sample.info,
             order(sample.info[, 'biology.batch'],
                   sample.info[, librarySize])), ]
      expr.data <- expr.data[, row.names(sample.info)]
      ps.ls <- lapply(
        selected.batches.ls, 
        function(x){
          index <- sample.info$biology.batch == x
          ls.data <- expr.data[ , index]
          low.ls <- Matrix::rowMeans(ls.data[ , 1:minSamplesForLibrarySizePS])
          high.ls <- Matrix::rowMeans(ls.data[ , c(ncol(ls.data)-(minSamplesForLibrarySizePS - 1)):ncol(ls.data) ])
          all.ls <- cbind(low.ls, high.ls)
          colnames(all.ls) <- rep(paste(x, 'LS', sep = '_'), 2)
          all.ls
        })
      ps.ls <- do.call(cbind, ps.ls)
      
    }else{
      message('error: there are not enough samples to create pseudo-samples for 
              removal of library size effects, you may want to lower minSamplesForLibrarySizePerBatch')
    }
  }else if (! include.ls){
    print('PRPS is not generated for librray size effects')
    ps.ls = list()
  }
  if(include.purity ){
    selected.biology.purity <- names(
      which(table(sample.info$biology) >= minSamplesForPurityPerBiology)
    ) 
    sample.info$purity <- sample.info[ , 'purity']
    if(length(selected.biology.purity) > 0){
      message('PRPS are generated for purity effects')
      sample.info <- sample.info[
        with(sample.info, order(biology, purity)), ]
      expr.data <- expr.data[, row.names(sample.info)]
      ps.purity <- lapply(
        selected.biology.purity,
        function(x) {
          index <- sample.info$biology == x
          purity.data <- expr.data[, index]
          low.pur <- Matrix::rowMeans(purity.data[, 1:minSamplesForPurityPS])
          high.pur <- Matrix::rowMeans(purity.data[, c(ncol(purity.data) - (minSamplesForPurityPS - 1)):ncol(purity.data)])
          all.purity <- cbind(low.pur, high.pur)
          colnames(all.purity) <- rep(paste(x, 'purity', sep = '_'), 2)
          all.purity
        })
      ps.purity <- do.call(cbind, ps.purity)
    }else{
      message('error: there are not enough samples to make pseudo-samples 
              for purity variation, you may want to lower minSamplesForPurityPerBiology')
    }
  } else if (!include.purity){
    print('PRPS is not generated for purity effects')
    ps.purity = list()
  }
  return(list(ps.batch = ps.batch, ps.ls = ps.ls, ps.purity = ps.purity))
}


library(BiocParallel)
RUV_III_PRPS_fast <- function(
  Y,
  M,
  ctl,
  k = NULL, 
  eta = NULL, 
  include.intercept = TRUE,
  average = FALSE, 
  fullalpha = NULL, 
  return.info = FALSE, 
  inputcheck = TRUE) 
{
  if (is.data.frame(Y) ) {
    Y <- data.matrix(Y)
  }
  m <- nrow(Y)
  n <- ncol(Y)
  M <- ruv::replicate.matrix(M)
  ctl <- tological(ctl, n)
  if (inputcheck) {
    if (m > n)
      warning(
        "m is greater than n!  This is not a problem itself, but may 
              indicate that you need to transpose your data matrix.  
              Please ensure that rows correspond to observations 
              (e.g. RNA-Seq assay) and columns correspond to features (e.g. genes).")
    if (sum(is.na(Y)) > 0)
      warning("Y contains missing values.  This is not supported.")
    if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) >
        0)
      warning("Y contains infinities.  This is not supported.")
  }
  Y <- ruv::RUV1(Y, eta, ctl, include.intercept = include.intercept)
  mu <- colMeans(Y)
  mu_mat <- rep(1, m) %*% t(mu)
  Y_stand <- Y - mu_mat
  if (ncol(M) >= m)
    newY <- Y
  else if (is.null(k)) {
    ycyctinv <- solve(Y[, ctl] %*% t(Y[, ctl]))
    newY <- (M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*% ycyctinv)) %*% Y
    fullalpha <- NULL
  } else if (k == 0) {
    newY <- Y
    fullalpha <- NULL
  } else {
    if (is.null(fullalpha) ) {
      Y0 <- fast.residop(Y, M)
      eig.vec <- svdObj <- BiocSingular::runSVD(
        x = Y0, 
        k = k, 
        BSPARAM = BiocSingular::bsparam()
        )
      fullalpha <- DelayedArray::t(svdObj$u[, 1:k, drop = FALSE]) %*% Y
    }
    alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
    ac <- alpha[, ctl, drop = FALSE]
    W <- Y_stand[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
    newY <- Y - W %*% alpha
  }
  if (average)
    newY <- ((1/apply(M, 2, sum)) * t(M)) %*% newY
  if (!return.info) {
    return(newY)
  } else {
    return(list(newY = newY, M = M, fullalpha = fullalpha,  W =  W))
  }
}

tological <- function(ctl, n) {
  ctl2 <- rep(FALSE, n)
  ctl2[ctl] <- TRUE
  return(ctl2)
}

fast.residop <- function(A, B){
  tBB = DelayedArray::t(B) %*% B
  tBB_inv = Matrix::solve(tBB)
  BtBB_inv = B %*% tBB_inv
  tBA = DelayedArray::t(B) %*% A
  BtBB_inv_tBA = BtBB_inv %*% tBA
  return(A - BtBB_inv_tBA)
}


library(BiocSingular)
.pca.fast <- function(data, nPcs, is.log) {
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

.silhouette.coeff <- function(
  pcs, 
  variable, 
  nPCs)
{
  d.matrix <- as.matrix(stats::dist(pcs[, seq_len(nPCs)]))
  summary(cluster::silhouette(
    as.numeric(as.factor(variable)), 
    d.matrix))$avg.width
}



library(ruv)
fastRUV_III_PRPS <- function (
  Y,
  Yrep,
  M, 
  ctl, 
  k = NULL, 
  eta = NULL, 
  include.intercept = TRUE,
  average = FALSE, 
  fullalpha = NULL, 
  return.info = FALSE, 
  inputcheck = TRUE)
{
  if (is.data.frame(Y))
    Y = data.matrix(Y)
  Yrep = data.matrix(Yrep)
  m = nrow(Yrep)
  m1 = nrow(Y)
  n = ncol(Y)
  M = replicate.matrix(M)
  ctl = tological(ctl, n)
  if (inputcheck) {
    if (m > n)
      warning("m is greater than n!  This is not a problem itself, but may
                indicate that you need to transpose your data matrix.
                Please ensure that rows correspond to observations (e.g. microarrays)
                and columns correspond to features (e.g. genes).")
    if (sum(is.na(Y)) > 0)
      warning("Y contains missing values.  This is not supported.")
    if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) >
        0)
      warning("Y contains infinities.  This is not supported.")
  }
  Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
  Yrep = RUV1(Yrep, eta, ctl, include.intercept = include.intercept)
  mu <- DelayedArray::colMeans(Y)
  mu_mat <- rep(1, m1) %*% t(mu)
  Y_stand <- Y - mu_mat
  if (ncol(M) >= m)
    newY = Y
  else if (is.null(k)) {
    ycyctinv = solve(Y[, ctl] %*% t(Y[, ctl]))
    newY = (M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*% ycyctinv)) %*% Y
    fullalpha = NULL
  }
  else if (k == 0) {
    newY = Y
    fullalpha = NULL
  }
  else {
    if (is.null(fullalpha)) {
      Y0 = residopFast(Yrep, M)
      k.eigVec = min(m - ncol(M), sum(ctl))
      eigVec = base::eigen(Y0 %*% DelayedArray::t(Y0), symmetric = TRUE)
      # eigVec = irlba::partial_eigen(
      #   Y0 %*% DelayedArray::t(Y0), 
      #   symmetric = TRUE, 
      #   n = k.eigVec
      #   )
      fullalpha = DelayedArray::t(eigVec$vectors[, seq_len(k.eigVec), drop = FALSE]) %*% Yrep
    }
    alpha = fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
    ac = alpha[, ctl, drop = FALSE]
    W = Y_stand[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
    newY = Y - W %*% alpha
  }
  if (average)
    newY = ((1/apply(M, 2, sum)) * t(M)) %*% newY
  if (!return.info)
    return(newY)
  else return(list(
    newY = newY,
    M = M,
    fullalpha = fullalpha,
    W = W,
    WA = W %*% alpha,
    alpha = alpha))
}


tological <- function (ctl, n)
{
  ctl2 = rep(FALSE, n)
  ctl2[ctl] = TRUE
  return(ctl2)
}


residopFast <- function(A, B){
  tBB = DelayedArray::t(B) %*% B
  tBB_inv = Matrix::solve(tBB)
  BtBB_inv = B %*% tBB_inv
  tBA = DelayedArray::t(B) %*% A
  BtBB_inv_tBA = BtBB_inv %*% tBA
  return(A - BtBB_inv_tBA)
}



# PCA plot with  with density ####
## PCS
## pc.var
## group.name
## group
## color
## strokeSize
## pointSize
## strokeColor
## alpha

.scatter.density.pc <- function(
  pcs, 
  pc.var, 
  group.name, 
  group, 
  color, 
  strokeSize, 
  pointSize, 
  strokeColor,
  alpha
){
  pair.pcs <- utils::combn(ncol(pcs), 2)
  pList <- list()
  for(i in 1:ncol(pair.pcs)){
    if(i == 1){
      x <- pair.pcs[1,i]
      y <- pair.pcs[2,i]
      p <- ggplot(mapping = aes(
        x = pcs[,x], 
        y = pcs[,y], 
        fill = group)) +
        xlab(paste0('PC', x, ' (', pc.var[x], '%)')) +
        ylab(paste0('PC', y, ' (', pc.var[y], '%)')) +
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
      p <- ggplot(mapping = aes(
        x = pcs[,x], 
        y = pcs[,y], 
        fill = group)) +
        xlab(paste0('PC', x, ' (',pc.var[x],  '%)')) +
        ylab(paste0('PC', y, ' (',pc.var[y], '%)')) +
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
    p <- p + theme(legend.position = "none")
    xdens <- cowplot::axis_canvas(p, axis = "x")+
      geom_density(
        mapping = aes(
          x = pcs[,x], 
          fill = group),
        alpha = 0.7, 
        size = 0.2
      ) +
      theme(legend.position = "none") +
      scale_fill_manual(values = color)
    
    ydens <- cowplot::axis_canvas(
      p, 
      axis = "y", 
      coord_flip = TRUE) +
      geom_density(
        mapping = aes(
          x = pcs[,y],
          fill = group),
        alpha = 0.7,
        size = 0.2) +
      theme(legend.position = "none") +
      scale_fill_manual(name = group.name, values = color) +
      coord_flip()
    
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


# RLE plot plot with  with density ####
.plot.rle <- function(
  rle.data, 
  variable, 
  rle.med , 
  zero.line.color, 
  rle.med.color, 
  show.tech.rep, 
  rle.med.point.size,
  rle.med.point.lwd,
  curve.lty,
  curve.lwd, 
  plot.lwd,
  curve.color, 
  ...
){
  if(!show.tech.rep){
    boxplot(
      rle.data,
      xaxt = 'n',
      yaxt = 'n' ,
      outline = FALSE,
      frame = FALSE,
      whisklty = 0,
      ylab = '',
      xlab = '',
      staplelty = 0,
      names = FALSE,
      ...
    )
    points(
      c(1:ncol(rle.data)),
      rle.med,
      bg = rle.med.color[factor(variable)],
      col = 'black',
      pch = 21,
      cex = rle.med.point.size,
      lwd = rle.med.point.lwd
    )
    abline(h = 0,
           col = zero.line.color,
           lwd = 3,
           lty = 2)
    axis(
      2,
      labels = T,
      mgp = c(3.5, .4, 0),
      lwd.ticks = 1,
      las = 1,
      cex.axis = 1
    )
    box(lwd = plot.lwd, bty = 'l')
    mtext('RLE', 2 , line = 2 , cex = 1)
    mtext('Samples', 1 , line = 2 , cex = 1)
    
  }
  if(show.tech.rep){
    boxplot(
      rle.data,
      xaxt = 'n',
      yaxt = 'n' ,
      ylab = '',
      xlab = '',
      outline = FALSE,
      frame = FALSE,
      whisklty = 0,
      staplelty = 0,
      names = FALSE,
      ...
    )
    points(
      c(1:ncol(rle.data)),
      rle.med,
      bg = rle.med.color[factor(variable)],
      col = 'black',
      pch = 21,
      cex = rle.med.point.size,
      lwd = rle.med.point.lwd
    )
    abline(h = 0,
           col = zero.line.color,
           lwd = 3,
           lty = 2)
    axis(
      2,
      mgp = c(3.5, .9, 0),
      lwd.ticks = 3,
      las = 1
    )
    box(lwd = plot.lwd, bty = 'l')
    mtext('RLE', 2 , line = 1 , cex = 1.5)
    mtext('Samples', 1 , line = 1 , cex = 1.5)
    tech.rep <- names(which(table(colnames(rle.data)) > 1))
    for(i in 1:length(tech.rep)){
      if(i %% 2 ==1){
        x.corr <- which(colnames(rle.data) == tech.rep[i] )
        y.corr <- rle.med[x.corr]
        for(j in 2:length(x.corr)){
          diagram::curvedarrow(
            from = c(x.corr[1], y.corr[1]), 
            to = c(x.corr[j], y.corr[j]),
            curve = -0.004, 
            arr.pos = 0.5, 
            lcol = curve.color,
            lty = curve.lty,
            lwd = curve.lwd,
            endhead = F,
            arr.length = 0,
            arr.adj = 1,
            arr.type = "triangle"
          )
        }
      }else{
        x.corr <- which(colnames(rle.data) == tech.rep[i] )
        y.corr <- rle.med[x.corr]
        for(j in 2:length(x.corr)){
          diagram::curvedarrow(
            from = c(x.corr[1], y.corr[1]), 
            to = c(x.corr[j], y.corr[j]),
            curve = 0.004, 
            arr.pos = 0.5, 
            lcol = curve.color,
            lty = curve.lty,
            lwd = curve.lwd,
            endhead = F,
            arr.length = 0,
            arr.adj = 1,
            arr.type = "triangle"
          )
        }
      }
      
    }
  }
}




#### color
dataSets.colors <- wesanderson::wes_palette(
  n = 4, 
  name = "GrandBudapest1")[c(1,2,4,3)]
names(dataSets.colors) <- c(
  'FPKM',
  'Quantile',
  'UQ',
  'RUV-III'
)


### pam50
### PAM50 colors
pam50.colors <- c(
  'gray',
  'red',
  'darkgreen',
  'navy',
  'cyan3', 
  'darkorange'
) 
names(pam50.colors) <- c(
  "Adjacent normal", 
  "Basal",
  "Her2",
  "LumA",
  "LumB", 
  "Normal like"
)


studies.colors <- c(
  'orange',
  'royalblue2',
  'plum2') 
names(studies.colors) <- c(
  'TCGA', 
  'GSE96058', 
  'GSE81538')

gg.theme <- theme(
  panel.background = element_blank(),
  axis.line = element_line(colour = 'black', size = .8),
  axis.title.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  axis.text.x = element_text(size = 8, angle = 25, hjust = 1),
  axis.text.y = element_text(size = 8),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 14),
  strip.text.x = element_text(size = 10)
)

