# Required R libraries and helper functions for the R shiny
# The helper functions are a subset of the Libraries_HelperFunctions.Rmd file

# R libraries
library(shiny)
library(SummarizedExperiment)
library(tidyverse)
library(BiocParallel)
library(BiocSingular)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(gridExtra)


# Principal component analysis using singular value decomposition
## data is gene expression matrix genes by samples
## is.log is logical factor indicating whether the data is log transformed or not
.pca <- function(data, nPCs = 10, is.log) {
  if (is.log) data <- data
  else data <- log2(data + 1)
  svd <- BiocSingular::runSVD(
    x = t(data),
    k = nPCs,
    BSPARAM = BiocSingular::bsparam(),
    center = TRUE,
    scale = FALSE
  )
  percent <- svd$d ^ 2 / sum(svd$d ^ 2) * 100
  percent <- sapply(seq_along(percent),
                    function(i) {
                      round(percent[i], 1)
                    })
  return(list(sing.val = svd,
              variation = percent))
}

# Spearman\Pearson correlation between individual genes and variables
## expr.data is gene expression matrix genes by samples
## is.log is logical factor indicating whether the data is log transformed or not
## variable is continuous variable such as library size 
## method whether Spearman or Pearson
## n.cores is number of cpus used for mclapply parallelization
## group a name for the variable such as library size or purity
.corr.gene.variable <- function(
  expr.data,
  is.log,
  variable,
  method,
  n.cores)
{
  if (is.log) expr.data <- expr.data
  else expr.data <- log2(expr.data + 1)
  rho <- parallel::mclapply(
    1:nrow(expr.data),
    function(x) {
      round(cor.test(
        x = expr.data[x,],
        y = variable,
        method = method)[[4]], 3)
    },
    mc.cores = n.cores)
  pval <- parallel::mclapply(
    1:nrow(expr.data),
    function(x) {
      cor.test(x = expr.data[x, ],
               y = variable,
               method = method)[[3]]
    },
    mc.cores = n.cores)
  
  results <- data.frame(
    genes = row.names(expr.data),
    rho = unlist(rho),
    pvalue = unlist(pval),
    adj.pvalue = p.adjust(unlist(pval), 'BH')
  )
  return(results)
}



# Analysis of variance between individual genes and variables 
## expr.data is gene expression matrix genes by samples
## is.log is logical factor indicating whether the data is log transformed or not
## n.cores is number of cpus used for mclapply parallelization
.Ftest <- function(data,
                   variable,
                   is.log,
                   n.cores)
{
  average.exp <- log2(rowMeans(data))
  if (is.log) data <- data
  else data <- log2(data + 1)
  f.test <- parallel::mclapply(
    1:nrow(data),
    function(x) {
      MASS::dropterm(lm(data[x , ] ~ variable), test = 'F')[c(5:6)]
    }
    , mc.cores = n.cores)
  f.test <- data.frame(
    Genes = row.names(data),
    FValue = round(unlist(lapply(f.test, function(x)
      x$`F Value`[2])), digits = 4) ,
    PValue = unlist(lapply(f.test, function(x)
      x$`Pr(F)`[2])),
    Adj.PValue = p.adjust(unlist(lapply(f.test, function(x)
      x$`Pr(F)`[2])), method = 'BH'),
    Mean = round(average.exp, digits = 2)
  )
  return(f.test)
}


# Scatter plot with  with density
### Function - here is my function to visualize PC
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
  pair.pcs <- utils::combn(colnames(pcs), 2)
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

# Silhouette coefficient analysis
## pcs is the matrix of principal components
## variable is a categorical variables such sample types
## nPcs is the number of principal components used to measure the distance
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


# Wilcoxon Rank Sum and Signed Rank Tests
## expr.data is gene expression matrix genes by samples
## is.log is logical factor indicating whether the data is log transformed or not
## variable is a categorical variables such sample types or batches
## n.cores is number of cpus used for mclapply parallelization
.wilcoxon.test <- function(
  expr.data, 
  is.log, 
  variable, 
  n.cores
){
  if(is.log){
    expr.data <- expr.data
  }
  else{
    expr.data <- log2(expr.data + 1)
  }
  pval <- parallel::mclapply(
    row.names(expr.data), 
    function(x) stats::wilcox.test(expr.data[x ,] ~ variable)[[3]], mc.cores = n.cores)
  results <- data.frame(
    genes = row.names(expr.data),
    pvalue = unlist(pval),
    ad.pvalue = p.adjust(p = unlist(pval), method = 'BH')
  )
  return(results)
}


# Pseudo-replicate of pseudo samples
## expr.data is gene expression matrix genes by samples
## sample.info is sample annotation file
## librarySize is a column in the sample.info the contains library size
## biology is columns in the sample.info the contains biological labels for samples
## batch is columns in the sample.info the contains batch information
## purity is columns in the sample.info the contains purity estimates
## include.ls indicates to consider PRPS for library size per batch 
## include.purity indicates to consider PRPS for purity per biology
## minSamplesPerBatchPS is the minimum number of samples per batch to creast PS
## minSamplesForPurityPerBiology is the minimum number of samples per biology to creast PS for purity
## minSamplesForPurityPS is the minimum number of samples to create PS for purity
## minSamplesForLibrarySizePerBatch is the minimum number of samples per batch to creast PS for library size
## minSamplesForLibrarySizePS is the minimum number of samples to create PS for library 

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
    message('error: there are not enough samples to create pseudo-samples for batch effects removal, you may want to lower minSamplesPerBatchPS')
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
          high.ls <- rowMeans(ls.data[ , c(ncol(ls.data)-(minSamplesForLibrarySizePS - 1)):ncol(ls.data) ])
          all <- cbind(low.ls, high.ls)
          colnames(all) <- rep(paste(x, 'LS', sep = '-'), 2)
          all
        })
      ps.ls <- do.call(cbind, ps.ls)
      
    }else{
      message('error: there are not enough samples to create pseudo-samples for removal of library size effects, you may want to lower minSamplesForLibrarySizePerBatch')
    }
  }else if (! include.ls){
    print('PRPS is not generated for librray size effects')
    ps.ls = list()
  }
  if(include.purity ){
    selected.biology.purity <- names(
      which(table(sample.info$biology) >= minSamplesForPurityPerBiology)
    ) 
    if(length(selected.biology.purity) > 0){
      message('PRPS are generated for purity effects')
      sample.info <- sample.info[
        with(sample.info,
             order(sample.info[, 'biology.batch'],
                   sample.info[, purity])),]
      expr.data <- expr.data[, row.names(sample.info)]
      ps.purity <- lapply(
        selected.biology.purity,
        function(x) {
          index <- sample.info$biology == x
          purity.data <- expr.data[, index]
          low.pur <- rowMeans(purity.data[, 1:minSamplesForPurityPS])
          high.pur <- rowMeans(purity.data[, c(ncol(purity.data) - (minSamplesForPurityPS - 1)):ncol(purity.data)])
          all <- cbind(low.pur, high.pur)
          colnames(all) <- rep(paste(x, 'purity', sep = '-'), 2)
          all
        })
      ps.purity <- do.call(cbind, ps.purity)
    }else{
      message('error: there are not enough samples to make pseudo-samples for purity variation, you may want to lower minSamplesForPurityPerBiology')
    }
  } else if (!include.purity){
    print('PRPS is not generated for purity effects')
    ps.purity = list()
  }
  return(list(ps.batch = ps.batch, ps.ls = ps.ls, ps.purity = ps.purity))
}


### Function
.fastRUVIII <- function(
  Y,
  M, ctl, 
  k = NULL, 
  eta = NULL,
  svd_k = 50, 
  include.intercept = TRUE, 
  average = FALSE,
  BPPARAM = SerialParam(), 
  BSPARAM = ExactParam(),
  fullalpha = NULL, 
  return.info = FALSE, 
  inputcheck = TRUE) 
{
  m <- nrow(Y)
  n <- ncol(Y)
  M <- ruv::replicate.matrix(M)
  ctl <- tological(ctl, n)
  
  ## Check the inputs
  if (inputcheck) {
    if (sum(is.na(Y)) > 0) {
      stop("Y contains missing values.  This is not supported.")
    }
    if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 
        0) {
      stop("Y contains infinities.  This is not supported.")
    }
  }
  ## RUV1 is a reprocessing step for RUVIII
  Y <- ruv::RUV1(Y, eta, ctl, include.intercept = include.intercept)
  if (class(BSPARAM) != "ExactParam") {
    svd_k <- min(m - ncol(M), sum(ctl), svd_k, na.rm = TRUE)
  } else {
    svd_k <- min(m - ncol(M), sum(ctl), na.rm = TRUE)
  }
  
  ## m represent the number of samples/observations ncol(M)
  ## represent the number of replicates If the replicate matrix
  ## is such that we have more replicates than samples, then
  ## RUV3 is not appropriate, thus, we return the Original input
  ## matrix
  if (ncol(M) >= m | k == 0) {
    newY <- Y
    fullalpha <- NULL
  } else {
    if (is.null(fullalpha)) 
    {
      ## The main RUVIII process Applies the residual operator of a
      ## matrix M to a matrix Y Y0 has the same dimensions as Y,
      ## i.e. m rows (observations) and n columns (genes).
      Y0 <- my_residop(Y, M)
      svdObj <- BiocSingular::runSVD(
        x = Y0, k = svd_k, BPPARAM = BPPARAM, BSPARAM = BSPARAM)
      
      fullalpha <- DelayedArray::t(svdObj$u[, seq_len(svd_k), drop = FALSE]) %*% Y
    }  ## End is.null(fullalpha)
    alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), , drop = FALSE]
    ac <- alpha[, ctl, drop = FALSE]
    W <- Y[, ctl] %*% DelayedArray::t(ac) %*% solve(ac %*% DelayedArray::t(ac))
    newY <- Y - W %*% alpha
  }  ## End else(ncol(M) >= m | k == 0)
  
  ## If the users want to get all the informations relating to
  ## the RUV, it can be done here.
  if (!return.info) {
    return(newY)
  } else {
    return(list(newY = newY, M = M, fullalpha = fullalpha, W = W))
  }
}
tological <- function(ctl, n) {
  ctl2 <- rep(FALSE, n)
  ctl2[ctl] <- TRUE
  return(ctl2)
}
my_residop <- function(A, B){
  tBB = DelayedArray::t(B) %*% B
  tBB_inv = Matrix::solve(tBB)
  BtBB_inv = B %*% tBB_inv
  tBA = DelayedArray::t(B) %*% A
  BtBB_inv_tBA = BtBB_inv %*% tBA
  return(A - BtBB_inv_tBA)
}
