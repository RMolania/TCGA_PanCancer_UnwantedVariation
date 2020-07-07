### Library #####################################
## tidy family
library(dplyr)
library(plyr)
library(reshape2)
library(data.table)
library(janitor)
library(ggpubr)
# library(tidyverse)
# library(tidyr)

## Parallelization
library(parallel)
library(doParallel)
library(foreach)
library(future.apply)
library(BiocParallel)

## Normalization and batch effects
library(scater)
library(Seurat)
library(glmpca)
library(scMerge)
library(scran)
library(sctransform)
library(SingleR)

## Statistical 
library(scales)
library(Matrix)
library(matrixStats)
library(quantreg)
library(MASS)

## Visualization
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(VennDiagram)
library(ComplexHeatmap)
library(plotrix)
library(vioplot)
library(ggbeeswarm)
library(cowplot)
library(gridExtra)
library(gridExtra)
library(viridis)

### Clustering
library(SC3)

## gene annotations
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GSA)
library(GSEABase)


## Cell type annotations
library(AUCell)
library(scry)

## class objects
library(SingleCellExperiment)


## TCGA database
# library(TCGAbiolinks)
library(HDF5Array)
library(SummarizedExperiment)

## reading files
library(readxl)


## sample scoring
library(singscore)




### Reading ecell sheet  ############################################
read.excel.allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(
    sheets, 
    function(y){
      readxl::read_excel(filename, sheet = y)
    }) 
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

### Tidy batch information  #########################################
tidy.bacth.info <- function(path, df) {
  temp.df <- read.delim(paste0(path, df))
  temp.df <- temp.df %>%
    dplyr::mutate(
      sample.id = substr(Sample, 1, 16),
      year = substr(ShipDate, 1, 4),
      month = substr(ShipDate, 6, 7),
      day = substr(ShipDate, 9, 10),
      center = substr(Sample, 27, 28),
      protion = substr(Sample, 18, 19),
      protion = substr(Sample, 18, 19),
      sample.id.a = substr(Sample, 1, 12),
      sample.id.b = substr(Sample, 1, 15),
      sample.id.c = substr(Sample, 1, 16),
      sample.id.d = substr(Sample, 1, 20),
      vial = substr(Sample, 16, 16)
    )
  temp.df
  
}
### Decoding TCGA barcode ###########################################
decode.tcga.barcode <- function(code, sep){
  temp <- sapply(
    code, 
    function(y) strsplit(
      x = y, 
      split = sep))
  temp <- as.data.frame(t(
    data.frame(temp))
  )
  colnames(temp) <-
    c('project.code',
      'tss.code',
      'participant.code',
      'tissue.code',
      'vial.code',
      'plate.code',
      'center.code')
  row.names(temp) <- gsub(
    '\\.',
    '-',
    row.names(temp)
  )
  temp
}


### PCA #####################################################
.pca <- function(data, is.log) {
  if(is.log == TRUE){
    data <- data
  }else{
    data <- log2(data + 1)
  }
  svd <- svd(apply(data, 1, function(x) scale(x, scale = FALSE, center = TRUE)))
  percent <- svd$d^2/sum(svd$d^2)*100
  percent <- sapply(seq_along(percent), function(i) {round(percent[i], 2)})
  return(list(sing.val = svd, var = percent))
}


### Median and IQR of RLE ####################################
rle.com <- function(data, is.log){
  if(is.log == TRUE){
    rle.df <- data - rowMedians(data)
    rle.med <- apply(rle.df, 2, median)
    rle.iqr <- apply(rle.df, 2, iqr)
  }else{
    data <- log2(data + 1)
    rle.df <- data - rowMedians(data)
    rle.med <- apply(rle.df, 2, median)
    rle.iqr <- apply(rle.df, 2, iqr)
    
  }
  return(list(rle.med = rle.med, rle.iqr = rle.iqr))
}



### PCA regression ###########################################
correlate.fun_two <- function(rot.data, batch, batch.levels) {
  #rot.data: some vector (numeric entries)
  #batch: some vector (categoric entries)
  a <- lm(rot.data ~ batch)
  result <- numeric(2)
  result[1] <- summary(a)$r.squared #coefficient of determination
  result[2] <- summary(a)$coefficients[2,4] #p-value (significance level)
  t.test.result <- t.test(rot.data[batch == batch.levels[1]],
                          rot.data[batch == batch.levels[2]], paired = FALSE)
  result[3] <- t.test.result$p.value
  result
}


correlate.fun_gen <- function(rot.data, batch){
  #rot.data: some vector (numeric covariate)
  #batch: some vector (categoric covariate)
  a <- lm(rot.data ~ batch)
  result <- numeric(2)
  result[1] <- summary(a)$r.squared #coefficient of determination
  F.test.result <- aov(rot.data ~ batch)
  F.test.summary <- summary(F.test.result)
  
  result[2] <- summary(a)$coefficients[2,4] #p-value (significance level)
  result[3] <- F.test.summary[[1]]$'Pr(>F)'[1] #p-value of the one-way anova test
  
  result
}

pcRegression <- function(pca.data, batch,n_top=50, tol=1e-16){
  batch.levels <- unique(batch)
  #make sure you do not try to assess more PCs than actually computed
  pca_rank = ncol(pca.data$x)
  max_comps <- min(pca_rank, n_top)
  if (length(pca.data$sdev) > pca_rank) {
    pca.data$sdev = pca.data$sdev[1:pca_rank]
  }
  
  if (length(batch.levels) == 2) {
    #r2.batch.raw <- r2.batch
    
    # for-loop replaced by correlate.fun and apply
    r2.batch <- apply(pca.data$x, 2, correlate.fun_two, batch, batch.levels)
    r2.batch <- t(r2.batch)
    colnames(r2.batch) <- c('R.squared', 'p.value.lm', 'p.value.t.test')
    
    r2.batch[r2.batch[,2] < tol, 2] <- tol
    r2.batch[r2.batch[,3] < tol, 3] <- tol
  } else {
    
    
    #r2.batch.raw <- r2.batch
    
    r2.batch <- apply(pca.data$x, 2, correlate.fun_gen, batch)
    
    r2.batch <- t(r2.batch)
    colnames(r2.batch) <- c('R.squared', 'p.value.lm', 'p.value.F.test')
    # for-loop replaced by correlate.fun and apply
    #for (k in 1:dim(r2.batch)[1]){
    #    a <- lm(pca.data$x[,k] ~ batch)
    #    r2.batch[k,1] <- summary(a)$r.squared #coefficient of determination
    #    r2.batch[k,2] <- summary(a)$coefficients['batch',4] #p-value (significance level)
    #}
    
    r2.batch[r2.batch[,2] < tol, 2] <- tol
    r2.batch[r2.batch[,3] < tol, 3] <- tol
  }
  
  argmin <- which(r2.batch[, 2] == min(r2.batch[, 2]))
  normal <- sum(pca.data$sdev^2)
  var <- round((pca.data$sdev)^2 / normal * 100,1)
  batch.var <- sum(r2.batch[,1]*var)/100
  setsignif <- p.adjust(r2.batch[1:max_comps,2], method = 'BH') < 0.05
  pcCorr <- sqrt(r2.batch[1:max_comps, 1])
  result <- list()
  result$maxVar <- var[argmin]
  result$PmaxVar <- r2.batch[argmin,2]
  result$pcNfrac <- mean(setsignif)
  result$pcRegscale <- sum(var[1:max_comps][setsignif])/sum(var[1:max_comps])
  result$maxCorr <- max(pcCorr)
  result$maxR2 <- max(r2.batch[1:max_comps, 1])
  result$msigPC   <- 1 - (min(c(which(setsignif), max_comps + 1)) - 1) / max_comps
  result$maxsigPC <- 1 - (min(c(which(pcCorr == max(pcCorr[setsignif])), max_comps + 1)) - 1) / max_comps
  result$R2Var <- batch.var
  result$ExplainedVar <- var
  result$r2 <- r2.batch
  result
}



### Ftest ##############################################
.Ftest <- function(data, var, is.log, ncores){
  if(is.log == TRUE){
    average.exp <- rowMeans(data)
    f.test <- mclapply(
      1:nrow(data), 
      function(x){
        dropterm (
          lm(data[x ,] ~ var),
          test = 'F')[c(5:6)] }
      , mc.cores = ncores)
    f.test <- data.frame(
      'FValue' = round(unlist(lapply(f.test, function(x) x$`F Value`[2])), digits = 2) , 
      'PValue' = unlist(lapply(f.test, function(x) x$`Pr(F)`[2])),
      'Mean' = round(average.exp, digits = 2))
  }else{
    data <- log2(data+1)
    average.exp <- rowMeans(data)
    f.test <- mclapply(
      1:nrow(data), 
      function(x){
        dropterm (
          lm(data[x ,] ~ var),
          test = 'F')[c(5:6)] }
      , mc.cores = ncores)
    f.test <- data.frame(
      'FValue' = round(unlist(lapply(f.test, function(x) x$`F Value`[2])), digits = 2) , 
      'PValue' = unlist(lapply(f.test, function(x) x$`Pr(F)`[2])),
      'Mean' = round(average.exp, digits = 2))
  }
  return(f.test)
}

### genes vs variable correlation ##############################################
.genes.var <- function(data, var, is.log){
  if(is.log == TRUE){
    data <- data
  }else{
    data <- log2(data + 1)
  }
  cor.coeff <- lapply(
    1:nrow(data),
    function(x){
      cor.test(data[x,] , var, method = 'spearman')[[4]]
    })
  return(round(unlist(cor.coeff), digits = 2))
  }


### scatter with density plot _ PCA ###########################################
.scatter.density.pc <- function(pcs, pc.var, group.name, group, color){
  pair.pcs <- utils::combn(ncol(pcs), 2)
  pList <- list()
  for(i in 1:ncol(pair.pcs)){
    if(i == 1){
      x <- pair.pcs[1,i]
      y <- pair.pcs[2,i]
      p <- ggplot(mapping = aes(
        x = pcs[,x], 
        y = pcs[,y], 
        color = group)) +
        xlab(paste0('PC', x, ' (',pc.var[x], ')')) +
        ylab(paste0('PC', y, ' (',pc.var[y], ')')) +
        geom_point() +
        theme(
          legend.position="right",
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.background = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.key = element_blank()
        ) +
        scale_color_manual(name = group.name, values = color)
      le <- get_legend(p)
    }else{
      x <- pair.pcs[1,i]
      y <- pair.pcs[2,i]
      p <- ggplot(mapping = aes(
        x = pcs[,x], 
        y = pcs[,y], 
        color = group)) +
        xlab(paste0('PC', x, ' (',pc.var[x], ')')) +
        ylab(paste0('PC', y, ' (',pc.var[y], ')')) +
        geom_point() +
        theme(
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position = "none") +
        scale_color_manual(
          values = color, 
          name = group.name
        )
      
    }
    p <- p + theme(legend.position = "none")
    xdens <- axis_canvas(p, axis = "x")+
      geom_density(
        mapping = aes(
          x = pcs[,x], 
          fill = group),
        alpha = 0.7, 
        size = 0.2
      ) +
      theme(legend.position = "none") +
      scale_fill_manual(values = color)
    
    ydens <- axis_canvas(
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


