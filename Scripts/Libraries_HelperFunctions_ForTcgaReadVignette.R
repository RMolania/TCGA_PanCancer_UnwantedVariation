# Required R libraries and helper functions for the vignette
# The helper functions are a subset of the Libraries_HelperFunctions.Rmd file

#=================== R libraries =================
library(ruv)
library(ggplot2)
library(cowplot)
library(mclust)
library(survival)
library(cowplot)
library(dplyr)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(ggfortify)
library(RColorBrewer)

#=================== PCA =================
# Principal component analysis using singular value decomposition (SVD)
## data: gene expression matrix genes by samples
## is.log:  logical factor indicating whether the data is log transformed or not
.pca <- function(data, is.log) {
  if (is.log)
    data <- data
  else
    data <- log2(data + 1)
  
  singValuDecomp <-
    base::svd(scale(
      x = t(data),
      center = TRUE,
      scale = FALSE
    ))
  percent <- singValuDecomp$d ^ 2 / sum(singValuDecomp$d ^ 2) * 100
  percent.var <-
    sapply(
      seq_along(percent),
      function(i) {
        round(percent[i], 1)
      })
  return(list(
    sing.val = singValuDecomp,
    variation = percent.var))
}
#=================== Correlation =================
# Spearman\Pearson correlation between individual genes and variables
## expr.data is gene expression matrix genes by samples
## is.log is logical factor indicating whether the data is log transformed or not
## variable is continuous variable such as library size 
## method whether Spearman or Pearson
## n.cores is number of cpus  for mclapply parallelization
## group a name for the variable such as library size or purity
.corr.gene.variable <- function(
  expr.data,
  is.log,
  variable,
  method,
  n.cores,
  group)
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
    pvalue = unlist(pval)
  )
  colnames(results) <- paste(
    group,
    colnames(results),
    sep = '_')
  return(results)
}

#=================== ANOVA =================
# Analysis of variance between individual genes and variables 
## expr.data is gene expression matrix genes by samples
## is.log is logical factor indicating whether the data is log transformed or not
## n.cores is number of cpus used for mclapply parallelization
.Ftest <- function(
  data,
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

#=================== Silhouette coefficient =================
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

#=================== Wilcoxon rank =================
# Wilcoxon Rank Sum and Signed Rank Test
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

#=================== PRPS =================
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
      message('error: there are not enough samples to create pseudo-samples for removal of library size effects, 
              you may want to lower minSamplesForLibrarySizePerBatch')
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
      message('error: there are not enough samples to make pseudo-samples for purity variation, 
              you may want to lower minSamplesForPurityPerBiology')
    }
  } else if (!include.purity){
    print('PRPS is not generated for purity effects')
    ps.purity = list()
  }
  return(list(ps.batch = ps.batch, ps.ls = ps.ls, ps.purity = ps.purity))
}

#=================== RUV-III normalization  =================
# Removing Unwanted Variation - III (RUV-III)
## Y raw gene expression matrix (samples by genes)
## M replicate matrix
## ctl negative control genes
## k the number of unwanted factors to use. 

RUV_III_PRPS <- function(
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
      Y0 <- ruv::residop(Y, M)
      fullalpha <- t(svd(Y0 %*% t(Y0))$u[, 1:min(m - ncol(M), sum(ctl)), drop = FALSE]) %*% Y
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


#=================== Survival analysis  =================
# This function depends on the below libraries 
# This function looks at the associations between different variables and survival outcome. These variables can be one of the below options that stratifies samples for survival analysis:

# expr : expression of a gene, will be split based on 33%-tile and 66%-tile (e.g. low, medium, high)\n
# score : score of a single signatre, will be split based on 33%-tile and 66%-tile (low, medium, high)\n
# covariate : A continouse covariate (e.g. age), will be split based on 33%-tile and 66%-tile (low, medium, high)\n
# score_expr : stratifies samples based on scores from a signature (high and low) and expression of a gene (high and low)\n
# covariate_expr : startifies samples according to covariate (age; high and low) and expression of a gene (high and low)\n
# score_covariate: stratifies samples according to scores from a single signature (high and low) and covariate (age; high and low)\n
# expr_expr : stratifies samples according to expression of two genes (high gene1/high gene2, high gene1/low gene2, etc)\n
# score_score : stratifies samples according to scores obtained from two signatures (high score1/high score2, high score1/low score2, etc)

survival_plot <- function(
  data = exprData,
  stratify = "score_score",
  annot = newAnnot,
  scoreCol =  c("TRM TGFb IL2 Sel Com", "Mes (Byers)"),
  gene = c("ITGAE", "ZNF683"),
  covariate = "age_at_initial_pathologic_diagnosis",
  isCategoricalCov = FALSE,
  timeCol = "OS.time",
  eventCol = "OS",
  nGroup = 2,
  confInt = F,
  mainTitle1,
  ylabel = "Survival",
  cols = c(brewer.pal(9, "Set1")[c(2, 3, 4, 5, 7, 8)],
           brewer.pal(8, "Dark2")[c(8, 1, 4, 6)]),
  nColLegend = 1,
  plotType = "autoplot") {
  annot[, timeCol] <- gsub("#N/A", NA, annot[, timeCol])
  comSamples <- intersect(rownames(annot), colnames(data))
  data <- data[, comSamples]
  annot <- annot[comSamples, ]

  data <- data[, complete.cases(annot[, timeCol])]
  annot <- annot[complete.cases(annot[, timeCol]),]
  
  annot[, timeCol] <- as.numeric(annot[, timeCol] )
  
  ##---------------------------------------- Scores
  if (!is.null(scoreCol)) {
    
    currentscoreCol <- scoreCol[1]
    median_score <- median(annot[, currentscoreCol])
    lowQ_score <- as.numeric(quantile(annot[, currentscoreCol], prob = 0.33))
    upQ_score <- as.numeric(quantile(annot[, currentscoreCol], prob = 0.66))
    
    # annot$scores_2status <-
    #   ifelse(
    #     annot[, currentscoreCol] >= median_score,
    #     paste("High ", currentscoreCol),
    #     paste("Low ", currentscoreCol)
    #   )
    
    annot$scores_2status[annot[, currentscoreCol] >= median_score] <- paste("High ", currentscoreCol)
    annot$scores_2status[annot[, currentscoreCol] < median_score] <- paste("Low ", currentscoreCol)
    
    annot$scores_3status[annot[, currentscoreCol] >= upQ_score] <-
      paste("High ", currentscoreCol)
    annot$scores_3status[annot[, currentscoreCol] <= lowQ_score] <-
      paste("Low ", currentscoreCol)
    annot$scores_3status[annot[, currentscoreCol] < upQ_score &
                           annot[, currentscoreCol] > lowQ_score] <-
      paste("Medium ", currentscoreCol)
    
    
    if (length(scoreCol) == 2) {
      currentscoreCol <- scoreCol[2]
      median_score <- median(annot[, currentscoreCol])
      lowQ_score <-
        quantile(annot[, currentscoreCol], prob = 0.33)
      upQ_score <-
        quantile(annot[, currentscoreCol], prob = 0.66)
      
      annot$scores2_2status <-
        ifelse(
          annot[, currentscoreCol] >= median_score,
          paste("High ", currentscoreCol),
          paste("Low ", currentscoreCol)
        )
      
      annot$scores2_3status[annot[, currentscoreCol] >= upQ_score] <-
        paste("High ", currentscoreCol)
      annot$scores2_3status[annot[, currentscoreCol] <= lowQ_score] <-
        paste("Low ", currentscoreCol)
      annot$scores2_3status[annot[, currentscoreCol] < upQ_score &
                              annot[, currentscoreCol] > lowQ_score] <-
        paste("Medium ", currentscoreCol)
    }
    if (length(scoreCol) > 2) {
      stop(paste0("You must specify maximum of 2 score columns at a time"))
    }
    ## save this new annotation data as sample annotation for the data
    # colData(data) <- newAnnot
  }
  
  ##-------------------------------------- Covariate
  
  if (!is.null(covariate)) {
    annot <- annot[complete.cases(annot[, covariate]), ]
    badcols <-
      c(
        "not reported",
        "NA",
        "Indeterminate",
        "[Not Applicable]",
        "[Not Available]",
        "[Discrepancy]",
        "[Unknown]",
        "Not Evaluable"
      )
    annot <- annot[ ! annot[, covariate] %in% badcols, ]
    
    comSamples <- intersect(row.names(annot), colnames(data))
    annot <- annot[comSamples, ]
    data <- data[, comSamples]
    
    # newAnnot <- colData(data)
    if(isCategoricalCov){
      annot[, covariate] <- as.factor(annot[, covariate])
    }
    else if(! isCategoricalCov) {
      annot[, covariate] <- as.numeric(annot[, covariate])
      median_cov <- median(annot[, covariate], na.rm = T)
      lowQ_cov <-
        as.numeric(quantile(annot[, covariate], prob = 0.33, na.rm = T))
      upQ_cov <-
        as.numeric(quantile(annot[, covariate], prob = 0.66, na.rm = T))
      
      annot$cov_2status <-
        ifelse(annot[, covariate] >= median_cov,
               "High covariate",
               "Low covariate")
      
      annot$cov_3status[annot[, covariate] >= upQ_cov] <-
        "High covariate"
      annot$cov_3status[annot[, covariate] <= lowQ_cov] <-
        "Low covariate"
      annot$cov_3status[annot[, covariate] < upQ_cov &
                          annot[, covariate] > lowQ_cov] <-
        "Medium covariate"
      
    }
  }
  
  currentData <- data
  
  ##------------------------------------- Gene expression
  if (!is.null(gene)) {
    if (sum(rownames(data) %in% gene) < 1) {
      stop(paste0(gene, " does not present in the row names of the expression data"))
    }
    if (length(gene) > 3) {
      stop(paste0("Please provide maximum of 3 genes at a time"))
    }
    currentGene <- gene[1]
    currentGeneIndx <- which(rownames(currentData) == currentGene)
    # newAnnot <- colData(currentData)
    
    ## calculate median and 33%-tile and 66%-tile of gene expression
    median_expr <- median(as.numeric(currentData[ currentGeneIndx, ]))
    lowQ_expr <- as.numeric(quantile(currentData[ currentGeneIndx, ], prob = 0.33))
    upQ_expr <- as.numeric(quantile(currentData[ currentGeneIndx, ], prob = 0.66))
    
    # annot$expr1_2status <-
    #   ifelse(
    #     currentData[ currentGeneIndx, ] >= median_expr,
    #     paste0("High ", currentGene),
    #     paste0("Low ", currentGene)
    #   )
    
    annot$expr1_2status[currentData[ currentGeneIndx, ] >= median_expr] <- paste0("High ", currentGene)
    annot$expr1_2status[currentData[ currentGeneIndx, ] < median_expr] <- paste0("Low ", currentGene)
    
    # annot$expr1_3status[as.numeric(currentData[ currentGeneIndx, ]) >= upQ_expr] <-
    #   paste0("High ", currentGene)
    
    annot$expr1_3status[currentData[ currentGeneIndx, ] >= upQ_expr] <-
      paste0("High ", currentGene)
    
    annot$expr1_3status[currentData[ currentGeneIndx, ] <= lowQ_expr] <-
      paste0("Low ", currentGene)
    annot$expr1_3status[currentData[ currentGeneIndx, ] < upQ_expr &
                          currentData[ currentGeneIndx, ] > lowQ_expr] <-
      paste0("Medium ", currentGene)
    
    ## save this new annotation as sample annotation for the data
    # colData(currentData) <- newAnnot
    
    if (length(gene) > 1) {
      currentGene <- gene[2]
      # newAnnot <- colData(currentData)
      currentGeneIndx <- which(rownames(currentData) == currentGene)
      
      ## calculate median and 33%-tile and 66%-tile of gene expression
      median_expr <-
        median(currentData[currentGeneIndx, ])
      lowQ_expr <-
        as.numeric(quantile(currentData[currentGeneIndx, ], prob = 0.33))
      upQ_expr <-
        as.numeric(quantile(currentData[currentGeneIndx, ], prob = 0.66))
      
      annot$expr2_2status <-
        ifelse(
          currentData[currentGeneIndx, ] >= median_expr,
          paste0("High ", currentGene),
          paste0("Low ", currentGene)
        )
      
      annot$expr2_3status[currentData[currentGeneIndx, ] >= upQ_expr] <-
        paste0("High ", currentGene)
      annot$expr2_3status[currentData[currentGeneIndx, ] <= lowQ_expr] <-
        paste0("Low ", currentGene)
      annot$expr2_3status[currentData[currentGeneIndx, ] < upQ_expr &
                            currentData[currentGeneIndx, ] > lowQ_expr] <-
        paste0("Medium ", currentGene)
      
      if (length(gene) == 3) {
        currentGene <- gene[3]
        currentGeneIndx <- which(rownames(currentData) == currentGene)
        
        ## calculate median and 33%-tile and 66%-tile of gene expression
        median_expr <-
          median(currentData[currentGeneIndx, ])
        lowQ_expr <-
          as.numeric(quantile(currentData[currentGeneIndx, ], prob = 0.33))
        upQ_expr <-
          as.numeric(quantile(currentData[currentGeneIndx, ], prob = 0.66))
        
        annot$expr3_2status <-
          ifelse(
            currentData[currentGeneIndx, ] >= median_expr,
            paste0("High ", currentGene),
            paste0("Low ", currentGene)
          )
        
        annot$expr3_3status[currentData[currentGeneIndx, ] >= upQ_expr] <-
          paste0("High ", currentGene)
        annot$expr3_3status[currentData[currentGeneIndx, ] <= lowQ_expr] <-
          paste0("Low ", currentGene)
        annot$expr3_3status[currentData[currentGeneIndx, ] < upQ_expr &
                              currentData[currentGeneIndx, ] > lowQ_expr] <-
          paste0("Medium ", currentGene)
        
      }
    }
    
    # colData(currentData) <- newAnnot
  }
  
  # currentData <-
  #   currentData[, complete.cases(annot[, timeCol])]
  
  ##------------------------------------- Check for stratification type
  if (stratify == "expr") {
    currentStrata <- paste0("expr1_", nGroup, "status")
    mainTitle <- gene[1]
  }
  if(stratify == "score"){
    currentStrata <- paste0("scores_", nGroup, "status")
    mainTitle <- scoreCol[1]
  }
  if(stratify == "covariate"){
    if(isCategoricalCov){
      currentStrata <- covariate
    }
    else if(!isCategoricalCov){
      currentStrata <- paste0("cov_", nGroup, "status")
    }
    mainTitle <- covariate
  }
  if (stratify == "score_expr") {
    currentSt_score <- paste0("scores_", nGroup, "status")
    currentSt_expr  <- paste0("expr1_", nGroup, "status")
    annot$score_expr <-
      paste0(annot[, currentSt_score],
             " / ",
             annot[, currentSt_expr])
    currentStrata <- "score_expr"
    mainTitle <- paste(scoreCol[1], " &\n", gene[1])
  }
  if (stratify == "covariate_expr") {
    if(is.null(covariate) | is.null(gene)){
      stop("Make sure both covriate and gene are provided")
    }
    
    if(isCategoricalCov){
      currentSt_cov <- covariate
    } else if(!isCategoricalCov){
      currentSt_cov   <- paste0("cov_", nGroup, "status")
    }
    
    #currentSt_cov  <- paste0("cov_", nGroup, "status")
    currentSt_expr <- paste0("expr1_", nGroup, "status")
    ## Remove samples with NA annotation as covariate:
    currentData <- currentData[ , complete.cases(annot[, currentSt_cov])]
    annot$cov_expr <-
      paste0(annot[, currentSt_cov],
             " / ",
             annot[, currentSt_expr])
    currentStrata <- "cov_expr"
    mainTitle <- paste(covariate, " &\n", gene[1])
  }
  if (stratify == "score_covariate") {
    
    if(isCategoricalCov){
      currentSt_cov <- covariate
    } else if(!isCategoricalCov){
      currentSt_cov   <- paste0("cov_", nGroup, "status")
    }
    
    currentSt_score <- paste0("scores_", nGroup, "status")
    
    annot$score_cov <-
      paste0(annot[, currentSt_score],
             " / ",
             annot[, currentSt_cov])
    currentStrata <- "score_cov"
    mainTitle <- paste(scoreCol[1], " &\n", covariate)
  }
  if (stratify == "expr_expr") {
    currentSt_expr1 <- paste0("expr1_", nGroup, "status")
    currentSt_expr2 <- paste0("expr2_", nGroup, "status")
    annot$expr_expr <-
      paste0(annot[, currentSt_expr1],
             " / ",
             annot[, currentSt_expr2])
    currentStrata <- "expr_expr"
    mainTitle <- paste(gene[1], " &\n", gene[2])
  }
  if (stratify == "score_score") {
    currentSt_score1 <- paste0("scores_", nGroup, "status")
    currentSt_score2 <- paste0("scores2_", nGroup, "status")
    annot$score_score <-
      paste0(
        annot[, currentSt_score1],
        " / ",
        annot[, currentSt_score2])
    currentStrata <- "score_score"
    mainTitle <- paste(scoreCol[1], " &\n", scoreCol[2])
  }
  
  if (stratify == "expr_expr_expr") {
    currentSt_expr1 <- paste0("expr1_", nGroup, "status")
    currentSt_expr2 <- paste0("expr2_", nGroup, "status")
    currentSt_expr3 <- paste0("expr3_", nGroup, "status")
    annot$expr_expr_expr <-
      paste0(annot[, currentSt_expr1],
             " / ",
             annot[, currentSt_expr2],
             " / ",
             annot[, currentSt_expr3])
    currentStrata <- "expr_expr_expr"
    mainTitle <- paste(gene[1], " &\n", gene[2], " & ", gene[3])
  }
  ##----------------------------- Fit survival model
  # annot2 <- annot
  # annot2[, timeCol] <- gsub("#N/A", NA, annot2[, timeCol])
  # annot2 <- annot2[complete.cases(annot2[, timeCol]), ]
  tt <- data.frame(table(annot[, currentStrata]))
  tt$Var1 <- as.character(tt$Var1)
  tt$Freq <- as.character(tt$Freq)
  for (i in 1:nrow(tt)) {
    annot$currentStrata_n[annot[, currentStrata] == tt$Var1[i]] <-
      paste0(tt$Var1[i], " (", tt$Freq[i], ")")
  }
  
  annot <- annot[, c("currentStrata_n", timeCol, eventCol)]
  annot[, timeCol] <- as.numeric(annot[, timeCol])
  annot[, eventCol] <- as.numeric(annot[, eventCol])
  
  fitValues <- survfit(Surv(time = annot[, timeCol],
                            event = annot[, eventCol]) ~
                         annot$currentStrata_n)
  
  ss <- survdiff(Surv(time =  annot[, timeCol],
                      event = annot[, eventCol]
  ) ~
    annot$currentStrata_n)
  # 
  #   ss <- survdiff(Surv(annot[, timeCol],
  #                       as.numeric(as.factor(
  #                         annot[, eventCol]
  #                       )) - 1) ~
  #                    annot$currentStrata_n)
  
  ##------------------------------- Calculate p-value
  ## This does not adjust for any covariates, unless the covariate option is included
  
  pval <-  ifelse (is.na(ss), next, (round(1 - pchisq(
    ss$chisq, length(ss$n) - 1
  ), 6)))[[1]]
  
  if(pval < 0.01){
    pval_add <- paste0("p-value (log-rank) < 0.01")
  }
  else{
    pval_add <- paste0("p-value (log-rank) = ", round(pval, 2))
  }
  
  ##------------------------------ Plot survival curve
  if(plotType == "autoplot"){
    p <-  autoplot(fitValues, surv.size = 1.5, conf.int = confInt) +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      ggtitle(paste0(
        mainTitle1, '\n' , 
        # " (Chisq = ", round(ss$chisq, 3),
        #"\n", 
        pval_add)) +
      ylab(ylabel) +
      xlab("Time (days)") +
      theme(
        panel.background = element_blank(),
        legend.position = 'bottom'
        
      ) +
      guides(
        color = guide_legend(ncol = nColLegend), 
        fill = guide_legend(ncol = nColLegend))
    
  }
  else if (plotType == "ggsurvplot"){
    p <- ggsurvplot (
      fitValues,
      data = annot,
      fun = "pct",
      pval = TRUE,
      # pval.method = TRUE,  ## Log Rank
      # test.for.trend = T,  ## when we have more than two groups
      conf.int = confInt,
      surv.median.line = "hv",
      # linetype = "strata",
      palette = cols,
      xlab = "Time",
      legend.title = mainTitle,
      # legend.labs = c("High score", "Low score"),
      legend = c(.2, .2),
      # break.time.by = 4,
      # risk.table = TRUE,
      # tables.height = 0.2,
      # tables.theme = theme_cleantable(),
      # risk.table.y.text.col = TRUE,
      # risk.table.y.text = TRUE
    )
  }
  p_pval <- list(plot = p, pval = pval)
  return(p_pval)
  
}


#=================== PCA plot with density  =================
## PCS: pca from 
## pc.var: pca variation output from .pca() function 
## group.legend name
## group: logical factor to color the PCA plot
## color: 
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

#=================== Visualization =================
# ggplot theme
ggplot.them <- theme(
  panel.background = element_blank(),
  axis.line = element_line(colour = 'black', size = 1),
  axis.title.x = element_text(size = 10),
  axis.title.y = element_text(size = 10),
  plot.title = element_text(size = 10),
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 8),
  strip.text.x = element_text(size = 8),
  legend.position = 'bottom'
)
gg.theme.2 <- theme(
  panel.background = element_blank(),
  axis.line = element_line(colour = 'black', size = 1),
  axis.title.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  axis.text.x = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  plot.title = element_text(size = 14),
  plot.margin = unit(c(.2,.2,.2,.2), "cm"),
  legend.position = 'none'
)

# Upper function for ggpair
upperfun <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_point(size = .8, alpha = 0.3) +
    geom_density2d() +
    geom_vline(xintercept = .0, color = 'red')+
    geom_hline(yintercept = .0, color = 'red')+
    scale_y_continuous(limits = c(-1, 1)) +
    scale_x_continuous(limits = c(-1, 1)) +
    geom_abline(col = 'cyan') +
    theme(
      axis.line = element_line(colour = 'black', size = 1),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 55),
      axis.title.y = element_text(size = 55),
    )}

# GeomSplitViolin 
GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <-
      transform(
        data,
        xminv = x - violinwidth * (x - xmin),
        xmaxv = x + violinwidth * (xmax - x)
      )
    grp <- data[1, "group"]
    newdata <-
      plyr::arrange(transform(data, x = if (grp %% 2 == 1)
        xminv
        else
          xmaxv), if (grp %% 2 == 1)
            y
        else-y)
    newdata <-
      rbind(newdata[1,], newdata, newdata[nrow(newdata),], newdata[1,])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <-
      round(newdata[1, "x"])
    
    if (length(draw_quantiles) > 0 &
        !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                1))
      quantiles <-
        ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <-
        data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <-
        rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <-
        GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin",
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
      ggplot2:::ggname("geom_split_violin",
                       GeomPolygon$draw_panel(newdata, ...))
    }
  })

geom_split_violin <- function(
  mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
  draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
  show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm,
      ...
    )
  )
}

# Selected colors
## Major time intervals
major.times.colors <- c(
  'turquoise4', 
  'yellow4'
  )
names(major.times.colors) <- c(
  '2010', 
  '2011:2014'
  )

## CMS subtypes
cms.colors <- c(
  'orange',
  'royalblue2',
  'plum2',
  'seagreen3',
  'gray70'
) 
names(cms.colors) <- c(
  'CMS1', 
  'CMS2', 
  'CMS3', 
  'CMS4',
  'Not classified'
)



cms.colors.so <- c(
  'orange',
  'royalblue2',
  'plum2',
  'seagreen3',
  'gray70',
  '#D9D9D9')
names(cms.colors.so) <- c(
  'CMS1',
  'CMS2',
  'CMS3',
  'CMS4',
  'Not classified',
  'Adjacent normal')


msi.colors.so <- c(
  'purple3',
  'hotpink2',
  'yellow2',
  'gray70',
  'red')
names(msi.colors.so) <- c(
  "MSS" ,
  "MSI-L",
  "MSI-H",
  "Indeterminate",
  "Adjacent normal")


## Datasets
dataSets.colors <- wesanderson::wes_palette(
  n = 4, 
  name = "GrandBudapest1")[c(1,2,4,3)]
names(dataSets.colors) <- c(
  'Raw counts',
  'FPKM',
  'FPKM.UQ',
  'RUV-III')

dataSets.colors.2 <- c(
  dataSets.colors,
  'darkgreen', 
  'blue', 
  'yellow3', 
  'cyan'
)
names(dataSets.colors.2)[5:8] <- c(
  'RUV-III-0.2', 
  'RUV-III-0.4', 
  'RUV-III-0.6', 
  'RUV-III-0.8')

dataSets.colors3 <- c(
  dataSets.colors, 
  'cyan')
names(dataSets.colors3)[5] <- 'RUV-III-P'

## Time- years
years.colors <- c(
  'purple4',
  'blue',
  'tan1',
  'darkgreen')
names(years.colors) <- c(
  '2010',
  '2011',
  '2013',
  '2014')
