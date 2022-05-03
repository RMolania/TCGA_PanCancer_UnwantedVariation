# RUV-III generation

#' @title Generate RUV-III Data Object
#'
#' @description This function is a part of the data analysis functionality of `tcgaCleaneR`. It captures both the uses PRPS values from library \code{tcgaCleaneR} combined with row counts and run SVD algorithm (\code{runSVD()}) from \code{BiocSingular} on the combined dataset. The function uses RUV-I algorithm from \code{ruv} as a pre-processing step to RUV-III.
#'
#' @param ruv.data S4 data object for RUV-III: A S4 data object with combined data including the row count from original filtered data using assay \code{HTseq_counts}, prps data for batch and library size. This data needs to be further converted to log scale and transposed.
#' @param ruv.rep S4 data matrix for RUV-III: A S4 data object that has been generated using \code{replicate.matrix} functionality from \code{ruv} package. This helps ruv to identify replicate samples.
#' @param ncg.set logical object: Set of Negative Controlled genes.
#' @param k Integer scalar specifying the number of unwanted factors to use. Default is NULL. Currently Value 1 represents the library size, 2 represents purity and 3 is time variation.
#' @param eta Gene-wise (as opposed to sample-wise) covariates. A matrix with n columns for \code{ruv::RUV1}. Default is NULL.
#' @param include.intercept Add an intercept term to eta if it does not include one already for \code{ruv::RUV1}. Default is True.
#' @param average Default is False.
#' @param fullalpha To perform RUV-III calculation. Default is NULL.
#' @param return.info logical: Do you want all the information related to RUV-III object. False gives all information whereas True gives only the
#' @param inputcheck logical: Check the inputs to identify if ruv.data contains missing values or infinite values.
#'
#' @return Based on the return.info we get either a S4 list will all information related to RUV-III object or just the RUV-III result.
#' @export
#'
#' @examples
#' \dontrun{
#' runRUV_III_PRPS(ruv.data = ruv.data, ruv.rep = ruv.rep, ncg.set = ncg.set, k=1, return.info = TRUE)
#' }

runRUV_III_PRPS <- function(ruv.data, ruv.rep, ncg.set, k = NULL, eta = NULL,
                    include.intercept = TRUE, average = FALSE,
                    fullalpha = NULL, return.info = FALSE, inputcheck = TRUE){

  if (is.data.frame(ruv.data) ) {
    ruv.data <- data.matrix(ruv.data)
  }
  Y <- ruv.data
  M <- ruv.rep

  m <- nrow(ruv.data)
  n <- ncol(ruv.data)
  M <- ruv::replicate.matrix(ruv.rep)

  tological <- function(ctl, n) {
    ctl2 <- rep(FALSE, n)
    ctl2[ctl] <- TRUE
    return(ctl2)
  }

  ctl <- tological(ncg.set, n)
  if (inputcheck) {
    if (n > m)
      warning(
        "ncol for ruv.data is greater than nrow!  This is not a problem itself, but may
              indicate that you need to transpose your data matrix.
              Please ensure that rows correspond to observations
              (e.g. RNA-Seq assay) and columns correspond to features (e.g. genes).")
    if (sum(is.na(ruv.data)) > 0)
      stop("ruv.data contains missing values.  This is not supported.")
    if (sum(ruv.data == Inf, na.rm = TRUE) + sum(ruv.data == -Inf, na.rm = TRUE) >
        0)
      stop("ruv.data contains infinities.  This is not supported.")
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
    return(list(new.ruv.data = newY, ruv.rep = M, fullalpha = fullalpha,  W =  W))
  }
}
