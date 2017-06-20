#' Compute variance of estimate from a "ranger" model
#' 
#' Computes variances for a prediction from a ranger-model, using the infinitesimal jackknife procedure
#' 
#' This function is a ranger-adaptation of the package \pkg{randomForestCI} of Wager et al. (2014). Their original can be found on github: \url{ https://github.com/swager/randomForestCI/}. All that the three functions below do differently is to allow using the ranger model object (which is computed much faster than the classical randomForest in \pkg{randomForest}). The functions are named differently to avoid confusion.
#' 
#' For classification problems, use \code{rangerInfJackMulticlass}; for regression problems use \code{rangerInfJack}. \code{rInfJack} is the internal workhorse called by the other two functions. It is probably of little relevance for the typical user.
#'
#' @aliases rangerInfJack
#' @aliases rangerInfJackMultiClass
#' @aliases rInfJack
#'
#' @param rf A ranger object, fitted using the function ranger in \pkg{ranger}, as a faster alternative to \pkg{randomForest}.
#'
#' @param pred A nrow(newdata) by no. of trees matrix which contains numeric predictions
#'        from a random forest trained with trees grown on bootstrap samples of the training data
#' @param inbag A number of obs. in the training data by no. of trees matrix giving the
#'        number of times the ith observation in the training data appeared in the bootstrap sample for the jth tree.
#' @param calibrate whether to apply calibration to mitigate Monte Carlo noise
#'        warning: if calibrate = FALSE, some variance estimates may be negative
#'                 due to Monte Carlo effects if the number of trees in rf is too small
#' @param used.trees set of trees to use for variance estimation; uses all tress if NULL
#'
#' @return For regressions, a two-column matrix is returned, with predictions in the first column and estimates of prediction variance in the second. For classification, a 2*k-column matrix is returned, with the probability estimates for the k classes in the first k columns, and their variance estimates in the next k columns.
#' 
#' @note The warning "No calibration with n<= 20" indicates that the variance estimates were not re-calibrated across the predictions. Please read the original paper to understand what that means (Wager et al. 2014).
#' 
#' @references Wager, Stefan, Trevor Hastie, and Bradley Efron. 2014. Confidence intervals for random forests: The jackknife and the infinitesimal jackknife. \emph{The Journal of Machine Learning Research} \bold{15} (1), 1625--1651.
#'
#' @author Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de>
#' 
#' @examples 
#' \dontrun{ # because the packages with the data are not in the dependencies of FReibier
#' # a binary classification example using the Titanic-data:
#' library(effects)
#' franger <- ranger(survived ~ age+sex+passengerClass, data=na.omit(TitanicSurvival), keep.inbag=T, replace=T)
#' rangerInfJackMulticlass(franger, newdata=TitanicSurvival[1:2,], calibrate=F)
#' 
#' # a regression example using the alfalfa data:
#' library(faraway)
#' fr <- ranger(yield ~ shade + irrigation + inoculum, data=alfalfa, keep.inbag=T, replace=T)
#' rangerInfJack(fr, newdata=alfalfa[1:2,], calibrate=F)
#' }

# The infinitesimal jackknife for random forests
#
# @param rf A random forest trained with replace = TRUE and keep.inbag = TRUE
# @param newdata A set of test points at which to evaluate standard errors
# @param calibrate whether to apply calibration to mitigate Monte Carlo noise
#        warning: if calibrate = FALSE, some variance estimates may be negative
#                 due to Monte Carlo effects if the number of trees in rf is too small
#' @rdname rangerInfJack
#' @export

rangerInfJack = function(rf, newdata, calibrate = TRUE) { 
  
  #
  # Setup
  #
  
  # produce a new entry of inbag as matrix:
  rf$inbag <- matrix(unlist(rf$inbag), nrow=length(rf$inbag[[1]]), ncol=length(rf$inbag))
  
  if (is.null(rf$inbag)) {
    stop("Random forest must be trained with keep.inbag = TRUE")
  }
  
  if (length(levels(factor(colSums(rf$inbag)))) > 1) {
    stop("The keep.inbag field must store the number of times each observation was used. Please make sure the version number of randomForest is 4.6-12 or higher.")
  }
  
  predictions = predict(rf, newdata, predict.all = TRUE)
  pred = predictions$predictions# predictions$individual
  # in case of classification, convert character labels to numeric (!)
  class(pred) = "numeric"
  
  results = rInfJack(pred, inbag=rf$inbag, calibrate, used.trees = NULL)
  return (results)
}

# The infinitesimal jackknife for random forests (multiclass target variable)
#
# @param rf A random forest trained with replace = TRUE and keep.inbag = TRUE
# @param newdata A set of test points at which to evaluate standard errors
# @param calibrate whether to apply calibration to mitigate Monte Carlo noise
#        warning: if calibrate = FALSE, some variance estimates may be negative
#                 due to Monte Carlo effects if the number of trees in rf is too small
#' @rdname rangerInfJack
#' @export

rangerInfJackMulticlass = function(rf, newdata, calibrate = TRUE) {
  
  #
  # Setup
  #
  
  if (is.null(rf$inbag)) {
    stop("Random forest must be trained with keep.inbag = TRUE")
  }
  
  # produce a new entry of inbag as matrix:
  rf$inbag <- matrix(unlist(rf$inbag), nrow=length(rf$inbag[[1]]), ncol=length(rf$inbag))
  
  if (length(levels(factor(colSums(rf$inbag)))) > 1) {
    stop("The keep.inbag field must store the number of times each observation was used. Please make sure the version number of randomForest is 4.6-12 or higher.")
  }
  
  predictions = predict(rf, newdata, predict.all = TRUE)
  pred = predictions$predictions #predictions$individual
  
  #number of classes
  rf$classes <- sort(unique(as.numeric(rf$predictions))) # cfd
  K <- length(rf$classes)
  
  #create empty dataframe to store results
  results.full <- data.frame(matrix(NA, nrow = nrow(newdata), ncol = K*2))
  
  #estimating variance for each class predictions
  for (k in seq(1, K)){
    #separate prediction matrix into K matrices following 1 vs. all logic
    pred.binary <- ifelse(pred == rf$classes[k], 1, 0) 
    #apply infinitessimal jackknife (ranger version)
    results <- rInfJack(pred.binary, rf$inbag, calibrate, used.trees = NULL)
    #extend column names of the final result with classes
    colnames(results.full)[c(k,k+K)] <- paste(colnames(results), rf$classes[k], sep = "_")
    #store results in the final dataframe
    results.full[,k] <- results[,1]
    results.full[,k+K] <- results[,2]
  }
  
  return (results.full)
}

#' @rdname rangerInfJack
#' @export

rInfJack = function(pred, inbag, calibrate = TRUE, used.trees = NULL) {
  
  # original: https://github.com/swager/randomForestCI/blob/master/R/infinitesimalJackknife.R
  
  if (is.null(used.trees)) {
    used.trees = 1:ncol(inbag)
  }
  
  pred = pred[, used.trees, drop=FALSE]
  
  # check if sampling without replacement
  no.replacement = (max(inbag) == 1)
  
  #
  # Extract tree-wise predictions and variable counts from random forest
  #
  
  B = length(used.trees)
  n = nrow(inbag)
  s = sum(inbag) / ncol(inbag)
  
  y.hat = rowMeans(pred)
  pred.centered = pred - rowMeans(pred)
  
  N = Matrix::Matrix(inbag[, used.trees], sparse = TRUE)
  N.avg = Matrix::rowMeans(N)
  
  #
  # Compute raw infinitesimal jackknife
  #
  
  if (B^2 > n * nrow(pred)) {
    
    C = Matrix::tcrossprod(N, pred.centered) -
      Matrix::Matrix(N.avg, nrow(N), 1) %*%
      Matrix::Matrix(rowSums(pred.centered), 1, nrow(pred.centered))
    raw.IJ = Matrix::colSums(C^2) / B^2
    
  } else {
    
    # Faster implementation when n is large. Uses the fact that
    # colSums((A - B)^2) = T1 - 2 * T2 + T3,
    # where T1 = diag(A'A), T2 = diag(B'A), and T3 = diag(B'B)
    
    NTN = Matrix::crossprod(N, N)
    NTNPT_T = Matrix::tcrossprod(pred.centered, NTN)
    T1 = Matrix::rowSums(pred.centered * NTNPT_T)
    
    RS = rowSums(pred.centered)
    NbarTN = Matrix::crossprod(N.avg, N)
    T2 = RS * Matrix::tcrossprod(NbarTN, pred.centered)
    
    T3 = sum(N.avg^2) * RS^2
    raw.IJ = as.numeric(T1 - 2 * T2 + T3) / B^2
    
  }
  
  #
  # Apply Monte Carlo bias correction
  #
  
  N.var = mean(Matrix::rowMeans(N^2) - Matrix::rowMeans(N)^2)
  boot.var = rowSums(pred.centered^2) / B
  bias.correction = n * N.var * boot.var / B
  vars = raw.IJ - bias.correction
  
  #
  # Finite sample correction
  #
  
  if (no.replacement) {
    
    variance.inflation = 1 / (1 - mean(inbag))^2
    vars = variance.inflation * vars
  }
  
  results = data.frame(y.hat=y.hat, var.hat=vars)
  
  if (nrow(results) <= 20) {
    calibrate = FALSE
    warning("No calibration with n <= 20")
  }
  
  #
  # If appropriate, calibrate variance estimates; this step in particular
  # ensures that all variance estimates wil be positive.
  #
  
  if (calibrate) {
    # Compute variance estimates using half the trees
    calibration.ratio = 2
    n.sample = ceiling(B / calibration.ratio)
    results.ss = rangerInfJack(pred, inbag, calibrate = FALSE, used.trees = sample(used.trees, n.sample))
    
    # Use this second set of variance estimates to estimate scale of Monte Carlo noise
    sigma2.ss = mean((results.ss$var.hat - results$var.hat)^2)
    delta = n.sample / B
    sigma2 = (delta^2 + (1 - delta)^2) / (2 * (1 - delta)^2) * sigma2.ss
    
    # Use Monte Carlo noise scale estimate for empirical Bayes calibration
    vars.calibrated = calibrateEB(vars, sigma2)
    results$var.hat = vars.calibrated
  }
  
  
  return(results)
}
