#' Compute AIC, AICc or BIC from R2 of a linear model
#'
#' Compute AIC-, AICc- or BIC-equivalent based on R2, sample size and number of model parameters (model rank)
#'
#' These functions are little helper functions if we need an AIC (or friends). The actual functions were taken from Burnham & Anderson (2002), but they also circulate the internet. Note that the actual value is different from the AIC() of a model, but the absolute differences between models should be comparable to those of differences between AICs.
#'
#' @aliases AICcFromR2 BICFromR2
#' @param R2 R2-value of the fitted linear model (not the adjusted R2).
#' @param n	 Sample size of the data used in the linear model.
#' @param p  Number of model parameters (i.e. the model's rank).
#'
#' @return AIC, AICc or BIC of the linear model, respectively; note that this value is supposed to be the AIC etc. up to a constant and will hence differ from the actual AIC of the model.
#'
#' @references Burnham, K.P. & Anderson, D.R. (2002) Model Selection and Multi-Model Inference: A Practical Information-Theoretical Approach. Springer, Berlin.
#'
#' @author Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de>
#'
#' @examples
#' set.seed(1)
#' x <- 1:100
#' y <- 2 + 2*x + rnorm(100, sd=5)
#' fm1 <- lm(y ~ x)
#' fm2 <- lm(y ~ poly(x, 3))
#' AICFromR2(summary(fm1)$r.squared, n=100, p=1)
#' AICFromR2(summary(fm1)$r.squared, n=100, p=3)
#' AIC(fm1, fm2)
#' AICcFromR2(summary(fm2)$r.squared, n=100, p=3)
#' BICFromR2(summary(fm1)$r.squared, n=100, p=3); BIC(fm1)
#' BICFromR2(summary(fm2)$r.squared, n=100, p=3); BIC(fm2)
#'
#' @export
AICFromR2 <- function(R2, n, p){
    # function to determine AICc for normally distributed data from their R2, sample size n and number of model parameters p
    # by Carsten F. Dormann
    (1-R2)*(n-1) + p*2 # according to Burnham & Anderson
}

#' @rdname AICFromR2
#' @export
AICcFromR2 <- function(R2, n, p){
    # function to determine AICc for normally distributed data from their R2, sample size n and number of model parameters p
    # by Carsten F. Dormann
    (1-R2)*(n-1) + p*2* (p+1)/(n-p-1) # according to Burnham & Anderson
}

#' @rdname AICFromR2
#' @export
BICFromR2 <- function(R2, n, p){
    # function to determine BIC for normally distributed data from their R2, sample size n and number of model parameters p
    # by Carsten F. Dormann
    (1-R2)*(n-1) + p*log(n)
}
