#' Range data between 0 and 1
#'
#' Re-scale the data to lie between 0 and 1
#'
#' A trivial function to re-scale data by subtracting the min and dividing by the difference between max and min. This function does not change the shape of the distribution.
#'
#' @param x Data vector to be ranged.
#' @param na.rm Logical; should NA-values be ignored when computing min and max? Defaults to TRUE:
#'
#' @return a vector of same length as x
#'
#' @author Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de>
#'
#' @seealso \code{scale}
#'
#' @examples
#' blubb <- rnorm(100, sd=15)
#' par(mfrow=c(1,2))
#' hist(blubb)
#' hist(rangeit(blubb))
#' par(mfrow=c(1,1))
#' @export
rangeit <- function(x, na.rm=TRUE){
  (x - min(x, na.rm=na.rm)) / (max(x, na.rm=na.rm) - min(x, na.rm=na.rm))
}