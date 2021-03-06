#' Construct a GLM-formula from a data.frame-object
#'
#' Constructs a formula-object from a \code{data.frame}, including quadratic terms and interactions
#'
#' This function reduces typing of lengthy variable names into a GLM-formula. Also, one doesn't have to think of how to nest \code{paste}-calls to achieve the same. And using the defaults it automatically puts in quadratic and first-order interaction terms, which probably are missing from many GLMs due to lazyness in typing them out.
#'
#' @param dataframe \code{data.frame} of the data to be used in the GLM, containing both the left- and right-hand side variables to be used. ALL columns will be used, thus you have to use subset or [,] to index the part of the data.frame you want represented in the analysis.
#' @param y.col	integer indicating the column in which to find the response variable (i.e. the left-hand side of the formula); defaults to 1.
#' @param quadratic  logical; TRUE (default) indicates that a quadratic term for each predictor will be added to the model formula; this makes only sense if all predictors are continuous.
#' @param interaction logical; TRUE (default) indicates that all first-order interactions will be added to the model formula
#'
#' @return a \code{formula} object to be used in a GLM
#'
#' @author Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de>
#'
#' @note This function could be expanded using ellipses to add terms ad lib; or to indicate which columns to use for X; or to indicate which are categorical.
#'
#' @examples
#' data(ChickWeight)
#' f <- formulaMaker(ChickWeight, quadratic=FALSE)
#' f
#' anova(lm(f, data=ChickWeight))
#'
#' @export
formulaMaker <- function(dataframe, y.col=1, quadratic=TRUE, interactions=TRUE){
  # makes a formula for GLM from dataframe column names, 
  # including quadratic effects and first-order interactions
  # by default, first column is taken to be the response (y); else, an integer giving the column with the response in "dataframe"
  # by Carsten F. Dormann
  if (quadratic && interactions) {
      f <- as.formula(paste(colnames(dataframe)[y.col], " ~ (", paste(colnames(dataframe[,-y.col]), collapse=" + ", sep=""), ")^2 + ", paste("I(", colnames(dataframe[,-y.col]), "^2)", collapse="+", sep="")))
  }
  
  if (quadratic & !interactions){
      f <- as.formula(paste(colnames(dataframe)[y.col], " ~ (", paste(colnames(dataframe[,-y.col]), collapse=" + ", sep=""), ") + ", paste("I(", colnames(dataframe[,-y.col]), "^2)", collapse="+", sep="")))
  }

  if (!quadratic & !interactions){
      f <- as.formula(paste(colnames(dataframe)[y.col], " ~ ", paste(colnames(dataframe[,-y.col]), collapse=" + ", sep="") ))
  }

  if (!quadratic & interactions){
      f <- as.formula(paste(colnames(dataframe)[y.col], " ~ (", paste(colnames(dataframe[,-y.col]), collapse=" + ", sep=""), ")^2"))
      # + ", paste("I(", colnames(dataframe[,-1]), "^2)", collapse="+", sep="")))
  }
  
  return(f)
}
