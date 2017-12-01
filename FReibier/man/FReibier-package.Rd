\encoding{UTF-8}
\name{FReibier-package}

\alias{FReibier-package}
\alias{FReibier}

\docType{package}

\title{
Miscellaneous R-functions of the Biometry group in Freiburg, Germany
}

\description{
This package contains an assortment of functions: AIC from R2, creating quadratic and interactions terms for a GLM-formula, error bars for ranger predictions, generating spatially autocorrelated data.}

\details{
This is a container for functions that otherwise would lie only on our computers. Some have a certain turnover and become deprecated as a better place is found for them (e.g. the rInfJack for ranger obviously sits better in \pkg{ranger}, as it now does).

\tabular{ll}{
Package: \tab Freibier \cr
Type: \tab Package\cr
Version: \tab 0.4 \cr
Date: \tab 2017-12-01\cr
License: \tab GPL \cr
}
}

\section{versionlog}{

\itemize{
  \item 0.4
  \describe{
    \item{Deprecated}{all functions related to computing CI for ranger (which are now available wihtin \pkg{ranger}, at least on github).}
    \item{I'm sure}{there will be more.}
  }
} % end of versionlog's itemize
} % end of section versionlog

\author{
Carsten F. Dormann

Maintainer: Carsten Dormann \email{carsten.dormann@biom.uni-freiburg.de}
}

\references{
}

\keyword{package}
