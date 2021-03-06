% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solow.R
\name{solow}
\alias{solow}
\alias{richardson}
\alias{richardson}
\title{Permute a community matrix to correct for uneven sampling}
\usage{
solow(data, return.richness = FALSE)

richardson(data, return.richness = FALSE)
}
\arguments{
\item{data}{an integer matrix with species in columns and sites/samples as rows (community data sensu vegan).}

\item{return.richness}{logical; if only the richness values are of interest, not the entire shuffled community matrix, this parameter can be set to TRUE; defaults to FALSE.}
}
\value{
a matrix with same dimension and row/column names as the original data; or, if \code{return.richness=TRUE}, a named vector of species richness for each site/sample.
}
\description{
Randomise a set of samples to estimate expected species richness according to Solow (1993) or Richardson & Richards (2008)
}
\details{
This function shuffles (= randomises = permutes) the entries of a community matrix (= set of samples). As a result, the new matrix has the same number of individuals in each sample, randomly drawn from all individuals across all samples. By doing so repeatedly, an expected number of species can be computed and compared to the observed. It avoids using rarefaction to the lowest sample size and thus wasting information. The algorithm was presented in Richardson & Richards (2008) but dates back to Solow (1993). Both names are available as function names. (Thanks to Caterina Penone for pointing out a stupid error of mine in an earlier version!)
}
\examples{
solow(spiders)
rowSums(richardson(spiders))
rowSums(spiders)
colSums(solow(spiders))
colSums(spiders)
# compare with original richness:
obs <- rowSums(spiders > 0)
rand <- replicate(20, richardson(spiders, return.richness=TRUE))
matplot(obs, jitter(rand), las=1, pch=1, col="black")
abline(0,1)
# above the line indicates observed data are clumped (= underdispersed)

}
\references{
Richardson, J. M. L. & Richards, M. H. (2008) A randomisation program to compare species-richness values. \emph{Insect Conservation And Diversity} \bold{1}, 135--141

Solow, A.R. (1993) A simple test for change in community structure. \emph{J. Animal Ecology} \bold{62}, 191-–193
}
\author{
Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de>
}
