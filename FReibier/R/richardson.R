#' Permute a community matrix to correct for uneven sampling
#'
#' Randomise a set of samples to estimate expected species richness according to Richardson & Richards (2008)
#'
#' This function shuffles (= randomises = permutes) the entries of a community matrix (= set of samples). As a result, the new matrix has the same number of individuals in each sample, randomly drawn from all individuals across all samples. By doing so repeatedly, an expected number of species can be computed and compared to the observed. It avoids using rarefaction to the lowest sample size and thus wasting information. The algorithm was presented in Richardson & Richards (2008). (Thanks to Caterina Penone for pointing out a stupid error of mine in an earlier version!)
#'
#' @param data an integer matrix with species in columns and sites/samples as rows (community data sensu vegan).
#' @param return.richness	logical; if only the richness values are of interest, not the entire shuffled community matrix, this parameter can be set to TRUE; defaults to FALSE.
#'
#' @return a matrix with same dimension and row/column names as the original data; or, if \code{return.richness=TRUE}, a named vector of species richness for each site/sample. 
#'
#' @author Carsten F. Dormann <carsten.dormann@biom.uni-freiburg.de>
#'
#' @references Richardson, J. M. L. & Richards, M. H. (2008) A randomisation program to compare species-richness values. \emph{Insect Conservation And Diversity} \bold{1}, 135--141.
#'
#' @examples
#' richardson(spiders)
#' rowSums(richardson(spiders))
#' rowSums(spiders)
#' colSums(richardson(spiders))
#' colSums(spiders)
#' # compare with original richness:
#' obs <- rowSums(spiders > 0)
#' rand <- replicate(20, richardson(spiders, return.richness=TRUE))
#' matplot(obs, jitter(rand), las=1, pch=1, col="black")
#' abline(0,1)
#' # above the line indicates observed data are clumped (= underdispersed)
#' 
#' @export
richardson <- function(data, return.richness=FALSE){
 # species names are the column names of the matrix
 list.of.individuals <- c()
 for (i in 1:nrow(data)){
   list.of.individuals  <- c(list.of.individuals, 
                             rep(colnames(data), times=data[i,])) # append the names to the vector
 }  

 species.per.site <- rowSums(data)
 inds <- 1:length(list.of.individuals)
 index <- list()

 for (j in 1:nrow(data)){
   index[[j]] <- sample(x=inds, size=species.per.site[j], replace=FALSE)
   inds <- inds[-which(inds %in% index[[j]])] # remove sampled individuals
 }

 #check:
 if (!all.equal(as.numeric(species.per.site), sapply(index, length))) stop("Unequal number of individuals: something when wrong here!")

 new.data <- data
 new.data[,] <- 0

 for (k in 1:nrow(new.data)){
   selection <- table(list.of.individuals[index[[k]]])
   new.data[k, match(names(selection), colnames(new.data))] <- selection
 }

 if (return.richness){
   out <- rowSums(new.data>0)
 } else {
   out <- new.data
 }

 return(out)

}