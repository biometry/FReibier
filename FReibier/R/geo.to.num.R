#' Transform geographic coordinates to numeric values
#'
#' This function transforms geographic coordinates, such as \code{"5N24E"} or \code{"7S27E"} to numeric values, such as \code{c(5, 24)} or \code{c(-7, 27)}.
#'
#' @param geo Geographic coordinate (Character).
#'
#' @return Returns the coordinates as numeric values. To be used for the creation of (\code{\link[raster]{extent}}) objects.
#'
#' @author Severin Hauenstein <severin.hauenstein@biom.uni-freiburg.de>
#'
#' @seealso Function is used in \code{\link{simData}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Transform "5N24E" or "7S27E" to numeric values
#' coords <- c("5N24E", "7S27E")
#' num.coords <- sapply(1:2, function(x) geo.to.num(coords[x]))
#' }
geo.to.num <- function(geo){
  # split geo into indiv. characters
  char <- strsplit(geo, split = "")[[1]]
  
  # extract North value 
  if("N" %in% char){ 
    w.SN <- grep("N", char)
    N <- as.numeric(Reduce(paste0, char[1:(w.SN-1)]))
  } else { # or if given as South (S) multiply with -1
    if ("S" %in% char){
      w.SN <- grep("S", char)
      N <- as.numeric(Reduce(paste0, char[1:(w.SN-1)])) * (-1)
    } else { # spelling error, if criteria are not met
      stop("Please define geographic coordinates, e.g.
           5N24E or 7S27E")
    }
  }
  # extract East value
  if (char[length(char)] == "E"){ # must be last Element
    E <- as.numeric(Reduce(paste0, char[(w.SN+1):(length(char)-1)]))
  } else {
    if (char[length(char)] == "W"){ # or if given as West (W) multiply with -1
      E <- as.numeric(Reduce(paste0, char[(w.SN+1):(length(char)-1)])) * (-1)
    } else { # spelling error, if criteria are not met
      stop("Please define geographic coordinates, e.g.
           5N24E or 7S27E")
    }
  }
  return(c(N,E))
  }