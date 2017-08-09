#' Wait for specific keypress
#'
#' This function allows to ask a question and prompts to a specific answer.
#' @param Q Question to be asked (Character).
#' @param A Vector of possible answers. Defaults to \code{c("y","n")}
#' @param Reminder If wrong answer was entered, reminds of the possible answers. Defaults to \code{"Please enter: y for yes, n for no."}
#' @return Returns the answer entered by the user (Character).
#'
#' @author Severin Hauenstein <severin.hauenstein@biom.uni-freiburg.de>
#'
#' @seealso Function is used in \code{\link{simData}}.
#' @export
#' @examples
#' # check with user if data really should be downloaded
#' \dontrun{
#' check.download <- keep.asking(Q = "Do you want to download data from http://www.worldclim.org?")
#' y
#' }
keep.asking <- function(Q, A = c("y","n"), Reminder = "Please enter: y for yes, n for no."){
  # print question
  cat(paste(Q, paste0(A, collapse = "/"), sep = "\t")) 
  check <- readline() # prompt answer
  while(!check %in% A){ # if answer not as requested
    cat(Reminder) 
    check <- readline() # prompt again
  }
  return(check)
}