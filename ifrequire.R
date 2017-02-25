ifrequire <- function(name){
	# little helper function to install package is not yet available
	if (name %in% installed.packages()[,1]){
		library(name, character.only=TRUE)	
	} else {
		install.packages(name)
		library(name, character.only=TRUE)
	}
}
