library(devtools)
setwd("~/dropbox/FReibier")
document() # process R-functions into .RD files, change namespace
setwd("..")
# build the .tar.gz file:
build("FReibier")  
# install on the computer:
install("FReibier")