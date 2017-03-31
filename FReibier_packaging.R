library(devtools)
setwd("~/Data/aktuell/Misc/R-stuff/FReibier")
document() # process R-functions into .RD files, change namespace
setwd("..")
# build the .tar.gz file:
build("FReibier", path="FReibier")  
# install on the computer:
install("FReibier")




# test things:

library(FReibier)
data(polis)
polis

colSums(spiders8)
