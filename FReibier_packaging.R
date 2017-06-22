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

richardson(spiders8)
solow(spiders8)

library(effects)
franger <- ranger(survived ~ age+sex+passengerClass, data=na.omit(TitanicSurvival), keep.inbag=T, replace=T)
rangerInfJackMulticlass(franger, newdata=TitanicSurvival[1, , drop=F], calibrate=F)
rangerInfJackMulticlass(franger, newdata=na.omit(TitanicSurvival[1:300,]), calibrate=T)

library(faraway)
fr <- ranger(yield ~ shade + irrigation + inoculum, data=alfalfa, keep.inbag=T, replace=T)
rangerInfJack(fr, newdata=alfalfa[1:2,], calibrate=F)
