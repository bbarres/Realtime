##############################################################################/
##############################################################################/
#Loading the data sets and packages needed for RealTime analyses
##############################################################################/
##############################################################################/

#loading the libraries
library(adegenet)
library(kinship2)

#loading the different necessary data set



##############################################################################/
#Writing info session for reproducibility####
##############################################################################/

sink("session_info.txt")
print(sessioninfo::session_info())
sink()
#inspired by an R gist of François Briatte: 
#https://gist.github.com/briatte/14e47fb0cfb8801f25c889edea3fcd9b


##############################################################################/
#END
##############################################################################/