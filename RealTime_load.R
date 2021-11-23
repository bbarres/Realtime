##############################################################################/
##############################################################################/
#Loading the data sets and packages needed for RealTime analyses
##############################################################################/
##############################################################################/

#loading the libraries
library(adegenet)
library(gdata)
library(kinship2)
library(RColorBrewer)
library(viridis)


#loading the different data sets

#information of the snp markers
#mrkrInf<-read.table()

#data by individuals, including all individuals (hd, families, parents and CC
RTdata<-read.table("data/datatot.txt",sep="\t",stringsAsFactors=FALSE,
                   header=TRUE)
#reordering the family levels
RTdata$family_simp<-factor(RTdata$family_simp,
                           levels=c("1","9","11","26","27","45","48",
                                    "51","70","71","72","73","74","75",
                                    "76","77","99","CC","PAR","hd"))

#limiting the data set to the experimental set up
dispo<-RTdata[RTdata$family_simp!="PAR",]
dispo<-drop.levels(dispo)

#data subset of individuals with microsatellite data. There is a total of 
#12 markers, but 3 were found to be of limited quality (too many missing
#data). Therefore the quality criteria to remove individuals is a maximum 
#of 2 microsatellite missing data for the 9 best microsatellite
micro.dat<-RTdata[RTdata$na.9micro<3 & !is.na(RTdata$na.9micro),]
micro.dat<-drop.levels(micro.dat)

#data subset of individuals with snp data. The quality of snp genotyping
#as been assessed beforehand and can be sorted with the 'Quality_SNPage' 
#column
snp.dat<-RTdata[RTdata$SNPage==1 & !is.na(RTdata$SNPage) 
                & RTdata$Quality_SNPage==1,]
snp.dat<-drop.levels(snp.dat)


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