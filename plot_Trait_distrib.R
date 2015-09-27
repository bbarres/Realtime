###############################################################################
###############################################################################
#Boxplot of the trait
###############################################################################
###############################################################################

#loading the libraries
library(gdata)

#set the working directory
setwd("~/work/Rfichiers/Githuber/Realtime_data")


###############################################################################
#Make boxplot for the size of the tree in 2015
###############################################################################

#first load the dataset
dataSNP<-read.table("Size2015.txt", header=TRUE, sep="\t")
#remove individuals with missing data
dataSNP<-dataSNP[!is.na(dataSNP$H_mai2015),]
#limit the dataset to the treatment 'nat' and 'renf'
dataSNP<-dataSNP[dataSNP$traitement!="lim",]
dataSNP<-drop.levels(dataSNP)

#extract the total number of individuals for each genotype
toteff<-table(dataSNP[,5])
#extract the total number of dead individuals for each genotype
dead<-table(dataSNP[dataSNP$H_mai2015==0,5])
#we plot the percentage of dead trees for each genotype
barplot(dead/toteff,main="% of dead trees",xlab="Genotype",
        names=paste(names(toteff)," (n=",toteff,")",sep=""),
        las=1,width=0.3,cex.names=0.7,horiz=FALSE)


#we now plot a boxplot of the size of tree still alive in 2015
dataSNP<-dataSNP[dataSNP$H_mai2015!=0,]
dd<-boxplot(dataSNP[,6]~dataSNP[,5],plot=F)
boxplot(dataSNP[,6]~dataSNP[,5],boxwex=0.3,
        names=paste(dd$names," (n=",dd$n,")",sep=""),
        main="SNP M97",xlab="Genotype",ylab="2015 Tree size (cm)")

#cleaning the R environment
rm(dataSNP,hh,dd,toteff,dead)

###############################################################################
#END
###############################################################################