##############################################################################/
##############################################################################/
#Loading the data sets and packages needed for RealTime analyses
##############################################################################/
##############################################################################/

##############################################################################/
#Loading the libraries####
##############################################################################/

library(adegenet)
library(ape)
library(FinePop2)
library(GAPIT3)
library(gdata)
library(genepop)
library(graph4lg)
library(gtools)
library(hierfstat)
library(kinship2)
library(plotrix)
library(PopGenReport)
library(poppr)
library(RColorBrewer)
library(visreg)
library(stringr)
library(tidyr)
library(vegan)
library(vioplot)
library(viridis)


##############################################################################/
#Loading and preparing the main data set####
##############################################################################/

#data by individuals, including all individuals (hd, families, parents and CC
RTdata<-read.table("data/datatot.txt",sep="\t",stringsAsFactors=FALSE,
                   header=TRUE)
#because some functions do not like "." within colnames, we replace them
colnames(RTdata)<-str_replace_all(colnames(RTdata),"[.]","_")
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
micro.dat<-RTdata[RTdata$na_9micro<3 & !is.na(RTdata$na_9micro),]
micro.dat<-drop.levels(micro.dat)

#data subset of individuals with snp data. The quality of snp genotyping
#as been assessed beforehand and can be sorted with the 'Quality_SNPage' 
#column
snp.dat<-RTdata[RTdata$SNPage==1 & !is.na(RTdata$SNPage) 
                & RTdata$Quality_SNPage==1,]
snp.dat<-drop.levels(snp.dat)


##############################################################################/
#Loading the necessary data sets for GWAS analyses####
##############################################################################/

AllTrait<-snp.dat[snp.dat$family_simp!="PAR" & snp.dat$family_simp!="CC",]
AllTrait<-AllTrait[,c("treat","Sample_ID","oid4_09","oid5_10",
                      "oid5_11","oid4_12","oid2_13","oid_moy",
                      "H12v","pgland","statut10")]
#rename the columns with more accurate names
colnames(AllTrait)<-c("treat","Taxa","oid4_09","oid5_10","oid5_11",
                      "oid4_12","oid2_13","Powdery mildew","Height",
                      "Acorn weight","Survival")
NatTrait<-AllTrait[AllTrait$treat=="exp",2:11]
LimTrait<-AllTrait[AllTrait$treat=="low",2:11]

# #to remove
# AllTrait<-read.table("data/pheno_final.txt",header=TRUE)
# colnames(AllTrait)[11:14]<-c("Powdery mildew","Height",
#                              "Acorn weight","Dead or Alive")
# #we remove the individuals without SNP data
# AllTrait<-AllTrait[AllTrait$SNPage==1 & !is.na(AllTrait$SNPage) &
#                      AllTrait$Quality_SNPage==1,]
# NatTrait<-AllTrait[AllTrait$treat=="exp",5:14]
# LimTrait<-AllTrait[AllTrait$treat=="low",5:14]
# #end remove

#loading the natural treatment genotype data
NatG<-read.delim("data/nat.hmp.txt",header=FALSE)

#loading the limited treatment genotype data
LimG<-read.delim("data/lim.hmp.txt",header=FALSE)


##############################################################################/
#Writing session information for reproducibility####
##############################################################################/

sink("session_info.txt")
print(sessioninfo::session_info())
sink()
#inspired by an R gist of François Briatte: 
#https://gist.github.com/briatte/14e47fb0cfb8801f25c889edea3fcd9b


##############################################################################/
#END
##############################################################################/