##############################################################################/
##############################################################################/
#Comparison of Fst between dead and survivor 
##############################################################################/
##############################################################################/

library(adegenet)
library(hierfstat)
library(poppr)

source("RealTime_load.R")


##############################################################################/
#using the dosage format for the snp data
##############################################################################/

compSNP<-read.snp("data/data.snp")
metaIND<-read.table("data/metadatIND.snp",sep="\t",header=T,
                    stringsAsFactors=TRUE,na.strings="..")
metaMARK<-read.table("data/metadatMARK.snp",sep="\t",header=T,
                     stringsAsFactors=TRUE,na.strings="..")
compSNP@loc.names<-as.character(metaMARK$loc.names)
#adding other information to the @other slot for individuals
indNames(compSNP)<-metaIND$nom_Ind
compSNP@other$famille<-metaIND$famille
compSNP@other$mother<-metaIND$mother
compSNP@other$father<-metaIND$father
compSNP@other$traitement<-metaIND$traitement
compSNP@other$traitsimp<-metaIND$trait_simp
compSNP@other$live_bin<-metaIND$live_bin


#limiting the data set to individuals of the experimental setting 
experSNP<-compSNP[compSNP$pop!="PAR" & compSNP$pop!="CC" 
                  & compSNP$pop!="A4" & compSNP$pop!="3P"]

#some individuals don't have the information of their death
experSNP@ind.names[experSNP@other$live_bin=="NA"]

#we remove these individuals until clarification
temp<-experSNP[experSNP@other$live_bin!="NA"]


##############################################################################/
#F statistics based on the microsatellite data
##############################################################################/

#the number of individuals with quality checked microsatellite data
dim(micro.dat)[1]

#preparing the microsatellite data for importation in adegenet
temp<-c("A11_all1","A11_all2","A15_all1","A15_all2","A3_all1","A3_all2",
        "AB_all1","AB_all2","AK_all1","AK_all2","AO_all1","AO_all2",
        "C_all1","C_all2","D20_all1","D20_all2","D31_all1","D31_all2",
        "F_all1","F_all2","G_all1","G_all2","S19_all1","S19_all2")
microName<-c("A11","A15","A3","AB","AK","AO","C","D20","D31","F","G","S19")
temp2<-data.frame(matrix(vector(), dim(micro.dat)[1], 12,
                         dimnames=list(c(),microName)),
                  stringsAsFactors=F)
for (i in 1:(length(temp)/2)) {
  temp2[,i]<-
  paste(micro.dat[,temp[c(2*i-1,2*i)]][,1],
        micro.dat[,temp[c(2*i-1,2*i)]][,2],
        sep="/")
}

microGen<-df2genind(temp2,ploidy=2,sep="/",
                    ind.names=micro.dat$Sample_ID,
                    pop=micro.dat$family_simp)

#this includes individual that are not directly related to the experiment, 
#such as parents and the controlled crosses and 2 individuals that were 
#probably "contaminant"
table(micro.dat$family_simp)



##############################################################################/
#F statistics based on the snp data
##############################################################################/

#list of the snp marker name
temp<-colnames(snp.dat)[48:(48+818)]

snpGen<-df2genind(snp.dat[,temp],ploidy=2,sep="/",
                  ind.names=snp.dat$Sample_ID,
                  pop=snp.dat$family_simp)
#some markers are not polymorphic or display only a very limited polymorphism
table(minorAllele(snpGen)>0.95)
table(minorAllele(snpGen)<0.05)

#we limit the data to markers with a minor allele frequency > 0.05
snpGen2<-snpGen[,loc=minorAllele(snpGen)<0.95 & minorAllele(snpGen)>0.05]


##############################################################################/
#END
##############################################################################/