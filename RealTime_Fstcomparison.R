##############################################################################/
##############################################################################/
#Comparison of Fst between dead and survivor 
##############################################################################/
##############################################################################/

library(genepop)
library(FinePop2)
library(PopGenReport)
library(graph4lg)

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
temp2<-data.frame(matrix(vector(),dim(micro.dat)[1],12,
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
#adding information in the 'other' slot
microGen@other$fam<-micro.dat$family_simp
microGen@other$treat<-micro.dat$exp
microGen@other$height<-micro.dat$Hdeb17
microGen@other$DoA<-micro.dat$live_bin
#combine experimental and dead or alive information
microGen@other$newPop<-paste(micro.dat$exp,micro.dat$live_bin,sep="")

#this includes individual that are not directly related to the experiment, 
#such as parents and the controlled crosses and 2 individuals that were 
#probably "contaminant"
table(micro.dat$family_simp)
#we limit the data set to the individuals belonging to the experimental
#set up
microGen<-microGen[(microGen@other$fam!="CC" & microGen@other$fam!="PAR")]
#we also remove individuals without dead or alive information
microGen<-microGen[!is.na(microGen@other$DoA)]
#checking the minor allele frequencies of the markers
minorAllele(microGen)

#set newPop as population
pop(microGen)<-microGen@other$newPop
table(pop(microGen))
pairwise.WCfst(genind2hierfstat(microGen))
genind_to_genepop(microGen,output="data/microGenAD.txt")
microGenAD<-"data/microGenAD.txt"
genedivFis(microGenAD,sizes=FALSE,"output/microGenAD.txt.Fis")
n.temp<-seppop(microGen)
Hobs<-do.call("c",lapply(n.temp,function(x) mean(summary(x)$Hobs)))
Hexp<-Hs(microGen)
Arich<-colMeans(allelic.richness(microGen,min.n=100)$Ar)
Ali.vs.Dea<-rbind(Hobs,Hexp,Arich)

#total vs surviving
temp<-repool(n.temp$exp1,n.temp$exp0)
pop(temp)<-rep("exp_init",times=nInd(temp))
temp<-repool(temp,n.temp$exp1)
temp2<-repool(n.temp$low1,n.temp$low0)
pop(temp2)<-rep("low_init",times=nInd(temp2))
temp<-repool(temp,temp2,n.temp$low1)
#number of individuals by populations
table(pop(temp))
pairwise.WCfst(genind2hierfstat(temp))
genind_to_genepop(temp,output="data/microGenDF.txt")
microGenAD<-"data/microGenDF.txt"
genedivFis(microGenAD,sizes=FALSE,"output/microGenDF.txt.Fis")
n.temp<-seppop(temp) 
Hobs<-do.call("c",lapply(n.temp,function(x) mean(summary(x)$Hobs)))
Hexp<-Hs(temp)
Arich<-colMeans(allelic.richness(temp,min.n=100)$Ar)
Deb.vs.Fin<-rbind(Hobs,Hexp,Arich)

#compute genetic distance between individuals
distMicro<-diss.dist(microGen,mat=FALSE)
treeMicro<-nj(distMicro)
plot(treeMicro,type="unr",show.tip=FALSE)
colovec<-c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(1,4,8)])
colovec<-c(colovec[1:3],"#ffffff",colovec[4:15],"#ffffff")
tiplabels(pch=20,col=colovec[as.numeric(microGen@other$fam)],cex=2)


##############################################################################/
#F statistics based on the snp data
##############################################################################/

#list of the snp marker names
temp<-colnames(snp.dat)[48:(48+818)]

snpGen<-df2genind(snp.dat[,temp],ploidy=2,sep="/",
                  ind.names=snp.dat$Sample_ID,
                  pop=snp.dat$family_simp)
#adding information in the 'other' slot
snpGen@other$fam<-snp.dat$family_simp
snpGen@other$treat<-snp.dat$exp
snpGen@other$height<-snp.dat$Hdeb17
snpGen@other$DoA<-snp.dat$live_bin
#combine experimental and dead or alive information
snpGen@other$newPop<-paste(snp.dat$exp,snp.dat$live_bin,sep="")

#this includes individual that are not directly related to the experiment, 
#such as parents and the controlled crosses and 2 individuals that were 
#probably "contaminant"
table(snp.dat$family_simp)
#we limit the data set to the individuals belonging to the experimental
#set up
snpGen<-snpGen[(snpGen@other$fam!="CC" & snpGen@other$fam!="PAR")]
#we also remove individuals without dead or alive information
snpGen<-snpGen[!is.na(snpGen@other$DoA)]

#some markers are not polymorphic or display only a very limited polymorphism
table(minorAllele(snpGen)>0.95)
table(minorAllele(snpGen)<0.05)

#we limit the data to markers with a minor allele frequency > 0.05
snpGen2<-snpGen[,loc=minorAllele(snpGen)<0.95 & minorAllele(snpGen)>0.05]

#set newPop as population
pop(snpGen2)<-snpGen2@other$newPop
table(pop(snpGen2))
pairwise.WCfst(genind2hierfstat(snpGen2))
n.temp<-seppop(snpGen2) 
Hobs<-do.call("c",lapply(n.temp,function(x) mean(summary(x)$Hobs)))
Hexp<-Hs(snpGen2)
Arich<-colMeans(allelic.richness(snpGen2,min.n=100)$Ar)
Ali.vs.Dea<-rbind(Hobs,Hexp,Arich)

#total vs surviving
temp<-repool(n.temp$exp1,n.temp$exp0)
pop(temp)<-rep("exp_init",times=nInd(temp))
temp<-repool(temp,n.temp$exp1)
temp2<-repool(n.temp$low1,n.temp$low0)
pop(temp2)<-rep("low_init",times=nInd(temp2))
temp<-repool(temp,temp2,n.temp$low1)
#number of individuals by populations
table(pop(temp))
pairwise.WCfst(genind2hierfstat(temp))
n.temp<-seppop(temp) 
Hobs<-do.call("c",lapply(n.temp,function(x) mean(summary(x)$Hobs)))
Hexp<-Hs(temp)
Arich<-colMeans(allelic.richness(temp,min.n=100)$Ar)
Deb.vs.Fin<-rbind(Hobs,Hexp,Arich)

#compute genetic distance between individuals
distSnp2<-diss.dist(snpGen2,mat=FALSE)
treeSnp2<-nj(distSnp2)
plot(treeSnp2,type="unr",show.tip=FALSE)
colovec<-c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(1,4,8)])
tiplabels(pch=20,col=colovec[as.numeric(snpGen2@other$fam)],cex=2)


##############################################################################/
#END
##############################################################################/