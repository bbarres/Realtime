##############################################################################/
##############################################################################/
#GWAS analyses using GAPIT3
##############################################################################/
##############################################################################/

# #at the time of writing this code, GAPIT3 was still under development on 
# #Github, so you have to run this to install the updated package
# devtools::install_github("jiabowang/GAPIT3",force=TRUE)

source("RealTime_load.R")
library(vioplot)

#to be removed when the code will be working properly
setwd("/myGAPIT")


##############################################################################/
#Model Blink/kinship for natural inoculum condition####
##############################################################################/

AllTrait<-read.table("pheno_final.txt",header=TRUE)
#we remove the individuals without SNP data
AllTrait<-AllTrait[AllTrait$SNPage==1 & is.na(AllTrait$SNPage)!=TRUE &
                     AllTrait$Quality_SNPage==1,]
NatTrait<-AllTrait[AllTrait$treat=="exp",5:14]
LimTrait<-AllTrait[AllTrait$treat=="low",5:14]

#loading genotype data
NatG<-read.delim("nat.hmp.txt",header=FALSE)
#alternative kinship matrix
NatLois<-read.table("nat_lois_mat.txt",header=FALSE)

#Blink method on powdery mildew phenotype
natGAPIT<-GAPIT(
  Y=NatTrait[,1:7],
  G=NatG,
  kinship.algorithm="Loiselle",
  #KI=NatLois,
  PCA.total=0,
  model="Blink"
  #,Random.model=TRUE
)

#Blink method on general phenotype
natGAPIT<-GAPIT(
  Y=NatTrait[,c(1,7:10)],
  G=NatG,
  kinship.algorithm="Loiselle",
  #KI=NatLois,
  PCA.total=0,
  model="Blink"
  #,Random.model=TRUE
)

#MLM method on powdery mildew trait
natGAPIT<-GAPIT(
  Y=NatTrait[,1:7],
  G=NatG,
  kinship.algorithm="Loiselle",
  #KI=NatLois,
  PCA.total=0,
  model="MLM"
)

colovec<-brewer.pal(12,"Paired")
plot(density(NatTrait$oid_moy,na.rm=TRUE))
plot(NatTrait$oid_moy[order(NatTrait$oid_moy)])
vioplot(NatTrait$oid_moy,plotCentre="points",col=colovec[5],border=colovec[5],
        lineCol="white",rectCol=colovec[3],names="oid_moy",las=1)
stripchart(NatTrait$oid_moy,method="jitter",col="red",
           vertical=TRUE,pch=19,cex=0.7,add=TRUE)
#export to .pdf 4 x 6 inches


##############################################################################/
#Model Blink/kinship for limited inoculum condition####
##############################################################################/

#loading genotype data
LimG<-read.delim("lim.hmp.txt",header=FALSE)
#alternative kinship matrix
LimLois<-read.table("lim_lois_mat.txt",header=FALSE)

#Blink method on powdery mildew phenotype
limGAPIT<-GAPIT(
  Y=LimTrait[,1:7],
  G=LimG,
  kinship.algorithm="Loiselle",
  #KI=LimLois,
  PCA.total=0,
  model="Blink"
  #,Random.model=TRUE
)

#Blink method on general phenotype
limGAPIT<-GAPIT(
  Y=LimTrait[,c(1,7:10)],
  G=LimG,
  kinship.algorithm="Loiselle",
  #KI=LimLois,
  PCA.total=0,
  model="Blink"
  #,Random.model=TRUE
)


##############################################################################/
#Model Blink/kinship for limited inoculum condition####
##############################################################################/



##############################################################################/
#END
##############################################################################/