##############################################################################/
##############################################################################/
#GWAS analyses using GAPIT3
##############################################################################/
##############################################################################/

# #at the time of writing this code, GAPIT3 was still under development on 
# #Github, so you have to run this to install the updated package
# devtools::install_github("jiabowang/GAPIT3",force=TRUE)

source("RealTime_load.R")

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
#END
##############################################################################/