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
#first attempt without Q and kinship specific matrices fot natural condition###
##############################################################################/

ToutTrait<-read.table("pheno_final.txt",header=TRUE)
NatG<-read.delim("nat.hmp.txt",header=FALSE)
NatLois<-read.table("nat_lois_mat.txt", head = FALSE)

#Blink method on oidium phenotype
natGAPIT<-GAPIT(
  Y=ToutTrait[,1:7],
  G=NatG,
  kinship.algorithm="Loiselle",
  #KI=NatLois,
  PCA.total=0,
  model="Blink"
  #,Random.model=TRUE
)

#Blink method on general phenotype
natGAPIT<-GAPIT(
  Y=ToutTrait[,c(1,7:10)],
  G=NatG,
  kinship.algorithm="Loiselle",
  #KI=NatLois,
  PCA.total=0,
  model="Blink"
  #,Random.model=TRUE
)

#MLM method on oidium trait
natGAPIT<-GAPIT(
  Y=ToutTrait[,1:7],
  G=NatG,
  kinship.algorithm="Loiselle",
  #KI=NatLois,
  PCA.total=0,
  model="MLM"
)

##############################################################################/
#first attempt without Q and kinship specific matrices fot natural condition###
##############################################################################/

ToutTrait<-read.table("pheno_final.txt",header=TRUE)
LimG<-read.delim("lim.hmp.txt",header=FALSE)
LimLois<-read.table("lim_lois_mat.txt", head = FALSE)

#Blink method on oidium phenotype
limGAPIT<-GAPIT(
  Y=ToutTrait[,1:7],
  G=LimG,
  kinship.algorithm="Loiselle",
  #KI=LimLois,
  PCA.total=0,
  model="Blink"
  #,Random.model=TRUE
)

#Blink method on general phenotype
limGAPIT<-GAPIT(
  Y=ToutTrait[,c(1,7:10)],
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