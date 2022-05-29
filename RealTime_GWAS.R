##############################################################################/
##############################################################################/
#GWAS analyses using GAPIT3
##############################################################################/
##############################################################################/

#at the time of writing this code, GAPIT3 was still under development on 
#Github
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)

#to be removed when the code will be working properly
setwd("/myGAPIT")


##############################################################################/
#first attempt without Q and kinship specific matrices fot natural condition###
##############################################################################/

ToutTrait<-read.table("pheno_final.txt",header=TRUE)
NatG<-read.delim("nat.hmp.txt",header=FALSE)
NatLois<-read.table("nat_lois_mat.txt", head = FALSE)

#Step 2: Run GAPIT
natGAPIT<-GAPIT(
  Y=ToutTrait,
  G=NatG,
  #KI=NatLois,
  PCA.total=2,
  model="Blink"
)


##############################################################################/
#first attempt without Q and kinship specific matrices fot natural condition###
##############################################################################/

ToutTrait<-read.table("pheno_final.txt",header=TRUE)
LimG<-read.delim("lim.hmp.txt",header=FALSE)
LimLois<-read.table("lim_lois_mat.txt", head = FALSE)

#Step 2: Run GAPIT
limGAPIT<-GAPIT(
  Y=ToutTrait,
  G=LimG,
  kinship.algorithm="Loiselle",
  Inter.Plot=TRUE,
  #KI=LimLois,
  PCA.total=0,
  model="Blink"
)


##############################################################################/
#END
##############################################################################/





#Tutorial 2: Using MLM 
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)
myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  G=myG,
  PCA.total=2,
  model="MLM"
)

#Tutorial 3: User defined Kinship and PCs
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)
myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)
myKI <- read.table("KSN.txt", head = FALSE)
myCV <- read.table("Copy of Q_First_Three_Principal_Components.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  G=myG,
  KI=myKI,
  CV=myCV,
)

#Tutorial 4: Genome Prediction
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)
myKI <- read.table("KSN.txt", head = FALSE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  G=myG,
  KI=myKI,
  PCA.total=3,
  model=c("gBLUP")
)

#Tutorial 5: Work with big data by spliting genotype Files
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  PCA.total=3,
  file.G="mdp_genotype_chr",
  file.Ext.G="hmp.txt",
  file.from=1,
  file.to=10
)

#Tutorial 6: Numeric Genotype Format
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)
myGD <- read.table("mdp_numeric.txt", head = TRUE)
myGM <- read.table("mdp_SNP_information.txt" , head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  GD=myGD,
  GM=myGM,
  PCA.total=3
)

#Tutorial 7: Numerical Multiple Files
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  PCA.total=3,
  file.GD="mdp_numeric",
  file.GM="mdp_SNP_information",
  file.Ext.GD="txt",
  file.Ext.GM="txt",
  file.from=1,
  file.to=3,
  
)


#Tutorial 8: Improve speed of computing PC and kinship by using fractional snps 
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  PCA.total=3,
  file.GD="mdp_numeric",
  file.GM="mdp_SNP_information",
  file.Ext.GD="txt",
  file.Ext.GM="txt",
  file.from=1,
  file.to=3,
  SNP.fraction=0.6
)

#Tutorial 9: Reduce memory usage by loading fragment of file, one at a time, by defini
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  PCA.total=3,
  file.GD="mdp_numeric",
  file.GM="mdp_SNP_information",
  file.Ext.GD="txt",
  file.Ext.GM="txt",
  file.from=1,
  file.to=3,
  SNP.fraction=0.6,
  file.fragment = 128
)

#Tutorial 10: Optimization for number of PCs based on BIC
#The result is saved in GAPIT.TraitName.BIC.Model.Selection.Results.csv
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)
myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  G=myG,
  PCA.total=3,
  Model.selection = TRUE
)

#Tutorial 11: SUPER GWAS method by Wang and et. al. (PLoS One, 2014)
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myCV <- read.table("Copy of Q_First_Three_Principal_Components.txt", head = TRUE)
myY  <- read.table("mdp_traits.txt", head = TRUE)
myG <- read.table("mdp_genotype_test.hmp.txt" , head = FALSE)

#Step 2: Run GAPIT
myGAPIT_SUPER <- GAPIT(
  Y=myY[,1:2],			
  G=myG,				
  CV=myCV,
  #PCA.total=3,				
  model="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1,
)


#Tutorial 12: Compare to Power against FDR for GLM,MLM,CMLM,ECMLM,SUPER(PLINK)
#Hint:Program runing time is more than 24 hours for repetition 100 times.
#Run description:Please refer to page 34,35 of the User manual on GAPIT website.
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myGD <-read.table("mdp_numeric.txt", head = TRUE)
myGM <-read.table("mdp_SNP_information.txt", head = TRUE)
myKI <- read.table("KSN.txt", head = FALSE)
myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)

#Step 2: Run GAPIT
#GAPIT.Power.compare.plink
GAPIT.Power.compare(
  GD=myGD,
  GM=myGM,
  WS=1000,
  rep=50,
  h2=0.75,
  Method=c("GLM","MLM","FarmCPU","Blink"),
  PCA.total=3,
  NQTN=20
)

#Tutorial 13: Genetic Prediction one time by cross validation
#-----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY<-read.table("mdp_traits.txt", head = TRUE)
myK<-read.table("KSN.txt", head = FALSE)

#Step 2: Run GAPIT
OnePre<-GAPIT.Prediction(
  myK=myK,
  y<-myY[,c(1,3)],
  ##y=y[,1:2],
  num=5
)
