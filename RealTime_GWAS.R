##############################################################################/
##############################################################################/
#GWAS analyses using GAPIT3
##############################################################################/
##############################################################################/

# #at the time of writing this code, GAPIT3 was still under development on 
# #Github, so you have to run this to install the updated package
# devtools::install_github("jiabowang/GAPIT3",force=TRUE)
source("RealTime_load.R")


##############################################################################/
#Loading the necessary data sets####
##############################################################################/

AllTrait<-read.table("data/pheno_final.txt",header=TRUE)
colnames(AllTrait)[11:14]<-c("Powdery mildew","Height",
                             "Acorn weight","Dead or Alive")
#we remove the individuals without SNP data
AllTrait<-AllTrait[AllTrait$SNPage==1 & is.na(AllTrait$SNPage)!=TRUE &
                     AllTrait$Quality_SNPage==1,]
NatTrait<-AllTrait[AllTrait$treat=="exp",5:14]
LimTrait<-AllTrait[AllTrait$treat=="low",5:14]

#loading the natural treatment genotype data
NatG<-read.delim("data/nat.hmp.txt",header=FALSE)

#loading the limited treatment genotype data
LimG<-read.delim("data/lim.hmp.txt",header=FALSE)


##############################################################################/
#Model Blink/kinship for exposed condition####
##############################################################################/

# #Blink method on powdery mildew phenotype
# natGAPIT<-GAPIT(
#   Y=NatTrait[,1:7],
#   G=NatG,
#   kinship.algorithm="Loiselle",
#   #KI=NatLois,
#   PCA.total=0,
#   model="Blink"
#   #,Random.model=TRUE
# )

nomTemp<-getwd()
dir.create(paste(nomTemp,"/output/natGWAS",sep=""))
setwd(paste(nomTemp,"/output/natGWAS",sep=""))
#Blink method on general phenotype
natGAPIT<-GAPIT(
  Y=NatTrait[,c(1,7:10)],
  G=NatG,
  SNP.effect="Add",
  PCA.total=0,
  model="Blink"
)
RezNatGAPIT<-read.table("GAPIT.Filter_GWAS_results.txt",header=TRUE,
                        sep=" ")
RezNatGAPIT$SNP<-str_replace_all(RezNatGAPIT$SNP,"~","_")
RezNatGAPIT$SNP<-str_replace_all(RezNatGAPIT$SNP,"-","_")
RezNatGAPIT$SNP<-str_replace_all(RezNatGAPIT$SNP,fixed("+"),"_")
setwd(nomTemp)
write.table(RezNatGAPIT,file="data/RezNatGAPIT.txt",sep="\t",
            quote=FALSE,row.names=FALSE)

# #MLM method on powdery mildew trait
# natGAPIT<-GAPIT(
#   Y=NatTrait[,c(1,7:10],
#   G=NatG,
#   kinship.algorithm="Loiselle",
#   #KI=NatLois,
#   PCA.total=0,
#   model="MLM"
# )


##############################################################################/
#Model Blink/kinship for limited inoculum condition####
##############################################################################/

# #Blink method on powdery mildew phenotype
# limGAPIT<-GAPIT(
#   Y=LimTrait[,1:7],
#   G=LimG,
#   kinship.algorithm="Loiselle",
#   #KI=LimLois,
#   PCA.total=0,
#   model="Blink"
#   #,Random.model=TRUE
# )

nomTemp<-getwd()
dir.create(paste(nomTemp,"/output/limGWAS",sep=""))
setwd(paste(nomTemp,"/output/limGWAS",sep=""))
#Blink method on general phenotype
limGAPIT<-GAPIT(
  Y=LimTrait[,c(1,7:10)],
  G=LimG,
  SNP.effect="Add",
  PCA.total=0,
  model="Blink"
)
RezLimGAPIT<-read.table("GAPIT.Filter_GWAS_results.txt",header=TRUE,
                        sep=" ")
RezLimGAPIT$SNP<-str_replace_all(RezLimGAPIT$SNP,"~","_")
RezLimGAPIT$SNP<-str_replace_all(RezLimGAPIT$SNP,"-","_")
RezLimGAPIT$SNP<-str_replace_all(RezLimGAPIT$SNP,fixed("+"),"_")
setwd(nomTemp)
write.table(RezLimGAPIT,file="data/RezLimGAPIT.txt",sep="\t",
            quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/