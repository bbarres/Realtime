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
#Model Blink/kinship for natural exposure####
##############################################################################/

#preparing the output folder tree structure
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
RezNatGAPIT<-read.table("GAPIT.Association.Filter_GWAS_results.csv",
                        header=TRUE,sep=",")
RezNatGAPIT$SNP<-str_replace_all(RezNatGAPIT$SNP,"~","_")
RezNatGAPIT$SNP<-str_replace_all(RezNatGAPIT$SNP,"-","_")
RezNatGAPIT$SNP<-str_replace_all(RezNatGAPIT$SNP,fixed("+"),"_")
#one SNP is not significant considering Pcor<0.01 criteria, so we remove it
#we also remove the first column which is unnecessary
RezNatGAPIT<-RezNatGAPIT[-3,-1]
setwd(nomTemp)
write.table(RezNatGAPIT,file="data/RezNatGAPIT.txt",sep="\t",
            quote=FALSE,row.names=FALSE)


##############################################################################/
#Model Blink/kinship for protected exposure####
##############################################################################/

#preparing the output folder tree structure
nomTemp<-getwd()
dir.create(paste(nomTemp,"/output/proGWAS",sep=""))
setwd(paste(nomTemp,"/output/proGWAS",sep=""))
#Blink method on general phenotype
limGAPIT<-GAPIT(
  Y=LimTrait[,c(1,7:10)],
  G=LimG,
  SNP.effect="Add",
  PCA.total=0,
  model="Blink"
)
RezLimGAPIT<-read.table("GAPIT.Association.Filter_GWAS_results.csv",
                        header=TRUE,sep=",")
RezLimGAPIT$SNP<-str_replace_all(RezLimGAPIT$SNP,"~","_")
RezLimGAPIT$SNP<-str_replace_all(RezLimGAPIT$SNP,"-","_")
RezLimGAPIT$SNP<-str_replace_all(RezLimGAPIT$SNP,fixed("+"),"_")
#one SNP is not significant considering Pcor<0.01 criteria, so we remove it
RezLimGAPIT<-RezLimGAPIT[-5,-1]
setwd(nomTemp)
write.table(RezLimGAPIT,file="data/RezProGAPIT.txt",sep="\t",
            quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/