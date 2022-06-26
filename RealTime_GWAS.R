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
#Model Blink/kinship for natural inoculum condition####
##############################################################################/

AllTrait<-read.table("data/pheno_final.txt",header=TRUE)
#we remove the individuals without SNP data
AllTrait<-AllTrait[AllTrait$SNPage==1 & is.na(AllTrait$SNPage)!=TRUE &
                     AllTrait$Quality_SNPage==1,]
NatTrait<-AllTrait[AllTrait$treat=="exp",5:14]
LimTrait<-AllTrait[AllTrait$treat=="low",5:14]

#loading genotype data
NatG<-read.delim("data/nat.hmp.txt",header=FALSE)

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
  kinship.algorithm="Loiselle",
  #KI=NatLois,
  PCA.total=0,
  model="Blink"
  #,Random.model=TRUE
)
setwd(nomTemp)

# #MLM method on powdery mildew trait
# natGAPIT<-GAPIT(
#   Y=NatTrait[,1:7],
#   G=NatG,
#   kinship.algorithm="Loiselle",
#   #KI=NatLois,
#   PCA.total=0,
#   model="MLM"
# )

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
LimG<-read.delim("data/lim.hmp.txt",header=FALSE)

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
  kinship.algorithm="Loiselle",
  #KI=LimLois,
  PCA.total=0,
  model="Blink"
  #,Random.model=TRUE
)
setwd(nomTemp)


##############################################################################/
#SNP genotype effect on the traits####
##############################################################################/

temp<-snp.dat[snp.dat$Sample_ID %in% NatTrait$Taxa,]
temp<-temp[,c(1,48:866)]
temp2<-temp[,colnames(temp)=="Sample_ID" | 
              colnames(temp)=="Entomo_CL7647CT8535_01_89"]
temp3<-merge(temp2,NatTrait,by.x="Sample_ID",by.y="Taxa")
boxplot(temp3$oid_moy~temp3$Entomo_CL7647CT8535_01_89,boxwex=0.3,las=1)
stripchart(temp3$oid_moy~temp3$Entomo_CL7647CT8535_01_89,
           method="jitter",col="red",vertical=TRUE,pch=19,cex=0.7,add=TRUE)
vioplot(temp3$oid_moy~temp3$Entomo_CL7647CT8535_01_89,boxwex=0.3,las=1)

table(temp3$Entomo_CL7647CT8535_01_89)


boxplot(temp3$statut~temp3$Entomo_CL7647CT8535_01_89,boxwex=0.3,las=1)
stripchart(temp3$statut~temp3$Entomo_CL7647CT8535_01_89,
           method="jitter",col="red",vertical=TRUE,pch=19,cex=0.7,add=TRUE)


##############################################################################/
#END
##############################################################################/