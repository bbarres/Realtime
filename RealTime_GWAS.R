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
#Model Blink/kinship for natural inoculum condition####
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
setwd(nomTemp)

# #MLM method on powdery mildew trait
# natGAPIT<-GAPIT(
#   Y=NatTrait[,c(1,7:10],
#   G=NatG,
#   kinship.algorithm="Loiselle",
#   #KI=NatLois,
#   PCA.total=0,
#   model="MLM"
# )




op<-par(mar=c(6, 5, 1, 1) + 0.1)
graf<-barplot(kdrDistr$n,ylim=c(0,1000),
              ylab="Number of individuals",cex.axis =1.3,cex.lab=2,
              las=1,xaxt="n",yaxt="n",bty="n",
              col=rep(thecol,4),
              border=NA,
              space=c(rep(0.1,3),1.4,rep(0.1,2),1.4,
                      rep(0.1,2),1.4,rep(0.1,2)),
              font.lab=2)
abline(h=c(50,100,200,400,600,800,1000),col=grey(0.8,0.8),lwd=2,lty=1)
barplot(kdrDistr$n,ylim=c(0,1100),
        ylab="Number of individuals",cex.axis =1.3,cex.lab=2,
        las=1,xaxt="n",yaxt="n",bty="n",
        col=rep(thecol,4),
        border=NA,
        space=c(rep(0.1,3),1.4,rep(0.1,2),1.4,
                rep(0.1,2),1.4,rep(0.1,2)),
        font.lab=2,add=TRUE)
axis(1,at=graf[c(2,5,8,11)],labels=FALSE, lwd=4)
axis(2,at=c(50,100,200,400,600,800,1000),
     labels=c(50,100,200,400,600,800,1000),lwd=4,las=1,font=2,cex.axis=1.1)
box(bty="l",lwd=4)
text(graf,kdrDistr$n+15,
     labels=as.character(kdrDistr$n),font=2)
mtext(levels(kdrDistr$species)[c(3,2,4,1)],at=graf[c(2,5,8,11)],
      line=1.5,cex=1.4,side=1,font=2)
mtext("Species", at=9.4,line=3,cex=2,side=1,
      font=2,padj=1)
legend(12,800,
       legend=c("R/R","R/S","S/S"),
       pch=15,col=thecol,bg=thecol,bty="n",cex=1.3,pt.cex=1.6,xpd=TRUE,
       ncol=1,x.intersp=1,y.intersp=0.8)
par(op)


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
setwd(nomTemp)


##############################################################################/
#SNP genotype effect on the traits####
##############################################################################/

temp<-snp.dat[snp.dat$Sample_ID %in% NatTrait$Taxa,]
temp<-temp[,c(1,48:866)]
temp2<-temp[,colnames(temp)=="Sample_ID" | 
              colnames(temp)==RezNatGAPIT$SNP[1]]
temp3<-merge(temp2,NatTrait,by.x="Sample_ID",by.y="Taxa")
barplot(as.data.frame(table(temp3$`Dead or Alive`,
                            temp3$Entomo_CL7647CT8535_01_89))$Freq,
        col=brewer.pal(8,"Dark2")[2:1],las=1,
        space=c(0.1,rep(c(0.1,1.4),2),0.1))


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