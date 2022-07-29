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
RezNatGAPIT$SNP<-str_replace_all(RezNatGAPIT$SNP,fixed("+"),"_")
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
axis(1,at=graf[c(2,5,8,11)],labels=FALSE,lwd=4)
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
RezLimGAPIT$SNP<-str_replace_all(RezLimGAPIT$SNP,"~","_")
RezLimGAPIT$SNP<-str_replace_all(RezLimGAPIT$SNP,"-","_")
RezLimGAPIT$SNP<-str_replace_all(RezLimGAPIT$SNP,fixed("+"),"_")
setwd(nomTemp)


##############################################################################/
#plotting significant SNP genotype effect on the traits####
##############################################################################/

#code for plotting all the significant SNP
temp<-snp.dat[snp.dat$Sample_ID %in% NatTrait$Taxa,]
temp<-temp[,c(1,48:866)]
tempb<-snp.dat[snp.dat$Sample_ID %in% LimTrait$Taxa,]
tempb<-tempb[,c(1,48:866)]

#for the Exposed treatment
for (i in 1:dim(RezNatGAPIT)[1]) {
  #data of the exposed treatment
  temp2<-temp[,colnames(temp)=="Sample_ID" | 
                colnames(temp)==RezNatGAPIT$SNP[i]]
  temp2<-merge(temp2,NatTrait,by.x="Sample_ID",by.y="Taxa")
  #temp2[,2]<-as.factor(temp2[,2])
  traiinterest<-substring(RezNatGAPIT[i,6],7)
  #data of the protected treatment
  temp3<-tempb[,colnames(tempb)=="Sample_ID" | 
                 colnames(tempb)==RezNatGAPIT$SNP[i]]
  temp3<-merge(temp3,LimTrait,by.x="Sample_ID",by.y="Taxa")
  #temp3[,2]<-as.factor(temp3[,2])
  
  pdf(file=paste("output/SNP_EXP_",traiinterest,"_",
                 RezNatGAPIT$SNP[i],".pdf",sep=""),
      width=12,height=4)
  op<-par(mfrow=c(1,4))
  vioplot(temp2$`Acorn weight`~temp2[,2],boxwex=0.3,las=1,
          col=brewer.pal(9,"Set1")[c(6,3,2)],
          xlab=RezNatGAPIT$SNP[i],
          ylab="Acorn weight (g)",main="Acorn weight")
  vioplot(temp2$`Height`~temp2[,2],boxwex=0.3,las=1,
          col=brewer.pal(9,"Set1")[c(6,3,2)],
          xlab=RezNatGAPIT$SNP[i],
          ylab="Height (cm)",main="Height")
  vioplot(temp2$`Powdery mildew`~temp2[,2],boxwex=0.3,las=1,
          col=brewer.pal(9,"Set1")[c(6,3,2)],
          xlab=RezNatGAPIT$SNP[i],
          ylab="Powdery mildew note",main="Powdery mildew")
  graf<-barplot(as.data.frame(table(temp2$`Dead or Alive`,
                                    temp2[,2]))$Freq,
                col=brewer.pal(12,"Set3")[6:7],las=1,
                space=c(0.1,rep(c(0.1,0.9),2),0.1),
                main="Dead or Alive",xlab=RezNatGAPIT$SNP[i])
  #abline(h=c(50,100,200,300),col=grey(0.8,0.8),lwd=2,lty=1)
  axis(1,at=(graf[c(1,3,5)]+graf[c(2,4,6)])/2,
       labels=colnames(table(temp2$`Dead or Alive`,
                             temp2[,2])))
  box(bty="o",lwd=1)
  par(op)
  dev.off()
  
  #second plot for the two treatment
  pdf(file=paste("output/SNP_EXP_DoA_",traiinterest,"_",
                 RezNatGAPIT$SNP[i],".pdf",sep=""),
      width=7,height=8)
  op<-par(mfcol=c(2,2))
  vioplot(temp2[,traiinterest]~temp2[,2],
          col=c("transparent"),sep=":",las=1,border="transparent",
          ylab=traiinterest,xlab="",frame.plot=FALSE,
          lineCol="transparent",rectCol="transparent",
          main=paste(RezNatGAPIT$SNP[i],"exposed"))
  vioplot(temp2[temp2[,11]==0,traiinterest]~temp2[temp2[,11]==0,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[6],sep=":",las=1,side="left",
          frame.plot=FALSE,add=TRUE)
  vioplot(temp2[temp2[,11]==1,traiinterest]~temp2[temp2[,11]==1,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[7],sep=":",las=1,side="right",
          frame.plot=FALSE,add=TRUE)
  box(bty="o",lwd=1)
  graf<-barplot(proportions(table(temp2$`Dead or Alive`,
                                  temp2[,2]),margin=2)[1,]*100,
                col=brewer.pal(12,"Set3")[6],las=1,space=1,
                ylim=c(0,100),names.arg="",
                main="Exposed death rate")
  axis(1,at=graf,
       labels=paste(colnames(table(temp2$`Dead or Alive`,temp2[,2])),
                    " (n=",table(temp2[,2]),")",
                    sep=""))
  box(bty="o",lwd=1)
  vioplot(temp3[,traiinterest]~temp3[,2],
          col=c("transparent"),sep=":",las=1,border="transparent",
          ylab=traiinterest,xlab="",frame.plot=FALSE,
          lineCol="transparent",rectCol="transparent",
          main=paste(RezNatGAPIT$SNP[i],"protected"))
  vioplot(temp3[temp3[,11]==0,traiinterest]~temp3[temp3[,11]==0,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[6],sep=":",las=1,side="left",
          frame.plot=FALSE,add=TRUE)
  vioplot(temp3[temp3[,11]==1,traiinterest]~temp3[temp3[,11]==1,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[7],sep=":",las=1,side="right",
          frame.plot=FALSE,add=TRUE)
  box(bty="o",lwd=1)
  grafb<-barplot(proportions(table(temp3$`Dead or Alive`,
                                   temp3[,2]),margin=2)[1,]*100,
                 col=brewer.pal(12,"Set3")[6],las=1,space=1,
                 ylim=c(0,100),names.arg="",
                 main="Protected death rate")
  axis(1,at=grafb,
       labels=paste(colnames(table(temp3$`Dead or Alive`,temp3[,2])),
                    " (n=",table(temp3[,2]),")",
                    sep=""))
  box(bty="o",lwd=1)
  par(op)
  dev.off()
  
}

#for the protected treatment
for (i in 1:dim(RezLimGAPIT)[1]) {
  #data of the exposed treatment
  temp2<-temp[,colnames(temp)=="Sample_ID" | 
                colnames(temp)==RezLimGAPIT$SNP[i]]
  temp2<-merge(temp2,NatTrait,by.x="Sample_ID",by.y="Taxa")
  #temp2[,2]<-as.factor(temp2[,2])
  traiinterest<-substring(RezLimGAPIT[i,6],7)
  #data of the protected treatment
  temp3<-tempb[,colnames(tempb)=="Sample_ID" | 
                 colnames(tempb)==RezLimGAPIT$SNP[i]]
  temp3<-merge(temp3,LimTrait,by.x="Sample_ID",by.y="Taxa")
  #temp3[,2]<-as.factor(temp3[,2])
  
  pdf(file=paste("output/SNP_PRO_",traiinterest,"_",
                 RezLimGAPIT$SNP[i],".pdf",sep=""),
      width=12,height=4)
  op<-par(mfrow=c(1,4))
  vioplot(temp3$`Acorn weight`~temp3[,2],boxwex=0.3,las=1,
          col=brewer.pal(9,"Set1")[c(6,3,2)],
          xlab=RezLimGAPIT$SNP[i],
          ylab="Acorn weight (g)",main="Acorn weight")
  vioplot(temp3$`Height`~temp3[,2],boxwex=0.3,las=1,
          col=brewer.pal(9,"Set1")[c(6,3,2)],
          xlab=RezLimGAPIT$SNP[i],
          ylab="Height (cm)",main="Height")
  vioplot(temp3$`Powdery mildew`~temp3[,2],boxwex=0.3,las=1,
          col=brewer.pal(9,"Set1")[c(6,3,2)],
          xlab=RezLimGAPIT$SNP[i],
          ylab="Powdery mildew note",main="Powdery mildew")
  graf<-barplot(as.data.frame(table(temp3$`Dead or Alive`,
                                    temp3[,2]))$Freq,
                col=brewer.pal(12,"Set3")[6:7],las=1,
                space=c(0.1,rep(c(0.1,0.9),2),0.1),
                main="Dead or Alive",xlab=RezLimGAPIT$SNP[i])
  #abline(h=c(50,100,200,300),col=grey(0.8,0.8),lwd=2,lty=1)
  axis(1,at=(graf[c(1,3,5)]+graf[c(2,4,6)])/2,
       labels=colnames(table(temp3$`Dead or Alive`,
                             temp3[,2])))
  box(bty="o",lwd=1)
  par(op)
  dev.off()
  
  #second plot for the two treatment
  pdf(file=paste("output/SNP_PRO_DoA_",traiinterest,"_",
                 RezLimGAPIT$SNP[i],".pdf",sep=""),
      width=7,height=8)
  op<-par(mfcol=c(2,2))
  vioplot(temp2[,traiinterest]~temp2[,2],
          col=c("transparent"),sep=":",las=1,border="transparent",
          ylab=traiinterest,xlab="",frame.plot=FALSE,
          lineCol="transparent",rectCol="transparent",
          main=paste(RezLimGAPIT$SNP[i],"exposed"))
  vioplot(temp2[temp2[,11]==0,traiinterest]~temp2[temp2[,11]==0,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[6],sep=":",las=1,side="left",
          frame.plot=FALSE,add=TRUE)
  vioplot(temp2[temp2[,11]==1,traiinterest]~temp2[temp2[,11]==1,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[7],sep=":",las=1,side="right",
          frame.plot=FALSE,add=TRUE)
  box(bty="o",lwd=1)
  graf<-barplot(proportions(table(temp2$`Dead or Alive`,
                                  temp2[,2]),margin=2)[1,]*100,
                col=brewer.pal(12,"Set3")[6],las=1,space=1,
                ylim=c(0,100),names.arg="",
                main="Exposed death rate")
  axis(1,at=graf,
       labels=paste(colnames(table(temp2$`Dead or Alive`,temp2[,2])),
                    " (n=",table(temp2[,2]),")",
                    sep=""))
  box(bty="o",lwd=1)
  vioplot(temp3[,traiinterest]~temp3[,2],
          col=c("transparent"),sep=":",las=1,border="transparent",
          ylab=traiinterest,xlab="",frame.plot=FALSE,
          lineCol="transparent",rectCol="transparent",
          main=paste(RezNatGAPIT$SNP[i],"protected"))
  vioplot(temp3[temp3[,11]==0,traiinterest]~temp3[temp3[,11]==0,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[6],sep=":",las=1,side="left",
          frame.plot=FALSE,add=TRUE)
  vioplot(temp3[temp3[,11]==1,traiinterest]~temp3[temp3[,11]==1,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[7],sep=":",las=1,side="right",
          frame.plot=FALSE,add=TRUE)
  box(bty="o",lwd=1)
  grafb<-barplot(proportions(table(temp3$`Dead or Alive`,
                                   temp3[,2]),margin=2)[1,]*100,
                 col=brewer.pal(12,"Set3")[6],las=1,space=1,
                 ylim=c(0,100),names.arg="",
                 main="Protected death rate")
  axis(1,at=grafb,
       labels=paste(colnames(table(temp3$`Dead or Alive`,temp3[,2])),
                    " (n=",table(temp3[,2]),")",
                    sep=""))
  box(bty="o",lwd=1)
  par(op)
  dev.off()
  
}


##############################################################################/
#plot for the SNP significant for Dead or Alive trait####
##############################################################################/

#code for plotting all the significant SNP
temp<-snp.dat[snp.dat$Sample_ID %in% NatTrait$Taxa,]
temp<-temp[,c(1,48:866)]
tempb<-snp.dat[snp.dat$Sample_ID %in% LimTrait$Taxa,]
tempb<-tempb[,c(1,48:866)]

#for the Exposed treatment
temp2<-temp[,colnames(temp)=="Sample_ID" | 
              colnames(temp)==RezNatGAPIT$SNP[14]]
temp2<-merge(temp2,NatTrait,by.x="Sample_ID",by.y="Taxa")
#temp2[,2]<-as.factor(temp2[,2])
traiinterest<-substring(RezNatGAPIT[14,6],7)
#data of the protected treatment
temp3<-tempb[,colnames(tempb)=="Sample_ID" | 
               colnames(tempb)==RezNatGAPIT$SNP[14]]
temp3<-merge(temp3,LimTrait,by.x="Sample_ID",by.y="Taxa")
#temp3[,2]<-as.factor(temp3[,2])

pdf(file=paste("output/DoA_SNP",traiinterest,"_",
               RezNatGAPIT$SNP[14],".pdf",sep=""),
    width=10,height=7)
op<-par(mfrow=c(2,3))
vioplot(temp2$`Acorn weight`~temp2[,2],boxwex=0.3,las=1,
        col=brewer.pal(9,"Set1")[c(6,3,2)],
        xlab=RezNatGAPIT$SNP[14],
        ylab="Acorn weight (g)",main="Acorn weight")
vioplot(temp2$`Height`~temp2[,2],boxwex=0.3,las=1,
        col=brewer.pal(9,"Set1")[c(6,3,2)],
        xlab=RezNatGAPIT$SNP[14],
        ylab="Height (cm)",main="Height")
vioplot(temp2$`Powdery mildew`~temp2[,2],boxwex=0.3,las=1,
        col=brewer.pal(9,"Set1")[c(6,3,2)],
        xlab=RezNatGAPIT$SNP[14],
        ylab="Powdery mildew note",main="Powdery mildew")
graf<-barplot(as.data.frame(table(temp2$`Dead or Alive`,
                                  temp2[,2]))$Freq,
              col=brewer.pal(12,"Set3")[6:7],las=1,
              space=c(0.1,rep(c(0.1,0.9),2),0.1),
              main="Dead or Alive",xlab=RezNatGAPIT$SNP[14])
#abline(h=c(50,100,200,300),col=grey(0.8,0.8),lwd=2,lty=1)
axis(1,at=(graf[c(1,3,5)]+graf[c(2,4,6)])/2,
     labels=colnames(table(temp2$`Dead or Alive`,
                           temp2[,2])))
box(bty="o",lwd=1)


graf<-barplot(proportions(table(temp2$`Dead or Alive`,
                                temp2[,2]),margin=2)[1,]*100,
              col=brewer.pal(12,"Set3")[6],las=1,space=1,
              ylim=c(0,100),names.arg="",
              main="Exposed death rate")
axis(1,at=graf,
     labels=paste(colnames(table(temp2$`Dead or Alive`,temp2[,2])),
                  " (n=",table(temp2[,2]),")",
                  sep=""))
box(bty="o",lwd=1)
grafb<-barplot(proportions(table(temp3$`Dead or Alive`,
                                 temp3[,2]),margin=2)[1,]*100,
               col=brewer.pal(12,"Set3")[6],las=1,space=1,
               ylim=c(0,100),names.arg="",
               main="Protected death rate")
axis(1,at=grafb,
     labels=paste(colnames(table(temp3$`Dead or Alive`,temp3[,2])),
                  " (n=",table(temp3[,2]),")",
                  sep=""))
box(bty="o",lwd=1)
par(op)
dev.off()


  
  
  

#for the protected treatment
  #data of the exposed treatment
  temp2<-temp[,colnames(temp)=="Sample_ID" | 
                colnames(temp)==RezLimGAPIT$SNP[6]]
  temp2<-merge(temp2,NatTrait,by.x="Sample_ID",by.y="Taxa")
  #temp2[,2]<-as.factor(temp2[,2])
  traiinterest<-substring(RezLimGAPIT[6,6],7)
  #data of the protected treatment
  temp3<-tempb[,colnames(tempb)=="Sample_ID" | 
                 colnames(tempb)==RezLimGAPIT$SNP[6]]
  temp3<-merge(temp3,LimTrait,by.x="Sample_ID",by.y="Taxa")
  #temp3[,2]<-as.factor(temp3[,2])
  
  pdf(file=paste("output/SNP_PRO_",traiinterest,"_",
                 RezLimGAPIT$SNP[6],".pdf",sep=""),
      width=12,height=4)
  op<-par(mfrow=c(1,4))
  vioplot(temp3$`Acorn weight`~temp3[,2],boxwex=0.3,las=1,
          col=brewer.pal(9,"Set1")[c(6,3,2)],
          xlab=RezLimGAPIT$SNP[6],
          ylab="Acorn weight (g)",main="Acorn weight")
  vioplot(temp3$`Height`~temp3[,2],boxwex=0.3,las=1,
          col=brewer.pal(9,"Set1")[c(6,3,2)],
          xlab=RezLimGAPIT$SNP[6],
          ylab="Height (cm)",main="Height")
  vioplot(temp3$`Powdery mildew`~temp3[,2],boxwex=0.3,las=1,
          col=brewer.pal(9,"Set1")[c(6,3,2)],
          xlab=RezLimGAPIT$SNP[6],
          ylab="Powdery mildew note",main="Powdery mildew")
  graf<-barplot(as.data.frame(table(temp3$`Dead or Alive`,
                                    temp3[,2]))$Freq,
                col=brewer.pal(12,"Set3")[6:7],las=1,
                space=c(0.1,rep(c(0.1,0.9),2),0.1),
                main="Dead or Alive",xlab=RezLimGAPIT$SNP[6])
  #abline(h=c(50,100,200,300),col=grey(0.8,0.8),lwd=2,lty=1)
  axis(1,at=(graf[c(1,3,5)]+graf[c(2,4,6)])/2,
       labels=colnames(table(temp3$`Dead or Alive`,
                             temp3[,2])))
  box(bty="o",lwd=1)
  par(op)
  dev.off()
  
  #second plot for the two treatment
  op<-par(mfcol=c(2,2))
  vioplot(temp2[,traiinterest]~temp2[,2],
          col=c("transparent"),sep=":",las=1,border="transparent",
          ylab=traiinterest,xlab="",frame.plot=FALSE,
          lineCol="transparent",rectCol="transparent",
          main=paste(RezLimGAPIT$SNP[6],"exposed"))
  vioplot(temp2[temp2[,11]==0,traiinterest]~temp2[temp2[,11]==0,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[6],sep=":",las=1,side="left",
          frame.plot=FALSE,add=TRUE)
  vioplot(temp2[temp2[,11]==1,traiinterest]~temp2[temp2[,11]==1,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[7],sep=":",las=1,side="right",
          frame.plot=FALSE,add=TRUE)
  box(bty="o",lwd=1)
  graf<-barplot(proportions(table(temp2$`Dead or Alive`,
                                  temp2[,2]),margin=2)[1,]*100,
                col=brewer.pal(12,"Set3")[6],las=1,space=1,
                ylim=c(0,100),names.arg="",
                main="Exposed death rate")
  axis(1,at=graf,
       labels=paste(colnames(table(temp2$`Dead or Alive`,temp2[,2])),
                    " (n=",table(temp2[,2]),")",
                    sep=""))
  box(bty="o",lwd=1)
  vioplot(temp3[,traiinterest]~temp3[,2],
          col=c("transparent"),sep=":",las=1,border="transparent",
          ylab=traiinterest,xlab="",frame.plot=FALSE,
          lineCol="transparent",rectCol="transparent",
          main=paste(RezNatGAPIT$SNP[6],"protected"))
  vioplot(temp3[temp3[,11]==0,traiinterest]~temp3[temp3[,11]==0,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[6],sep=":",las=1,side="left",
          frame.plot=FALSE,add=TRUE)
  vioplot(temp3[temp3[,11]==1,traiinterest]~temp3[temp3[,11]==1,2],
          plotCentre="line",
          col=brewer.pal(12,"Set3")[7],sep=":",las=1,side="right",
          frame.plot=FALSE,add=TRUE)
  box(bty="o",lwd=1)
  grafb<-barplot(proportions(table(temp3$`Dead or Alive`,
                                   temp3[,2]),margin=2)[1,]*100,
                 col=brewer.pal(12,"Set3")[6],las=1,space=1,
                 ylim=c(0,100),names.arg="",
                 main="Protected death rate")
  axis(1,at=grafb,
       labels=paste(colnames(table(temp3$`Dead or Alive`,temp3[,2])),
                    " (n=",table(temp3[,2]),")",
                    sep=""))
  box(bty="o",lwd=1)
  par(op)
  dev.off()


##############################################################################/
#END
##############################################################################/