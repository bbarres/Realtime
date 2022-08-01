##############################################################################/
##############################################################################/
#Plot related to the GWAS results
##############################################################################/
##############################################################################/

source("RealTime_load.R")
#in order to load the results of the GWAS analyses, you need to run the 
#RealTime_GWAS.R script before this script


##############################################################################/
#Loading the necessary data sets####
##############################################################################/

RezNatGAPIT<-read.table("data/RezNatGAPIT.txt",sep="\t",header=TRUE)
RezLimGAPIT<-read.table("data/RezLimGAPIT.txt",sep="\t",header=TRUE)


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