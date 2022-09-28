##############################################################################/
##############################################################################/
#Plot related to the GWAS results
##############################################################################/
##############################################################################/

source("RealTime_load.R")


##############################################################################/
#Loading the necessary data sets####
##############################################################################/

#in order to load the results of the GWAS analyses, you need to run the 
#RealTime_GWAS.R script before this script
RezNatGAPIT<-read.table("data/RezNatGAPIT.txt",sep="\t",header=TRUE)
RezNatGAPIT$treat<-"Natural"
RezLimGAPIT<-read.table("data/RezLimGAPIT.txt",sep="\t",header=TRUE)
RezLimGAPIT$treat<-"Protected"
Rezcomb<-rbind(RezNatGAPIT,RezLimGAPIT)
#defining colors vector
colnat<-brewer.pal(11,"PRGn")[7:10]
colprt<-brewer.pal(11,"RdBu")[7:10]


##############################################################################/
#plotting significant SNP genotype effect on ALL the traits####
##############################################################################/

#code for plotting all the significant SNP
temp<-snp.dat[snp.dat$Sample_ID %in% NatTrait$Taxa,]
temp<-temp[,c(1,48:866)]
tempb<-snp.dat[snp.dat$Sample_ID %in% LimTrait$Taxa,]
tempb<-tempb[,c(1,48:866)]

pdf(file=paste("output/SNP_interest",".pdf",sep=""),
    width=15,height=10)

nf<-layout(matrix(c(1,2,2,2,3,3,3,4,4,4,5,5,5,
                    6,9,9,9,10,10,10,11,11,11,12,12,12,
                    6,9,9,9,10,10,10,11,11,11,12,12,12,
                    6,9,9,9,10,10,10,11,11,11,12,12,12,
                    7,13,13,13,14,14,14,15,15,15,16,16,16,
                    7,13,13,13,14,14,14,15,15,15,16,16,16,
                    7,13,13,13,14,14,14,15,15,15,16,16,16,
                    8,17,17,17,18,18,18,19,19,19,20,20,20,
                    8,17,17,17,18,18,18,19,19,19,20,20,20,
                    8,17,17,17,18,18,18,19,19,19,20,20,20),
                  10,13,byrow=TRUE))
op<-par(mar=c(0,0,0,0),lwd=2)
plot.new()
plot.new()
text(0.5,0.5,"Acorn weight (g)",font=2,cex=2.5)
plot.new()
text(0.5,0.5,"Height (cm)",font=2,cex=2.5)
plot.new()
text(0.5,0.5,"Powdery mildew\ninfection",font=2,cex=2.5)
plot.new()
text(0.5,0.7,"Mortality rate (%)",font=2,cex=2.5)
legend(0.25,0.55,legend=c("Natural treatment","Protected treatment"),
       fill=c(colnat[3],colprt[3]),border=c(colnat[3],colprt[3]),
       density=c(-1,20),angle=c(0,60),bty="n",
       cex=1.5,y.intersp=1.3,x.intersp=0.5,ncol=1,xpd=TRUE)
plot.new()
text(0.5,0.5,Rezcomb[1,1],srt=90,font=2,cex=1.5)
plot.new()
text(0.5,0.5,Rezcomb[11,1],srt=90,font=2,cex=1.5)
plot.new()
text(0.5,0.5,Rezcomb[16,1],srt=90,font=2,cex=1.5)

par(mar=c(3,3,1,1))
for (i in c(1,11,16)) {
  #data of the exposed treatment
  temp2<-temp[,colnames(temp)=="Sample_ID" | 
                colnames(temp)==Rezcomb$SNP[i]]
  temp2<-merge(temp2,NatTrait,by.x="Sample_ID",by.y="Taxa")
  #temp2[,2]<-as.factor(temp2[,2])
  traiinterest<-substring(Rezcomb[i,6],7)
  treatme<-Rezcomb[i,7]
  #data of the protected treatment
  temp3<-tempb[,colnames(tempb)=="Sample_ID" | 
                 colnames(tempb)==Rezcomb$SNP[i]]
  temp3<-merge(temp3,LimTrait,by.x="Sample_ID",by.y="Taxa")
  #temp3[,2]<-as.factor(temp3[,2])
  if(treatme=="Natural") {
    tempo<-temp2
    colo<-colnat
  } else {
    tempo<-temp3
    colo<-colprt
  }
  interme<-proportions(table(temp2$`Survival`,
                             temp2[,2]),margin=2)[1,]*100
  interme<-as.data.frame(interme)
  intermeb<-proportions(table(temp3$`Survival`,
                              temp3[,2]),margin=2)[1,]*100
  intermeb<-as.data.frame(intermeb)
  interme<-cbind(interme,intermeb)
  colnames(interme)<-c("Natural","Protected")
  interme$nNatur<-table(temp2[,2])
  interme$nProte<-table(temp3[,2])
  propsurv<-pivot_longer(interme[,c(1:2)],cols=1:2,
                         names_to="treat",values_to="Mortality Rate")
  totaleff<-pivot_longer(interme[,c(3:4)],cols=1:2,
                         names_to="treat",values_to="Total number")
  
  vioplot(tempo$`Acorn weight`~tempo[,2],boxwex=0.3,las=1,
          border=colo[2:4],lineCol=colo[2:4],rectCol=colo[2:4],
          col=colo[1],ann=FALSE,ylim=c(0,13.5),
          xlab=Rezcomb$SNP[i],lty=1,lwd=3,axes=0,yaxt="n",
          ylab="Acorn weight (g)",main="Acorn weight")
  axis(1,labels=paste(colnames(table(tempo$`Acorn weight`,tempo[,2])),
                      " (n=",table(tempo[,2]),")",
                      sep=""),at=c(1,2,3),cex.axis=1.3,lwd=2)
  axis(2,at=seq(from=1,to=13,by=2),las=1,cex.axis=1.3,lwd=2)
  box(bty="o")
  vioplot(tempo$`Height`~tempo[,2],boxwex=0.3,las=1,
          border=colo[2:4],lineCol=colo[2:4],rectCol=colo[2:4],
          col=colo[1],ann=FALSE,ylim=c(0,110),
          xlab=Rezcomb$SNP[i],lty=1,lwd=3,axes=0,yaxt="n",
          ylab="Height (cm)",main="Height")
  axis(1,labels=paste(colnames(table(tempo$`Height`,tempo[,2])),
                      " (n=",table(tempo[,2]),")",
                      sep=""),at=c(1,2,3),cex.axis=1.3,lwd=2)
  axis(2,at=seq(from=5,to=115,by=15),las=1,cex.axis=1.3,lwd=2)
  box(bty="o")
  vioplot(tempo$`Powdery mildew`~tempo[,2],boxwex=0.3,las=1,
          border=colo[2:4],lineCol=colo[2:4],rectCol=colo[2:4],
          col=colo[1],ann=FALSE,ylim=c(0,70),
          xlab=Rezcomb$SNP[i],lty=1,lwd=3,axes=0,yaxt="n",
          ylab="Powdery mildew note",main="Powdery mildew")
  axis(1,labels=paste(colnames(table(tempo$`Powdery mildew`,tempo[,2])),
                      " (n=",table(tempo[,2]),")",
                      sep=""),at=c(1,2,3),cex.axis=1.3,lwd=2)
  axis(2,at=seq(from=0,to=70,by=10),las=1,cex.axis=1.3,lwd=2)
  box(bty="o")
  graf<-barplot(propsurv$`Mortality Rate`,density=c(-1,20),angle=c(0,60),
                col=c(colnat[2],colprt[2],colnat[3],colprt[3],
                         colnat[4],colprt[4]),
                border=c(colnat[2],colprt[2],colnat[3],colprt[3],
                      colnat[4],colprt[4]),
                ylim=c(-1,110),ann=FALSE,
                ylab="Mortality Rate",lwd=2,
                cex.axis=1,cex.lab=1.5,las=1,yaxt="n",bty="n",
                space=c(0,0.1,rep(c(0.8,0.1),times=2)),font.lab=2)
  text(labels=paste("n=",totaleff$`Total number`,sep=""),
       x=graf,y=propsurv$`Mortality Rate`+5,cex=1.1,font=3,srt=0)
  axis(1,at=graf[seq(1,6,2)]+0.55,cex.axis=1.3,lwd=2,
       colnames(table(tempo$`Survival`,tempo[,2])))
  axis(2,at=seq(from=0,to=100,by=20),las=1,cex.axis=1.3,lwd=2)
  box(bty="o")

}

par(op)
#export to .pdf 15 x 10 inches  
dev.off()


##############################################################################/
#plotting significant SNP genotype effect on ALL the traits (old version)####
##############################################################################/

#for the SNP significant in the Exposed treatment
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
  graf<-barplot(as.data.frame(table(temp2$`Survival`,
                                    temp2[,2]))$Freq,
                col=brewer.pal(12,"Set3")[6:7],las=1,
                space=c(0.1,rep(c(0.1,0.9),2),0.1),
                main="Survival",xlab=RezNatGAPIT$SNP[i])
  #abline(h=c(50,100,200,300),col=grey(0.8,0.8),lwd=2,lty=1)
  axis(1,at=(graf[c(1,3,5)]+graf[c(2,4,6)])/2,
       labels=colnames(table(temp2$`Survival`,
                             temp2[,2])))
  box(bty="o",lwd=1)
  par(op)
  dev.off()
  
  #second plot for the distribution in the two treatments
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
  graf<-barplot(proportions(table(temp2$`Survival`,
                                  temp2[,2]),margin=2)[1,]*100,
                col=brewer.pal(12,"Set3")[6],las=1,space=1,
                ylim=c(0,100),names.arg="",
                main="Exposed death rate")
  axis(1,at=graf,
       labels=paste(colnames(table(temp2$`Survival`,temp2[,2])),
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
  grafb<-barplot(proportions(table(temp3$`Survival`,
                                   temp3[,2]),margin=2)[1,]*100,
                 col=brewer.pal(12,"Set3")[6],las=1,space=1,
                 ylim=c(0,100),names.arg="",
                 main="Protected death rate")
  axis(1,at=grafb,
       labels=paste(colnames(table(temp3$`Survival`,temp3[,2])),
                    " (n=",table(temp3[,2]),")",
                    sep=""))
  box(bty="o",lwd=1)
  par(op)
  dev.off()
  
}


#for the SNP significant for the protected treatment
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
  graf<-barplot(as.data.frame(table(temp3$`Survival`,
                                    temp3[,2]))$Freq,
                col=brewer.pal(12,"Set3")[6:7],las=1,
                space=c(0.1,rep(c(0.1,0.9),2),0.1),
                main="Survival",xlab=RezLimGAPIT$SNP[i])
  #abline(h=c(50,100,200,300),col=grey(0.8,0.8),lwd=2,lty=1)
  axis(1,at=(graf[c(1,3,5)]+graf[c(2,4,6)])/2,
       labels=colnames(table(temp3$`Survival`,
                             temp3[,2])))
  box(bty="o",lwd=1)
  par(op)
  dev.off()
  
  #second plot for the distribution in the two treatments
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
  graf<-barplot(proportions(table(temp2$`Survival`,
                                  temp2[,2]),margin=2)[1,]*100,
                col=brewer.pal(12,"Set3")[6],las=1,space=1,
                ylim=c(0,100),names.arg="",
                main="Exposed death rate")
  axis(1,at=graf,
       labels=paste(colnames(table(temp2$`Survival`,temp2[,2])),
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
  grafb<-barplot(proportions(table(temp3$`Survival`,
                                   temp3[,2]),margin=2)[1,]*100,
                 col=brewer.pal(12,"Set3")[6],las=1,space=1,
                 ylim=c(0,100),names.arg="",
                 main="Protected death rate")
  axis(1,at=grafb,
       labels=paste(colnames(table(temp3$`Survival`,temp3[,2])),
                    " (n=",table(temp3[,2]),")",
                    sep=""))
  box(bty="o",lwd=1)
  par(op)
  dev.off()
  
}


##############################################################################/
#plot for the SNP significant for Survival and Powdery mildew trait####
##############################################################################/

#code for plotting all the significant SNP
temp<-snp.dat[snp.dat$Sample_ID %in% NatTrait$Taxa,]
temp<-temp[,c(1,48:866)]
tempb<-snp.dat[snp.dat$Sample_ID %in% LimTrait$Taxa,]
tempb<-tempb[,c(1,48:866)]

#for the Exposed treatment
temp2<-temp[,colnames(temp)=="Sample_ID" | 
              colnames(temp)==RezNatGAPIT$SNP[12]]
temp2<-merge(temp2,NatTrait,by.x="Sample_ID",by.y="Taxa")
#temp2[,2]<-as.factor(temp2[,2])
traiinterest<-substring(RezNatGAPIT[12,6],7)
#data of the protected treatment
temp3<-tempb[,colnames(tempb)=="Sample_ID" | 
               colnames(tempb)==RezNatGAPIT$SNP[12]]
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
graf<-barplot(as.data.frame(table(temp2$`Survival`,
                                  temp2[,2]))$Freq,
              col=brewer.pal(12,"Set3")[6:7],las=1,
              space=c(0.1,rep(c(0.1,0.9),2),0.1),
              main="Survival",xlab=RezNatGAPIT$SNP[14])
#abline(h=c(50,100,200,300),col=grey(0.8,0.8),lwd=2,lty=1)
axis(1,at=(graf[c(1,3,5)]+graf[c(2,4,6)])/2,
     labels=colnames(table(temp2$`Survival`,
                           temp2[,2])))
box(bty="o",lwd=1)

graf<-barplot(proportions(table(temp2$`Survival`,
                                temp2[,2]),margin=2)[1,]*100,
              col=brewer.pal(12,"Set3")[6],las=1,space=1,
              ylim=c(0,100),names.arg="",
              main="Exposed death rate")
axis(1,at=graf,
     labels=paste(colnames(table(temp2$`Survival`,temp2[,2])),
                  " (n=",table(temp2[,2]),")",
                  sep=""))
box(bty="o",lwd=1)
grafb<-barplot(proportions(table(temp3$`Survival`,
                                 temp3[,2]),margin=2)[1,]*100,
               col=brewer.pal(12,"Set3")[6],las=1,space=1,
               ylim=c(0,100),names.arg="",
               main="Protected death rate")
axis(1,at=grafb,
     labels=paste(colnames(table(temp3$`Survival`,temp3[,2])),
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
graf<-barplot(as.data.frame(table(temp3$`Survival`,
                                  temp3[,2]))$Freq,
              col=brewer.pal(12,"Set3")[6:7],las=1,
              space=c(0.1,rep(c(0.1,0.9),2),0.1),
              main="Survival",xlab=RezLimGAPIT$SNP[6])
#abline(h=c(50,100,200,300),col=grey(0.8,0.8),lwd=2,lty=1)
axis(1,at=(graf[c(1,3,5)]+graf[c(2,4,6)])/2,
     labels=colnames(table(temp3$`Survival`,
                           temp3[,2])))
box(bty="o",lwd=1)
par(op)
dev.off()

#second plot for the distribution for the two treatments
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
graf<-barplot(proportions(table(temp2$`Survival`,
                                temp2[,2]),margin=2)[1,]*100,
              col=brewer.pal(12,"Set3")[6],las=1,space=1,
              ylim=c(0,100),names.arg="",
              main="Exposed death rate")
axis(1,at=graf,
     labels=paste(colnames(table(temp2$`Survival`,temp2[,2])),
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
grafb<-barplot(proportions(table(temp3$`Survival`,
                                 temp3[,2]),margin=2)[1,]*100,
               col=brewer.pal(12,"Set3")[6],las=1,space=1,
               ylim=c(0,100),names.arg="",
               main="Protected death rate")
axis(1,at=grafb,
     labels=paste(colnames(table(temp3$`Survival`,temp3[,2])),
                  " (n=",table(temp3[,2]),")",
                  sep=""))
box(bty="o",lwd=1)
par(op)
dev.off()


##############################################################################/
#END
##############################################################################/