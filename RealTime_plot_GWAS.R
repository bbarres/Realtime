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
#Figure 10: plotting the 3 significant SNP effect on ALL the traits####
##############################################################################/

#code for plotting all the significant SNP
temp<-snp.dat[snp.dat$Sample_ID %in% NatTrait$Taxa,]
temp<-temp[,c(1,48:866)]
tempb<-snp.dat[snp.dat$Sample_ID %in% LimTrait$Taxa,]
tempb<-tempb[,c(1,48:866)]

pdf("output/Figure_10_SNPinterest.pdf",width=15,height=10)
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
text(0.5,0.5,"Height\nin 2012 (cm)",font=2,cex=2.5)
plot.new()
text(0.5,0.5,"Mean infection\n(2009-2013)",font=2,cex=2.5)
plot.new()
text(0.5,0.7,"Mortality (%)",font=2,cex=2.5)
legend(0.25,0.55,legend=c("Natural exposure","Protected exposure"),
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
  #data of the exposed exposure
  temp2<-temp[,colnames(temp)=="Sample_ID" | 
                colnames(temp)==Rezcomb$SNP[i]]
  temp2<-merge(temp2,NatTrait,by.x="Sample_ID",by.y="Taxa")
  #temp2[,2]<-as.factor(temp2[,2])
  traiinterest<-substring(Rezcomb[i,6],7)
  treatme<-Rezcomb[i,7]
  #data of the protected exposure
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
          xaxt="n",lty=1,lwd=3,axes=0,yaxt="n",
          ylab="Acorn weight (g)",main="Acorn weight")
  axis(1,labels=paste(colnames(table(tempo$`Acorn weight`,tempo[,2])),
                      " (n=",table(tempo[,2]),")",
                      sep=""),at=c(1,2,3),cex.axis=1.3,lwd=2)
  axis(2,at=seq(from=1,to=13,by=2),las=1,cex.axis=1.3,lwd=2)
  box(bty="o")
  vioplot(tempo$`Height`~tempo[,2],boxwex=0.3,las=1,
          border=colo[2:4],lineCol=colo[2:4],rectCol=colo[2:4],
          col=colo[1],ann=FALSE,ylim=c(0,110),
          xaxt="n",lty=1,lwd=3,axes=0,yaxt="n",
          ylab="Height (cm)",main="Height")
  axis(1,labels=paste(colnames(table(tempo$`Height`,tempo[,2])),
                      " (n=",table(tempo[,2]),")",
                      sep=""),at=c(1,2,3),cex.axis=1.3,lwd=2)
  axis(2,at=seq(from=5,to=115,by=15),las=1,cex.axis=1.3,lwd=2)
  box(bty="o")
  vioplot(tempo$`Powdery mildew`~tempo[,2],boxwex=0.3,las=1,
          border=colo[2:4],lineCol=colo[2:4],rectCol=colo[2:4],
          col=colo[1],ann=FALSE,ylim=c(0,70),
          xaxt="n",lty=1,lwd=3,axes=0,yaxt="n",
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
#Figure SXX: plotting significant SNP genotype effect on ALL the traits####
##############################################################################/

#code for plotting all the significant SNP from the GWAS analysis on 
#multiple files
temp<-snp.dat[snp.dat$Sample_ID %in% NatTrait$Taxa,]
temp<-temp[,c(1,48:866)]
tempb<-snp.dat[snp.dat$Sample_ID %in% LimTrait$Taxa,]
tempb<-tempb[,c(1,48:866)]

for (j in 1:(dim(Rezcomb)[1]/4)) {
  pdf(paste("output/Figure_SNPinterest_ALL",j,".pdf",sep=""),
      width=15,height=13)
  nf<-layout(matrix(c(1,2,2,2,3,3,3,4,4,4,5,5,5,
                      6,10,10,10,11,11,11,12,12,12,13,13,13,
                      6,10,10,10,11,11,11,12,12,12,13,13,13,
                      6,10,10,10,11,11,11,12,12,12,13,13,13,
                      7,14,14,14,15,15,15,16,16,16,17,17,17,
                      7,14,14,14,15,15,15,16,16,16,17,17,17,
                      7,14,14,14,15,15,15,16,16,16,17,17,17,
                      8,18,18,18,19,19,19,20,20,20,21,21,21,
                      8,18,18,18,19,19,19,20,20,20,21,21,21,
                      8,18,18,18,19,19,19,20,20,20,21,21,21,
                      9,22,22,22,23,23,23,24,24,24,25,25,25,
                      9,22,22,22,23,23,23,24,24,24,25,25,25,
                      9,22,22,22,23,23,23,24,24,24,25,25,25),
                    13,13,byrow=TRUE))
  op<-par(mar=c(0,0,0,0),lwd=2)
  plot.new()
  plot.new()
  text(0.5,0.5,"Acorn weight (g)",font=2,cex=2.5)
  plot.new()
  text(0.5,0.5,"Height\nin 2012 (cm)",font=2,cex=2.5)
  plot.new()
  text(0.5,0.5,"Mean infection\n(2009-2013)",font=2,cex=2.5)
  plot.new()
  text(0.5,0.7,"Mortality (%)",font=2,cex=2.5)
  legend(0.25,0.55,legend=c("Natural exposure","Protected exposure"),
         fill=c(colnat[3],colprt[3]),border=c(colnat[3],colprt[3]),
         density=c(-1,20),angle=c(0,60),bty="n",
         cex=1.5,y.intersp=1.3,x.intersp=0.5,ncol=1,xpd=TRUE)
  
  for (i in (1+4*(j-1)):(4*j)) {
    plot.new()
    text(0.5,0.5,Rezcomb[i,1],srt=90,font=2,cex=1.5)
  }
  
  par(mar=c(3,3,1,1))
  for (i in (1+4*(j-1)):(4*j)) {
    #data of the exposed exposure
    temp2<-temp[,colnames(temp)=="Sample_ID" | 
                  colnames(temp)==Rezcomb$SNP[i]]
    temp2<-merge(temp2,NatTrait,by.x="Sample_ID",by.y="Taxa")
    #temp2[,2]<-as.factor(temp2[,2])
    traiinterest<-substring(Rezcomb[i,6],7)
    treatme<-Rezcomb[i,7]
    #data of the protected exposure
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
            xaxt="n",lty=1,lwd=3,axes=0,yaxt="n",
            ylab="Acorn weight (g)",main="Acorn weight")
    axis(1,labels=paste(colnames(table(tempo$`Acorn weight`,tempo[,2])),
                        " (n=",table(tempo[,2]),")",
                        sep=""),at=c(1,2,3),cex.axis=1.3,lwd=2)
    axis(2,at=seq(from=1,to=13,by=2),las=1,cex.axis=1.3,lwd=2)
    box(bty="o")
    vioplot(tempo$`Height`~tempo[,2],boxwex=0.3,las=1,
            border=colo[2:4],lineCol=colo[2:4],rectCol=colo[2:4],
            col=colo[1],ann=FALSE,ylim=c(0,110),
            xaxt="n",lty=1,lwd=3,axes=0,yaxt="n",
            ylab="Height (cm)",main="Height")
    axis(1,labels=paste(colnames(table(tempo$`Height`,tempo[,2])),
                        " (n=",table(tempo[,2]),")",
                        sep=""),at=c(1,2,3),cex.axis=1.3,lwd=2)
    axis(2,at=seq(from=5,to=115,by=15),las=1,cex.axis=1.3,lwd=2)
    box(bty="o")
    vioplot(tempo$`Powdery mildew`~tempo[,2],boxwex=0.3,las=1,
            border=colo[2:4],lineCol=colo[2:4],rectCol=colo[2:4],
            col=colo[1],ann=FALSE,ylim=c(0,70),
            xaxt="n",lty=1,lwd=3,axes=0,yaxt="n",
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
  dev.off()
  
}



#code for plotting all the significant SNP from the GWAS analysis on 
#a unique file
temp<-snp.dat[snp.dat$Sample_ID %in% NatTrait$Taxa,]
temp<-temp[,c(1,48:866)]
tempb<-snp.dat[snp.dat$Sample_ID %in% LimTrait$Taxa,]
tempb<-tempb[,c(1,48:866)]

pdf("output/Figure_SNPinterest_ALL.pdf",width=15,height=49)
nf<-layout(matrix(c(1,2,2,2,3,3,3,4,4,4,5,5,5,
                    6,22,22,22,23,23,23,24,24,24,25,25,25,
                    6,22,22,22,23,23,23,24,24,24,25,25,25,
                    6,22,22,22,23,23,23,24,24,24,25,25,25,
                    7,26,26,26,27,27,27,28,28,28,29,29,29,
                    7,26,26,26,27,27,27,28,28,28,29,29,29,
                    7,26,26,26,27,27,27,28,28,28,29,29,29,
                    8,30,30,30,31,31,31,32,32,32,33,33,33,
                    8,30,30,30,31,31,31,32,32,32,33,33,33,
                    8,30,30,30,31,31,31,32,32,32,33,33,33,
                    9,34,34,34,35,35,35,36,36,36,37,37,37,
                    9,34,34,34,35,35,35,36,36,36,37,37,37,
                    9,34,34,34,35,35,35,36,36,36,37,37,37,
                    10,38,38,38,39,39,39,40,40,40,41,41,41,
                    10,38,38,38,39,39,39,40,40,40,41,41,41,
                    10,38,38,38,39,39,39,40,40,40,41,41,41,
                    11,42,42,42,43,43,43,44,44,44,45,45,45,
                    11,42,42,42,43,43,43,44,44,44,45,45,45,
                    11,42,42,42,43,43,43,44,44,44,45,45,45,
                    12,46,46,46,47,47,47,48,48,48,49,49,49,
                    12,46,46,46,47,47,47,48,48,48,49,49,49,
                    12,46,46,46,47,47,47,48,48,48,49,49,49,
                    13,50,50,50,51,51,51,52,52,52,53,53,53,
                    13,50,50,50,51,51,51,52,52,52,53,53,53,
                    13,50,50,50,51,51,51,52,52,52,53,53,53,
                    14,54,54,54,55,55,55,56,56,56,57,57,57,
                    14,54,54,54,55,55,55,56,56,56,57,57,57,
                    14,54,54,54,55,55,55,56,56,56,57,57,57,
                    15,58,58,58,59,59,59,60,60,60,61,61,61,
                    15,58,58,58,59,59,59,60,60,60,61,61,61,
                    15,58,58,58,59,59,59,60,60,60,61,61,61,
                    16,62,62,62,63,63,63,64,64,64,65,65,65,
                    16,62,62,62,63,63,63,64,64,64,65,65,65,
                    16,62,62,62,63,63,63,64,64,64,65,65,65,
                    17,66,66,66,67,67,67,68,68,68,69,69,69,
                    17,66,66,66,67,67,67,68,68,68,69,69,69,
                    17,66,66,66,67,67,67,68,68,68,69,69,69,
                    18,70,70,70,71,71,71,72,72,72,73,73,73,
                    18,70,70,70,71,71,71,72,72,72,73,73,73,
                    18,70,70,70,71,71,71,72,72,72,73,73,73,
                    19,74,74,74,75,75,75,76,76,76,77,77,77,
                    19,74,74,74,75,75,75,76,76,76,77,77,77,
                    19,74,74,74,75,75,75,76,76,76,77,77,77,
                    20,78,78,78,79,79,79,80,80,80,81,81,81,
                    20,78,78,78,79,79,79,80,80,80,81,81,81,
                    20,78,78,78,79,79,79,80,80,80,81,81,81,
                    21,82,82,82,83,83,83,84,84,84,85,85,85,
                    21,82,82,82,83,83,83,84,84,84,85,85,85,
                    21,82,82,82,83,83,83,84,84,84,85,85,85),
                  49,13,byrow=TRUE))
op<-par(mar=c(0,0,0,0),lwd=2)
plot.new()
plot.new()
text(0.5,0.5,"Acorn weight (g)",font=2,cex=2.5)
plot.new()
text(0.5,0.5,"Height\nin 2012 (cm)",font=2,cex=2.5)
plot.new()
text(0.5,0.5,"Mean infection\n(2009-2013)",font=2,cex=2.5)
plot.new()
text(0.5,0.7,"Mortality (%)",font=2,cex=2.5)
legend(0.25,0.55,legend=c("Natural exposure","Protected exposure"),
       fill=c(colnat[3],colprt[3]),border=c(colnat[3],colprt[3]),
       density=c(-1,20),angle=c(0,60),bty="n",
       cex=1.5,y.intersp=1.3,x.intersp=0.5,ncol=1,xpd=TRUE)

for (i in 1:dim(Rezcomb)[1]) {
  plot.new()
  text(0.5,0.5,Rezcomb[i,1],srt=90,font=2,cex=1.5)
}

par(mar=c(3,3,1,1))
for (i in 1:dim(Rezcomb)[1]) {
  #data of the exposed exposure
  temp2<-temp[,colnames(temp)=="Sample_ID" | 
                colnames(temp)==Rezcomb$SNP[i]]
  temp2<-merge(temp2,NatTrait,by.x="Sample_ID",by.y="Taxa")
  #temp2[,2]<-as.factor(temp2[,2])
  traiinterest<-substring(Rezcomb[i,6],7)
  treatme<-Rezcomb[i,7]
  #data of the protected exposure
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
          xaxt="n",lty=1,lwd=3,axes=0,yaxt="n",
          ylab="Acorn weight (g)",main="Acorn weight")
  axis(1,labels=paste(colnames(table(tempo$`Acorn weight`,tempo[,2])),
                      " (n=",table(tempo[,2]),")",
                      sep=""),at=c(1,2,3),cex.axis=1.3,lwd=2)
  axis(2,at=seq(from=1,to=13,by=2),las=1,cex.axis=1.3,lwd=2)
  box(bty="o")
  vioplot(tempo$`Height`~tempo[,2],boxwex=0.3,las=1,
          border=colo[2:4],lineCol=colo[2:4],rectCol=colo[2:4],
          col=colo[1],ann=FALSE,ylim=c(0,110),
          xaxt="n",lty=1,lwd=3,axes=0,yaxt="n",
          ylab="Height (cm)",main="Height")
  axis(1,labels=paste(colnames(table(tempo$`Height`,tempo[,2])),
                      " (n=",table(tempo[,2]),")",
                      sep=""),at=c(1,2,3),cex.axis=1.3,lwd=2)
  axis(2,at=seq(from=5,to=115,by=15),las=1,cex.axis=1.3,lwd=2)
  box(bty="o")
  vioplot(tempo$`Powdery mildew`~tempo[,2],boxwex=0.3,las=1,
          border=colo[2:4],lineCol=colo[2:4],rectCol=colo[2:4],
          col=colo[1],ann=FALSE,ylim=c(0,70),
          xaxt="n",lty=1,lwd=3,axes=0,yaxt="n",
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
dev.off()


##############################################################################/
#END
##############################################################################/