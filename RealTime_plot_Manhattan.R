##############################################################################/
##############################################################################/
#Manhattan plot of the GWAS analyses results
##############################################################################/
##############################################################################/

source("RealTime_load.R")
#in order to load the results of the GWAS analyses, you need to run the 
#RealTime_GWAS.R script before this script


##############################################################################/
#function to import and prepare the output files from GAPIT analysis####
##############################################################################/

#defining a function to import and prepare the data set from GAPIT
#output file
readManDa<-function(pathtoGWASrez) {
  datMan<-read.csv(file=pathtoGWASrez)
  #we change the chromosome information into a factor with convenient
  #levels order
  datMan$Chromosome<-factor(datMan$Chromosome,
                               levels=c("1","2","3","4","5","6","7",
                                        "8","9","10","11","12","13",
                                        "MT","Pltd"))
  #there is only 12 chromosome in Q. robur, the 13 stand for an 
  #unlocated SNP on the genome at the time of the analysis. The level
  #"13" is therefore turn into "NA"
  levels(datMan$Chromosome)[13]<-"Na"
  levels(datMan$Chromosome)[14]<-"Mt"
  levels(datMan$Chromosome)[15]<-"Cp"
  #we turn Position from integer to numeric
  datMan$Position<-as.numeric(datMan$Position)
  datMan<-datMan[order(datMan$Chromosome,datMan$Position),]
  datMan$logp<- -log10(datMan$P.value)
  datMan$grpNum<-as.numeric(datMan$Chromosome)
  datMan$decal<-rep(c(0,cumsum(tapply(datMan$Position,
                                      datMan$Chromosome,max))[1:14]),
                    times=table(datMan$Chromosome))
  datMan$signiCat<-cut(datMan$FDR_Adjusted_P.values,
                       breaks=c(0,0.0001,0.001,0.01,1),
                       labels=c("<0.0001","<0.001","<0.01","ns"))
  return(datMan)
}


##############################################################################/
#Loading the data using the defined function####
##############################################################################/

AcNat<-readManDa("output/natGWAS/GAPIT.Blink.Acorn weight.GWAS.Results.csv")
AcLim<-readManDa("output/limGWAS/GAPIT.Blink.Acorn weight.GWAS.Results.csv")
HeNat<-readManDa("output/natGWAS/GAPIT.Blink.Height.GWAS.Results.csv")
HeLim<-readManDa("output/limGWAS/GAPIT.Blink.Height.GWAS.Results.csv")
PMNat<-readManDa("output/natGWAS/GAPIT.Blink.Powdery mildew.GWAS.Results.csv")
PMLim<-readManDa("output/limGWAS/GAPIT.Blink.Powdery mildew.GWAS.Results.csv")
SuNat<-readManDa("output/natGWAS/GAPIT.Blink.Survival.GWAS.Results.csv")
SuLim<-readManDa("output/limGWAS/GAPIT.Blink.Survival.GWAS.Results.csv")


##############################################################################/
#function to plot a Manhattan plot####
##############################################################################/

ManhaPlot<-function(datMan,colovec,colosign="red",
                    decalCHR=0,desiXax=1,ylimi=c(0,10)) {
  datMan$posabsol<-datMan$decal+
    datMan$Position+
    (datMan$grpNum-1)*decalCHR
  plot(datMan$posabsol,datMan$logp,pch=19,las=1,font.axis=2,
       ylab=expression(-log[10](pval)),cex.lab=1.5,cex=1.0,
       xaxt="n",yaxt="n",bty="n",xlab="",ylim=ylimi,
       col=colovec[as.numeric(even(as.numeric(datMan$Chromosome)))+1])
  points(datMan[datMan$signiCat=="<0.01",]$posabsol,
         datMan[datMan$signiCat=="<0.01",]$logp,
         pch=21,bg=colosign[1],cex=1.5)
  points(datMan[datMan$signiCat=="<0.001",]$posabsol,
         datMan[datMan$signiCat=="<0.001",]$logp,
         pch=21,bg=colosign[2],cex=1.5)
  points(datMan[datMan$signiCat=="<0.0001",]$posabsol,
         datMan[datMan$signiCat=="<0.0001",]$logp,
         pch=21,bg=colosign[3],cex=1.5)
  axis(2,lwd=2,las=1)
  if(desiXax==1){
    axis(1,lwd=2,cex.axis=1,
         at=c(c(0,cumsum(tapply(datMan$Position,datMan$Chromosome,max))[1:14])+
                tapply(datMan$Position,datMan$Chromosome,max)/2+
                (0:14)*decalCHR),
         lab=levels(datMan$Chromosome),las=1,hadj=0.5,padj=-0.5,font=2)
  } else {
    axis(3,lwd=2,
         at=c(c(0,cumsum(tapply(datMan$Position,datMan$Chromosome,max))[1:14])+
                tapply(datMan$Position,datMan$Chromosome,max)/2+
                (0:14)*decalCHR),
         lab=NA)
  }
  box(lwd=2)
}


##############################################################################/
#Ploting using the function####
##############################################################################/

#Figure for one trait
colovec<-c(brewer.pal(12,"Paired")[1:2],brewer.pal(12,"Paired")[3:4])
colosign<-brewer.pal(9,"YlOrRd")[c(4,6,8)]
op<-par(mfrow=c(2,1),mar=c(1,5,3,0.5))
#Acorn weight
ManhaPlot(AcNat,colovec[1:2],colosign,
          decalCHR=40000000,desiXax=1,ylimi=c(0,8))
title(main="Acorn weight",font=2,cex.main=3)
legend(-50000000,8.7,legend=c("Exposed","Protected"),xpd=TRUE,
       pch=22,pt.cex=4,pt.lwd=4,cex=1.3,
       pt.bg=colovec[c(1,3)],col=colovec[c(2,4)],bty="n",
       x.intersp=0.5,y.intersp=1.6)
legend(150000000,8.5,legend=c("P<0.01","P<0.001","P<0.0001"),xpd=TRUE,
       pch=21,pt.cex=2,pt.lwd=1,cex=1.2,
       pt.bg=colosign,col="black",bty="n",
       x.intersp=0.4,y.intersp=1.1)
par(mar=c(2.5,5,1.5,0.5))
ManhaPlot(AcLim,colovec[3:4],colosign,
          decalCHR=40000000,desiXax=0,ylimi=c(8,0))
par(op)
#export to .pdf 15 x 7 inches

#Figure for the article, combining the 4 traits GWAS results
colovec<-c(brewer.pal(12,"Paired")[1:2],brewer.pal(12,"Paired")[3:4])
colosign<-brewer.pal(9,"YlOrRd")[c(4,6,8)]
pdf(file="output/Figure_Manhat.pdf",width=10,height=13)
op<-par(mfrow=c(8,1),mar=c(1,5,2.5,0.5))
#survival
ManhaPlot(SuNat,colovec[3:4],colosign,
          decalCHR=40000000,desiXax=1,ylimi=c(0,11))
title(main="Survival",font=2,cex.main=3)
legend(-35000000,11.5,legend=c("Natural","Protected"),xpd=TRUE,
       pch=22,pt.cex=4,pt.lwd=4,cex=1.4,
       pt.bg=colovec[c(3,1)],col=colovec[c(4,2)],bty="n",
       x.intersp=1.1,y.intersp=1.4)
legend(150000000,11.5,legend=c("Pcor<0.01","Pcor<0.001","Pcor<0.0001"),
       xpd=TRUE,pch=21,pt.cex=2,pt.lwd=1,cex=1.4,
       pt.bg=colosign,col="black",bty="n",
       x.intersp=0.8,y.intersp=0.9)
par(mar=c(2,5,1.5,0.5))
ManhaPlot(SuLim,colovec[1:2],colosign,
          decalCHR=40000000,desiXax=0,ylimi=c(11,0))
#Powdery mildew
par(mar=c(1,5,2.5,0.5))
ManhaPlot(PMNat,colovec[3:4],colosign,
          decalCHR=40000000,desiXax=1,ylimi=c(0,8))
title(main="Powdery mildew",font=2,cex.main=3)
par(mar=c(2,5,1.5,0.5))
ManhaPlot(PMLim,colovec[1:2],colosign,
          decalCHR=40000000,desiXax=0,ylimi=c(8,0))
#Height
par(mar=c(1,5,2.5,0.5))
ManhaPlot(HeNat,colovec[3:4],colosign,
          decalCHR=40000000,desiXax=1,ylimi=c(0,8))
title(main="Height",font=2,cex.main=3)
par(mar=c(2,5,1.5,0.5))
ManhaPlot(HeLim,colovec[1:2],colosign,
          decalCHR=40000000,desiXax=0,ylimi=c(8,0))
#Acorn weight
par(mar=c(1,5,2.5,0.5))
ManhaPlot(AcNat,colovec[3:4],colosign,
          decalCHR=40000000,desiXax=1,ylimi=c(0,9))
title(main="Acorn weight",font=2,cex.main=3)
par(mar=c(2,5,1.5,0.5))
ManhaPlot(AcLim,colovec[1:2],colosign,
          decalCHR=40000000,desiXax=0,ylimi=c(9,0))
par(op)
#export 10 x 13 inches in .pdf
dev.off()


##############################################################################/
#END
##############################################################################/