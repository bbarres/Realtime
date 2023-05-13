##############################################################################/
##############################################################################/
#Individual genetic diversity analyses
##############################################################################/
##############################################################################/

#loading the necessary packages and data sets
source("RealTime_load.R")

#loading and preparing the dataset
#you first need to run the 'RealTime_GENHET.R' script
HetVivMort<-read.table(file="data/Hetind_DoA.txt",sep="\t",
                       header=TRUE,colClasses="character")

HetVivMortExp<-HetVivMort[HetVivMort$moda=="exp",]
#removing global
HetVivMortExp<-HetVivMortExp[HetVivMortExp$fam!="Global",]
HetVivMortExp$fam[HetVivMortExp$fam=="1"]<-"01"
HetVivMortExp$fam[HetVivMortExp$fam=="9"]<-"09"

HetVivMortLim<-HetVivMort[HetVivMort$moda=="low",]
#removing global
HetVivMortLim<-HetVivMortLim[HetVivMortLim$fam!="Global",]
HetVivMortLim$fam[HetVivMortLim$fam=="1"]<-"01"
HetVivMortLim$fam[HetVivMortLim$fam=="9"]<-"09"


##############################################################################/
#Figure 8: comparing Dead or Alive heterozygosities indices by family####
##############################################################################/

#plot of the PHt index for both natural and protected treatment by family
pdf(file="output/Figure_8_PHt.pdf",width=10,height=10)
op<-par(mfrow=c(2,1),mar=c(2,4,4,1))
#semi violin plot for the natural treatment
HetAli<-HetVivMortExp[HetVivMortExp$vivmor==1,]
HetDea<-HetVivMortExp[HetVivMortExp$vivmor==0,]
colovec<-c(brewer.pal(12,"Set3")[6:7],
           brewer.pal(9,"Set1")[1:2])
vioplot(as.numeric(HetVivMortExp$PHt)~HetVivMortExp$fam,xaxt="n",yaxt="n",
        col = c("transparent"),sep=":",las=1,border="transparent",
        ylab="PHt",xlab="",ylim=c(0.07,0.35),frame.plot=FALSE,cex.main=2,
        lineCol="transparent",rectCol="transparent",main="Natural exposure")
axis(1,lwd=2,at=c(1:15),labels=levels(as.factor(HetVivMortExp$fam)),font=2)
axis(2,lwd=2,las=1,font=2)
vioplot(as.numeric(HetDea$PHt)~HetDea$fam,plotCentre="line",
        col=colovec[1],sep=":",las=1,side="left",frame.plot=FALSE,
        add=TRUE)
stripchart(as.numeric(HetDea$PHt)~HetDea$fam,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col=colovec[3],
           at=c(1:15)-0.2,cex=0.8)
vioplot(as.numeric(HetAli$PHt)~HetAli$fam,plotCentre="line",
        col=colovec[2],sep=":",las=1,side="right",frame.plot=FALSE,
        add=TRUE)
stripchart(as.numeric(HetAli$PHt)~HetAli$fam,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col=colovec[4],
           at=c(1:15)+0.2,cex=0.8)
segments(c(1:15)-0.2,aggregate(as.numeric(HetDea$PHt),
                               list(HetDea$fam),FUN=mean)[,2],
         c(1:15)+0.2,aggregate(as.numeric(HetAli$PHt),
                               list(HetAli$fam),FUN=mean)[,2],
         col=grey(0.95,0.9),lwd=6)
points(x=c(1:15)-0.2,y=aggregate(as.numeric(HetDea$PHt),
                                 list(HetDea$fam),FUN=mean)[,2],
       pch=21,bg=colovec[3],cex=1.2)
points(x=c(1:15)+0.2,y=aggregate(as.numeric(HetAli$PHt),
                                 list(HetAli$fam),FUN=mean)[,2],
       pch=21,bg=colovec[4],cex=1.2)
box(bty="l",lwd=2)
legend(0.1,0.15,c("dead","alive"),fill=colovec[1:2],cex=1.5,
       bty="n",x.intersp=0.7,y.intersp=0.9,xpd=TRUE)
#semi violin plot for the protected exposure
HetAli<-HetVivMortLim[HetVivMortLim$vivmor==1,]
HetDea<-HetVivMortLim[HetVivMortLim$vivmor==0,]
colovec<-c(brewer.pal(12,"Set3")[6:7],
           brewer.pal(9,"Set1")[1:2])
vioplot(as.numeric(HetVivMortLim$PHt)~HetVivMortLim$fam,xaxt="n",yaxt="n",
        col = c("transparent"),sep=":",las=1,border="transparent",
        ylab="PHt",xlab="",ylim=c(0.07,0.35),frame.plot=FALSE,
        lineCol="transparent",rectCol="transparent",cex.main=2,
        main="Protected exposure")
axis(1,lwd=2,at=c(1:15),labels=levels(as.factor(HetVivMortExp$fam)),font=2)
axis(2,lwd=2,las=1,font=2)
vioplot(as.numeric(HetDea$PHt)~HetDea$fam,plotCentre="line",
        col=colovec[1],sep=":",las=1,side="left",frame.plot=FALSE,
        add=TRUE)
stripchart(as.numeric(HetDea$PHt)~HetDea$fam,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col=colovec[3],
           at=c(1:15)-0.2,cex=0.8)
vioplot(as.numeric(HetAli$PHt)~HetAli$fam,plotCentre="line",
        col=colovec[2],sep=":",las=1,side="right",frame.plot=FALSE,
        add=TRUE)
stripchart(as.numeric(HetAli$PHt)~HetAli$fam,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col=colovec[4],
           at=c(1:15)+0.2,cex=0.8)
segments(c(1:15)-0.2,aggregate(as.numeric(HetDea$PHt),
                               list(HetDea$fam),FUN=mean)[,2],
         c(1:15)+0.2,aggregate(as.numeric(HetAli$PHt),
                               list(HetAli$fam),FUN=mean)[,2],
         col=grey(0.95,0.9),lwd=6)
points(x=c(1:15)-0.2,y=aggregate(as.numeric(HetDea$PHt),
                                 list(HetDea$fam),FUN=mean)[,2],
       pch=21,bg=colovec[3],cex=1.2)
points(x=c(1:15)+0.2,y=aggregate(as.numeric(HetAli$PHt),
                                 list(HetAli$fam),FUN=mean)[,2],
       pch=21,bg=colovec[4],cex=1.2)
box(bty="l",lwd=2)
par(op)
#export to .pdf 10 x 10 inches
dev.off()


##############################################################################/
#Figure S7: Correlation between individual heterozygosity indices####
##############################################################################/

#a function to compute the absolute correlation between pairs of variables
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#loading and preparing the dataset
#you first need to run the 'RealTime_GENHET.R' script
HetVivMort<-read.table(file="data/Hetind_DoA.txt",sep="\t",
                       header=TRUE)
HetVivMortGlob<-HetVivMort[HetVivMort$fam=="Global",]

#plot the graph
pdf(file="output/Figure_S7_corHeteroz.pdf",width=10,height=10)
pairs(HetVivMortGlob[,c(2:6)],las=1,cex.main=2,
      main="Correlation between heterozygosity indices",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#export to pdf 10 x 7 inches
dev.off()


##############################################################################/
#Figure S16: Plot of the distribution of mortality by PHt classes####
##############################################################################/

#you first need to run the 'RealTime_GENHET.R' script
pdf(file="output/Figure_SX_mortaPHtclass.pdf",width=5,height=8)
colovec<-c(brewer.pal(12,"Set3")[6:7],
           brewer.pal(9,"Set1")[1:2])
HetVivMort<-read.table(file="data/Hetind_DoA.txt",sep="\t",
                       header=TRUE)
HetVivMortFam<-HetVivMort[HetVivMort$fam!="Global",]
HetVivMortFam$catHet<-cut(HetVivMortFam$PHt,
                          breaks=c(0.07,0.22,0.24,0.26,0.28,0.30,0.32,0.36))
effectif<-colSums(table(HetVivMortFam$vivmor,HetVivMortFam$catHet))
freqMor<-proportions(table(HetVivMortFam$vivmor,HetVivMortFam$catHet),
                     margin=2)*100
temp<-barplot(freqMor,las=1,main="Mortality rate by PHt classes",
              col=colovec[1:2],axes=FALSE,axisnames=FALSE,space=0.5)
axis(1,at=temp,labels=FALSE,lwd=3,font=2)
text(temp+0.3,par("usr")[1]-10,labels=names(effectif),srt=-60,
     xpd=TRUE,cex=1,font=2)
axis(2,lwd=3,font=2,cex.axis=1.2,las=1)
box(bty="l",lwd=3)
text(temp,102,paste("n=",effectif,sep=""),font=3,cex=0.9,xpd=TRUE)
legend(-2,115,c("dead","alive"),fill=colovec[1:2],cex=1.3,
       bty="n",x.intersp=0.5,y.intersp=0.7,xpd=TRUE)
#export to .pdf 5 x 8 inches
dev.off()


##############################################################################/
#Additional Figure: global data comparing Exposed and Limited treatments####
##############################################################################/

HetVivMort<-read.table(file="data/Hetind_DoA.txt",sep="\t",
                       header=TRUE)
#preparing the data set
HetVivMortGlob<-HetVivMort[HetVivMort$fam=="Global",]
#semi violin plot for the exposed treatment
HetAli<-HetVivMortGlob[HetVivMortGlob$vivmor==1,]
HetDea<-HetVivMortGlob[HetVivMortGlob$vivmor==0,]
colovec<-c(brewer.pal(12,"Set3")[6:7],
           brewer.pal(9,"Set1")[1:2])
vioplot(as.numeric(HetVivMortGlob$PHt)~HetVivMortGlob$moda,
        col = c("transparent"),sep=":",las=1,border="transparent",
        ylab="PHt",xlab="",ylim=c(0.07,0.35),frame.plot=FALSE,
        lineCol="transparent",rectCol="transparent",
        main="Comparison between treatments")
vioplot(as.numeric(HetDea$PHt)~HetDea$moda,plotCentre="line",
        col=colovec[1],sep=":",las=1,side="left",frame.plot=FALSE,
        add=TRUE)
stripchart(as.numeric(HetDea$PHt)~HetDea$moda,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col=colovec[3],
           at=c(1:2)-0.2)
vioplot(as.numeric(HetAli$PHt)~HetAli$moda,plotCentre="line",
        col=colovec[2],sep=":",las=1,side="right",frame.plot=FALSE,
        add=TRUE)
stripchart(as.numeric(HetAli$PHt)~HetAli$moda,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col=colovec[4],
           at=c(1:2)+0.2)
segments(c(1:2)-0.2,aggregate(as.numeric(HetDea$PHt),
                              list(HetDea$moda),FUN=mean)[,2],
         c(1:2)+0.2,aggregate(as.numeric(HetAli$PHt),
                              list(HetAli$moda),FUN=mean)[,2],
         col=grey(0.95,0.9),lwd=6)
points(x=c(1:2)-0.2,y=aggregate(as.numeric(HetDea$PHt),
                                list(HetDea$moda),FUN=mean)[,2],
       pch=21,bg=colovec[3],col="black")
points(x=c(1:2)+0.2,y=aggregate(as.numeric(HetAli$PHt),
                                list(HetAli$moda),FUN=mean)[,2],
       pch=21,bg=colovec[4],col="black")
box(bty="l")
legend(0.2,0.4,c("dead","alive"),fill=colovec[1:2],cex=1.3,
       bty="n",x.intersp=0.5,y.intersp=0.7,xpd=TRUE)
#export to .pdf 6 x 8 inches


##############################################################################/
#Additional Figure: Heterozygosity computation: total vs surviving####
##############################################################################/

#reshaping the data set for computation
n.temp<-seppop(snpGen2,treatOther=TRUE)
n.temp.other<-lapply(n.temp,as.data.frame(other))
temp<-repool(n.temp$exp1,n.temp$exp0)
pop(temp)<-rep("expinit",times=nInd(temp))
temp<-repool(temp,n.temp$exp1)
temp2<-repool(n.temp$low1,n.temp$low0)
pop(temp2)<-rep("lowinit",times=nInd(temp2))
snpGen3<-repool(temp,temp2,n.temp$low1)
snpGen3@other<-rbind(n.temp.other$exp1,n.temp.other$exp0,
                     n.temp.other$exp1,n.temp.other$low1,
                     n.temp.other$low0,n.temp.other$low1)
colnames(snpGen3@other)<-c("fam","treat","height","DoA","newPop")
snpGen3@other$newPop<-pop(snpGen3)

nomFam<-popNames(snpGen3)
temp<-genind2df(snpGen3,oneColPerAll=TRUE)
temp[temp=="A"]<-10
temp[temp=="T"]<-20
temp[temp=="C"]<-30
temp[temp=="G"]<-40
tempop<-temp[,1]
sampleid<-row.names(temp)
temp<-temp[,-1]
temp<-data.frame(lapply(temp,as.numeric))
temp$sampleid<-sampleid
temp<-temp[,c(1639,1:1638)]
temp2<-GENHET(dat=temp,estimfreq="T",locname=nomSNP$simpleNames)
temp2<-as.data.frame(temp2)
temp2$pop<-tempop
temp2$moda<-stringr::str_match(tempop, "(...)(.*)")[,2]
temp2$vivmor<-stringr::str_match(tempop, "(...)(.*)")[,3]
temp2$fam<-"Global"
HetStarSto<-temp2

#same thing by family
pop(snpGen3)<-snpGen3@other$fam
snpGen3fam<-seppop(snpGen3)
nomFam<-popNames(snpGen3)
for (i in 1:length(nomFam)) {
  pop(snpGen3fam[[i]])<-snpGen3fam[[i]]@other$newPop
  temp<-genind2df(snpGen3fam[[i]],oneColPerAll=TRUE)
  temp[temp=="A"]<-10
  temp[temp=="T"]<-20
  temp[temp=="C"]<-30
  temp[temp=="G"]<-40
  tempop<-temp[,1]
  sampleid<-row.names(temp)
  temp<-temp[,-1]
  temp<-data.frame(lapply(temp,as.numeric))
  temp$sampleid<-sampleid
  temp<-temp[,c(1639,1:1638)]
  temp2<-GENHET(dat=temp,estimfreq="T",locname=nomSNP$simpleNames)
  temp2<-as.data.frame(temp2)
  temp2$pop<-tempop
  temp2$moda<-stringr::str_match(tempop, "(...)(.*)")[,2]
  temp2$vivmor<-stringr::str_match(tempop, "(...)(.*)")[,3]
  temp2$fam<-nomFam[i]
  HetStarSto<-rbind(HetStarSto,temp2)
}

#in order to have the good order of the categories, we change 
#the vivmor catefories
HetStarSto[HetStarSto$vivmor=="init","vivmor"]<-"beg"
HetStarSto[HetStarSto$vivmor=="1","vivmor"]<-"end"

HetStarStoExp<-HetStarSto[HetStarSto$moda=="exp",]
HetStarStoLim<-HetStarSto[HetStarSto$moda=="low",]

#the plot
op<-par(mfrow=c(5,1))
vioplot(as.numeric(HetStarStoExp$PHt)~HetStarStoExp$vivmor:HetStarStoExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Exposed / PHt")
vioplot(as.numeric(HetStarStoExp$Hs_obs)~
          HetStarStoExp$vivmor:HetStarStoExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Exposed / Hs_obs")
vioplot(as.numeric(HetStarStoExp$Hs_exp)~
          HetStarStoExp$vivmor:HetStarStoExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Exposed / Hs_exp")
vioplot(as.numeric(HetStarStoExp$IR)~HetStarStoExp$vivmor:HetStarStoExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Exposed / IR")
vioplot(as.numeric(HetStarStoExp$HL)~HetStarStoExp$vivmor:HetStarStoExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Exposed / HL")
par(op)
#export to .pdf 20 x 20 inches

op<-par(mfrow=c(5,1))
vioplot(as.numeric(HetStarStoLim$PHt)~HetStarStoLim$vivmor:HetStarStoLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Limited / PHt")
vioplot(as.numeric(HetStarStoLim$Hs_obs)~
          HetStarStoLim$vivmor:HetStarStoLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Limited / Hs_obs")
vioplot(as.numeric(HetStarStoLim$Hs_exp)~
          HetStarStoLim$vivmor:HetStarStoLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Limited / Hs_exp")
vioplot(as.numeric(HetStarStoLim$IR)~HetStarStoLim$vivmor:HetStarStoLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Limited / IR")
vioplot(as.numeric(HetStarStoLim$HL)~HetStarStoLim$vivmor:HetStarStoLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Limited / HL")
par(op)
#export to .pdf 20 x 20 inches


##############################################################################/
#END
##############################################################################/