##############################################################################/
##############################################################################/
#Plot the map of the RealTime experiment
##############################################################################/
##############################################################################/

source("RealTime_load.R")


##############################################################################/
#yearly evolution of the survival rate####
##############################################################################/

#limiting the data set to useful data
evolsurv<-dispo[dispo$family_simp!="CC" & dispo$family_simp!="hd" & 
                   !is.na(dispo$an_mort) & dispo$an_mort!="g",
                 c("Sample_ID","bloc","PU","treat","an_mort","family_simp")]
evolsurv$grpbloc<-paste(evolsurv$bloc,evolsurv$PU,evolsurv$treat,sep="")
#defining a color vector
colovec<-c(brewer.pal(12,"Paired")[4],brewer.pal(12,"Paired")[2])

survEv<-t(apply(table(evolsurv$grpbloc,evolsurv$an_mort),1,cumsum))
survEv<-100*(survEv[,10]-survEv)/survEv[,10]
survEv<-survEv[,-10]
survEv.bloc<-factor(substr(rownames(survEv),1,2))
survEv.treat<-factor(substr(rownames(survEv),3,5))
matplot(t(survEv),type="b",ylim=c(0,100),las=1,
        col=colovec[as.numeric(survEv.treat)],lwd=2,
        lty=as.numeric(survEv.treat),
        pch=c(0,15,1,2,3,16,17,5,6),
        xaxt="n",yaxt="n",bty="n",
        ylab="Percentage of survival",xlab="Year",
        font.lab=2,cex.lab=1.5)
axis(1,lwd=2,cex.axis=1,at=c(1:dim(survEv)[2]),
     lab=colnames(survEv),las=1,font=2)
axis(2,lwd=2,las=1,font=2)
box(lwd=2)
legend(1,55,legend=rep(NA,6),title="Natural\ntreatment",title.adj=0,
       col=colovec[1],lty=1,pch=c(22,21,24,3,23,25),pt.bg="white",
       bty="n",cex=1.3,title.cex=1,y.intersp=1.5,x.intersp=0.5,lwd=2)
legend(3,55,legend=rep(NA,3),title="Protected\ntreatment",title.adj=0,
       col=colovec[2],lty=2,pch=c(15,16,17),pt.bg="white",
       bty="n",cex=1.3,title.cex=1,y.intersp=1.5,x.intersp=0.5,lwd=2)
#export to .pdf 6 x 6 inches


##############################################################################/
#yearly level of infection####
##############################################################################/

#limiting the data set to useful data
evolMild<-dispo[dispo$family_simp!="CC" & dispo$family_simp!="hd" & 
                  !is.na(dispo$an_mort) & dispo$an_mort!="g",
                c("Sample_ID","bloc","PU","treat","family_simp",
                  "oid4_09","oid5_10","oid5_11","oid4_12","oid2_13",
                  "oid_16","oid_17","oid_moy")]
evolMild$oid_14<-0
evolMild$oid_15<-0
evolMild<-evolMild[,c("Sample_ID","bloc","PU","treat","family_simp",
                      "oid4_09","oid5_10","oid5_11","oid4_12","oid2_13",
                      "oid_14","oid_15","oid_16","oid_17","oid_moy")]
evolMild$grpbloc<-paste(evolMild$bloc,evolMild$PU,evolMild$treat,sep="")
colnames(evolMild)<-c("Sample_ID","bloc","PU","treat","family_simp",
                      "m2009","m2010","m2011","m2012","m2013",
                      "m2014","m2015","m2016","m2017","oid_moy","grpbloc")
#defining a color vector
colovec<-c(brewer.pal(12,"Paired")[4],brewer.pal(12,"Paired")[2])

infectM<-data.frame("expmean"=as.numeric(),"lowmean"=as.numeric(),
                    "expsd"=as.numeric(),"lowsd"=as.numeric())
for (i in 6:14) {
  infectM<-rbind(infectM,
                 c(tapply(evolMild[,i],evolMild$treat,FUN=mean,na.rm=TRUE),
                   tapply(evolMild[,i],evolMild$treat,FUN=sd,na.rm=TRUE)))
}
colnames(infectM)<-c("expmean","lowmean","expsd","lowsd")
rownames(infectM)<-c(2009:2017)

temp1<-pivot_longer(infectM,cols=1:2,names_to="treat",values_to="mean")
graf<-barplot(temp1$mean,density=c(-1,20),angle=c(0,60),border=NA,
              col=colovec,ylim=c(-1,105),ylab="Percentage of infection",
              cex.axis=1,cex.lab=1.5,las=1,yaxt="n",bty="n",
              space=c(0,0.25,rep(c(1,0.25),times=8)),font.lab=2)
axis(1,at=graf[seq(1,17,2)]+0.625,labels=rownames(infectM),
     lwd=2,font=2,cex.axis=0.9)
axis(2,at=seq(0,100,20),labels=seq(0,100,20),lwd=2,las=1,
     font=2,cex.axis=1)
box(bty="o",lwd=2)
#adding standard deviation
temp2<-pivot_longer(infectM,cols=3:4,names_to="treat",values_to="sd")
plotCI(x=graf,y=temp1$mean,
       ui=temp1$mean+temp2$sd,
       li=temp1$mean,
       pch=NA,lwd=1.5,add=TRUE)
#for cosmetic reason, we superimpose again the barplot
par(lwd=2)
barplot(temp1$mean,density=c(-1,20),angle=c(0,60),border=colovec,
        col=colovec,ylim=c(-1,105),ylab="Percentage of infection",
        cex.axis=1,cex.lab=1.5,las=1,yaxt="n",bty="n",
        space=c(0,0.25,rep(c(1,0.25),times=8)),font.lab=2,add=TRUE)
legend(-1,103,legend=c("Natural treatment","Protected treatment"),
       fill=colovec,density=c(-1,20),angle=c(0,60),bty="n",border=colovec,
       cex=1.3,y.intersp=1.5,x.intersp=0.5)
par(lwd=1)
#export to .pdf 6 x 6 inches


##############################################################################/
#yearly evolution of seedlings height####
##############################################################################/

#limiting the data set to useful data
evolHeight<-dispo[dispo$family_simp!="CC" & dispo$family_simp!="hd" & 
                    !is.na(dispo$an_mort) & dispo$an_mort!="g",
                  c("Sample_ID","bloc","PU","treat","family_simp","H09v",
                    "H10v","H11v","H12v","H14v","H15v","H16v","H17v")]
evolHeight$grpbloc<-paste(evolHeight$bloc,evolHeight$PU,
                          evolHeight$treat,sep="")
#defining a color vector
colovec<-c(brewer.pal(12,"Paired")[4],brewer.pal(12,"Paired")[2])

HeighTr<-data.frame("expmean"=as.numeric(),"lowmean"=as.numeric(),
                    "expsd"=as.numeric(),"lowsd"=as.numeric())
for (i in 6:13) {
  HeighTr<-rbind(HeighTr,
                 c(tapply(evolHeight[,i],evolHeight$treat,
                          FUN=mean,na.rm=TRUE),
                   tapply(evolHeight[,i],evolHeight$treat,
                          FUN=sd,na.rm=TRUE)))
}
colnames(HeighTr)<-c("expmean","lowmean","expsd","lowsd")
HeighTr$year<-c(2009:2012,2014:2017)

plot(HeighTr$expmean~HeighTr$year,ylim=c(0,70),type="b",las=1,
     xaxt="n",yaxt="n",bty="n",col=colovec[1],
     ylab="Seedling height (cm)",xlab="Year",
     font.lab=2,cex.lab=1.5,lwd=2,pch=22,bg="white")
plotCI(x=HeighTr$year,y=HeighTr$expmean,
       ui=HeighTr$expmean,
       li=HeighTr$expmean-HeighTr$expsd,
       pch=NA,lwd=1.5,add=TRUE)
plotCI(x=HeighTr$year,y=HeighTr$lowmean,
       ui=HeighTr$lowmean+HeighTr$lowsd,
       li=HeighTr$lowmean,
       pch=NA,lwd=1.5,add=TRUE)
points(HeighTr$expmean~HeighTr$year,type="b",
       col=colovec[1],lwd=2.5,pch=22,bg="white",cex=1.5)
points(HeighTr$lowmean~HeighTr$year,type="b",lty=2,
       col=colovec[2],lwd=2.5,pch=19,cex=1.5)
axis(1,lwd=2,cex.axis=1,at=c(2009:2017),
     lab=colnames(survEv),las=1,font=2)
axis(2,lwd=2,las=1,font=2)
box(lwd=2)
legend(2009.5,10,legend=c("Natural treatment","Protected treatment"),
       col=colovec,lty=c(1,2),pch=c(22,19),pt.bg="white",
       bty="n",cex=1.1,title.cex=1.5,y.intersp=1.8,x.intersp=0.5,lwd=2)
#export to .pdf 6 x 6 inches


##############################################################################/
#Figure combining the 3 plot of phenotypic traits evolution: portrait####
##############################################################################/


op<-par(mfrow=c(3,1),mar=c(4.1,5.1,2.1,1.1))
#Figure A
matplot(t(survEv),type="b",ylim=c(0,100),las=1,
        col=colovec[as.numeric(survEv.treat)],lwd=2,
        lty=as.numeric(survEv.treat),bg="white",
        pch=c(22,19,22,22,22,19,19,22,22),
        xaxt="n",yaxt="n",bty="n",
        ylab="Percentage of survival",xlab="Year",
        font.lab=2,cex.lab=2,cex=1.8)
axis(1,lwd=2,cex.axis=1.5,at=c(1:dim(survEv)[2]),
     lab=colnames(survEv),las=1,font=2)
axis(2,lwd=2,las=1,font=2,cex.axis=1.5)
box(lwd=2)
legend(1.5,45,legend=c("Natural treatment","Protected treatment"),
       col=colovec,lty=c(1,2),pch=c(22,19),pt.bg="white",
       bty="n",cex=1.8,y.intersp=0.65,x.intersp=0.5,lwd=2)


#Figure B
graf<-barplot(temp1$mean,density=c(-1,20),angle=c(0,60),border=NA,
              col=colovec,ylim=c(-1,105),ylab="Percentage of infection",
              xlab="Year",cex.axis=1,cex.lab=2,las=1,yaxt="n",bty="n",
              space=c(0,0.25,rep(c(1,0.25),times=8)),font.lab=2)
axis(1,at=graf[seq(1,17,2)]+0.625,labels=rownames(infectM),
     lwd=2,font=2,cex.axis=1.5)
axis(2,at=seq(0,100,20),labels=seq(0,100,20),lwd=2,las=1,
     font=2,cex.axis=1.5)
box(bty="o",lwd=2)
#adding standard deviation
temp2<-pivot_longer(infectM,cols=3:4,names_to="treat",values_to="sd")
plotCI(x=graf,y=temp1$mean,
       ui=temp1$mean+temp2$sd,
       li=temp1$mean,
       pch=NA,lwd=1.5,add=TRUE)
#for cosmetic reason, we superimpose again the barplot
par(lwd=2)
barplot(temp1$mean,density=c(-1,20),angle=c(0,60),border=colovec,
        col=colovec,ylim=c(-1,105),ylab="Percentage of infection",
        cex.axis=1,cex.lab=1.5,las=1,yaxt="n",bty="n",
        space=c(0,0.25,rep(c(1,0.25),times=8)),font.lab=2,add=TRUE)
legend(-1,105,legend=c("Natural treatment","Protected treatment"),
       fill=colovec,density=c(-1,20),angle=c(0,60),bty="n",border=colovec,
       cex=1.5,y.intersp=0.8,x.intersp=0.5)
par(lwd=1)

#Figure C
plot(HeighTr$expmean~HeighTr$year,ylim=c(0,70),type="b",las=1,
     xaxt="n",yaxt="n",bty="n",col=colovec[1],
     ylab="Seedling height (cm)",xlab="Year",
     font.lab=2,cex.lab=2,lwd=2,pch=22,bg="white")
plotCI(x=HeighTr$year,y=HeighTr$expmean,
       ui=HeighTr$expmean,
       li=HeighTr$expmean-HeighTr$expsd,
       pch=NA,lwd=1.5,add=TRUE)
plotCI(x=HeighTr$year,y=HeighTr$lowmean,
       ui=HeighTr$lowmean+HeighTr$lowsd,
       li=HeighTr$lowmean,
       pch=NA,lwd=1.5,add=TRUE)
points(HeighTr$expmean~HeighTr$year,type="b",
       col=colovec[1],lwd=2.5,pch=22,bg="white",cex=2)
points(HeighTr$lowmean~HeighTr$year,type="b",lty=2,
       col=colovec[2],lwd=2.5,pch=19,cex=2)
axis(1,lwd=2,cex.axis=1.5,at=c(2009:2017),
     lab=colnames(survEv),las=1,font=2)
axis(2,lwd=2,las=1,font=2,cex.axis=1.5)
box(lwd=2)
legend(2009.8,18,legend=c("Natural treatment","Protected treatment"),
       col=colovec,lty=c(1,2),pch=c(22,19),pt.bg="white",
       bty="n",cex=1.8,y.intersp=0.65,x.intersp=0.5,lwd=2)
par(op)
#export to .pdf 6 x 14 portrait


##############################################################################/
#Figure combining the 3 plot of phenotypic traits evolution: landscape####
##############################################################################/

pdf(file="output/truc.pdf",width=18,height=6)
op<-par(mfrow=c(1,3),mar=c(4.1,5.1,2.1,1.1))
#Figure A
matplot(t(survEv),type="b",ylim=c(0,100),las=1,
        col=colovec[as.numeric(survEv.treat)],lwd=2,
        lty=as.numeric(survEv.treat),bg="white",
        pch=c(22,19,22,22,22,19,19,22,22),
        xaxt="n",yaxt="n",bty="n",
        ylab="Percentage of survival",xlab="Year",
        font.lab=2,cex.lab=2,cex=1.8)
axis(1,lwd=2,cex.axis=1.5,at=c(1:dim(survEv)[2]),
     lab=colnames(survEv),las=1,font=2)
axis(2,lwd=2,las=1,font=2,cex.axis=1.5)
box(lwd=2)
legend(1.5,30,legend=c("Natural treatment","Protected treatment"),
       col=colovec,lty=c(1,2),pch=c(22,19),pt.bg="white",
       bty="n",cex=2,y.intersp=1.3,x.intersp=0.8,lwd=2,seg.len=2)

#Figure B
graf<-barplot(temp1$mean,density=c(-1,20),angle=c(0,60),border=NA,
              col=colovec,ylim=c(-1,105),ylab="Percentage of infection",
              xlab="Year",cex.axis=1,cex.lab=2,las=1,yaxt="n",bty="n",
              space=c(0,0.25,rep(c(1,0.25),times=8)),font.lab=2)
axis(1,at=graf[seq(1,17,2)]+0.625,labels=rownames(infectM),
     lwd=2,font=2,cex.axis=1.5)
axis(2,at=seq(0,100,20),labels=seq(0,100,20),lwd=2,las=1,
     font=2,cex.axis=1.5)
box(bty="o",lwd=2)
#adding standard deviation
temp2<-pivot_longer(infectM,cols=3:4,names_to="treat",values_to="sd")
plotCI(x=graf,y=temp1$mean,
       ui=temp1$mean+temp2$sd,
       li=temp1$mean,
       pch=NA,lwd=1.5,add=TRUE)
#for cosmetic reason, we superimpose again the barplot
par(lwd=2)
barplot(temp1$mean,density=c(-1,20),angle=c(0,60),border=colovec,
        col=colovec,ylim=c(-1,105),ylab="Percentage of infection",
        cex.axis=1,cex.lab=1.5,las=1,yaxt="n",bty="n",
        space=c(0,0.25,rep(c(1,0.25),times=8)),font.lab=2,add=TRUE)
legend(0,98,legend=c("Natural treatment","Protected treatment"),
       fill=colovec,density=c(-1,20),angle=c(0,60),bty="n",border=colovec,
       cex=2,y.intersp=1.3,x.intersp=0.8)
par(lwd=1)

#Figure C
plot(HeighTr$expmean~HeighTr$year,ylim=c(0,70),type="b",las=1,
     xaxt="n",yaxt="n",bty="n",col=colovec[1],
     ylab="Seedling height (cm)",xlab="Year",
     font.lab=2,cex.lab=2,lwd=2,pch=22,bg="white")
plotCI(x=HeighTr$year,y=HeighTr$expmean,
       ui=HeighTr$expmean,
       li=HeighTr$expmean-HeighTr$expsd,
       pch=NA,lwd=1.5,add=TRUE)
plotCI(x=HeighTr$year,y=HeighTr$lowmean,
       ui=HeighTr$lowmean+HeighTr$lowsd,
       li=HeighTr$lowmean,
       pch=NA,lwd=1.5,add=TRUE)
points(HeighTr$expmean~HeighTr$year,type="b",
       col=colovec[1],lwd=2.5,pch=22,bg="white",cex=2)
points(HeighTr$lowmean~HeighTr$year,type="b",lty=2,
       col=colovec[2],lwd=2.5,pch=19,cex=2)
axis(1,lwd=2,cex.axis=1.5,at=c(2009:2017),
     lab=colnames(survEv),las=1,font=2)
axis(2,lwd=2,las=1,font=2,cex.axis=1.5)
box(lwd=2)
legend(2009.5,11,legend=c("Natural treatment","Protected treatment"),
       col=colovec,lty=c(1,2),pch=c(22,19),pt.bg="white",
       bty="n",cex=2,y.intersp=1.3,x.intersp=0.8,lwd=2,seg.len=2)
par(op)
#export to .pdf 6 x 18 landscape
dev.off()


##############################################################################/
#END
##############################################################################/