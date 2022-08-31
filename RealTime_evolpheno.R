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
        pch=as.numeric(survEv.bloc)+14,
        xaxt="n",yaxt="n",bty="n",
        ylab="Percentage of survival",xlab="Year",
        font.lab=2,cex.lab=1.5)
axis(1,lwd=2,cex.axis=1,at=c(1:dim(survEv)[2]),
     lab=colnames(survEv),las=1,font=2)
axis(2,lwd=2,las=1,font=2)
box(lwd=2)
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
graf<-barplot(temp1$mean,
              col=colovec,ylim=c(-1,105),ylab="Percentage of infection",
              cex.axis=1,cex.lab=1.5,las=1,yaxt="n",bty="n",border=NA,
              space=c(0.1,0.1,rep(c(1,0.1),times=8)),font.lab=2)
axis(1,at=graf[seq(1,17,2)]+0.55,labels=rownames(infectM),
     lwd=2,font=2,cex.axis=0.9)
axis(2,at=seq(0,100,20),labels=seq(0,100,20),lwd=2,las=1,
     font=2,cex.axis=1)
box(bty="o",lwd=2)
#adding standard deviation
temp2<-pivot_longer(infectM,cols=3:4,names_to="treat",values_to="sd")
plotCI(x=graf,y=temp1$mean,
       ui=temp1$mean+temp2$sd,
       li=temp1$mean,
       pch=NA,lwd=2,add=TRUE)
#for cosmetic reason, we superimpose again the barplot
barplot(temp1$mean,
        col=colovec,ylim=c(-1,105),ylab="Percentage of infection",
        cex.axis=1,cex.lab=1.5,las=1,yaxt="n",bty="n",border=NA,
        space=c(0.1,0.1,rep(c(1,0.1),times=8)),font.lab=2,add=TRUE)

legend(-1,103,legend=c("Natural treatment",
                     "Protected treatment"),
       pch=15,col=colovec,bg=colovec,bty="n",cex=1.3,
       pt.cex=1.4,xpd=TRUE,y.intersp=1.5,x.intersp=0.5)


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



##############################################################################/
#END
##############################################################################/