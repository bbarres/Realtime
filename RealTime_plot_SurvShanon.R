##############################################################################/
##############################################################################/
#Evolution of Shannon index by families in each plot
##############################################################################/
##############################################################################/

source("RealTime_load.R")


##############################################################################/
#Computation of the Shannon index####
##############################################################################/

#limiting the data set to useful data
evolShan<-dispo[dispo$family_simp!="CC" & dispo$family_simp!="hd" & 
                  !is.na(dispo$an_mort) & dispo$an_mort!="g" & 
                  dispo$family_simp!="26",
                c("Sample_ID","bloc","PU","treat","an_mort","family_simp")]
evolShan$grpbloc<-paste(evolShan$bloc,evolShan$PU,evolShan$treat,sep="")
levels(evolShan$family_simp)[1:2]<-c("01","09")
evolShan<-drop.levels(evolShan)

#computing the Shannon index
temp<-table(evolShan$family_simp,evolShan$grpbloc)
Shanplot<-apply(temp,2,diversity,base=2)
temp<-table(evolShan[evolShan$an_mort>2009,]$family_simp,
            evolShan[evolShan$an_mort>2009,]$grpbloc)
Shanplot<-rbind(Shanplot,apply(temp,2,diversity,base=2))
temp<-table(evolShan[evolShan$an_mort>2010,]$family_simp,
            evolShan[evolShan$an_mort>2010,]$grpbloc)
Shanplot<-rbind(Shanplot,apply(temp,2,diversity,base=2))
temp<-table(evolShan[evolShan$an_mort>2011,]$family_simp,
            evolShan[evolShan$an_mort>2011,]$grpbloc)
Shanplot<-rbind(Shanplot,apply(temp,2,diversity,base=2))
temp<-table(evolShan[evolShan$an_mort>2012,]$family_simp,
            evolShan[evolShan$an_mort>2012,]$grpbloc)
Shanplot<-rbind(Shanplot,apply(temp,2,diversity,base=2))
temp<-table(evolShan[evolShan$an_mort>2013,]$family_simp,
            evolShan[evolShan$an_mort>2013,]$grpbloc)
Shanplot<-rbind(Shanplot,apply(temp,2,diversity,base=2))
temp<-table(evolShan[evolShan$an_mort>2014,]$family_simp,
            evolShan[evolShan$an_mort>2014,]$grpbloc)
Shanplot<-rbind(Shanplot,apply(temp,2,diversity,base=2))
temp<-table(evolShan[evolShan$an_mort>2015,]$family_simp,
            evolShan[evolShan$an_mort>2015,]$grpbloc)
Shanplot<-rbind(Shanplot,apply(temp,2,diversity,base=2))
temp<-table(evolShan[evolShan$an_mort>2016,]$family_simp,
            evolShan[evolShan$an_mort>2016,]$grpbloc)
Shanplot<-rbind(Shanplot,apply(temp,2,diversity,base=2))
temp<-table(evolShan[evolShan$an_mort>2017,]$family_simp,
            evolShan[evolShan$an_mort>2017,]$grpbloc)
Shanplot<-rbind(Shanplot,apply(temp,2,diversity,base=2))
row.names(Shanplot)<-c("Emerging","2009","2010","2011",
                       "2012","2013","2014","2015","2016",
                       "2017")


##############################################################################/
#Figure 7: yearly evolution of the shannon diversity index by bloc####
##############################################################################/

#defining a color vector
colovec<-c(brewer.pal(12,"Paired")[4],brewer.pal(12,"Paired")[2])

#code to generate the figure
pdf(file="output/Figure_7_ShannonEvol.pdf",width=7,height=6)
op<-par(mar=c(5.1,4.3,1.1,1.1))
matplot(Shanplot,type="b",ylim=c(3.25,3.85),las=1,
        col=colovec[c(1,2,1,1,1,2,2,1,1)],lwd=2,
        lty=c(1,2,1,1,1,2,2,1,1),
        pch=c(22,19,22,22,22,19,19,22,22),
        xaxt="n",yaxt="n",bty="n",
        ylab="Shannon index",xlab="Year",
        font.lab=2,cex.lab=1.5)
axis(1,lwd=2,cex.axis=0.9,at=c(1:10),
     lab=c("Emerg.","2009","2010","2011",
           "2012","2013","2014","2015","2016",
           "2017"),las=1,font=2)
axis(2,lwd=2,las=1,font=2)
box(lwd=2)
legend(1,3.4,legend="Natural exposure",
       col=colovec[1],lty=1,pch=22,pt.bg="white",
       bty="n",cex=1.3,title.cex=1,y.intersp=1.5,x.intersp=0.5,lwd=2)
legend(1,3.35,legend="Protected exposure",
       col=colovec[2],lty=2,pch=19,pt.bg="white",
       bty="n",cex=1.3,title.cex=1,y.intersp=1.5,x.intersp=0.5,lwd=2)
par(op)
#export to .pdf 7 x 6 inches
dev.off()


##############################################################################/
#END
##############################################################################/