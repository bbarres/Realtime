##############################################################################/
##############################################################################/
#Evolution of diversity and mortality rate by families
##############################################################################/
##############################################################################/

source("RealTime_load.R")


##############################################################################/
#yearly evolution of the shannon diversity index by bloc####
##############################################################################/

#limiting the data set to useful data
evolShan<-dispo[dispo$family_simp!="CC" & dispo$family_simp!="hd" & 
                  !is.na(dispo$an_mort) & dispo$an_mort!="g" & 
                dispo$family_simp!="26",
                c("Sample_ID","bloc","PU","treat","an_mort","family_simp")]
evolShan$grpbloc<-paste(evolShan$bloc,evolShan$PU,evolShan$treat,sep="")
levels(evolShan$family_simp)[1:2]<-c("01","09")
evolShan<-drop.levels(evolShan)
#defining a color vector
colovec<-c(brewer.pal(12,"Paired")[4],brewer.pal(12,"Paired")[2])

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

pdf(file="output/Figure_ShannonEvol.pdf",width=7,height=6)
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
legend(1,3.4,legend="Natural treatment",
       col=colovec[1],lty=1,pch=22,pt.bg="white",
       bty="n",cex=1.3,title.cex=1,y.intersp=1.5,x.intersp=0.5,lwd=2)
legend(1,3.35,legend="Protected treatment",
       col=colovec[2],lty=2,pch=19,pt.bg="white",
       bty="n",cex=1.3,title.cex=1,y.intersp=1.5,x.intersp=0.5,lwd=2)
par(op)
#export to .pdf 7 x 6 inches
dev.off()


##############################################################################/
#progeny survival by treatment and powdery mildew infection####
##############################################################################/

#preparing the data set
progSurv<-dispo[dispo$family_simp!="CC" & dispo$family_simp!="hd" & 
                  !is.na(dispo$an_mort) & dispo$an_mort!="g" & 
                  dispo$family_simp!="26",
                c("Sample_ID","bloc","PU","family_simp",
                  "treat","an_mort","oid_moy")]
levels(progSurv$family_simp)[1:2]<-c("01","09")
progSurv<-drop.levels(progSurv)
progSurv$DoA<-progSurv$an_mort
progSurv[progSurv$DoA!="vivant",]$DoA<-"dead"
progSurv[progSurv$DoA=="vivant",]$DoA<-"alive"
progSurv$DoA<-as.factor(progSurv$DoA)
progSurv$treat<-as.factor(progSurv$treat)

#computing the survival rate by family by treatment
SurvProp<-as.data.frame(prop.table(table(progSurv$treat,
                                         progSurv$DoA,
                                         progSurv$family_simp),
                     margin=c(1,3))*100)
SurvProp<-SurvProp[SurvProp$Var2=="alive",]
SurvProp<-pivot_wider(SurvProp[,-c(2)],names_from=Var1,values_from=Freq)
colnames(SurvProp)<-c("Progeny","Natural","Protected")
SurvProp$PMinf<-tapply(progSurv$oid_moy,INDEX=progSurv$family_simp,
                       FUN=mean,na.rm=TRUE)

#ploting the results (version 1)
pdf(file="output/Figure5_SurvByProg.pdf",width=7,height=7)
op<-par(mar=c(5.1,4.3,1,1))
coloor<-brewer.pal(9,"YlOrRd") #with yellow to red gradient
nbbrak<-c(20,21,22,23,24,25,26,27,28)
plot(SurvProp[,c(3,2)],xlim=c(0,100),ylim=c(0,100),bty="n",
     las=1,pch=21,cex=3,col="black",xaxt="n",yaxt="n",
     bg=coloor[as.numeric(cut(SurvProp$PMinf,breaks=nbbrak))],
     ylab="Natural treatment survival rate (%)",
     xlab="Protected treatment survival rate (%)",
     font.lab=2,cex.lab=1.5)
text(SurvProp[,c(3,2)],labels=SurvProp$Progeny)
abline(0,1,lty=5,col=grey(0.6,1),lwd=3)
legend(8,98,legend=c("]20-21]","]21-22]","]22-23]","]23-24]",
                    "]24-25]","]25-26]","]26-27]","]27-28]"),
       bty="n",fill=coloor,title="Mean\nPM infection")
axis(1,lwd=2,cex.axis=1,las=1,font=2)
axis(2,lwd=2,las=1,font=2)
box(lwd=2)
par(op)
dev.off()

#ploting the results (version 2)
pdf(file="output/Figure5_SurvByProg2.pdf",width=7,height=7)
op<-par(mar=c(5.1,4.3,1,1))
coloor<-brewer.pal(11,"RdYlGn")[11:1] #with green to red gradient
nbbrak<-11 #depend on the coloor chosen
as.numeric(cut(SurvProp$PMinf,breaks=nbbrak))
plot(SurvProp[,c(3,2)],xlim=c(0,100),ylim=c(0,100),bty="n",
     las=1,pch=21,cex=3,col="black",xaxt="n",yaxt="n",
     bg=coloor[as.numeric(cut(SurvProp$PMinf,breaks=nbbrak))],
     ylab="Natural treatment survival rate (%)",
     xlab="Protected treatment survival rate (%)",
     font.lab=2,cex.lab=1.5)
text(SurvProp[,c(3,2)],labels=SurvProp$Progeny)
abline(0,1,lty=5,col=grey(0.6,1),lwd=3)
legend(8,98,legend=levels(cut(SurvProp$PMinf,breaks=nbbrak)),
       bty="n",fill=coloor,title="Mean\nPM infection")
axis(1,lwd=2,cex.axis=1,las=1,font=2)
axis(2,lwd=2,las=1,font=2)
box(lwd=2)
par(op)
dev.off()


##############################################################################/
#END
##############################################################################/