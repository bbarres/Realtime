##############################################################################/
##############################################################################/
#Progeny survival in Protected vs Natural plot and PM infection level
##############################################################################/
##############################################################################/

source("RealTime_load.R")


##############################################################################/
#Figure 5: progeny survival by exposure and powdery mildew infection####
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

#computing the survival rate by family by exposure
SurvProp<-as.data.frame(prop.table(table(progSurv$treat,
                                         progSurv$DoA,
                                         progSurv$family_simp),
                                   margin=c(1,3))*100)
SurvProp<-SurvProp[SurvProp$Var2=="alive",]
SurvProp<-pivot_wider(SurvProp[,-c(2)],names_from=Var1,values_from=Freq)
colnames(SurvProp)<-c("Progeny","Natural","Protected")
SurvProp$PMinf<-tapply(progSurv$oid_moy,INDEX=progSurv$family_simp,
                       FUN=mean,na.rm=TRUE)

#ploting the results
pdf(file="output/Figure_5_SurvByTreat.pdf",width=7,height=7)
op<-par(mar=c(5.1,4.3,1,1))
coloor<-brewer.pal(11,"RdYlGn")[c(11:9,7,5,3:1)] #with green to red gradient
#another color gradient brewer.pal(9,"YlOrRd") #with yellow to red gradient
nbbrak<-c(20,21,22,23,24,25,26,27,28)
plot(SurvProp[,c(3,2)],xlim=c(0,100),ylim=c(0,100),bty="n",
     las=1,pch=21,cex=3,col="black",xaxt="n",yaxt="n",
     bg=coloor[as.numeric(cut(SurvProp$PMinf,breaks=nbbrak))],
     ylab="Natural exposure survival rate (%)",
     xlab="Protected exposure survival rate (%)",
     font.lab=2,cex.lab=1.5)
text(SurvProp[,c(3,2)],labels=SurvProp$Progeny)
abline(0,1,lty=5,col=grey(0.6,1),lwd=3)
legend(8,98,legend=c("]20-21]","]21-22]","]22-23]","]23-24]",
                     "]24-25]","]25-26]","]26-27]","]27-28]"),
       bty="n",fill=coloor,title="Mean infection\n(2009-2013)")
axis(1,lwd=2,cex.axis=1,las=1,font=2)
axis(2,lwd=2,las=1,font=2)
box(lwd=2)
par(op)
dev.off()


##############################################################################/
#END
##############################################################################/