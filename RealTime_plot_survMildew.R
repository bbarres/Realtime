##############################################################################/
##############################################################################/
#Figure 3: modelisation of the survival as a function of mildew infection
##############################################################################/
##############################################################################/

#loading the necessary packages and data sets
source("RealTime_load.R")


##############################################################################/
#modelisation of the survival####
##############################################################################/

temp<-AllTrait
colnames(temp)[8]<-"Powdery_mildew"
colnames(temp)[10]<-"Acorn_weight"
survOid<-glm(Survival~Acorn_weight+Powdery_mildew,data=temp,family=binomial)


##############################################################################/
#plot using visreg package####
##############################################################################/

colovec<-c(brewer.pal(12,"Set3")[6:7],
           brewer.pal(9,"Set1")[1:2])
pdf(file="output/Figure_3_GLMmildew.pdf",width=7,height=5)
op<-par(mar=c(5.1,4.1,1.1,1.1))
visreg(survOid,"Powdery_mildew",scale="response",rug=2,ylim=c(0,1),
       xlab="Mean infection (2009-2013) (%)",
       ylab="Survival (2017)")
rug(temp[temp$Survival==0,]$Powdery_mildew,side=1,col=colovec[1])
rug(temp[temp$Survival==1,]$Powdery_mildew,side=3,col=colovec[2])
box()
par(op)
#export to pdf 7 x 5 inches
dev.off()


##############################################################################/
#END
##############################################################################/