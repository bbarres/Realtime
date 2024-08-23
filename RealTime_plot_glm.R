##############################################################################/
##############################################################################/
#Figure 3: modelisation of the survival as a function of mildew infection
##############################################################################/
##############################################################################/

#loading the necessary packages and data sets
source("RealTime_load.R")

#preparing the dataset
temp<-RTdata[RTdata$family_simp!="PAR" & RTdata$family_simp!="CC" &
               RTdata$family_simp!="hd" & !is.na(RTdata$statut10),]
temp<-temp[,c("treat","Sample_ID","family_simp","oid4_09","oid5_10",
              "oid5_11","oid4_12","oid2_13","oid_moy",
              "H09v","H12v","pgland","gel_13","statut10")]
#rename the columns with more accurate names
colnames(temp)<-c("Treatment","Taxa","family_simp","oid4_09","oid5_10",
                  "oid5_11","oid4_12","oid2_13","Powdery_mildew","Height",
                  "H09v","Acorn_weight","Frost_damage","Survival")
temp$Treatment<-as.factor(temp$Treatment)
levels(temp$Treatment)<-c("Natural","Protected")


##############################################################################/
#Model 1: effect of acorn weight and exposure####
##############################################################################/

survTreatAndAcorn<-glm(Survival~Treatment+Acorn_weight,
                       data=temp,family=binomial)
summary(survTreatAndAcorn)


##############################################################################/
#Figure S10/A20: plot using visreg for acorn weight effect visualization####
##############################################################################/

pdf(file="output/Figure_S10_GLMacorn.pdf",width=7,height=5)
op<-par(mar=c(4.3,4.1,1.9,1.1))
#defining a color vector
colovec<-c(brewer.pal(12,"Paired")[1:2],brewer.pal(12,"Paired")[3:4],
           brewer.pal(11,"RdYlGn")[1])
visreg(survTreatAndAcorn,"Acorn_weight",scale="response",
       rug=FALSE,ylim=c(0,1),
       by="Treatment",overlay=TRUE,
       xlab="Acorn weight (g)",
       ylab="Survival (2017)",
       line=list(col=colovec[c(4,2)]),
       fill=list(col=colovec[c(3,1)]))
#defining a color vector
colovec<-c(brewer.pal(12,"Set3")[6:7],
           brewer.pal(9,"Set1")[1:2])
rug(jitter(temp[temp$Survival==0,]$Acorn_weight),side=1,col=colovec[1])
rug(jitter(temp[temp$Survival==1,]$Acorn_weight),side=3,col=colovec[2])
box()
par(op)
#export to pdf 7 x 5 inches
dev.off()


##############################################################################/
#Model 2: effect of Powdery mildew mean infection and acorn weight####
##############################################################################/

survOid<-glm(Survival~Acorn_weight+Powdery_mildew,
             data=temp,family=binomial)
summary(survOid)


##############################################################################/
#Figure 3: plot using visreg for mildew infection effect visualization####
##############################################################################/

pdf(file="output/Figure_3_GLMmildew.pdf",width=7,height=5)
op<-par(mar=c(5.1,4.1,1.1,1.1))
#defining a color vector
colovec<-c(brewer.pal(12,"Set3")[6:7],
           brewer.pal(9,"Set1")[1:2])
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
#Model 3 for plot: effect of Height in 2009 and exposure####
##############################################################################/

survH09<-glm(Survival~H09v*Treatment,
             data=temp,family=binomial)
summary(survH09)


##############################################################################/
#Figure S11/A21: plot using visreg Height in 2009 effect visualization ####
##############################################################################/

pdf(file="output/Figure_S11_GLMH09.pdf",width=7,height=5)
op<-par(mar=c(4.3,4.1,1.9,1.1))
#defining a color vector
colovec<-c(brewer.pal(12,"Paired")[1:2],brewer.pal(12,"Paired")[3:4],
           brewer.pal(11,"RdYlGn")[1])
visreg(survH09,"H09v",scale="response",xlim=c(0,110),
       rug=2,ylim=c(0,1),by="Treatment",overlay=TRUE,
       xlab="Height in 2009 (cm)",
       ylab="Survival (2017)",
       line=list(col=colovec[c(4,2)]),
       fill=list(col=colovec[c(3,1)]))
colovec<-c(brewer.pal(12,"Set3")[6:7],
           brewer.pal(9,"Set1")[1:2])
rug(temp[temp$Survival==0,]$H09v,side=1,col=colovec[1])
rug(temp[temp$Survival==1,]$H09v,side=3,col=colovec[2])
box()
par(op)
#export to pdf 7 x 5 inches
dev.off()


##############################################################################/
#Model 4: effect of Frost damage####
##############################################################################/

survFrost<-glm(Survival~Acorn_weight+Powdery_mildew+Frost_damage,
               data=temp,family=binomial)
summary(survFrost)


##############################################################################/
#Figure S12/A22: plot using visreg package for Frost damage visualization####
##############################################################################/

pdf(file="output/Figure_S12_GLMFrost.pdf",width=7,height=5)
op<-par(mar=c(4.3,4.1,1.9,1.1))
#defining a color vector
colovec<-c(brewer.pal(12,"Paired")[1:2],brewer.pal(12,"Paired")[3:4],
           brewer.pal(11,"RdYlGn")[1])
visreg(survFrost,"Powdery_mildew",scale="response",
       rug=2,ylim=c(0,1),by="Frost_damage",overlay=TRUE,
       xlab="Mean infection (2009-2013) (%)",
       ylab="Survival (2017)")
colovec<-c(brewer.pal(12,"Set3")[6:7],
           brewer.pal(9,"Set1")[1:2])
rug(temp[temp$Survival==0,]$Powdery_mildew,side=1,col=colovec[1])
rug(temp[temp$Survival==1,]$Powdery_mildew,side=3,col=colovec[2])
box()
par(op)
#export to pdf 7 x 5 inches
dev.off()


##############################################################################/
#Model 5: full model####
##############################################################################/

survFull<-glm(Survival~Acorn_weight+Treatment*family_simp+Frost_damage,
              data=temp,family=binomial)
summary(survFull)
anova(survFull,test="Chisq")


##############################################################################/
#END
##############################################################################/