##############################################################################/
##############################################################################/
#Figure S4: Sorting individuals by general genotyping quality
##############################################################################/
##############################################################################/

library(RColorBrewer)
#loading and preparing the data set
genoQual<-read.table("data/dataSup/RT_s_pre_ST.txt",
                     sep="\t",header=T,dec=".")
summary(genoQual)
colnames(genoQual)


##############################################################################/
#Preparing the data set####
##############################################################################/

#defining the threshold
seuilp50min<-(mean(genoQual$p50.GC,na.rm=T)-0.01)
seuilp10min<-(mean(genoQual$p10.GC,na.rm=T)-0.015)
seuilCRatemin<-(mean(genoQual$Call.Rate,na.rm=T)-0.01)
xrange<-c(min(genoQual$Call.Rate[genoQual$Call.Rate!=0],
              na.rm=T)-0.01,
          max(genoQual$Call.Rate[genoQual$Call.Rate!=0],
              na.rm=T)+0.01)
yrange10<-c(min(genoQual$p10.GC,na.rm=T)-0.01,
          max(genoQual$p10.GC,na.rm=T)+0.01)
yrange50<-c(min(genoQual$p50.GC,na.rm=T)-0.01,
          max(genoQual$p50.GC,na.rm=T)+0.01)
#defining a color vector
coVec<-c(rgb(0,0,0,max=255,alpha=30),
         rgb(228,26,26,max=255,alpha=150),
         rgb(55,126,184,max=255,alpha=150),
         rgb(77,175,74,max=255,alpha=150))


##############################################################################/
#plotting Figure S4####
##############################################################################/

#this code will produce a pdf file in the output folder
pdf(file="output/Figure_S4.pdf",width=11,height=6)
op<-par(mfrow=c(1,2))
plot(p10.GC~Call.Rate,xlim=xrange,ylim=yrange10,
     data=genoQual[genoQual$p50.GC>=seuilp50min & 
                     genoQual$p10.GC>=seuilp10min &
                     genoQual$Call.Rate>=seuilCRatemin,],
     col=coVec[1],pch=1,cex=0.6)
abline(h=c(seuilp10min),lty=2,col=brewer.pal(9,"Set1")[2],lwd=3)
abline(v=c(seuilCRatemin),lty=2,col=brewer.pal(9,"Set1")[3],lwd=3)
points((p10.GC[p10.GC<seuilp10min])~(Call.Rate[p10.GC<seuilp10min]),
       col=coVec[3],pch=1,cex=0.8,data=genoQual)
points((p10.GC[Call.Rate<seuilCRatemin & p10.GC>seuilp10min])~
         (Call.Rate[Call.Rate<seuilCRatemin & p10.GC>seuilp10min]),
       col=coVec[4],pch=1,cex=0.8,data=genoQual)

plot(p50.GC~Call.Rate,xlim=xrange,ylim=yrange50,
     data=genoQual[genoQual$p50.GC>=seuilp50min & 
                     genoQual$p10.GC>=seuilp10min &
                     genoQual$Call.Rate>=seuilCRatemin,],
     col=coVec[1],pch=1,cex=0.6)
abline(h=c(seuilp50min),lty=2,col=brewer.pal(9,"Set1")[1],lwd=3)
abline(v=c(seuilCRatemin),lty=2,col=brewer.pal(9,"Set1")[3],lwd=3)
points((p50.GC[p50.GC<seuilp50min])~(Call.Rate[p50.GC<seuilp50min]),
       col=coVec[2],pch=1,cex=0.8,data=genoQual)
points((p50.GC[Call.Rate<seuilCRatemin & p50.GC>seuilp50min])~
         (Call.Rate[Call.Rate<seuilCRatemin & p50.GC>seuilp50min]),
       col=coVec[4],pch=1,cex=0.8,data=genoQual)
par(op)
dev.off()

#list of individuals excluded because of poor global genotyping quality
excl<-(genoQual[genoQual$p50.GC<seuilp50min|
                  genoQual$p10.GC<seuilp10min|
                  genoQual$Call.Rate<seuilCRatemin,
             c("Sample.ID","Array.Info.Sentrix.ID",
               "Array.Info.Sentrix.Position")])


##############################################################################/
#END
##############################################################################/