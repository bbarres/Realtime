##############################################################################/
##############################################################################/
#Figure S4: Sorting individuals by general genotyping quality
##############################################################################/
##############################################################################/


#loading and preparing the additionnal dataset
genoQual<-read.table("data/dataSup/Ind_Geno_Qual_ST.txt",
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
#adding a column with different notes for different quality of genotyping
genoQual$qnote<-0
genoQual$qnote<-ifelse(genoQual$Call.Rate<seuilCRatemin,
                       genoQual$qnote+1,
                       genoQual$qnote+0)
genoQual$qnote<-ifelse(genoQual$p10.GC<seuilp10min,
                       genoQual$qnote+1,
                       genoQual$qnote+0)
genoQual$qnote<-ifelse(genoQual$p50.GC<seuilp50min,
                       genoQual$qnote+1,
                       genoQual$qnote+0)

#defining a color vector
coVec<-c(rgb(0,0,0,max=255,alpha=40),
         rgb(255,255,51,max=255,alpha=200),
         rgb(255,127,0,max=255,alpha=200),
         rgb(228,26,28,max=255,alpha=200))


##############################################################################/
#Figure S4: plotting the quality of individuals genotyped with SNP####
##############################################################################/

#this code will produce a pdf file in the output folder
pdf(file="output/Figure_S4_QualInd.pdf",width=11,height=6)
op<-par(mfrow=c(1,2))
plot(p10.GC~Call.Rate,xlim=xrange,ylim=yrange10,
     data=genoQual[genoQual$p50.GC>=seuilp50min & 
                     genoQual$p10.GC>=seuilp10min &
                     genoQual$Call.Rate>=seuilCRatemin,],
     xlab="Call Rate",ylab="10% GenCall",
     col=coVec[1],pch=1,cex=0.6,las=1)
abline(h=c(seuilp10min),lty=2,col=brewer.pal(9,"Set1")[2],lwd=3)
abline(v=c(seuilCRatemin),lty=2,col=brewer.pal(9,"Set1")[3],lwd=3)
points(p10.GC~Call.Rate,pch=21,cex=0.7,
       bg=coVec[qnote+1],col=rgb(0,0,0,max=255,alpha=200),
       data=genoQual[genoQual$qnote!=0,])
plot(p50.GC~Call.Rate,xlim=xrange,ylim=yrange50,
     data=genoQual[genoQual$p50.GC>=seuilp50min & 
                     genoQual$p10.GC>=seuilp10min &
                     genoQual$Call.Rate>=seuilCRatemin,],
     xlab="Call Rate",ylab="50% GenCall",
     col=coVec[1],pch=1,cex=0.6,las=1)
abline(h=c(seuilp50min),lty=2,col=brewer.pal(9,"Set1")[4],lwd=3)
abline(v=c(seuilCRatemin),lty=2,col=brewer.pal(9,"Set1")[3],lwd=3)
points(p50.GC~Call.Rate,pch=21,cex=0.7,
       bg=coVec[qnote+1],col=rgb(0,0,0,max=255,alpha=200),
       data=genoQual[genoQual$qnote!=0,])
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