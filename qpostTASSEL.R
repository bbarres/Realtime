###############################################################################
###############################################################################
#Post treatment of the TASSEL output file
###############################################################################
###############################################################################

#loading the libraries
library(qvalue)

#set the working directory
setwd("~/work/Rfichiers/Githuber/Realtime_data")


###############################################################################
#turn the p-value into q-value
###############################################################################

#here is a function to compute q-value. It takes the output file of TASSEL
#as an input. The second argument of the function is the q threshold chosen

qpostassel<-function(stat,qalpha) 
{
  stat<-stat[order(stat$Trait),]
  statNA<-stat[is.na(stat$p)==FALSE,]
  traitID<-levels(statNA$Trait)
  nbT<-as.numeric(length(traitID))
  colSup<-c()
  for (i in 1:nbT) {
    TR<-statNA[statNA$Trait==traitID[i],]
    qval<-qvalue(TR$p,fdr.level=qalpha)
    novCol<-(rep(traitID[i],dim(as.data.frame(qval$qvalues))[1]))
    novCol<-cbind(novCol,TR$Marker,as.data.frame(qval$qvalues)
                  ,as.data.frame(qval$significant))
    colSup<-rbind(colSup,novCol)
  }
  statFin<-cbind(statNA,colSup)
  names(statFin)[(dim(statFin)[2]-2):(dim(statFin)[2])]<-
    c("Marker","qval","qsignif")
  return(statFin)
}

NATRENF_stat<-read.table("rez_NATRENF_stats.txt",sep="\t",header=T)
q_NAT_stat<-qpostassel(NATRENF_stat,0.05)


###############################################################################
#code for a qq-plot
###############################################################################

unif<-runif(table(q_NAT_stat$Trait)[[1]],0,1)

pdf(file="plot_NAT_Trait_QQ.pdf",width=40, height=10)
op<-par(mfrow=c(2,8))
for (i in 1:length(levels(q_NAT_stat$Trait))) {
  range<-q_NAT_stat[q_NAT_stat$Trait==(levels(q_NAT_stat$Trait))[i],]
  plot(pty="s",type="o",
       sort(-log(unif[1:dim(q_NAT_stat[q_NAT_stat$Trait
                                    ==(levels(q_NAT_stat$Trait))[i],])[1]])),
       sort(-log(range$p)),
       col=as.factor(range$qsignif)[order(-log(range$p))],
       ylim=c(0,max(-log(range$p))),
       xaxt="s",yaxt="s", xlab="", ylab="",pch=21,cex=2,lwd=3,
       main=levels(q_NAT_stat$Trait)[i],cex.main=2.5)
  abline(0,1,lty=3,lwd=2)
}
par(op)
dev.off()


###############################################################################
#END
###############################################################################