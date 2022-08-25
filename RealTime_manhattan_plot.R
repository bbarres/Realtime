##############################################################################/
##############################################################################/
#Manhattan plot of the GWAS analyses results
##############################################################################/
##############################################################################/

source("RealTime_load.R")
#in order to load the results of the GWAS analyses, you need to run the 
#RealTime_GWAS.R script before this script


##############################################################################/
#function to import and prepare the output files from GAPIT analysis####
##############################################################################/

#defining a function to import and prepare the data set from GAPIT
#output file
readManDa<-function(pathtoGWASrez) {
  datMan<-read.csv(file=pathtoGWASrez)
  #we change the chromosome information into a factor with convenient
  #levels order
  datMan$Chromosome<-factor(datMan$Chromosome,
                               levels=c("1","2","3","4","5","6","7",
                                        "8","9","10","11","12","13",
                                        "MT","Pltd"))
  #there is only 12 chromosome in Q. robur, the 13 stand for an 
  #unlocated SNP on the genome at the time of the analysis. The level
  #"13" is therefore turn into "NA"
  levels(datMan$Chromosome)[13]<-"Na"
  levels(datMan$Chromosome)[14]<-"Mt"
  levels(datMan$Chromosome)[15]<-"Cp"
  #we turn Position from integer to numeric
  datMan$Position<-as.numeric(datMan$Position)
  datMan<-datMan[order(datMan$Chromosome,datMan$Position),]
  datMan$logp<- -log10(datMan$P.value)
  datMan$grpNum<-as.numeric(datMan$Chromosome)
  datMan$decal<-rep(c(0,cumsum(tapply(datMan$Position,
                                      datMan$Chromosome,max))[1:14]),
                    times=table(datMan$Chromosome))
  datMan$signiCat<-cut(datMan$FDR_Adjusted_P.values,
                       breaks=c(0,0.0001,0.001,0.01,1),
                       labels=c("<0.0001","<0.001","<0.01","ns"))
  return(datMan)
}


##############################################################################/
#function to plot a Manhattan plot####
##############################################################################/

ManhaPlot<-function(datMan,colovec,colosign="red",
                    decalCHR=0,desiXax=1,ylimi=c(0,10)) {
  datMan$posabsol<-datMan$decal+
    datMan$Position+
    (datMan$grpNum-1)*decalCHR
  plot(datMan$posabsol,datMan$logp,pch=19,las=1,font.axis=2,
       ylab=expression(-log[10](pval)),cex.lab=1.5,cex=1.0,
       xaxt="n",yaxt="n",bty="n",xlab="",ylim=ylimi,
       col=colovec[as.numeric(even(as.numeric(datMan$Chromosome)))+1])
  points(datMan[datMan$signiCat=="<0.01",]$posabsol,
         datMan[datMan$signiCat=="<0.01",]$logp,
         pch=21,bg=colosign[1],cex=1.5)
  points(datMan[datMan$signiCat=="<0.001",]$posabsol,
         datMan[datMan$signiCat=="<0.001",]$logp,
         pch=21,bg=colosign[2],cex=1.5)
  points(datMan[datMan$signiCat=="<0.0001",]$posabsol,
         datMan[datMan$signiCat=="<0.0001",]$logp,
         pch=21,bg=colosign[3],cex=1.5)
  axis(2,lwd=2,las=1)
  if(desiXax==1){
    axis(1,lwd=2,cex.axis=1,
         at=c(c(0,cumsum(tapply(datMan$Position,datMan$Chromosome,max))[1:14])+
                tapply(datMan$Position,datMan$Chromosome,max)/2+
                (0:14)*decalCHR),
         lab=levels(datMan$Chromosome),las=1,hadj=0.5,padj=-0.5,font=2)
  } else {
    axis(1,tcl=0.5,lwd=2,
         at=c(c(0,cumsum(tapply(datMan$Position,datMan$Chromosome,max))[1:14])+
                tapply(datMan$Position,datMan$Chromosome,max)/2+
                (0:14)*decalCHR),
         lab=NA)
  }

}


##############################################################################/
#running the functions to obtain a beautiful Manhattan plot####
##############################################################################/

AcNat<-readManDa("output/natGWAS/GAPIT.Blink.Acorn weight.GWAS.Results.csv")
AcLim<-readManDa("output/limGWAS/GAPIT.Blink.Acorn weight.GWAS.Results.csv")
HeNat<-readManDa("output/natGWAS/GAPIT.Blink.Height.GWAS.Results.csv")
HeLim<-readManDa("output/limGWAS/GAPIT.Blink.Height.GWAS.Results.csv")
PMNat<-readManDa("output/natGWAS/GAPIT.Blink.Powdery mildew.GWAS.Results.csv")
PMLim<-readManDa("output/limGWAS/GAPIT.Blink.Powdery mildew.GWAS.Results.csv")
SuNat<-readManDa("output/natGWAS/GAPIT.Blink.Survival.GWAS.Results.csv")
SuLim<-readManDa("output/limGWAS/GAPIT.Blink.Survival.GWAS.Results.csv")

colovec<-c(brewer.pal(12,"Paired")[1:2],brewer.pal(12,"Paired")[3:4])
colosign<-brewer.pal(9,"YlOrRd")[c(4,6,8)]
op<-par(mfrow=c(8,1),mar=c(1,5,3,0))
#Acorn weight
ManhaPlot(AcNat,colovec[1:2],colosign,
          decalCHR=40000000,desiXax=1,ylimi=c(0,8))
title(main="Acorn weight",font=2,cex.main=3)
legend(-50000000,9,legend=c("Exposed","Protected"),xpd=TRUE,
       pch=22,pt.cex=3,pt.lwd=4,cex=1.2,
       pt.bg=colovec[c(1,3)],col=colovec[c(2,4)],bty="n",
       x.intersp=0.4,y.intersp=1.6)
legend(150000000,9,legend=c("P<0.01","P<0.001","P<0.0001"),xpd=TRUE,
       pch=21,pt.cex=2,pt.lwd=1,cex=1.2,
       pt.bg=colosign,col="black",bty="n",
       x.intersp=0.4,y.intersp=1.2)
par(mar=c(3,5,1,0))
ManhaPlot(AcLim,colovec[3:4],colosign,
          decalCHR=40000000,desiXax=0,ylimi=c(8,0))
#Height
par(mar=c(1,5,3,0))
ManhaPlot(HeNat,colovec[1:2],colosign,
          decalCHR=40000000,desiXax=1,ylimi=c(0,8))
title(main="Height",font=2,cex.main=3)
par(mar=c(3,5,1,0))
ManhaPlot(HeLim,colovec[3:4],colosign,
          decalCHR=40000000,desiXax=0,ylimi=c(8,0))
#Powdery mildew
par(mar=c(1,5,3,0))
ManhaPlot(PMNat,colovec[1:2],colosign,
          decalCHR=40000000,desiXax=1,ylimi=c(0,8))
title(main="Powdery mildew",font=2,cex.main=3)
par(mar=c(3,5,1,0))
ManhaPlot(PMLim,colovec[3:4],colosign,
          decalCHR=40000000,desiXax=0,ylimi=c(8,0))
par(mar=c(1,5,3,0))
ManhaPlot(SuNat,colovec[1:2],colosign,
          decalCHR=40000000,desiXax=1,ylimi=c(0,8))
title(main="Survival",font=2,cex.main=3)
par(mar=c(3,5,1,0))
ManhaPlot(SuLim,colovec[3:4],colosign,
          decalCHR=40000000,desiXax=0,ylimi=c(8,0))
par(op)

#export 13.5 x 7 inches in .pdf




library(qqman)
op<-par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
manhattan(gwasResults,ylim=c(0,10),cex=2.2,cex.lab=2.5,
          font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=4)
par(mar=c(0,5,3,3))
manhattan(gwasResults,ylim=c(10,0),cex=2.2,cex.lab=2.5,
          font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=4,xlab="",xaxt="n")
par(op)





#this is an adaptation of the manhattan plot fonction coded by stephen Turner, 
#only very few things were changed to fit my need

# Stephen Turner
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html

# Last updated: Tuesday, April19, 2011
# R code for making manhattan plots and QQ plots from plink output files. 
# manhattan() with GWAS data this can take a lot of memory, recommended for use on 64bit machines only, for now. 
# Altnernatively, use bmanhattan() , i.e., base manhattan. uses base graphics. way faster.


## This is for testing purposes.
# set.seed(42)
# nchr=23
# nsnps=1000
# d=data.frame(
#     SNP=sapply(1:(nchr*nsnps), function(x) paste("rs",x,sep='')),
#     CHR=rep(1:nchr,each=nsnps), 
#     BP=rep(1:nsnps,nchr), 
#   P=runif(nchr*nsnps)
# )
# annotatesnps <- d$SNP[7550:7750]

# manhattan plot using base graphics
manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, ...) {
  
  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
  d$logp = -log10(d$P)
  d$pos=NA
  ticks=NULL
  lastbase=0
  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  if (ymax=="max") ymax<-ceiling(max(d$logp))
  #if (ymax<8) ymax<-8 #I comment this line because I don't really see the point 
  #to have a parameter for ymax if it's set to 8 afterwards
  
  numchroms=length(unique(d$CHR))
  if (numchroms==1) {
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
  } else {
    for (i in unique(d$CHR)) {
      if (i==1) {
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)+50000 #increasing the distance between chromosomes by 50000 bp
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
      }
      ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    }
  }
  
  if (numchroms==1) {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  }	else {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
    axis(1, at=ticks, lab=unique(d$CHR), ...)
    icol=1
    for (i in unique(d$CHR)) {
      with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
      icol=icol+1
    }
  }
  
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, logp, col="red", ...)) #color modification
  }
  
  if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (genomewideline) abline(h=genomewideline, col="red")
}


## Make a pretty QQ plot of p-values
qq = function(pvector, ...) {
  if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")
  pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( ppoints(length(pvector) ))
  plot(e,o,pch=19,cex=1, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
  abline(0,1,col="red")
}




#and now we are using it to plot

library(RColorBrewer)
display.brewer.all()
mypalette<-brewer.pal(8,"Dark2")


dd<-read.table("IndPrecMo.txt",header=T)
manhattan(dd,colors=mypalette[1:2],pch=16,limitchromosomes=1:13,cex=1.2, 
          suggestiveline=F,genomewideline=F, main="Manhattan plot",
          annotate=dd$SNP[dd$QVAL=="TRUE"])

manhattan(dd,pch=16,limitchromosomes=1:13,cex=1.2, 
          suggestiveline=F,genomewideline=F, main="Manhattan plot IndPreMo",
          annotate=dd$SNP[dd$QVAL=="TRUE"])

pdf(file="ddddd.pdf",width=25,height=5)
manhattan(dd,pch=16,limitchromosomes=1:13,cex=0.5, ymax=5,
          suggestiveline=F,genomewideline=F, main="Manhattan plot IndPreMo",
          annotate=dd$SNP[dd$QVAL=="TRUE"])
dev.off()

