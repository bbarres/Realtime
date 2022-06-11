##############################################################################/
##############################################################################/
#GWAS analyses using GAPIT3
##############################################################################/
##############################################################################/

# #at the time of writing this code, GAPIT3 was still under development on 
# #Github, so you have to run this to install the updated package
# devtools::install_github("jiabowang/GAPIT3",force=TRUE)

source("RealTime_load.R")


##############################################################################/
#function GENHET (Aurélie Coulon)####
##############################################################################/

#this code is available here: 
#http://aureliecoulon.net/0bf520f8_e6b3_4e12_8fc9_4867ce9ebad2.html

"GENHET"<-
  function(dat,estimfreq,locname,alfuser){
    
    nbloc=(ncol(dat)-1)/2
    nbind=nrow(dat)
    
    #estimation of allele frequencies (only if estimfreq=T)
    
    if(estimfreq=="T")
      
    {
      
      #creation of the list of alleles
      datv=vector(length=nbind*nbloc*2)
      for (i in 2:ncol(dat)) {
        datv[(nrow(dat)*(i-2)+1):(nrow(dat)*(i-1))]=dat[,i]
      }
      al=sort(na.omit(unique(datv)))
      
      #count of the number of times each allele appears + nb of missing data
      alcount=matrix(nrow=(length(al)+1),ncol=(nbloc+1))
      alcount[,1]=c(al,NA)
      for(j in 1:(nrow(alcount)-1))
        for(k in 1:(ncol(alcount)-1))
          alcount[j,(k+1)]=sum(dat[,(k*2):(k*2+1)]==alcount[j,1],na.rm=T)
      for(l in 2:ncol(alcount))
        alcount[nrow(alcount),l]=(2*nbind-sum(alcount[1:(nrow(alcount)-1),l]))
      
      
      #creation of the table of allele frequencies
      alfreq=matrix(nrow=length(al),ncol=(nbloc+1))
      colnames(alfreq)=c("Allele",locname)
      alfreq[,1]=al
      for(m in (1:nrow(alfreq)))
        for (n in 2:ncol(alfreq)){
          alfreq[m,n]=alcount[m,n]/(nbind*2-alcount[nrow(alcount),n])
        }
    }
    
    else alfreq=alfuser
    
    
    dat=as.data.frame(dat)
    library(gtools)
    res=matrix(nrow=nrow(dat),ncol=6)
    colnames(res)=c("sampleid","PHt","Hs_obs","Hs_exp","IR","HL")
    res[,1]=as.character(dat[,1])
    
    
    
    #estimation of E per locus (for HL and Hs_exp)
    E=vector(length=nbloc)
    alfreq2=alfreq[,2:ncol(alfreq)]*alfreq[,2:ncol(alfreq)]
    for(k in 1:ncol(alfreq2)) E[k]=1-sum(alfreq2[,k],na.rm=T)
    
    
    #estimation of the mean heterozygosity per locus
    mHtl=vector(length=nbloc)
    ctNAl=0
    ctHtl=0
    for(l in 1:ncol(dat))
    {if (even(l)==T)
    {
      for (m in 1:nrow(dat))
      { if (is.na(dat[m,l])==T) ctNAl=(ctNAl+1)
      else if (is.na(dat[m,(l+1)])==T) ctNAl=(ctNAl+1)
      else if (dat[m,l]!=dat[m,(l+1)]) ctHtl=(ctHtl+1)
      }
      mHtl[l/2]=ctHtl/(nrow(dat)-ctNAl)
      ctNAl=0
      ctHtl=0
    }
    }
    
    #the program in itself
    
    ctHt=0
    ctNA=0
    
    ctHm=0
    smHtl=0
    mmHtl=0
    sE=0
    mE=0
    
    sfl=0
    
    sEh=0
    sEj=0
    
    for(i in 1:nrow(dat))
    { for (j in 2:(nbloc*2))
    { if (even(j)==T)
    {
      if (is.na(dat[i,j])==T) ctNA=(ctNA+1)
      else if (is.na(dat[i,(j+1)])==T) ctNA=(ctNA+1)
      else {
        if (dat[i,j]!=dat[i,(j+1)])
        {
          ctHt=(ctHt+1)
          sEj=sEj+E[j/2]
        }
        else sEh=sEh+E[j/2]
        smHtl=smHtl+mHtl[j/2]
        sE=sE+E[j/2]
        sfl=sfl+alfreq[alfreq[,1]==as.numeric(dat[i,j]),
                       (j/2+1)]+alfreq[alfreq[,1]==as.numeric(dat[i,j+1]),
                                       (j/2+1)]
      }
    }
    }
      res[i,2]=ctHt/(nbloc-ctNA)
      mmHtl=smHtl/(nbloc-ctNA)
      res[i,3]=(ctHt/(nbloc-ctNA))/mmHtl
      mE=sE/(nbloc-ctNA)
      res[i,4]=(ctHt/(nbloc-ctNA))/mE
      ctHm=nbloc-ctHt-ctNA
      res[i,5]=(2*ctHm-sfl)/(2*(nbloc-ctNA)-sfl)
      res[i,6]=sEh/(sEh+sEj)
      ctHt=0
      ctNA=0
      ctHm=0
      smHtl=0
      mmHtl=0
      sE=0
      mE=0
      sfl=0
      sEh=0
      sEj=0
    }
    return(res)
  }


##############################################################################/
#preparing the data set####
##############################################################################/

#preparing the data set
#list of the snp marker names
temp<-colnames(snp.dat)[48:(48+818)]
#creating a dataframe for the correspondence between true 
#names and simple names
nomSNP<-as.data.frame(temp)
colnames(nomSNP)<-"realNames"
nomSNP$simpleNames<-paste("snp",1:819,sep="")

snpGen<-df2genind(snp.dat[,temp],ploidy=2,sep="/",
                  ind.names=snp.dat$Sample_ID,
                  pop=snp.dat$family_simp)
#adding information in the 'other' slot
snpGen@other$fam<-snp.dat$family_simp
snpGen@other$treat<-snp.dat$exp
snpGen@other$height<-snp.dat$Hdeb17
snpGen@other$DoA<-snp.dat$live_bin
#combine experimental and dead or alive information
snpGen@other$newPop<-paste(snp.dat$exp,snp.dat$live_bin,sep="")
locNames(snpGen)<-paste("snp",1:819,sep="")

#we limit the data set to the individuals belonging to the experimental
#set up
snpGen<-snpGen[(snpGen@other$fam!="CC" & snpGen@other$fam!="PAR")]
#we also remove individuals without dead or alive information
snpGen<-snpGen[!is.na(snpGen@other$DoA)]


##############################################################################/
#Heterozygosity dead vs alive####
##############################################################################/

#set newPop as population
snpGen2<-snpGen
pop(snpGen2)<-snpGen2@other$newPop
nomFam<-popNames(snpGen2)
temp<-genind2df(snpGen2,oneColPerAll=TRUE)
temp[temp=="A"]<-10
temp[temp=="T"]<-20
temp[temp=="C"]<-30
temp[temp=="G"]<-40
tempop<-temp[,1]
sampleid<-row.names(temp)
temp<-temp[,-1]
temp<-data.frame(lapply(temp,as.numeric))
temp$sampleid<-sampleid
temp<-temp[,c(1639,1:1638)]
temp2<-GENHET(dat=temp,estimfreq="T",locname=nomSNP$simpleNames)
temp2<-as.data.frame(temp2)
temp2$pop<-tempop
temp2$moda<-stringr::str_match(tempop, "(...)(.*)")[,2]
temp2$vivmor<-stringr::str_match(tempop, "(...)(.*)")[,3]
temp2$fam<-"Global"
HetVivMort<-temp2

#same thing by family
snpGenfam<-seppop(snpGen)
nomFam<-popNames(snpGen)
for (i in 1:length(nomFam)) {
  pop(snpGenfam[[i]])<-snpGenfam[[i]]@other$newPop
  temp<-genind2df(snpGenfam[[i]],oneColPerAll=TRUE)
  temp[temp=="A"]<-10
  temp[temp=="T"]<-20
  temp[temp=="C"]<-30
  temp[temp=="G"]<-40
  tempop<-temp[,1]
  sampleid<-row.names(temp)
  temp<-temp[,-1]
  temp<-data.frame(lapply(temp,as.numeric))
  temp$sampleid<-sampleid
  temp<-temp[,c(1639,1:1638)]
  temp2<-GENHET(dat=temp,estimfreq="T",locname=nomSNP$simpleNames)
  temp2<-as.data.frame(temp2)
  temp2$pop<-tempop
  temp2$moda<-stringr::str_match(tempop, "(...)(.*)")[,2]
  temp2$vivmor<-stringr::str_match(tempop, "(...)(.*)")[,3]
  temp2$fam<-nomFam[i]
  HetVivMort<-rbind(HetVivMort,temp2)
}

HetVivMortExp<-HetVivMort[HetVivMort$moda=="exp",]
#removing global
HetVivMortExp<-HetVivMortExp[HetVivMortExp$fam!="Global",]
HetVivMortExp$fam[HetVivMortExp$fam=="1"]<-"01"
HetVivMortExp$fam[HetVivMortExp$fam=="9"]<-"09"

HetVivMortLim<-HetVivMort[HetVivMort$moda=="low",]
#removing global
HetVivMortLim<-HetVivMortLim[HetVivMortLim$fam!="Global",]
HetVivMortLim$fam[HetVivMortLim$fam=="1"]<-"01"
HetVivMortLim$fam[HetVivMortLim$fam=="9"]<-"09"

op<-par(mfrow=c(5,1))
vioplot(as.numeric(HetVivMortExp$PHt)~HetVivMortExp$vivmor:HetVivMortExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Dead or Alive: Exposed / PHt")
vioplot(as.numeric(HetVivMortExp$Hs_obs)~
          HetVivMortExp$vivmor:HetVivMortExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Dead or Alive: Exposed / Hs_obs")
vioplot(as.numeric(HetVivMortExp$Hs_exp)~
          HetVivMortExp$vivmor:HetVivMortExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Dead or Alive: Exposed / Hs_exp")
vioplot(as.numeric(HetVivMortExp$IR)~HetVivMortExp$vivmor:HetVivMortExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Dead or Alive: Exposed / IR")
vioplot(as.numeric(HetVivMortExp$HL)~HetVivMortExp$vivmor:HetVivMortExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Dead or Alive: Exposed / HL")
par(op)
#export to .pdf 20 x 20 inches

op<-par(mfrow=c(5,1))
vioplot(as.numeric(HetVivMortLim$PHt)~HetVivMortLim$vivmor:HetVivMortLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Dead or Alive: Limited / PHt")
vioplot(as.numeric(HetVivMortLim$Hs_obs)~
          HetVivMortLim$vivmor:HetVivMortLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Dead or Alive: Limited / Hs_obs")
vioplot(as.numeric(HetVivMortLim$Hs_exp)~
          HetVivMortLim$vivmor:HetVivMortLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Dead or Alive: Limited / Hs_exp")
vioplot(as.numeric(HetVivMortLim$IR)~HetVivMortLim$vivmor:HetVivMortLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Dead or Alive: Limited / IR")
vioplot(as.numeric(HetVivMortLim$HL)~HetVivMortLim$vivmor:HetVivMortLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Dead or Alive: Limited / HL")
par(op)
#export to .pdf 20 x 20 inches

op<-par(mfrow=c(2,1))
#semi violin plot for the exposed treatment
HetAli<-HetVivMortExp[HetVivMortExp$vivmor==1,]
HetDea<-HetVivMortExp[HetVivMortExp$vivmor==0,]
colovec<-c(brewer.pal(12,"Set3")[6:7],
           brewer.pal(9,"Set1")[1:2])
vioplot(as.numeric(HetVivMortExp$PHt)~HetVivMortExp$fam,
        col = c("transparent"),sep=":",las=1,border="transparent",
        ylab="PHt",xlab="",ylim=c(0.07,0.35),frame.plot=FALSE,
        lineCol="transparent",rectCol="transparent",main="Exposed treatment")
vioplot(as.numeric(HetDea$PHt)~HetDea$fam,plotCentre="line",
        col=colovec[1],sep=":",las=1,side="left",frame.plot=FALSE,
        add=TRUE)
stripchart(as.numeric(HetDea$PHt)~HetDea$fam,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col=colovec[3],
           at=c(1:15)-0.2)
vioplot(as.numeric(HetAli$PHt)~HetAli$fam,plotCentre="line",
        col=colovec[2],sep=":",las=1,side="right",frame.plot=FALSE,
        add=TRUE)
stripchart(as.numeric(HetAli$PHt)~HetAli$fam,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col=colovec[4],
           at=c(1:15)+0.2)
segments(c(1:15)-0.2,aggregate(as.numeric(HetDea$PHt),
                               list(HetDea$fam),FUN=mean)[,2],
         c(1:15)+0.2,aggregate(as.numeric(HetAli$PHt),
                               list(HetAli$fam),FUN=mean)[,2],
         col=grey(0.95,0.9),lwd=6)
points(x=c(1:15)-0.2,y=aggregate(as.numeric(HetDea$PHt),
                                  list(HetDea$fam),FUN=mean)[,2],
       pch=19,col=colovec[3])
points(x=c(1:15)+0.2,y=aggregate(as.numeric(HetAli$PHt),
                                  list(HetAli$fam),FUN=mean)[,2],
       pch=19,col=colovec[4])
box(bty="l")
legend(-0.5,0.43,c("dead","alive"),fill=colovec[1:2],cex=1.3,
       bty="n",x.intersp=0.2,y.intersp=0.7,xpd=TRUE)

#semi violin plot for the exposed treatment
HetAli<-HetVivMortLim[HetVivMortLim$vivmor==1,]
HetDea<-HetVivMortLim[HetVivMortLim$vivmor==0,]
colovec<-c(brewer.pal(12,"Set3")[6:7],
           brewer.pal(9,"Set1")[1:2])
vioplot(as.numeric(HetVivMortLim$PHt)~HetVivMortLim$fam,
        col = c("transparent"),sep=":",las=1,border="transparent",
        ylab="PHt",xlab="",ylim=c(0.07,0.35),frame.plot=FALSE,
        lineCol="transparent",rectCol="transparent",main="Limited treatment")
vioplot(as.numeric(HetDea$PHt)~HetDea$fam,plotCentre="line",
        col=colovec[1],sep=":",las=1,side="left",frame.plot=FALSE,
        add=TRUE)
stripchart(as.numeric(HetDea$PHt)~HetDea$fam,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col=colovec[3],
           at=c(1:15)-0.2)
vioplot(as.numeric(HetAli$PHt)~HetAli$fam,plotCentre="line",
        col=colovec[2],sep=":",las=1,side="right",frame.plot=FALSE,
        add=TRUE)
stripchart(as.numeric(HetAli$PHt)~HetAli$fam,vertical=TRUE,
           method="jitter",pch=21,add=TRUE,col=colovec[4],
           at=c(1:15)+0.2)
segments(c(1:15)-0.2,aggregate(as.numeric(HetDea$PHt),
                               list(HetDea$fam),FUN=mean)[,2],
         c(1:15)+0.2,aggregate(as.numeric(HetAli$PHt),
                               list(HetAli$fam),FUN=mean)[,2],
         col=grey(0.95,0.9),lwd=6)
points(x=c(1:15)-0.2,y=aggregate(as.numeric(HetDea$PHt),
                                 list(HetDea$fam),FUN=mean)[,2],
       pch=19,col=colovec[3])
points(x=c(1:15)+0.2,y=aggregate(as.numeric(HetAli$PHt),
                                 list(HetAli$fam),FUN=mean)[,2],
       pch=19,col=colovec[4])
box(bty="l")
par(op)

#export to .pdf 10 x 12 inches

var.test(as.numeric(HetVivMortExp$PHt)~HetVivMortExp$vivmor:HetVivMortExp$fam)



##############################################################################/
#Heterozygosity total vs surviving####
##############################################################################/

#reshaping the data set for computation
n.temp<-seppop(snpGen2,treatOther=TRUE)
n.temp.other<-lapply(n.temp,as.data.frame(other))
temp<-repool(n.temp$exp1,n.temp$exp0)
pop(temp)<-rep("expinit",times=nInd(temp))
temp<-repool(temp,n.temp$exp1)
temp2<-repool(n.temp$low1,n.temp$low0)
pop(temp2)<-rep("lowinit",times=nInd(temp2))
snpGen3<-repool(temp,temp2,n.temp$low1)
snpGen3@other<-rbind(n.temp.other$exp1,n.temp.other$exp0,
                  n.temp.other$exp1,n.temp.other$low1,
                  n.temp.other$low0,n.temp.other$low1)
colnames(snpGen3@other)<-c("fam","treat","height","DoA","newPop")
snpGen3@other$newPop<-pop(snpGen3)

nomFam<-popNames(snpGen3)
temp<-genind2df(snpGen3,oneColPerAll=TRUE)
temp[temp=="A"]<-10
temp[temp=="T"]<-20
temp[temp=="C"]<-30
temp[temp=="G"]<-40
tempop<-temp[,1]
sampleid<-row.names(temp)
temp<-temp[,-1]
temp<-data.frame(lapply(temp,as.numeric))
temp$sampleid<-sampleid
temp<-temp[,c(1639,1:1638)]
temp2<-GENHET(dat=temp,estimfreq="T",locname=nomSNP$simpleNames)
temp2<-as.data.frame(temp2)
temp2$pop<-tempop
temp2$moda<-stringr::str_match(tempop, "(...)(.*)")[,2]
temp2$vivmor<-stringr::str_match(tempop, "(...)(.*)")[,3]
temp2$fam<-"Global"
HetStarSto<-temp2

#same thing by family
pop(snpGen3)<-snpGen3@other$fam
snpGen3fam<-seppop(snpGen3)
nomFam<-popNames(snpGen3)
for (i in 1:length(nomFam)) {
  pop(snpGen3fam[[i]])<-snpGen3fam[[i]]@other$newPop
  temp<-genind2df(snpGen3fam[[i]],oneColPerAll=TRUE)
  temp[temp=="A"]<-10
  temp[temp=="T"]<-20
  temp[temp=="C"]<-30
  temp[temp=="G"]<-40
  tempop<-temp[,1]
  sampleid<-row.names(temp)
  temp<-temp[,-1]
  temp<-data.frame(lapply(temp,as.numeric))
  temp$sampleid<-sampleid
  temp<-temp[,c(1639,1:1638)]
  temp2<-GENHET(dat=temp,estimfreq="T",locname=nomSNP$simpleNames)
  temp2<-as.data.frame(temp2)
  temp2$pop<-tempop
  temp2$moda<-stringr::str_match(tempop, "(...)(.*)")[,2]
  temp2$vivmor<-stringr::str_match(tempop, "(...)(.*)")[,3]
  temp2$fam<-nomFam[i]
  HetStarSto<-rbind(HetStarSto,temp2)
}

#in order to have the good order of the categories, we change 
#the vivmor catefories
HetStarSto[HetStarSto$vivmor=="init","vivmor"]<-"beg"
HetStarSto[HetStarSto$vivmor=="1","vivmor"]<-"end"

HetStarStoExp<-HetStarSto[HetStarSto$moda=="exp",]
HetStarStoLim<-HetStarSto[HetStarSto$moda=="low",]


op<-par(mfrow=c(5,1))
vioplot(as.numeric(HetStarStoExp$PHt)~HetStarStoExp$vivmor:HetStarStoExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Exposed / PHt")
vioplot(as.numeric(HetStarStoExp$Hs_obs)~
          HetStarStoExp$vivmor:HetStarStoExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Exposed / Hs_obs")
vioplot(as.numeric(HetStarStoExp$Hs_exp)~
          HetStarStoExp$vivmor:HetStarStoExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Exposed / Hs_exp")
vioplot(as.numeric(HetStarStoExp$IR)~HetStarStoExp$vivmor:HetStarStoExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Exposed / IR")
vioplot(as.numeric(HetStarStoExp$HL)~HetStarStoExp$vivmor:HetStarStoExp$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Exposed / HL")
par(op)
#export to .pdf 20 x 20 inches

op<-par(mfrow=c(5,1))
vioplot(as.numeric(HetStarStoLim$PHt)~HetStarStoLim$vivmor:HetStarStoLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Limited / PHt")
vioplot(as.numeric(HetStarStoLim$Hs_obs)~
          HetStarStoLim$vivmor:HetStarStoLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Limited / Hs_obs")
vioplot(as.numeric(HetStarStoLim$Hs_exp)~
          HetStarStoLim$vivmor:HetStarStoLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Limited / Hs_exp")
vioplot(as.numeric(HetStarStoLim$IR)~HetStarStoLim$vivmor:HetStarStoLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Limited / IR")
vioplot(as.numeric(HetStarStoLim$HL)~HetStarStoLim$vivmor:HetStarStoLim$fam,
        col = c("orange","yellow"),sep=":",las=1,
        ylab="Value",xlab="Pheno:Family",
        main="Begin vs end: Limited / HL")
par(op)
#export to .pdf 20 x 20 inches

##############################################################################/
#END
##############################################################################/