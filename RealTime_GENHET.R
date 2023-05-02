##############################################################################/
##############################################################################/
#Individual genetic diversity analyses
##############################################################################/
##############################################################################/

#loading the necessary packages and data sets
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
snpGen@other$treat<-snp.dat$treat
snpGen@other$height<-snp.dat$Hdeb17
snpGen@other$DoA<-snp.dat$statut10
#combine experimental and dead or alive information
snpGen@other$newPop<-paste(snp.dat$treat,snp.dat$statut10,sep="")
locNames(snpGen)<-paste("snp",1:819,sep="")

#we limit the data set to the individuals belonging to the experimental
#set up
snpGen<-snpGen[(snpGen@other$fam!="CC" & snpGen@other$fam!="PAR")]
#we also remove individuals without dead or alive information
snpGen<-snpGen[!is.na(snpGen@other$DoA)]


##############################################################################/
#Heterozygosity dead vs alive computation####
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

#export of the final result file
write.table(HetVivMort,file="data/Hetind_DoA.txt",sep="\t",
            quote=FALSE,row.names=FALSE)


##############################################################################/
#Comparison of mean and variance of PHt####
##############################################################################/

HetVivMort<-read.table(file="data/Hetind_DoA.txt",sep="\t",
                       header=TRUE)
HetVivMort$vivmor<-as.factor(HetVivMort$vivmor)

VivMortGlob<-HetVivMort[HetVivMort$fam=="Global",]
DoANat<-VivMortGlob[VivMortGlob$moda=="exp",]
DoAPrt<-VivMortGlob[VivMortGlob$moda=="low",]

#looking at the distribution
plot(density(DoANat$PHt),main="Natural treatment")
lines(density(DoANat[DoANat$vivmor=="0",]$PHt),col="red")
lines(density(DoANat[DoANat$vivmor=="1",]$PHt),col="blue")
qqnorm(DoANat$PHt)
qqline(DoANat$PHt)
#tests for normality
shapiro.test(DoANat$PHt)
ks.test(DoANat$PHt,'pnorm')
#comparison of mean PHt between Dead and Alive
wilcox.test(PHt~vivmor,data=DoANat,exact=FALSE)
#comparison of the variance between Dead and Alive
fligner.test(PHt~vivmor,data=DoANat)
mean(DoANat[DoANat$vivmor=="0",]$PHt)
sqrt(var(DoANat[DoANat$vivmor=="0",]$PHt))
var(DoANat[DoANat$vivmor=="0",]$PHt)
mean(DoANat[DoANat$vivmor=="1",]$PHt)
sqrt(var(DoANat[DoANat$vivmor=="1",]$PHt))
var(DoANat[DoANat$vivmor=="1",]$PHt)

#looking at the distribution
plot(density(DoAPrt$PHt),main="Protected treatment")
lines(density(DoAPrt[DoAPrt$vivmor=="0",]$PHt),col="red")
lines(density(DoAPrt[DoAPrt$vivmor=="1",]$PHt),col="blue")
qqnorm(DoAPrt$PHt)
qqline(DoAPrt$PHt)
shapiro.test(DoAPrt$PHt)
ks.test(DoAPrt$PHt,'pnorm')
#comparison of mean PHt between Dead and Alive
wilcox.test(PHt~vivmor,data=DoAPrt,exact=FALSE)
#comparison of the variance between Dead and Alive
fligner.test(PHt~vivmor,data=DoAPrt)
mean(DoAPrt[DoAPrt$vivmor=="0",]$PHt)
sqrt(var(DoAPrt[DoAPrt$vivmor=="0",]$PHt))
var(DoAPrt[DoAPrt$vivmor=="0",]$PHt)
mean(DoAPrt[DoAPrt$vivmor=="1",]$PHt)
sqrt(var(DoAPrt[DoAPrt$vivmor=="1",]$PHt))
var(DoAPrt[DoAPrt$vivmor=="1",]$PHt)


##############################################################################/
#END
##############################################################################/