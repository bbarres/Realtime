##############################################################################/
##############################################################################/
#SNP plot functions
##############################################################################/
##############################################################################/

#loading the necessary packages and data sets
source("RealTime_load.R")


##############################################################################/
#function for plotting ellipse####
##############################################################################/

#The following code has been retrieved from a function written by 
#Peter D. M. Macdonald of McMaster University
ellipse<-function(hlaxa=1,hlaxb=1,theta=0,xc=0,yc=0,newplot=F,npoints=100,...)
{
  a<-seq(0,2*pi,length=npoints+1)
  x<-hlaxa*cos(a)
  y<-hlaxb*sin(a)
  alpha<-angle(x,y)
  rad<-sqrt(x^2+y^2)
  xp<-rad*cos(alpha+theta)+xc
  yp<-rad*sin(alpha+theta)+yc
  if (newplot)
    plot(xp,yp,type="l",...)
  else lines(xp,yp,...)
  invisible()
}

angle<-function(x,y)
{
  angle2<-function(xy) {
    x<-xy[1]
    y<-xy[2]
    if(x>0) {
      atan(y/x)
    }
    else {
      if (x<0&y!=0) {
        atan(y/x)+sign(y)*pi
      }
      else {
        if (x<0&y==0) {
          pi
        }
        else {
          if(y!=0) {
            (sign(y)*pi)/2
          }
          else {
            NA
          }
        }
      }
    }
  }
  apply(cbind(x,y),1,angle2)
  
}


##############################################################################/
#function for creating table list for each SNP####
##############################################################################/

#this function take a FDT output file from the GenomeStudio software. The 
#second variable 'nbSNP' is the number of SNP markers to process. It returns 
#a list of tables, one for each SNP processed
preplotSNP<-function(matnom,nbSNP)
{
  nbIND<-(dim(matnom)[2]-10)/4
  matnom1<-matnom
  nomcol<-c("GType","Score","Theta","R")
  nomcolmat<-c(colnames(matnom[1:10]),rep(nomcol,nbIND))
  names(matnom)<-nomcolmat
  listeSNP<-vector(mode="list", length=nbSNP)
  for (i in 1:nbSNP) {
    SNP<-data.frame(check.names=F)
    NomInd<-c()
    for (j in 1:nbIND) {
      SNP<-rbind(SNP,matnom[i,(4*(j-1)+11):(4*(j-1)+14)])
      NomInd<-c(NomInd,substr(colnames(matnom1[(4*(j-1)+14)]),1,7)) 
    }
    SNP<-cbind(SNP,NomInd)[,c(5,1,2,3,4)]
    listeSNP[[i]]<-SNP	
  }
  return(listeSNP)
}


##############################################################################/
#function for plotting each SNP####
##############################################################################/

#this function take as variable 'singleSNP', a list produce by 
#the 'preplot' function (see above), 'SNPstat' a file with the 
#characteristics of the clusters for each SNP imported from a SNP_T output 
#file of the GenomeStudio software with minor modifications to the column 
#names (see an example in the script below) and 'FileName' a chain of 
#character for the file name
SNPlot<-function(singleSNP,SNPstat,FileName) #one plot per page
{
  nbpag<-(length(singleSNP))
  pdf(file=paste("output\\",FileName,".pdf",sep=""))
  for (i in 1:nbpag) {
    levels(singleSNP[[i]][,2])<-list("AA"="AA","AB"="AB","BB"="BB","NC"="NC")
    attach(singleSNP[[i]])
    plot(Theta,R,col=c("red", "purple", "blue", "black")[as.numeric(GType)],
         cex=0.8,xlim=c(0,1),las=1,xlab="Normalized Theta",ylab="Normalized R",
         ylim=c(0,ifelse((summary(R)[6])>1,(summary(R)[6]),1)),
         main=paste(SNPstat$Name[i]," // score ",SNPstat$Aux[i]))
    detach(singleSNP[[i]])
    attach(SNPstat)
    ellipse(AA_T_Dev[i],AA_R_Dev[i],0,AA_T_Mean[i],AA_R_Mean[i],
            col="red",lwd=3)
    ellipse(AB_T_Dev[i],AB_R_Dev[i],0,AB_T_Mean[i],AB_R_Mean[i],
            col="purple",lwd=3)
    ellipse(BB_T_Dev[i],BB_R_Dev[i],0,BB_T_Mean[i],BB_R_Mean[i],
            col="blue",lwd=3)
    detach(SNPstat)
  }
  dev.off()
}


##############################################################################/
#importing the supplementary datasets####
##############################################################################/

#supplementary data importation for SNP results plotting
coordIndSNP<-read.table("data/dataSup/SNP_Ind_coord_FDT.txt",
                        header=T,sep="\t",dec=".",stringsAsFactors=TRUE)
#a glimpse to the data structure
coordIndSNP[1:10,1:14]

#supplementary data importation of the cluster statistic for each SNP
statTable<-read.table("data/dataSup/SNP_Stat_Tab_T.txt",
                      header=T,sep="\t",dec=".",stringsAsFactors=TRUE)
names(statTable)
statTable<-statTable[,-c(13,23)]
colnames(statTable)<-c("Index","Name","Chr","Position","ChiTest100",
                       "Het_Excess","AA_Freq","AB_Freq","BB_Freq",
                       "Call_Freq","Minor_Freq","Aux","P-C_Errors",
                       "P-P-C_Errors","Rep_Errors","10_GC","50_GC",
                       "SNP","Calls","no_calls","Plus_Minus_Strand",
                       "Address","GenTrain_Score","Orig_Score",
                       "Edited","Cluster_Sep","AA_T_Mean","AA_T_Dev",
                       "AB_T_Mean","AB_T_Dev","BB_T_Mean","BB_T_Dev",
                       "AA_R_Mean","AA_R_Dev","AB_R_Mean","AB_R_Dev",
                       "BB_R_Mean","BB_R_Dev","Address2","Norm_ID")

#list of individuals to plot (excluding bad quality individuals, individuals 
#that belong to CC and PAR categories and individuals that have not been
#genotyped with the SNP chip)
listIndBarCode<-RTdata[RTdata$Quality_SNPage==1 & 
                         !is.na(RTdata$Quality_SNPage) &
                         RTdata$pb_robot_SNPage==0 & 
                         RTdata$family_simp!="CC" &
                         RTdata$family_simp!="PAR","barcode"]

#addition to make to the name of individuals
subcolname<-c(".GType",".Score",".Theta",".R")

#list of columns names to keep
listcolretain<-c(colnames(coordIndSNP)[1:10],
                 paste(rep(listIndBarCode,each=length(subcolname)),
                       rep(subcolname,length(listIndBarCode)),
                       sep=""))

#limiting the dataset to the individuals included in the study
coordIndSNP<-coordIndSNP[,listcolretain]


##############################################################################/
#plotting all the 1536 SNP markers (take some time to run)####
##############################################################################/

#first we prepare the list for the 1536 SNP on the chip
mylist<-preplotSNP(coordIndSNP,1536)
#here is a glimpse to the structure of the dataset
mylist[[5]][1:20,1:5]
#now we plot the SNP and produce the Figure file
SNPlot(mylist,statTable,"Raw_1536SNP")


##############################################################################/
#plotting all the SNP markers excluded after visual inspection####
##############################################################################/

#limiting the dataset to excluded SNP after visual inspection
statTableExcl<-statTable[statTable$Aux=="31",]
coordIndSNPExcl<-coordIndSNP[statTableExcl$Index,]
#ploting the excluded SNP
mylist<-preplotSNP(coordIndSNPExcl,dim(statTableExcl)[1])
SNPlot(mylist,statTableExcl,"Raw_ExcludedSNP")


##############################################################################/
#plotting all the SNP markers saved after visual inspection####
##############################################################################/

#limiting the dataset to saved SNP after visual inspection
statTableSave<-statTable[statTable$Aux=="13",]
coordIndSNPSave<-coordIndSNP[statTableSave$Index,]
#ploting the saved SNP
mylist<-preplotSNP(coordIndSNPSave,dim(statTableSave)[1])
SNPlot(mylist,statTableSave,"Raw_SavedSNP")


##############################################################################/
#Figure S5/A15: plotting SNP markers example for score 1 (good) and 3 (bad)####
##############################################################################/

statTableGooBad<-statTable[c(77,81,103,11,2,13,15,45,93),]
coordIndSNPGooBad<-coordIndSNP[c(77,81,103,11,2,13,15,45,93),]
mylist<-preplotSNP(coordIndSNPGooBad,dim(statTableGooBad)[1])

#ploting the Figure S5/A15
pdf(file="output/Figure_S5_GoodBadSNP.pdf",width=10,height=10)
colovec<-brewer.pal(9,"Set1")[c(3,1)]
op<-par(mfrow=c(3,3),pty="s")
for (j in (1:9)) {
  levels(mylist[[j]][, 2])<-list("AA"="AA","AB"="AB","BB"="BB","NC"="NC")
  attach(mylist[[j]])
  plot(Theta,R,col=c("red", "purple", "blue", "black")[as.numeric(GType)],
       cex=0.8,xlim=c(0,1),las=1,
       ylim=c(0,ifelse((summary(R)[6])>1,(summary(R)[6]),1)),
       xlab="Normalized Theta",ylab="Normalized R",
       main=paste(statTableGooBad$Name[j]," // score ",
                  statTableGooBad$Aux[j]),
       cex.lab=2)
  box(col=colovec[as.numeric(as.factor(statTableGooBad$Aux))[j]],lwd=2)
  detach(mylist[[j]])
  attach(statTableGooBad)
  ellipse(AA_T_Dev[j],AA_R_Dev[j],0,AA_T_Mean[j],AA_R_Mean[j],
          col="red",lwd=3)
  ellipse(AB_T_Dev[j],AB_R_Dev[j],0,AB_T_Mean[j],AB_R_Mean[j],
          col="purple",lwd=3)
  ellipse(BB_T_Dev[j],BB_R_Dev[j],0,BB_T_Mean[j],BB_R_Mean[j],
          col="blue",lwd=3)
  detach(statTableGooBad)
}
par(op)
dev.off()


##############################################################################/
#Figure S6/A16: plotting SNP markers example for saved and excluded SNP####
##############################################################################/

statTableSavExcl<-rbind(statTableSave[c(21,26,2),],
                        statTableExcl[c(1,5,20),])
coordIndSNPSavExcl<-rbind(coordIndSNPSave[c(21,26,2),],
                          coordIndSNPExcl[c(1,5,20),])
mylist<-preplotSNP(coordIndSNPSavExcl,dim(statTableSavExcl)[1])

#ploting the Figure S6
pdf(file="output/Figure_S6_SavExclSNP.pdf",width=10,height=7)
colovec<-brewer.pal(9,"Set1")[c(3,1)]
op<-par(mfrow=c(2,3),pty="s")
for (j in (1:6)) {
  levels(mylist[[j]][, 2])<-list("AA"="AA","AB"="AB","BB"="BB","NC"="NC")
  attach(mylist[[j]])
  plot(Theta,R,col=c("red", "purple", "blue", "black")[as.numeric(GType)],
       cex=0.8,xlim=c(0,1),las=1,
       xlab="Normalized Theta",ylab="Normalized R",
       ylim=c(0,ifelse((summary(R)[6])>1,(summary(R)[6]),1)),
       main=paste(statTableSavExcl$Name[j]," // score ",
                  statTableSavExcl$Aux[j]),
       cex.lab=2)
  box(col=colovec[as.numeric(as.factor(statTableSavExcl$Aux))[j]],lwd=2)
  detach(mylist[[j]])
  attach(statTableSavExcl)
  ellipse(AA_T_Dev[j],AA_R_Dev[j],0,AA_T_Mean[j],AA_R_Mean[j],
          col="red",lwd=3)
  ellipse(AB_T_Dev[j],AB_R_Dev[j],0,AB_T_Mean[j],AB_R_Mean[j],
          col="purple",lwd=3)
  ellipse(BB_T_Dev[j],BB_R_Dev[j],0,BB_T_Mean[j],BB_R_Mean[j],
          col="blue",lwd=3)
  detach(statTableSavExcl)
}
par(op)
dev.off()


##############################################################################/
#END
##############################################################################/