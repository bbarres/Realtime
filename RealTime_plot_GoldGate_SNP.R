##############################################################################/
##############################################################################/
#SNP plot functions
##############################################################################/
##############################################################################/

#loading the necessary packages and data sets
source("RealTime_load.R")


##############################################################################/
#peripheral functions####
##############################################################################/

#Functions to draw an ellipse. The following code has been retrieved from a 
#function written by Peter D. M. Macdonald of McMaster University
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
#Data formatting function####
##############################################################################/

#supplementary data importation for SNP results plotting
coordIndSNP<-read.table("data/dataSup/RT_comp_11P_FDT.txt",
                        header=T,sep="\t",dec=".")
#a glimpse to the data structure
coordIndSNP[1:10,1:14]

#supplementary data importation of the cluster statistic for each SNP
statTable<-read.table("data/dataSup/RT_comp_11P_SNP_T.txt",
                      header=T,sep="\t",dec=".")
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

#list of individuals to plot (excluding bad quality individual, individuals 
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


#a function to create a table for each SNP
preplotSNP<-function(matnom,nbSNP)
{
	matnom<-matnom[order(matnom$Index),]
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


mylist<-preplotSNP(coordIndSNP,1061)
#un petit aper?u d'un des fichiers contenus dans la list
mylist[[5]][1:20,1:5]


#une fois qu'on a cr?? la liste des SNP comme on la voulait, on fait une fonction de 
#plotage des SNPs

#cette fonction demande en entr?e un fichier de type 'singleSNP'produit par la fonction 
#'preplotSNP'et un fichier de type 'SNPstat' qui correspond ? l'importation du fichier 
#des caract?ristiques des SNP (fichier produit par GenomeStudio) dans lequel on a 
#notamment les coordonn?es des ellipses. Il faut aussi noter qu'il est important de 
#charger les fonction 'ellipse' et 'angle' en d?but de ce fichier pour pouvoir faire 
#tourner les fonctions ? l'int?rieur de cette fonction


SNPlot<-function(singleSNP,SNPstat) #avec 1 figure par page
{
	SNPstat<-SNPstat[order(SNPstat$Index),] 
	nbpag<-(length(singleSNP))
	pdf(file=paste("figSNP",".pdf"))
	for (i in 1:nbpag) {
		levels(singleSNP[[i]][,2])<-list("AA"="AA","AB"="AB","BB"="BB","NC"="NC")
		attach(singleSNP[[i]])
		plot(Theta,R,col=c("red", "purple", "blue", "black")[as.numeric(GType)],
		cex=0.8,xlim=c(0,1),ylim=c(0,ifelse((summary(R)[6])>1,(summary(R)[6]),1)),
		main=paste(SNPstat$Name[i]," // note ",SNPstat$Aux[i]))
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

SNPlot(mylist,statTable)

#faisons le pour tous les SNP
mylist<-preplotSNP(coordIndSNP,1536) #attention c'est tr?s long !
SNPlot(mylist,statTable)


#avec "n" figures par page, ici "n=6", les figures viennent d'un m?me fichier

SNPlot6<-function(singleSNP,SNPstat) 
{
	SNPstat<-SNPstat[order(SNPstat$Index),] 
	nbpag<-ceiling((length(singleSNP))/6)
	for (i in 1:nbpag) {
		pdf(file=paste((6*(i-1)+1),"_a_",6*i,".pdf"))
		op<-par(mfrow=c(2,3),pty="s")
		for (j in ((6*(i-1))+1):(i*6)) {
			levels(singleSNP[[j]][, 2])<-list("AA"="AA","AB"="AB","BB"="BB","NC"="NC")
			attach(singleSNP[[j]])
			plot(Theta,R, col = c("red", "purple", "blue", "black")[as.numeric(GType)],
			cex=0.5,xlim=c(0,1),ylim=c(0,ifelse((summary(R)[6])>1,(summary(R)[6]),1)),
			main=paste(SNPstat$Name[i]," // note ",SNPstat$Aux[i]))
			detach(singleSNP[[j]])
			attach(SNPstat)
			ellipse(AA_T_Dev[j],AA_R_Dev[j],0,AA_T_Mean[j],AA_R_Mean[j],
			        col="red",lwd=3)
			ellipse(AB_T_Dev[j],AB_R_Dev[j],0,AB_T_Mean[j],AB_R_Mean[j],
			        col="purple",lwd=3)
			ellipse(BB_T_Dev[j],BB_R_Dev[j],0,BB_T_Mean[j],BB_R_Mean[j],
			        col="blue",lwd=3)
			detach(SNPstat)
		}
		par(op)
		dev.off()
	}
}

SNPlot6(mylist,statTable)


#figure c?t? ? c?te pour pouvoir comparer deux annotations diff?rentes par exemples, il faut 
#que les fichiers de donn?es contiennent exactement les m?me SNP et les m?mes individus, sinon 
#?a donne n'importe quoi

biSNPlot<-function(singleSNP1,SNPstat1,singleSNP2,SNPstat2) #avec 1 figure par page
{
	SNPstat1<-SNPstat1[order(SNPstat1$Index),]
	SNPstat2<-SNPstat2[order(SNPstat2$Index),]
	nbpag<-(length(singleSNP1))
	pdf(file=paste("figSNPbi2",".pdf"),width=14)
	op<-par(mfrow=c(1,2),pty="s")
	for (i in 1:nbpag) {
		levels(singleSNP1[[i]][, 2])<-list("AA"="AA","AB"="AB","BB"="BB","NC"="NC")
		attach(singleSNP1[[i]])
		plot(Theta,R, col = c("red", "purple", "blue", "black")[as.numeric(GType)],
		cex=0.8,xlim=c(0,1),ylim=c(0,ifelse((summary(R)[6])>1,(summary(R)[6]),1)),
		main=paste("jeux A ",SNPstat1$Name[i]," // note ",SNPstat1$Aux[i]))
		detach(singleSNP1[[i]])
		attach(SNPstat1)
		ellipse(AA_T_Dev[i],AA_R_Dev[i],0,AA_T_Mean[i],AA_R_Mean[i],col="red",lwd=3)
		ellipse(AB_T_Dev[i],AB_R_Dev[i],0,AB_T_Mean[i],AB_R_Mean[i],col="purple",lwd=3)
		ellipse(BB_T_Dev[i],BB_R_Dev[i],0,BB_T_Mean[i],BB_R_Mean[i],col="blue",lwd=3)
		detach(SNPstat1)
		levels(singleSNP2[[i]][, 2])<-list("AA"="AA","AB"="AB","BB"="BB","NC"="NC")
		attach(singleSNP2[[i]])
		plot(Theta,R, col = c("red", "purple", "blue", "black")[as.numeric(GType)],
		cex=0.8,xlim=c(0,1),ylim=c(0,ifelse((summary(R)[6])>1,(summary(R)[6]),1)),
		main=paste("jeux B ",SNPstat2$Name[i]," // note ",SNPstat2$Aux[i]))
		detach(singleSNP2[[i]])
		attach(SNPstat2)
		ellipse(AA_T_Dev[i],AA_R_Dev[i],0,AA_T_Mean[i],AA_R_Mean[i],col="red",lwd=3)
		ellipse(AB_T_Dev[i],AB_R_Dev[i],0,AB_T_Mean[i],AB_R_Mean[i],col="purple",lwd=3)
		ellipse(BB_T_Dev[i],BB_R_Dev[i],0,BB_T_Mean[i],BB_R_Mean[i],col="blue",lwd=3)
		detach(SNPstat2)
	}
	par(op)
	dev.off()
}

#run d'exemple de biplot, ne pas oublier d'enlever les di?ses dans les fichier de type 
#'SNP_Table'
coordIndSNP1<-read.table("rs_brute_FDT.txt", header=T, sep="\t", dec=".")
statTable1<-read.table("rs_brute_SNP_T.txt", header=T, sep="\t", dec=".")
colnames(statTable1)<-c("Index","Name","Chr","Position","ChiTest100","Het_Excess","AA_Freq",
"AB_Freq","BB_Freq","Call_Freq","Minor_Freq","Aux","P-C_Errors","P-P-C_Errors","Rep_Errors",
"p10_GC","p50_GC","SNP","Calls","no_calls","Plus/Minus_Strand","Address","GenTrain_Score",
"Orig_Score","Edited","Cluster_Sep","AA_T_Mean","AA_T_Dev","AB_T_Mean","AB_T_Dev","BB_T_Mean",
"BB_T_Dev","AA_R_Mean","AA_R_Dev","AB_R_Mean","AB_R_Dev","BB_R_Mean","BB_R_Dev","Address2",
"Norm_ID")
coordIndSNP2<-read.table("rs_FDT.txt", header=T, sep="\t", dec=".")
statTable2<-read.table("rs_SNP_T.txt", header=T, sep="\t", dec=".")
colnames(statTable2)<-c("Index","Name","Chr","Position","ChiTest100","Het_Excess","AA_Freq",
"AB_Freq","BB_Freq","Call_Freq","Minor_Freq","Aux","P-C_Errors","P-P-C_Errors","Rep_Errors",
"10_GC","50_GC","SNP","Calls","no_calls","Plus/Minus_Strand","Address","GenTrain_Score",
"Orig_Score","Edited","Cluster_Sep","AA_T_Mean","AA_T_Dev","AB_T_Mean","AB_T_Dev","BB_T_Mean",
"BB_T_Dev","AA_R_Mean","AA_R_Dev","AB_R_Mean","AB_R_Dev","BB_R_Mean","BB_R_Dev","Address2",
"Norm_ID")
mylist1<-preplotSNP(coordIndSNP1,1536)
mylist2<-preplotSNP(coordIndSNP2,1536)
biSNPlot(mylist1,statTable1,mylist2,statTable2)



statSNP<-statSNP[order(statSNP$Index),] 
attach(statSNP)
ellipse(AA_T_Dev[6],AA_R_Dev[6],0,AA_T_Mean[6],AA_R_Mean[6],col="red",lwd=3)
ellipse(AB_T_Dev[6],AB_R_Dev[6],0,AB_T_Mean[6],AB_R_Mean[6],col="purple",lwd=3)
ellipse(BB_T_Dev[6],BB_R_Dev[6],0,BB_T_Mean[6],BB_R_Mean[6],col="blue",lwd=3)
detach(statSNP)
par(op)

******************************************************************************************
#il faut utiliser les g?notypes pour les couleurs, on utilise cette astuce pour que 
#tous les SNP aient les 4 types de levels, pour ?viter d'avoir des probl?mes de couleur 
#par la suite
levels(mylist[[6]][, 1])<-list("AA"="AA","AB"="AB","BB"="BB","NC"="NC")
levels(mylist[[6]][, 1])
#?a marche pour changer les noms aussi
##levels(mylist[[6]][, 1])<-list("red"="AA","purple"="AB","blue"="BB","black"="NC")

attach(mylist[[6]])
plot(Theta,R, col = c("red", "purple", "blue", "black")[as.numeric(GType)],cex=0.5,
xlim=c(0,1),ylim=c(0,1.2),main=statSNP$Name[6])
detach(mylist[[6]])
statSNP<-statSNP[order(statSNP$Index),] 
attach(statSNP)
ellipse(AA_T_Dev[6],AA_R_Dev[6],0,AA_T_Mean[6],AA_R_Mean[6],col="red",lwd=3)
ellipse(AB_T_Dev[6],AB_R_Dev[6],0,AB_T_Mean[6],AB_R_Mean[6],col="purple",lwd=3)
ellipse(BB_T_Dev[6],BB_R_Dev[6],0,BB_T_Mean[6],BB_R_Mean[6],col="blue",lwd=3)
detach(statSNP)

*****************************************************************************************


#ensuite pour faire tourner une fonction sur la liste des SNP
lapply()


ls()

rm(list=ls())

j<-1
matnom<-coordIndSNP
colnames(matnom[(4*(j-1)+14)])


m <- cbind(1, 1:7) # the '1' (= shorter vector) is recycled
m
m <- cbind(m, 8:14)[, c(1, 3, 2)] # insert a column
m
cbind(1:7, diag(3))# vector is subset -> warning

cbind(0, rbind(1, 1:3))
cbind(I=0, X=rbind(a=1, b=1:3))  # use some names
xx <- data.frame(I=rep(0,2))
cbind(xx, X=rbind(a=1, b=1:3))   # named differently

cbind(0, matrix(1, nrow=0, ncol=4))#> Warning (making sense)
dim(cbind(0, matrix(1, nrow=2, ncol=0)))#-> 2 x 1

## deparse.level
dd <- 10
rbind(1:4, c=2, "a++" = 10, dd, deparse.level=0)# middle 2 rownames
rbind(1:4, c=2, "a++" = 10, dd, deparse.level=1)# 3 rownames (default)
rbind(1:4, c=2, "a++" = 10, dd, deparse.level=2)# 4 rownames


##############################################################################/
#END
##############################################################################/