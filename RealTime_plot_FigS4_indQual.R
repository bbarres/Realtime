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




#############################################################################
#S?lection des SNP
#############################################################################


#on va essayer de se faire une id?e de quels crit?res retenir pour faire de 
#la s?lection automatique des bons et des mauvais SNP

SNPTable<-read.table("cc_SNP_T.txt", header=T, sep="\t", dec=".")
colnames(SNPTable)<-c("Index","Name","Chr","Position","ChiTest100",
                      "Het_Excess","AA_Freq","AB_Freq","BB_Freq","Call_Freq",
                      "Minor_Freq","Aux","P-C_Errors","P-P-C_Errors",
                      "Rep_Errors","p10_GC","p50_GC","SNP","Calls",
                      "no_calls","Plus/Minus_Strand","Address",
                      "GenTrain_Score","Orig_Score","Edited","Cluster_Sep",
                      "AA_T_Mean","AA_T_Dev","AB_T_Mean","AB_T_Dev",
                      "BB_T_Mean","BB_T_Dev","AA_R_Mean","AA_R_Dev",
                      "AB_R_Mean","AB_R_Dev","BB_R_Mean","BB_R_Dev",
                      "Address2","Norm_ID")


#premi?re chose, si on regarde les densit? de distribution du centre moyen 
#des 3 clusters, on observe cette distribution
plot(density(SNPTable$BB_T_Mean), xlim=c(0,1),col="blue",lwd=2)
lines(density(SNPTable$AA_T_Mean),col="red",lwd=2)
lines(density(SNPTable$AB_T_Mean),col="purple",lwd=2)

#on peut plotter facilement les valeurs individuelles avec la fonction 'rug'
rug((SNPTable$AB_T_Mean),col="black",lwd=0.2)

#essayons de faire la m?me chose pour les diff?rentes cat?gories de SNP 
#1=bon, 2=passable,3=pas bon, 4=bizarre et -1=monomorphe

summary(SNPTable)
op<-par(mfrow=c(2,2))
SNPTable1<-SNPTable[SNPTable$Aux==1,]
plot(density(SNPTable1$BB_T_Mean), xlim=c(0,1),col="blue",lwd=2)
lines(density(SNPTable1$AA_T_Mean),col="red",lwd=2)
lines(density(SNPTable1$AB_T_Mean),col="purple",lwd=2)
SNPTable2<-rbind(SNPTable[SNPTable$Aux==2,],SNPTable[SNPTable$Aux==23,])
plot(density(SNPTable2$AA_T_Mean), xlim=c(0,1),col="red",lwd=2)
lines(density(SNPTable2$BB_T_Mean),col="blue",lwd=2)
lines(density(SNPTable2$AB_T_Mean),col="purple",lwd=2)
SNPTable3<-SNPTable[SNPTable$Aux==3,]
plot(density(SNPTable3$BB_T_Mean), xlim=c(0,1),col="blue",lwd=2)
lines(density(SNPTable3$AA_T_Mean),col="red",lwd=2)
lines(density(SNPTable3$AB_T_Mean),col="purple",lwd=2)
SNPTable4<-SNPTable[SNPTable$Aux==4,]
plot(density(SNPTable4$BB_T_Mean), xlim=c(0,1),col="blue",lwd=2)
lines(density(SNPTable4$AA_T_Mean),col="red",lwd=2)
lines(density(SNPTable4$AB_T_Mean),col="purple",lwd=2)
par(op)

#?videmment tout ?a ne donne qu'une vague id?e des crit?res de s?lection 
#des SNP qui seraient bons. On remarque quand m?me qu'on perd un double pic 
#de densit? des 2 homozygotes quand on ne retient que les bons SNP. 

#on peut aussi plotter les individus apr?s les avoir class? selon un crit?re 
#et les color?s selon que ce sont des bons ou des mauvais SNP
SNPTable<-read.table("cc_SNP_T.txt", header=T, sep="\t", dec=".")
SNPTable[1:15,]
colnames(SNPTable)
SNPTable<-SNPTable[order(SNPTable$Call.Freq),]
attach(SNPTable)
plot(Call.Freq,col=as.factor(Aux))
#on voit tr?s bien qu'au d?but il y a quasiment que des bleus (qui sont 
#les -1), donc sous un certain seuil il faudra faire sauter les SNP.
plot(Call.Freq,col=as.factor(Aux),xlim=c(80,300))
detach(SNPTable)

#on voit tr?s bien qu'au d?but il y a quasiment que des bleus (qui sont les
#-1), donc sous un certain seuil de 'CallFreq' il faudrait faire sauter les 
#SNP.
plot(Call.Freq,col=as.factor(Aux),xlim=c(80,300))

#faisons le avec un autre jeux de donn?es pour voir
SNPTable<-read.table("rw_SNP_T.txt", header=T, sep="\t", dec=".")
SNPTable[1:15,]
colnames(SNPTable)
SNPTable<-SNPTable[order(SNPTable$Call.Freq),]
attach(SNPTable)
plot(Call.Freq,col=as.factor(Aux))
detach(SNPTable)


#On peut bien sur r?aliser ces genres de graphique pour d'autres crit?res, 
#mais essayons plut?t de faire par une approche de type classification 
#qui sera multivari? et qui prendra directement en comptes les diff?rents 
#crit?res

#D'abord simplifions la matrice en enlevant tout ce qui n'est pas de type 
#num?rique quantitatif.

ccSNPTable<-read.table("cc_SNP_T.txt", header=T, sep="\t", dec=".")
colnames(ccSNPTable)<-c("Index","Name","Chr","Position","ChiTest100",
                      "Het_Excess","AA_Freq","AB_Freq","BB_Freq","Call_Freq",
                      "Minor_Freq","Aux","P-C_Errors","P-P-C_Errors",
                      "Rep_Errors","p10_GC","p50_GC","SNP","Calls",
                      "no_calls","Plus/Minus_Strand","Address",
                      "GenTrain_Score","Orig_Score","Edited","Cluster_Sep",
                      "AA_T_Mean","AA_T_Dev","AB_T_Mean","AB_T_Dev",
                      "BB_T_Mean","BB_T_Dev","AA_R_Mean","AA_R_Dev",
                      "AB_R_Mean","AB_R_Dev","BB_R_Mean","BB_R_Dev",
                      "Address2","Norm_ID")

ccSNPTable<-ccSNPTable[,-c(3,4,5,6,13,14,15,18,21,22,25,39,40)]

#observons les diff?rentes variables quantitatives pour "humer" un peu les 
#les donn?es

pairs(ccSNPTable[,c(9:15)])
pairs(ccSNPTable[,c(16:27)])

#encore un peu de nettoyage pour que la matrice soit analysable
ccSNPTable$Aux<-(as.factor(ccSNPTable$Aux))
tccSNPTable<-ccSNPTable[order(ccSNPTable$Index),]
tccSNPTable<-tccSNPTable[,-c(1,2)]

#on charge la librairie pour faire les analyse de type CART
#ensuite essayons une premi?re m?thode de type donc Classification and 
#Regression Tree
library(rpart)
tete<-rpart(Aux~.,data=tccSNPTable,minsplit=1,xval=20)
library(tree)
tete2<-rpart(Aux~ .,data=tccSNPTable)
plot(tete)
text(tete)
plot(tete2)
text(tete2)

plot(tete,compress=T, margin=0.05,branch=0.8)
text(tete,use.n=T,cex=0.5)

#une petite fonction pour transformer un dataframe SNP_T pour pouvoir faire 
#des analyses statistiques dessus ensuites,?a consiste essentiellement ? 
#renommer les colonnes pour que les noms soit plus court et ne pose pas de 
#probl?me, et en plus on enl?ve quelques colonnes qui ne servent ? rien 
#pour les analyses. De plus on reclasse la matrice avec la colonne "Index" 
#pour pas qu'il n'y ait de probl?me d'ordre si jamais on veut comparer 
#plusieurs SNP_T avec les m?mes SNP mais sur des individus diff?rents. A 
#noter ?galement que cette fonction me sert ici ? la pr?paration d'une  
#analyse de type CART pour d?finir des r?gles de d?cision sur la qualit? 
#des SNP. Pour cela une note de qualit? a ?t? attribu? ? chaque SNP apr?s 
#v?rification "? l'oeil" et plac?e dans la colonne "Aux" du SNP_T

simpleSNPT<-function(SNPTable)
{
  colnames(SNPTable)<-c("Index","Name","Chr","Position","ChiTest100",
                      "Het_Excess","AA_Freq","AB_Freq","BB_Freq","Call_Freq",
                      "Minor_Freq","Aux","P-C_Errors","P-P-C_Errors",
                      "Rep_Errors","p10_GC","p50_GC","SNP","Calls",
                      "no_calls","Plus/Minus_Strand","Custom_cluster",
                      "Address","GenTrain_Score","Orig_Score","Edited",
                      "Cluster_Sep","AA_T_Mean","AA_T_Dev","AB_T_Mean",
                      "AB_T_Dev","BB_T_Mean","BB_T_Dev","AA_R_Mean",
                      "AA_R_Dev","AB_R_Mean","AB_R_Dev","BB_R_Mean",
                      "BB_R_Dev","Address2","Norm_ID")
  SNPTable<-SNPTable[,-c(3,4,5,6,13,14,15,18,21,22,23,25,40,41)]
  SNPTable$Aux<-(as.factor(SNPTable$Aux))
  tSNPTable<-SNPTable[order(SNPTable$Index),]
  return(tSNPTable)
}

#exemple d'utilisation de la cette nouvelle fonction avec 'rpart'
rsSNPTable<-read.table("rs_SNP_T.txt", header=T, sep="\t", dec=".")
trs_SNPT<-simpleSNPT(rsSNPTable)
pairs(trs_SNPT[,c(9:15)],col=trs_SNPT$Aux)
pairs(trs_SNPT[,c(16:27)],col=trs_SNPT$Aux)
library(rpart)
tete<-rpart(Aux~.,data=trs_SNPT,minsplit=1,xval=100)
summary(tete)
plotcp(tete)
plot(tete,compress=T, margin=0.05,branch=1)
text(tete,use.n=T,cex=0.7)

#apr?s avoir cr?? une liste des SNP comme dans 'plotterSNP', on va choisir
#des exemples dans les cas de chacunes des feuilles du CART

exempGRA<-cbind(tete$where,rsSNPTable$Aux)
colnames(exempGRA)<-c("leave","groupe")
exempGRA<-as.data.frame(exempGRA)
table(exempGRA)

#histoire de d?tecter les num?ros des SNP qui sont ? prendre pour exemple 
#on fait une petite fonction, elle prend pour chaque feuille le premier 
#SNP dans chaque cat?gorie pr?sente dans la feuille. Fonction ? perfectionner 
#parce que l? ?a n'est pas tr?s g?n?rique

creaExleaf<-function(num)
{
  exleaf<-c()
  exleaf<-rbind(exleaf,rownames(exempGRA[exempGRA$leave==num & 
    exempGRA$groupe==-1,][1,1:2]))
  exleaf<-rbind(exleaf,rownames(exempGRA[exempGRA$leave==num & 
    exempGRA$groupe==1,][1,1:2]))
  exleaf<-rbind(exleaf,rownames(exempGRA[exempGRA$leave==num & 
    exempGRA$groupe==2,][1,1:2]))
  exleaf<-rbind(exleaf,rownames(exempGRA[exempGRA$leave==num & 
    exempGRA$groupe==3,][1,1:2]))
  return(exleaf)
}

#et voil? des exemples d'utilisation. Apr?s avoir sortie la liste des feuilles 
#on fait tourner la fonction pour chaque feuille

exempGRA<-cbind(tete$where,rsSNPTable$Aux)
colnames(exempGRA)<-c("leave","groupe")
exempGRA<-as.data.frame(exempGRA)
table(exempGRA)

exleaf5<-creaExleaf(5)
exleaf6<-creaExleaf(6)
exleaf9<-creaExleaf(9)
exleaf10<-creaExleaf(10)
exleaf11<-creaExleaf(11)
exleaf13<-creaExleaf(13)
exleaf14<-creaExleaf(14)
exleaf17<-creaExleaf(17)
exleaf19<-creaExleaf(19)
exleaf20<-creaExleaf(20)
exleaf21<-creaExleaf(21)


#ensuite on fait des graphes avec la fonction cr??e dans 'plotterSNP' et 
#avec la fonction lapply et SNPlot modifi?, d'abord on modifie la fonction 
#de SNPlot

SNPlotmod<-function(singleSNP)#,SNPstat) #avec 1 figure par page
{
		levels(singleSNP[, 2])<-list("AA"="AA","AB"="AB","BB"="BB","NC"="NC")
		attach(singleSNP)
		plot(Theta,R, col = c("red", "purple", "blue", "black")
         [as.numeric(GType)],cex=0.8,xlim=c(0,1),ylim=c(0,
         ifelse((summary(R)[6])>1,(summary(R)[6]),1)),space=0,
         ann=F,xaxt="n",yaxt="n",asp=1) #,
		#main=paste(SNPstat$Name," // note ",SNPstat$Aux))
		detach(singleSNP)
		#attach(SNPstat)
		#ellipse(AA_T_Dev,AA_R_Dev,0,AA_T_Mean,AA_R_Mean,col="red",lwd=3)
		#ellipse(AB_T_Dev,AB_R_Dev,0,AB_T_Mean,AB_R_Mean,col="purple",lwd=3)
		#ellipse(BB_T_Dev,BB_R_Dev,0,BB_T_Mean,BB_R_Mean,col="blue",lwd=3)
		#detach(SNPstat)
	
}

exempList<-list(mylist1[[6]],mylist1[[10]],mylist1[[7]]) 
#,statTable1[statTable1$Index==6,])
op<-par(mfrow=c(1,4),pty="s")
lapply(exempList,SNPlotmod)
par(op)


#tout ce qui pr?c?de est assez brouillon, donc on va faire une fonction 
#pour cr?er une sous liste de la liste compl?te et qui sera ranger comme 
#il faut, puis on utilisera une des fonction plotSNP existante en la 
#modifiant un petit peu

#pour faire focntionner cette fonction, il nous faut un objet r?sultat 
#d'une analyse faite par 'rpart', et aussi la table correspondante de 
#type "_SNP_T" simplifi?e avec 'simpleSNPT', et aussi la liste des SNP 
#cr??e directement ? partir des donn?es Full Data Table sorties de 
#GenomeStudio


listeExCART<-function(classif,SNP_T,listSNP)
{
  #newwhere<-classif$where[order(as.numeric(attr(classif$where,"names")))]
  temp<-as.data.frame(cbind(classif$where,SNP_T$Aux))
  rownames(temp)<-c(1:(dim(temp)[1]))
  colnames(temp)<-c("leave","groupe")
  recap<-table(temp)
  nbleaf<-dim(recap)[1]
  nbgrp<-dim(recap)[2]
  listSNPex<-c()
  for (i in 1:nbleaf) {
    num<-as.numeric(rownames(recap)[i])
    exx<-c()
    for (j in 1:nbgrp) {
      exx<-c(exx,as.numeric(rownames(temp
      [temp$leave==num & temp$groupe==as.numeric
       (colnames(recap)[j]),][1,1:2])))
    }
    listSNPex<-c(listSNPex,exx)
  }
  selec<-vector(mode="list",length=nbleaf*nbgrp)
  for (i in 1:(nbleaf*nbgrp)) {
    Name<-c()
    ttrr<-listSNPex[[i]]
    selec[[i]]<-listSNP[[ttrr]]
    qest<-is.data.frame(listSNP[[ttrr]])
    ifelse(qest==TRUE,attr(selec[[i]],"Name")<-paste(SNP_T$Name[ttrr],
                      ' class ',SNP_T$Aux[ttrr]),Name<-'Missing Class')
  }
  #return(listSNPex)
  return(selec)
}


test<-listeExCART(tete,rsSNPTable,mylist2)

#et maintenant on modifie un peu (beaucoup) la fonction SNPlot, 
#pour que l'on puisse faire une figure avec un exemple de chaque classe 
#Aux, dans chaque feuille. Elle prend en compte une 'classif' de type CART, 
#une liste des SNP-exemples (sortie de la fonction 'listeExCART') et une 
#table des statistiques des SNP simplifi?e avec la fonction 'simpleSNPT'. 
#En sortie on a une figure avec une colonne par feuille de l'arbre CART, 
#et dans chacune de ces colonnes un graphe par classe Aux (quand il y en a). 

SNPlotCART<-function(classif,exempSNP,SNP_T) 
{
  #newwhere<-classif$where[order(as.numeric(attr(classif$where,"names")))]
  temp<-as.data.frame(cbind(classif$where,SNP_T$Aux))
  colnames(temp)<-c("leave","groupe")
  recap<-table(temp)
  nbleaf<-dim(recap)[1]
  nbgrp<-dim(recap)[2]
  donneManq<-data.frame(rbind(cbind("BBBBBBB","NC",0,0,0),
                              cbind("BBBBBBB","NC",0,1,1)))
  colnames(donneManq)<-c("NomInd","GType","Score","Theta","R")
  donneManq$Theta<-as.numeric(donneManq$Theta)
  donneManq$R<-as.numeric(donneManq$R)
  donneManq$Theta<-(donneManq$Theta-1)
  donneManq$R<-(donneManq$R-1)
  pdf(file=paste("ExempGrp",".pdf"),width=(7*nbleaf),height=(7*nbgrp))
  op<-par(mfcol=c(nbgrp,nbleaf),pty="s")
  exemp<-vector(mode="list", length=nbleaf*nbgrp)
  for (i in 1:(nbleaf*nbgrp)) {
    qest<-is.data.frame(exempSNP[[i]])
    ifelse(qest==TRUE,(exemp[[i]]<-exempSNP[[i]]),(exemp[[i]]<-donneManq))
    levels(exemp[[i]][,2])<-list("AA"="AA","AB"="AB","BB"="BB","NC"="NC")
    attach(exemp[[i]])
    yymax<-ifelse((summary(R)[6])>=1,(summary(R)[6]),1)
		plot(Theta,R, col = c("red", "purple", "blue", "black")
         [as.numeric(GType)],cex=2,xlim=c(0,1),ylim=c(0,yymax),ann=T,
         xaxt="s",yaxt="s",main=attr(exemp[[i]],"Name"),pch=19)
    detach(exemp[[i]])
	}
	par(op)
  dev.off()
}

SNPlotCART(tete,test,rsSNPTable)


#pour voir si tout fonctionne et donner une id?e de la d?marche ? suivre 
#pour faire une analyse de type CART et une figure voil? le code qui utilise 
#les fonctions d?finies ci-dessus. L'expemple concerne le jeux de donn?e 
#"cc"
cc_SNP_T<-read.table("CC_SNP_T.txt", header=T, sep="\t", dec=".",
                     comment.char="")
cc_FDT<-read.table("CC_FDT.txt", header=T, sep="\t", dec=".",comment.char="")
cc_list<-preplotSNP(cc_FDT,1536)
library(rpart)
T_cc_SNP_T<-simpleSNPT(cc_SNP_T)
arbreDeci<-rpart(Aux~.,data=T_cc_SNP_T[,-c(1:2)],minsplit=1,xval=100)
plotcp(arbreDeci)
plot(arbreDeci,compress=T, margin=0.05,branch=1)
text(arbreDeci,use.n=T,cex=0.7)
choixExem<-listeExCART(arbreDeci,T_cc_SNP_T,cc_list)
summary(choixExem)
SNPlotCART(arbreDeci,choixExem,T_cc_SNP_T) #produit un pdf dans WorkingDirect 




