##############################################################################/
##############################################################################/
#Plot the map of the RealTime experiment
##############################################################################/
##############################################################################/

source("RealTime_load.R")


##############################################################################/
#evolution of the phenotypic traits####
##############################################################################/

#limiting the dataset to usefull data
evolPheno<-dispo[dispo$family_simp!="CC" & dispo$family_simp!="hd" & 
                   !is.na(dispo$an_mort) & dispo$an_mort!="g",
                 c("Sample_ID","bloc","PU","treat","an_mort","family_simp")]
evolPheno$grpbloc<-paste(evolPheno$bloc,evolPheno$PU,evolPheno$treat,sep="")
#defining a color vector
colovec<-c(brewer.pal(12,"Paired")[4],brewer.pal(12,"Paired")[2])

temp<-t(apply(table(evolPheno$grpbloc,evolPheno$an_mort),1,cumsum))
temp<-100*(temp[,10]-temp)/temp[,10]
temp<-temp[,-10]
temp.bloc<-factor(substr(rownames(temp),1,2))
temp.treat<-factor(substr(rownames(temp),3,5))
matplot(t(temp),type="b",ylim=c(0,100),las=1,
        col=colovec[as.numeric(temp.treat)],lwd=2,
        lty=as.numeric(temp.treat),
        pch=as.numeric(temp.bloc)+14,
        xaxt="n",yaxt="n",bty="n",
        ylab="Percentage of survival",xlab="Year",
        font.lab=2,cex.lab=1.5)
axis(1,lwd=2,cex.axis=1,at=c(1:dim(temp)[2]),
     lab=colnames(temp),las=1,font=2)
axis(2,lwd=2,las=1,font=2)
box(lwd=2)
#export to .pdf 6 x 6 inches


##############################################################################/
#END
##############################################################################/