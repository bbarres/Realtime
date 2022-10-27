##############################################################################/
##############################################################################/
#Plot the map of the RealTime experiment
##############################################################################/
##############################################################################/

source("RealTime_load.R")


##############################################################################/
#map of the experimental set up####
##############################################################################/

#how many individuals in each treatment for each family
table(dispo$family_simp,dispo$trait)

colovec<-c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(1,4,8)])
colovec<-c(colovec[1:3],"#ffffff",colovec[4:15],"#ffffff","#000000",
           "#ffffff","#ffffff")
pdf(file="output/Figure_S2_map.pdf",width=20,height=12)
op<-par(mar=c(0,0,15,0))
plot(dispo$coord_X,dispo$coord_Y,bty="n",ann=FALSE,axes=FALSE,
     pch=c(rep(21,17),24,22)[as.numeric(dispo$family_simp)],
     bg=colovec[as.numeric(dispo$family_simp)])
text(x=c(20,620,1220,20,620,1220,20,620,1220),
     y=c(760,760,760,455,455,455,145,145,145),
     labels=c("A1 / natural","A2 / protected","A3 / natural",
              "B1 / natural","B2 / natural","B3 / protected",
              "C1 / protected","C2 / natural","C3 / natural"),
     xpd=TRUE,cex=2,font=2,adj=c(0,0))
legend(x=1430,y=1020,horiz=FALSE,x.intersp=1,xpd=TRUE,pt.cex=2,
       text.width=50,bty="n",ncol=4,
       legend=as.character(levels(dispo$family_simp))[c(1:3,5:16,18,19)],
       pch=c(rep(21,15),24,22),y.intersp=1.8,
       pt.bg=colovec[-c(4,17)],title="Families",title.cex=2)
title(main="Map of the experimental setup",cex.main=4,line=10)
par(op)
#export to .pdf 20 x 12 inches
dev.off()


##############################################################################/
#map of the dead tree####
##############################################################################/

colovec<-c(brewer.pal(12,"Paired")[c(rep(6,9),9,4)])
op<-par(mar=c(0,0,15,0))
plot(dispo$coord_X,dispo$coord_Y,bty="n",ann=FALSE,axes=FALSE,
     pch=c(rep(21,17),24,21,22)[as.numeric(dispo$family_simp)],
     bg=colovec[as.numeric(as.factor(dispo$an_mort))])
text(x=c(20,620,1220,20,620,1220,20,620,1220),
     y=c(760,760,760,455,455,455,145,145,145),
     labels=c("A1 / natural","A2 / protected","A3 / natural",
              "B1 / natural","B2 / natural","B3 / protected",
              "C1 / protected","C2 / natural","C3 / natural"),
     xpd=TRUE,cex=2,font=2,adj=c(0,0))
legend(x=230,y=950,horiz=FALSE,x.intersp=0.5,y.intersp=0.8,
       xpd=TRUE,pt.cex=3,bty="n",text.font=2,
       legend=c("Dead","Acorn","Alive"),
       pch=c(rep(15,3)),
       col=colovec[c(1,10,11)])
legend(x=60,y=950,horiz=FALSE,x.intersp=0.5,y.intersp=0.8,
       xpd=TRUE,pt.cex=3,bty="n",text.font=2,
       legend=c("Experiment","CC","hd"),
       pch=c(21,24,22))
title(main="Dead or alive map",cex.main=4,line=10)
par(op)
#export to .pdf 20 x 12 inches


#some individuals from the experimental setup are not included in the 
#phenotypic data set
dispo[dispo$family_simp!="CC" & dispo$family_simp!="hd" & 
              is.na(dispo$an_mort),"Sample_ID"]


##############################################################################/
#map of the tree height####
##############################################################################/

#defining a vector to chose the columns with tree height information
temp<-c("Hfin09","Hfin10","Hfin11","Hfin12",
        "Hdeb14","Hdeb15","Hdeb16","Hdeb17")
min(dispo[,temp])
cut(as.numeric(as.matrix(dispo[,temp])),49,na.rm=TRUE,levels=FALSE)

colovec<-viridis(2)
op<-par(mar=c(0,0,15,0))
plot(dispo$coord_X,dispo$coord_Y,bty="n",ann=FALSE,axes=FALSE,
     pch=c(rep(21,17),24,21,22)[as.numeric(dispo$family_simp)],
     bg=colovec[as.numeric(dispo$an_mort)])
text(x=c(20,620,1220,20,620,1220,20,620,1220),
     y=c(760,760,760,455,455,455,145,145,145),
     labels=c("A1 / natural","A2 / protected","A3 / natural",
              "B1 / natural","B2 / natural","B3 / protected",
              "C1 / protected","C2 / natural","C3 / natural"),
     xpd=TRUE,cex=2,font=2,adj=c(0,0))
legend(x=300,y=950,horiz=FALSE,x.intersp=1,y.intersp=0.5,
       xpd=TRUE,pt.cex=3,bty="n",text.font=2,
       legend=c("Dead","Acorn","Alive"),
       pch=c(rep(15,3)),
       col=colovec[c(1,10,11)])
legend(x=60,y=950,horiz=FALSE,x.intersp=1,y.intersp=0.5,
       xpd=TRUE,pt.cex=3,bty="n",text.font=2,
       legend=c("Experiment","CC","hd"),
       pch=c(21,24,22))
title(main="Dead or alive map",cex.main=4,line=10)
par(op)
#export to .pdf 20 x 12 inches


##############################################################################/
#END
##############################################################################/