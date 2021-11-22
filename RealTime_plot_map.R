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
op<-par(mar=c(0,0,15,0))
plot(dispo$coord_X,dispo$coord_Y,bty="n",ann=FALSE,axes=FALSE,
     pch=c(rep(21,17),24,21,22)[as.numeric(dispo$family_simp)],
     bg=colovec[as.numeric(dispo$family_simp)])
text(x=c(20,620,1220,20,620,1220,20,620,1220),
     y=c(760,760,760,455,455,455,145,145,145),
     labels=c("A1 / nat","A2 / lim","A3 / renf",
              "B1 / renf","B2 / nat","B3 / lim",
              "C1 / lim", "C2 / renf","C3 / nat"),
     xpd=TRUE,cex=2,font=2,adj=c(0,0))
legend(x=300,y=900,horiz=TRUE,x.intersp=0.2,xpd=TRUE,pt.cex=2,
       text.width=seq(from=0,to=1,by=0.06),bty="n",
       legend=as.character(levels(dispo$family_simp))[c(1:3,5:16,18,20)],
       pch=c(rep(21,15),24,22),
       pt.bg=colovec[-c(4,17)])
title(main="Map of the experimental setup",cex.main=4,line=10)
par(op)
#export to .pdf 20 x 12 inches


##############################################################################/
#map of the dead tree####
##############################################################################/

colovec<-c(brewer.pal(12,"Paired")[c(rep(6,9),9,4)])
op<-par(mar=c(0,0,15,0))
plot(dispo$coord_X,dispo$coord_Y,bty="n",ann=FALSE,axes=FALSE,
     pch=c(rep(21,17),24,21,22)[as.numeric(dispo$family_simp)],
     bg=colovec[as.numeric(dispo$an_mort)])
text(x=c(20,620,1220,20,620,1220,20,620,1220),
     y=c(760,760,760,455,455,455,145,145,145),
     labels=c("A1 / nat","A2 / lim","A3 / renf",
              "B1 / renf","B2 / nat","B3 / lim",
              "C1 / lim", "C2 / renf","C3 / nat"),
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
#map of the tree height####
##############################################################################/

#defining a vector to chose the columns with tree height information
temp<-c("Hfin09","Hfin10","Hfin11","Hfin12",
        "Hdeb14","Hdeb15","Hdeb16","Hdeb17")
min(dispo[,temp])


colovec<-viridis(2)
op<-par(mar=c(0,0,15,0))
plot(dispo$coord_X,dispo$coord_Y,bty="n",ann=FALSE,axes=FALSE,
     pch=c(rep(21,17),24,21,22)[as.numeric(dispo$family_simp)],
     bg=colovec[as.numeric(dispo$an_mort)])
text(x=c(20,620,1220,20,620,1220,20,620,1220),
     y=c(760,760,760,455,455,455,145,145,145),
     labels=c("A1 / nat","A2 / lim","A3 / renf",
              "B1 / renf","B2 / nat","B3 / lim",
              "C1 / lim", "C2 / renf","C3 / nat"),
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


#map of the different families
plan<-read.table("data/tout.txt")
erreur<-read.table("data/etiq_dif_microsat.txt")
bord<-read.table("data/bordure_nonechant.txt")
plot(plan, pch=19, bg="blue")
points(erreur, pch=19, bg="red", col="red")
points(bord, pch=19, bg="green4", col="green4")

#map of the flagshoot for years 2010 and 2011
plan<-read.table("data/tout.txt")
bord<-read.table("data/bordures.txt")
drap2010<-read.table("data/drap2010.txt")
drap2011<-read.table("data/drap2011.txt")
plot(plan, pch=21, col="black")
points(bord, pch=21, bg="green4", col="black")
points(drap2010, pch=20, bg="blue", col="red2")
points(drap2011, pch=20, bg="orange", col="orange")

#map of the individuals of the 3P family
plan<-read.table("data/tout.txt")
bord<-read.table("data/bordures.txt")
fam_3P<-read.table("data/3P_fam.txt")
plot(plan, pch=21, col="black")
points(bord, pch=21, bg="green4", col="black")
points(fam_3P, pch=20, bg="blue", col="orange")

#export to .pdf 17 x 10 inches