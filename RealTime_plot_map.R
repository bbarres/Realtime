###############################################################################
###############################################################################
#Plot the map of the RealTime experiment
###############################################################################
###############################################################################

#set the working directory
setwd("~/work/Rfichiers/Githuber/Realtime_data")


###############################################################################
#load the data
###############################################################################


#ce premier plan a servit ? identifier s'il y avait des soucis sur le dispositif 
#en v?rifiant que la famille suppos?e correspondait ? ce que l'on trouvait en 
#microsatellites

plan<-read.table("tout.txt")
erreur<-read.table("etiq_dif_microsat.txt")
bord<-read.table("bordure_nonechant.txt")
plot(plan, pch=19, bg="blue")
points(erreur, pch=19, bg="red", col="red")
points(bord, pch=19, bg="green4", col="green4")


#un plan pour la localisation des drapeaux sur les ann?es 2010 et 2011

plan<-read.table("tout.txt")
bord<-read.table("bordures.txt")
drap2010<-read.table("drap2010.txt")
drap2011<-read.table("drap2011.txt")
plot(plan, pch=21, col="black")
points(bord, pch=21, bg="green4", col="black")
points(drap2010, pch=20, bg="blue", col="red2")
points(drap2011, pch=20, bg="orange", col="orange")


#un plan pour la localisation des individus de la famille 3P
plan<-read.table("tout.txt")
bord<-read.table("bordures.txt")
fam_3P<-read.table("3P_fam.txt")
plot(plan, pch=21, col="black")
points(bord, pch=21, bg="green4", col="black")
points(fam_3P, pch=20, bg="blue", col="orange")


###############################################################################
#END
###############################################################################