##############################################################################/
##############################################################################/
#Plot the map of the RealTime experiment
##############################################################################/
##############################################################################/

##############################################################################/
#load the data####
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


##############################################################################/
#END
##############################################################################/