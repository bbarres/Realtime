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

(rowSums(table(evolPheno$grpbloc,evolPheno$an_mort))-table(evolPheno$grpbloc,evolPheno$an_mort))/rowSums(table(evolPheno$grpbloc,evolPheno$an_mort))

cumsum(table(evolPheno$grpbloc,evolPheno$an_mort)[1,])



##############################################################################/
#END
##############################################################################/