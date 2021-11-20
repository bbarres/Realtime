##############################################################################/
##############################################################################/
#Comparison of Fst between dead and survivor 
##############################################################################/
##############################################################################/

library(adegenet)
library(hierfstat)
library(poppr)


##############################################################################/
#using the dosage format for the snp data
##############################################################################/

compSNP<-read.snp("data/data.snp")
metaIND<-read.table("data/metadatIND.snp",sep="\t",header=T,
                    stringsAsFactors=TRUE,na.strings="..")
metaMARK<-read.table("data/metadatMARK.snp",sep="\t",header=T,
                     stringsAsFactors=TRUE,na.strings="..")
compSNP@loc.names<-as.character(metaMARK$loc.names)
#adding other information to the @other slot for individuals
indNames(compSNP)<-metaIND$nom_Ind
compSNP@other$famille<-metaIND$famille
compSNP@other$mother<-metaIND$mother
compSNP@other$father<-metaIND$father
compSNP@other$traitement<-metaIND$traitement
compSNP@other$traitsimp<-metaIND$trait_simp
compSNP@other$live_bin<-metaIND$live_bin


#limiting the data set to individuals of the experimental setting 
experSNP<-compSNP[compSNP$pop!="PAR" & compSNP$pop!="CC" 
                  & compSNP$pop!="A4" & compSNP$pop!="3P"]

#some individuals don't have the information of their death
experSNP@ind.names[experSNP@other$live_bin=="NA"]

#we remove these individuals until clarification
temp<-experSNP[experSNP@other$live_bin!="NA"]





##############################################################################/
#END
##############################################################################/