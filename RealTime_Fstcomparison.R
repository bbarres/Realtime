##############################################################################/
##############################################################################/
#Comparison of Fst between dead and survivor 
##############################################################################/
##############################################################################/


library(adegenet)
first<-read.snp("data/data.snp")

meta<-read.table("data/metadatIND.snp",sep="\t",header=T,
                 stringsAsFactors=TRUE,na.strings="..")

indNames(first)<-meta$nom_Ind
first@pop<-meta$famille
first@other$mother<-meta$mother
first@other$father<-meta$Father
first@other$traitement<-meta$traitement
first@loc.names<-as.character(meta_mark$loc.names)
second<-first[first$pop!="PAR" & first$pop!="CC" & first$pop!="A4"]


first$pop


meta




##############################################################################/
#END
##############################################################################/