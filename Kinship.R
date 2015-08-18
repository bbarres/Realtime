###############################################################################
###############################################################################
#Preparing Kinship matrices for the Mixed linear model analysis
###############################################################################
###############################################################################

#loading the libraries
library(kinship2)

#set the working directory
setwd("~/work/Rfichiers/Githuber/Realtime_data")


###############################################################################
#Kinship matrix based on the pedigree
###############################################################################

pedifile<-read.table("pedfile.txt",header=T)
kinMORG<-2*kinship(id=pedifile$id,mom=pedifile$mother,dad=pedifile$father)
write.table(kinMORG,file="kinMORG.txt",sep='\t')

#you can also compute matrices for subsets of the entire dataset
pedifile_nat<-read.table("pedfile_nat.txt",header=T)
kinMORG_nat<-2*kinship(id=pedifile_nat$id,mom=pedifile_nat$mother,
                       dad=pedifile_nat$father)
write.table(kinMORG_nat,file="kinMORG_nat.txt",sep='\t')

pedifile_lim<-read.table("pedfile_lim.txt",header=T)
kinMORG_lim<-2*kinship(id=pedifile_lim$id,mom=pedifile_lim$mother,
                       dad=pedifile_lim$father)
write.table(kinMORG_lim,file="kinMORG_lim.txt",sep='\t')

pedifile_renf<-read.table("pedfile_renf.txt",header=T)
kinMORG_renf<-2*kinship(id=pedifile_renf$id,mom=pedifile_renf$mother,
                        dad=pedifile_renf$father)
write.table(kinMORG_renf,file="kinMORG_renf.txt",sep='\t')


###############################################################################
#Kinship matrix based on the pedigree
###############################################################################

#Here are some codes to transform the output of the Spagedi 1.5a software into 
#a matrix suitable for the analysis in TASSEL. For that you just need to run 
#Spagedi choosing the relatedness index you are interested in and then extract 
#from the output file the matrix of pairwise kinship coefficient and the list 
#of inbreeding coefficient of each individual. We then combine this two files, 
#putting in the diagonal of the matrix 0.5*(1+inbreeding coefficient)

kinLOISlim<-read.table("kin_snp_LOIS_nat.txt",header=T,sep="\t")

#il faut aussi récupérer le coefficient d'inbreeding dans le même fichier 
#de sortie SpaGeDi. D'abord faire le calcul 0.5*(1+inbreed) et intituler 
#la nouvelle colonne 'intrakin'
kinintra_nat<-read.table("kin_snp_intraLOIS_nat.txt",header=T,sep="\t")

#voilà les transfo qu'il faut faire
dim(kinLOISnat)[1]
kinintra_nat$intrakin[kinintra_nat$Name.i==rownames(kinLOISnat)[2]]
for (i in 1:dim(kinLOISnat)[1]) {
  kinLOISnat[i,i]<-kinintra_nat$intrakin[kinintra_nat$Name.i==rownames(kinLOISnat)[i]]
}
ifelse(min(kinLOISnat)<0,kinLOISnat<-kinLOISnat+abs(min(kinLOISnat)),kinLOISnat)

#maintenant voici la fonction qui va mettre la valeur intra individu dans 
#la diagonale. En plus cette fonction va éliminer les valeurs nulles de la 
#matrice en recherchant et en ajoutant la valeur absolue de la valeur min 
#de la matrice à tous les cellules de la matrice. En entrée il faut les 
#deux fichiers décrits ci-dessus qui sont arrangés à partir de la sortie 
#SpaGeDi
insertinbreed<-function(kin,inbreed)
{
  for (i in 1:dim(kin)[1]) {
    kin[i,i]<-inbreed$intrakin[inbreed$Name.i==rownames(kin)[i]]
  }
  ifelse(min(kin)<0,kin<-kin+abs(min(kin)),kin)
  return(kin)
}  

#démonstration de l'utilisation de cette fonction, d'abord SNP NAT
kinLOISnat<-read.table("kin_snp_LOIS_nat.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_snp_intraLOIS_nat.txt",header=T,sep="\t")
transLOISnat<-2*insertinbreed(kinLOISnat,kinintra_nat)
write.table(transLOISnat,file="snp_NAT_kinLOIS.txt")

kinRITLnat<-read.table("kin_snp_RITL_nat.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_snp_intraRITL_nat.txt",header=T,sep="\t")
transRITLnat<-2*insertinbreed(kinRITLnat,kinintra_nat)
write.table(transRITLnat,file="snp_NAT_kinRITL.txt")


#SNP LIM
kinLOISlim<-read.table("kin_snp_LOIS_lim.txt",header=T,sep="\t")
kinintra_lim<-read.table("kin_snp_intraLOIS_lim.txt",header=T,sep="\t")
transLOISlim<-2*insertinbreed(kinLOISlim,kinintra_lim)
write.table(transLOISlim,file="snp_LIM_kinLOIS.txt")

kinRITLlim<-read.table("kin_snp_RITL_lim.txt",header=T,sep="\t")
kinintra_lim<-read.table("kin_snp_intraRITL_lim.txt",header=T,sep="\t")
transRITLlim<-2*insertinbreed(kinRITLlim,kinintra_lim)
write.table(transRITLlim,file="snp_LIM_kinRITL.txt")

#SNP RENF
kinLOISrenf<-read.table("kin_snp_LOIS_renf.txt",header=T,sep="\t")
kinintra_renf<-read.table("kin_snp_intraLOIS_renf.txt",header=T,sep="\t")
transLOISrenf<-2*insertinbreed(kinLOISrenf,kinintra_renf)
write.table(transLOISrenf,file="snp_RENF_kinLOIS.txt")

kinRITLrenf<-read.table("kin_snp_RITL_renf.txt",header=T,sep="\t")
kinintra_renf<-read.table("kin_snp_intraRITL_renf.txt",header=T,sep="\t")
transRITLrenf<-2*insertinbreed(kinRITLrenf,kinintra_renf)
write.table(transRITLrenf,file="snp_RENF_kinRITL.txt")

#SAT NAT, vu les résultats sur les SNP on ne garde que LOISEL, même si le 
#RITLAND est quasi aussi bon. Quand il y a un individu avec des données 
#complétement manquantes, on complète par des 0 aussi bien dans la matrice 
#que pour l'inbreed, donc 0.5 pour le kinship intra (0.5*(1+inb))

kinLOISnat<-read.table("kin_sat_LOIS_nat.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_sat_intraLOIS_nat.txt",header=T,sep="\t")
transLOISnat<-2*insertinbreed(kinLOISnat,kinintra_nat)
write.table(transLOISnat,file="sat_NAT_kinLOIS.txt")

#calcul pour voir si on prend un sous ensemble de SNP référence vs candidat, 
#est-ce que ça a un impact sur la correction ? Pour simplifier on ne fait 
#ça que sur les données nat pour l'instant et qu'avec l'apparentement de 
#LOISEL

#d'abord gènes 'références'
kinLOISnat<-read.table("kin_snp_LOIS_natREF.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_snp_intraLOIS_natREF.txt",header=T,sep="\t")
transLOISnatREF<-2*insertinbreed(kinLOISnat,kinintra_nat)
write.table(transLOISnatREF,file="snp_NAT_REF_kinLOIS.txt")

#puis les gènes 'candidats'
kinLOISnat<-read.table("kin_snp_LOIS_natPATHO.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_snp_intraLOIS_natPATHO.txt",header=T,sep="\t")
transLOISnatPATHO<-2*insertinbreed(kinLOISnat,kinintra_nat)
write.table(transLOISnatPATHO,file="snp_NAT_PATHO_kinLOIS.txt")


#rebelotte pour les jeux de données criblés avec les MAF

kinRITLnat<-read.table("kin_snp_RITL_nat.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_snp_intraRITL_nat.txt",header=T,sep="\t")
transRITLnat<-2*insertinbreed(kinRITLnat,kinintra_nat)
write.table(transRITLnat,file="snp_NAT_kinRITL.txt")

kinLOISnat<-read.table("kin_snp_LOIS_nat.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_snp_intraLOIS_nat.txt",header=T,sep="\t")
transLOISnat<-2*insertinbreed(kinLOISnat,kinintra_nat)
write.table(transLOISnat,file="snp_NAT_kinLOIS.txt")

kinRITLnat<-read.table("kin_sat_RITL_nat.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_sat_intraRITL_nat.txt",header=T,sep="\t")
transRITLnat<-2*insertinbreed(kinRITLnat,kinintra_nat)
write.table(transRITLnat,file="sat_NAT_kinRITL.txt")

kinRITLnat<-read.table("kin_snp_PATHO_RITL_nat.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_snp_intraPATHO_RITL_nat.txt",header=T,sep="\t")
transRITLnat<-2*insertinbreed(kinRITLnat,kinintra_nat)
write.table(transRITLnat,file="snp_NAT_kinPATHO_RITL.txt")

kinRITLnat<-read.table("kin_snp_REF_RITL_nat.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_snp_intraREF_RITL_nat.txt",header=T,sep="\t")
transRITLnat<-2*insertinbreed(kinRITLnat,kinintra_nat)
write.table(transRITLnat,file="snp_NAT_kinREF_RITL.txt")

kinLOISnat<-read.table("kin_snp_PATHO_LOIS_nat.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_snp_intraPATHO_LOIS_nat.txt",header=T,sep="\t")
transLOISnat<-2*insertinbreed(kinLOISnat,kinintra_nat)
write.table(transLOISnat,file="snp_NAT_kinPATHO_LOIS.txt")

kinLOISnat<-read.table("kin_snp_REF_LOIS_nat.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_snp_intraREF_LOIS_nat.txt",header=T,sep="\t")
transLOISnat<-2*insertinbreed(kinLOISnat,kinintra_nat)
write.table(transLOISnat,file="snp_NAT_kinREF_LOIS.txt")

#rerebelotte pour les jeux de données criblés avec les MAF + LD, attention 
#les noms de fichiers sont parfois identiques à ce que l'on trouve au-dessus, 
#il faut donc faire attention de se trouver dans le bon dossier

kinRITLnat<-read.table("kin_snp_RITL_nat.txt",header=T,sep="\t")
kinintra_nat<-read.table("kin_snp_intraRITL_nat.txt",header=T,sep="\t")
transRITLnat<-2*insertinbreed(kinRITLnat,kinintra_nat)
write.table(transRITLnat,file="snp_NAT_kinRITL.txt")

kinRITLlim<-read.table("kin_snp_RITL_lim.txt",header=T,sep="\t")
kinintra_lim<-read.table("kin_snp_intraRITL_lim.txt",header=T,sep="\t")
transRITLlim<-2*insertinbreed(kinRITLlim,kinintra_lim)
write.table(transRITLlim,file="snp_LIM_kinRITL.txt")

kinRITLrenf<-read.table("kin_snp_RITL_renf.txt",header=T,sep="\t")
kinintra_renf<-read.table("kin_snp_intraRITL_renf.txt",header=T,sep="\t")
transRITLrenf<-2*insertinbreed(kinRITLrenf,kinintra_renf)
write.table(transRITLrenf,file="snp_RENF_kinRITL.txt")
