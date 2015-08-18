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

#here is the function to insert the inbreeding coefficient in the kinship 
#matrix
insertinbreed<-function(kin,inbreed)
{
  for (i in 1:dim(kin)[1]) {
    kin[i,i]<-inbreed$intrakin[inbreed$Name.i==rownames(kin)[i]]
  }
  ifelse(min(kin)<0,kin<-kin+abs(min(kin)),kin)
  return(kin)
}


#We process the file for the 'lim' treatment, using a kinship computed with 
#610 SNPs with a LD lower than 0.8 and a MAF>0.05
#first we load the kinship matrix
kinLOISlim<-read.table("lim_ldmaf_LOIS.txt",header=TRUE,row.names=1,sep="\t")
#then we load the inbreeding coefficient and turn it into 0.5*(1+inbreed coeff)
kinintra_lim<-read.table("lim_ldmaf_intraLOIS.txt",header=T,sep="\t")
kinintra_lim<-cbind(kinintra_lim,"intrakin"=0.5*(1+kinintra_lim$ALL.LOCI))
#finally we used the function to include the inbreeding coefficient value in 
#the matrix and export the dataset to be used in TASSEL
transLOISlim<-2*insertinbreed(kinLOISlim,kinintra_lim)
write.table(transLOISlim,file="snp_LIM_kinLOIS.txt")

#We process the file for the 'lim' treatment, using a kinship computed with 
#610 SNPs with a LD lower than 0.8 and a MAF>0.05
#first we load the kinship matrix
kinLOISnatrenf<-read.table("natrenf_ldmaf_LOIS.txt",header=TRUE,
                           row.names=1,sep="\t")
#then we load the inbreeding coefficient and turn it into 0.5*(1+inbreed coeff)
kinintra_natrenf<-read.table("natrenf_ldmaf_intraLOIS.txt",header=T,sep="\t")
kinintra_natrenf<-cbind(kinintra_natrenf,
                        "intrakin"=0.5*(1+kinintra_natrenf$ALL.LOCI))
#finally we used the function to include the inbreeding coefficient value in 
#the matrix and export the dataset to be used in TASSEL
transLOISnatrenf<-2*insertinbreed(kinLOISnatrenf,kinintra_natrenf)
write.table(transLOISnatrenf,file="snp_NATRENF_kinLOIS.txt")


###############################################################################
#END
###############################################################################