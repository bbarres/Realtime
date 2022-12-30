[![DOI](https://zenodo.org/badge/41293576.svg)](https://zenodo.org/badge/latestdoi/41293576) 
# Significant demographic and genetic impacts of powdery mildew in a young oak cohort 
*This repository contains the R code and SAS scripts used for the data analyses and production of the figures of the related article*

![alt text](https://am3pap005files.storage.live.com/y4mXUw4rIo7I6I-UHrfoPc32YcaZYoOIM_R2-WK8ZLDnMrfurYfJ5FV2WjlTIh_idbaCJMyDFBGPAnA3lwmFOaY_6M5ra2vfzfOiobz1ENwdeA1QbGCTFYXnkZznKUDZXRVNRsXyB7KCZzHJLIFI1B8rqivCd0_12NtVbUs-X7a5FWMEHdFaMUfnwWvsiDjm8JU?width=1584&height=588&cropmode=none)


## Context



## R Datasets

In this section, you will find the list of the data sets used in this study for R analyses and Figures production. The data files can be found in the "data" folder. For the data tables, the name of the different variables are listed and explained as well. There are 5 data sets used in this study.

- **datatot.txt:** the main data set which contains all the . Each line correspond to one individuals and the following information for each individuals can be found in this table:
    - *Sample_ID*: unique identification number for each sample
    - *family*: which family the individual belongs to. Only the 15 principal families with a number as a family ID ("1", "9", "11", "27", "45", "48", "51", "70", "71", "72", "73", "74", "75", "76", "77") are of interest for this study. The family has been checked using a parentage analysis to confirm the mother tree when genetic data were available. If no genetic data were available, this is simply the original family ID of the acorn. Using genotyping tools, a few individuals were identified as not belonging to one of the 15 original families ("26" and "99"). Individuals used to surround the individuals included in the experimental design to limit boundary effects are labeled "hd". Individuals labelled "3P", "3P tv", "A4", "cc" and "cc-cast" belong to a controlled crossing not investigated in this study. Individuals labeled "asc" are potential father trees. 
    - *family_simp*: same information as the "family" variable, but with simplified coding (all the controlled crossing individuals are labeled "CC" and the "asc" individuals are labeled "PAR")
    - *parent_id*: ID of the identified mother for the individual
    - *bloc*: block ID in the experimental design
    - *PU*: unitary plot ID
    - *rg*: rank number within the unitary plot. There are 12 rank per block
    - *n*: individual number on the rank. There are a maximum of 48 position on a rank
    - *coord_X*: X coordinates in the general experimental setting, used for producing the map of the experimental setting in the supplementary material
    - *coord_Y*: Y coordinates in the general experimental setting, used for producing the map of the experimental setting in the supplementary material
    - *na.9micro*: number of missing data for an individual for the 9 best microsatellites
    - *na.12micro*: number of missing data for an individual for the 12 microsatellites
    - *barcode*: barcode ID for the GoldenGate Assay
    - *A11_all1* to *S19_all2* (24 columns): each pairs of column give the genotype for a microsatellite marker. The number represent the allele size in base pair
    - *nb_conta_12*: number of contaminations for the 12 microsatellite markers. A contamination is counted when 3 or more alleles have been detected for a locus
    - *nb_conta_09*: number of contaminations for the 9 best microsatellite markers. A contamination is counted when 3 or more alleles have been detected for a locus
    - *SNPage*: is the individual have been included in the GoldenGate Assay (1=included; 0=not included)
    - *Quality_SNPage*: for the GoldenGate genotyped individuals, was the quality of genotyping good or not (1=good quality; 0=poor quality)
    - *pb_robot_SNPage*: for some individuals, the volume of DNA added to the PCR reaction was reduced due to a problem with the pipetting robot (1=pipetting problem; 0=no pipetting problem)
    - *mother.snp*: the mother ID based on a parentage analysis on the 819 SNP information. If one the parental ID match the original ID of the acorn, the mother ID is kept. If the parentage analysis identify another mother tree, the ID is changed to the genetic match. If the parentage analysis doesn't identify one of the original mother tree, the ID is set to missing (="NA")
    - *father.snp*: the father ID based on a parentage analysis on the 819 SNP information. It can only be one of the mother tree or one of the adult tree from the vicinity of the mother trees that were sampled while collecting the acorn
    - *live_bin*: is the juvenile still alive in 2017 ? ("1"=yes; "0"=no)
    - *live_year*: year of the death of the juvenile (from "2009"" to "2017"; "vivant" if the juvenile was still alive in 2017)
    - *na.snp*: total number of missing information for all 819 selected SNP
    - *REF\~CL371CT472_02-163* to *REF\~CL9715CT16278_02-801* (819 columns): the allele combination for each of the 819 SNPs successfully genotyped. Oak being diploids, the genotype is indicated by 2 letters coding for the nucleotide (one for each copy of the genome) separated by a slash
    - *trait_init*:
    - *phenP2_09*:
    - *nbp09*:
    - *oid1_09*:
    - *oid2_09*:
    - *oid3_09*:
    - *oid4_09*:
    - *nec1_09*:
    - *nec2_09*:
    - *herbi1_09*:
    - *herbi2_09*:
    - *herbi3_09*:
    - *herbi4_09*:
    - *phen1_10*:
    - *phen2_10*:
    - *phen3_10*:
    - *phen4_10*:
    - *phen5_10*:
    - *phenP2_10*:
    - *nbp10*:
    - *mortapex10*:
    - *oid1_10*:
    - *oid2_10*:
    - *oid3_10*:
    - *oid3P1_10*:
    - *oid3P2_10*:
    - *oid4_10*:
    - *oid5_10*:
    - *herbi1_10*:
    - *herbi2_10*:
    - *nbp11*:
    - *phen1_11*:
    - *phen2_11*:
    - *phen3_11*:
    - *phen4_11*:
    - *phenP2_11*:
    - *phenP3_11*:
    - *oid1_11*:
    - *oid2_11*:
    - *oid3_11*:
    - *oid4_11*:
    - *oid4P2_11*:
    - *oid5_11*:
    - *herbi_11*:
    - *phen1_12*:
    - *phen2_12*:
    - *phen3_12*:
    - *phen4_12*:
    - *phen5_12*:
    - *oid1_12*:
    - *oid2_12*:
    - *oid3_12*:
    - *oid4_12*:
    - *nbp_12*:
    - *phen1_13*:
    - *gel_13*:
    - *oid1_13*:
    - *oid2_13*:
    - *oid3_13*:
    - *oid1_14*:
    - *oid2_14*:
    - *oid_16*:
    - *oid_17*:
    - *drap10*:
    - *drap11*:
    - *drap12*:
    - *drap13*:
    - *drap15*:
    - *freq_drap*:
    - *pgland*: mass of the acorn (grams)
    - *date_em*: date of acorn raising
    - *Hfin09*:
    - *Hfin10*:
    - *Hfin11*:
    - *Hfin12*:
    - *Hdeb14*:
    - *Hdeb15*:
    - *Hdeb16*:
    - *Hdeb17*:
    - *diam16*:
    - *an_mort*:
    - *H09v*:
    - *H10v*:
    - *H11v*:
    - *H12v*:
    - *H14v*:
    - *H15v*:
    - *H16v*:
    - *H17v*:
    - *cr10*:
    - *cr11*:
    - *cr12*:
    - *cr13*:
    - *cr14*:
    - *cr15*:
    - *cr16*:
    - *esp*:
    - *statut*:
    - *treat*:
    - *oid_moy*:
    - *statut10*:
    
- **lim.hmp.txt.txt** and **nat.hmp.txt:** two data sets formatted for the GWAS analyses, for the protected and natural treatment respectively
    - *indiv_ID*:


## R scripts

In this section, you will find the list of the different scripts used in the article with a brief description of their purpose.

- **RealTime_load.R:** the script to load the different data sets, functions and packages that are necessary for the data analyses and representation in the R environment.
- **RealTime_evolpheno.R:** script to analyze and plot the evolutionary patterns of the main phenotypic traits over the years.
- **RealTime_plot_SurvShannon.R:** script to compute the Shannon index evolution and to plot the Figure of the survival and Shannon index evolution.
- **RealTime_Fstcomparison.R:** script to compute F-statistics and compare them across the different treatment and time. The output files are stored in a Genepop folder created within the output folder.
- **RealTime_GENHET.R:** script to compute intra-individual heterozygosity indices. The code to produce the related Figure is also included.
- **RealTime_GWAS.R:** script to perform the Genome-Wide Association analyses for both treatment. This script will produce two folders each containing results from the GWAS analyses. These output files are necessary for plotting the manhattan and GWAS plot.
- **RealTime_plot_manhattan.R:** script to plot the manhattan plot Figure.
- **RealTime_plot_GWAS.R:** script to plot the GWAS results for the SNP of interest.


## R session info

For reproducibility purpose, you will find all the information about the versions of R, Rstudio, OS etc., as well as the list and version number of the packages used at the time of publishing this script in the **session_info.txt** file.


## SAS Datasets

In this section, you will find the list of the data sets used in this study for SAS analyses and Figures production. The data files can be found in the "dataSAS" folder. For the data tables, the name of the different variables are listed and explained as well. There are 5 data sets used in this study.


## SAS scripts

In this section, you will find the list of the different scripts used in the article with a brief description of their purpose.


## Citation

You will be soon (hopefully) able to cite the related study as follow: 
+ Barrès B., Saint-Jean G., Lepoittevin C., Burban C., Garnier-Géré P., Dutech C. and Desprez-Loustau M.-L. [Significant demographic and genetic impacts of powdery mildew in a young oak cohort. *journal name*.](https://)

If you want to use (some of) the code found on this page or if you want to cite this repository: 
+ Benoit Barrès and Marie-Laure Desprez-Loustau. [Supporting data and code for: Significant demographic and genetic impacts of powdery mildew in a young oak cohort. Zenodo.](https://zenodo.org/badge/latestdoi/sss)
