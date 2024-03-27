[![DOI](https://zenodo.org/badge/33980368.svg)](https://zenodo.org/badge/latestdoi/33980368)
# Significant demographic and genetic impacts of powdery mildew in a young oak cohort
*This repository contains the R code and SAS scripts used for the data analyses and production of the figures of the related article*

![alt text](https://am3pap005files.storage.live.com/y4mXUw4rIo7I6I-UHrfoPc32YcaZYoOIM_R2-WK8ZLDnMrfurYfJ5FV2WjlTIh_idbaCJMyDFBGPAnA3lwmFOaY_6M5ra2vfzfOiobz1ENwdeA1QbGCTFYXnkZznKUDZXRVNRsXyB7KCZzHJLIFI1B8rqivCd0_12NtVbUs-X7a5FWMEHdFaMUfnwWvsiDjm8JU?width=1584&height=588&cropmode=none)


## Context
The impact of pathogens on their host populations in the wild can be manifold, but remains a relatively unexplored subject. Studying this type of impact in a wild system is complex, not least because of the number of factors (in addition to the pathogen studied) that can influence the demography of the host species: herbivory, competition with other species, other pathogens... In this study, we used a semi-controlled experimental set-up to mimic a regeneration cohort of young oaks (*Quercus robur*), in order to study the impact of oak powdery mildew (*Eysiphe alphitoides* and *E. quercicola*) on the evolution of the composition of this cohort. The system included two treatments: natural infection and infection controlled by a fungicide treatment. Acorns were planted and emerging oaks were monitored for 9 years. Numerous phenotypic traits were recorded annually, and a random sub-sample of planted individuals were genotyped for 819 SNPs identified on genes of interest. This enabled us to assess the impact of oak powdery mildew pressure on cohort demography and diversity over time. 


## R Data sets (/data)
In this section, you will find the list of the data sets used in this study for R analyses and Figures production. The data files can be found in the "data" folder. For the data tables, the name of the different variables are listed and explained as well. There are 3 R data sets used in this study.

- **datatot.txt:** the main data set which contains all the phenotypic and genotypic information for all the individuals of the experimental setup. Each line correspond to one individuals and the following information for each individuals can be found in this table:
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
    - *live_bin*: is the seedling still alive in 2017 ? ("1"=yes; "0"=no)
    - *live_year*: year of the death of the seedlings (from "2009"" to "2017"; "vivant" if the seedling was still alive in 2017)
    - *na.snp*: total number of missing information for all 819 selected SNP
    - *REF\~CL371CT472_02-163* to *REF\~CL9715CT16278_02-801* (819 columns): the allele combination for each of the 819 SNPs successfully genotyped. Oak being diploids, the genotype is indicated by 2 letters coding for the nucleotide (one for each copy of the genome) separated by a slash
    - *trait_init*: initial coding of the treatment of the plot ("nat" for natural treatment, "renf" for treatment with additional inoculum and "lim" for plot treated with a fungicide)
    - *phenP2_09*: scoring of the phenology of the second shoot of the individual on a scale from 0 to 4 for 2009 (NA for missing information)
    - *nbp09*: number of shoots in 2009, ranging from 1 to 4 (NA for missing information)
    - *oid1_09*: first powdery mildew infection note for 2009 (NA for missing information)
    - *oid2_09*: second powdery mildew infection note for 2009 (NA for missing information)
    - *oid3_09*: third powdery mildew infection note for 2009 (NA for missing information)
    - *oid4_09*: fourth powdery mildew infection note for 2009 (NA for missing information)
    - *nec1_09*: first note of leaf necrosis for the year 2009 (NA for missing information)
    - *nec2_09*: second note of leaf necrosis for the year 2009 (NA for missing information)
    - *herbi1_09*: first herbivore damage note for 2009 (NA for missing information)
    - *herbi2_09*: second herbivore damage note for 2009 (NA for missing information)
    - *herbi3_09*: third herbivore damage note for 2009 (NA for missing information)
    - *herbi4_09*: fourth herbivore damage note for 2009 (NA for missing information)
    - *phen1_10*: first scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2010 (NA for missing information)
    - *phen2_10*: second scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2010 (NA for missing information)
    - *phen3_10*: third scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2010 (NA for missing information)
    - *phen4_10*: fourth scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2010 (NA for missing information)
    - *phen5_10*: fifth scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2010 (NA for missing information)
    - *phenP2_10*: scoring of the phenology of the second shoot of the individual on a scale from 0 to 4 for 2010 (NA for missing information)
    - *nbp10*: number of shoots in 2010, ranging from 1 to 5 (NA for missing information)
    - *mortapex10*: observed apical mortality in 2010 (NA for missing information)
    - *oid1_10*: first powdery mildew infection note for 2010 (NA for missing information)
    - *oid2_10*: second powdery mildew infection note for 2010 (NA for missing information)
    - *oid3_10*: third powdery mildew infection note for 2010 (NA for missing information)
    - *oid3P1_10*: third powdery mildew infection note for 2010 on the first shoot only (NA for missing information)
    - *oid3P2_10*: third powdery mildew infection note for 2010 on the second shoot only (NA for missing information)
    - *oid4_10*: fourth powdery mildew infection note for 2010 (NA for missing information)
    - *oid5_10*: fifth powdery mildew infection note for 2010 (NA for missing information)
    - *herbi1_10*: first herbivore damage note for 2010 (NA for missing information)
    - *herbi2_10*: second herbivore damage note for 2010 (NA for missing information)
    - *nbp11*: number of shoots in 2011, ranging from 1 to 4 (NA for missing information)
    - *phen1_11*: first scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2011 (NA for missing information)
    - *phen2_11*: second scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2011 (NA for missing information)
    - *phen3_11*: third scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2011 (NA for missing information)
    - *phen4_11*: fourth scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2011 (NA for missing information)
    - *phenP2_11*: scoring of the phenology of the second shoot of the individual on a scale from 0 to 4 for 2011 (NA for missing information)
    - *phenP3_11*: scoring of the phenology of the third shoot of the individual on a scale from 0 to 4 for 2011 (NA for missing information)
    - *oid1_11*: first powdery mildew infection note for 2011 (NA for missing information)
    - *oid2_11*: second powdery mildew infection note for 2011 (NA for missing information)
    - *oid3_11*: third powdery mildew infection note for 2011 (NA for missing information)
    - *oid4_11*: fourth powdery mildew infection note for 2011 (NA for missing information)
    - *oid4P2_11*: fourth powdery mildew infection note for 2011 on the second shoot only (NA for missing information)
    - *oid5_11*: fifth powdery mildew infection note for 2011 (NA for missing information)
    - *herbi_11*: herbivore damage note for 2011 (NA for missing information)
    - *phen1_12*: first scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2012 (NA for missing information)
    - *phen2_12*: second scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2012 (NA for missing information)
    - *phen3_12*: third scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2012 (NA for missing information)
    - *phen4_12*: fourth scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2012 (NA for missing information)
    - *phen5_12*: fifth scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2012 (NA for missing information)
    - *oid1_12*: first powdery mildew infection note for 2012 (NA for missing information)
    - *oid2_12*: second powdery mildew infection note for 2012 (NA for missing information)
    - *oid3_12*: third powdery mildew infection note for 2012 (NA for missing information)
    - *oid4_12*: fourth powdery mildew infection note for 2012 (NA for missing information)
    - *nbp_12*: number of shoots in 2012, ranging from 1 to 5 (NA for missing information)
    - *phen1_13*: first scoring of the phenology of the first shoot of the individual on a scale from 0 to 4 for 2013 (NA for missing information)
    - *gel_13*: individual affected by frost in 2013 ("0" = no, "1"= yes, "NA" = missing information)
    - *oid1_13*: first powdery mildew infection note for 2013 (NA for missing information)
    - *oid2_13*: second powdery mildew infection note for 2013 (NA for missing information)
    - *oid3_13*: third powdery mildew infection note for 2013 (NA for missing information)
    - *oid1_14*: first powdery mildew infection note for 2014 (NA for missing information)
    - *oid2_14*: second powdery mildew infection note for 2014 (NA for missing information)
    - *oid_16*: powdery mildew infection note for 2016 (NA for missing information)
    - *oid_17*: powdery mildew infection note for 2017 (NA for missing information)
    - *drap10*: flagshoot symptoms detected in 2010 ("0" = no, "1"= yes, "NA" = missing information)
    - *drap11*: flagshoot symptoms detected in 2011 ("0" = no, "1"= yes, "NA" = missing information)
    - *drap12*: flagshoot symptoms detected in 2012 ("0" = no, "1"= yes, "NA" = missing information)
    - *drap13*: flagshoot symptoms detected in 2013 ("0" = no, "1"= yes, "NA" = missing information)
    - *drap15*: flagshoot symptoms detected in 2015 ("0" = no, "1"= yes, "NA" = missing information)
    - *freq_drap*: frequency of flagshoot symptoms detection between 2010 and 2015 (from 0 to 3, "NA" = missing information)
    - *pgland*: weight of the acorn (grams)
    - *date_em*: date of acorn raising
    - *Hfin09*: height measured at the end of 2009 (cm)
    - *Hfin10*: height measured at the end of 2010 (cm)
    - *Hfin11*: height measured at the end of 2011 (cm)
    - *Hfin12*: height measured at the end of 2012 (cm)
    - *Hdeb14*: height measured at the beginning of 2014 (cm)
    - *Hdeb15*: height measured at the beginning of 2015 (cm)
    - *Hdeb16*: height measured at the beginning of 2016 (cm)
    - *Hdeb17*: height measured at the beginning of 2017 (cm)
    - *diam16*: diameter measured in 2016 (cm)
    - *an_mort*: year of death of the seedling oak
    - *H09v*: size at the end of the year 2009 in cm (measured from the highest living shoot of the individual)
    - *H10v*: size at the end of the year 2010 in cm (measured from the highest living shoot of the individual)
    - *H11v*: size at the end of the year 2011 in cm (measured from the highest living shoot of the individual)
    - *H12v*: size at the end of the year 2012 in cm (measured from the highest living shoot of the individual)
    - *H14v*: size at the end of the year 2014 in cm (measured from the highest living shoot of the individual)
    - *H15v*: size at the end of the year 2015 in cm (measured from the highest living shoot of the individual)
    - *H16v*: size at the end of the year 2016 in cm (measured from the highest living shoot of the individual)
    - *H17v*: size at the end of the year 2017 in cm (measured from the highest living shoot of the individual)
    - *cr10*: growth in cm for the year 2010 (can be negative)
    - *cr11*: growth in cm for the year 2011 (can be negative)
    - *cr12*: growth in cm for the year 2012 (can be negative)
    - *cr13*: growth in cm for the year 2013 (can be negative)
    - *cr14*: growth in cm for the year 2014 (can be negative)
    - *cr15*: growth in cm for the year 2015 (can be negative)
    - *cr16*: growth in cm for the year 2016 (can be negative)
    - *esp*:
    - *statut*: dead or alive status of the plant in 2017 coded as "viva" for plant still alive and "mort" for dead plants (NA for missing data)
    - *treat*: final exposure category "exp" for exposed to powdery mildew and "low" for plot with powdery mildew infection limited by the use of fungicide
    - *oid_moy*: mean powdery mildew infection score between 2009 and 2013
    - *statut10*: dead or alive status of the plant in 2017 coded as a binary variable (0=dead; 1=alive)
    
- **pro.hmp.txt** and **nat.hmp.txt:** two data sets formatted for the GWAS analyses using the *hapmap* format, for the protected and natural treatment, respectively
    - *rs#*: SNP ID
    - *alleles*: possible nucleotide allele for each SNP
    - *chrom*: ID number of the chromosome the SNP is located on. Range from 1 to 12 for the 12 chromosome of *Quercus robur*. The *13* code is used for unmapped SNP. *MT* is for SNP located on mitochondrial DNA and *Plst* for SNP located on plastide
    - *pos*: the position of the SNP on the chromosome
    - *strand*: not relevant for this analysis (filled with "+")
    - *assembly#*: from which genomic assembly the chromosome and position on the chromosome of the SNP has been established
    - *center*: the laboratory from which the sample were processed
    - *protLSID*: not relevant for this analysis (filled with "NA")
    - *assayLSID*: not relevant for this analysis (filled with "NA")
    - *panelLSID*: ID of the panel from which the SNP has been selected from. For these dataset, all the SNP were designed for this experiment
    - *QCcode*: not relevant for this analysis (filled with "NA")
    - *A2B18* to *C1K14*: each of this column names stands for a genotyped individual. In each column, you will find the genotypes for the successfully genotyped SNP. The value is either the biallelic code for the genotype or "NA" when the genotyping failed


## R scripts
In this section, you will find the list of the different scripts used in the article with a brief description of their purpose.

- **RealTime_Fstcomparison.R:** script to compute F-statistics and compare them across the different treatment and time. The output files are stored in a Genepop folder created within the output folder.
- **RealTime_GENHET.R:** script to compute intra-individual heterozygosity indices. 
- **RealTime_GWAS.R:** script to perform the Genome-Wide Association analyses for both treatment. This script will produce two folders each containing results from the GWAS analyses. These output files are necessary for plotting the manhattan and GWAS plot (see scripts below).
- **RealTime_load.R:** the script to load the different data sets, functions and packages that are necessary for the data analyses and representation in the R environment.
- **RealTime_plot_evolpheno.R:** script to analyze and plot the temporal evolution of the main phenotypic traits over the years (Figure 2). The code for the distribution of tree height by family and by powdery mildew exposures (Figure 6) is also included in this script. 
- **RealTime_plot_FigS2_map.R:** script to plot the map (Figure S2) of the experimental setup.
- **RealTime_plot_FigS4_indQual.R:** script to plot the scatterplot used to determined the individuals genotyped with the SNP markers with sufficient quality.
- **RealTime_plot_FigsIndHeter.R:** script to produce the Figures related to intra-individual heterozygosity (Figure 8, Figure S7, Figure S14 and Figure S16). The ***RealTime_GENHET.R*** must be run prior to this script in order to produce the necessary entry files. 
- **RealTime_plot_glm.R:** script to plot the Figures (Figure 3, Figure S10, Figure S11 and Figure S12) related to the GLM analyses performed in SAS.
- **RealTime_plot_GoldGate_SNP.R:** script to plot the scatterplot of the genotypes examples of SNP markers of variable quality. Figure S5 and Figure S6 are produce using this script.
- **RealTime_plot_GWAS.R:** script to plot the distribution of phenotype traits by genotypes for SNP of interest (Figure 10).
- **RealTime_plot_manhattan.R:** script to plot the manhattan plot Figure 9 for the four phenotypic traits of interest and for the two exposures (Natural and Protected).
- **RealTime_plot_SurvProtecVsNat.R:** script to plot the progeny survival percentages across families in Protected experimental plots versus Natural powdery mildew infected plots (Figure 5).
- **RealTime_plot_SurvShannon.R:** script to compute the Shannon index for each yera and to plot the Shannon index evolution over years (Figure 7).


## R session info
For reproducibility purpose, you will find all the information about the versions of R, Rstudio, OS etc., as well as the list and version number of the packages used at the time of publishing this script in the **session_info.txt** file.


## SAS Data sets (/dataSAS)
In this section, you will find the list of the data sets used in this study for SAS analyses. The data files can be found in the "dataSAS" folder. For the data tables, the name of the different variables are listed and explained as well. There are 2 data sets used in this study for analyses in SAS.
- **RT_pop.csv:** the main data set which includes the information for all individuals:
    - *Ind:* unique identification number for each sample
    - *fameff:* which family the individual belongs to ("1", "9", "11", "27", "45", "48", "51", "70", "71", "72", "73", "74", "75", "76", "77") 
    - *bloc:*	block ID in the experimental design
    - *PU:* unitary plot ID
    - *rg*: rank number within the unitary plot. There are 12 rank per block
    - *n*: individual number on the rank. There are a maximum of 48 position on a rank
    - *fam:* original family ID before correction using the genotyping information
    - *trait:* initial coding of the exposure of the plot ("nat" and "renf" for natural exposure, "lim" for protected exposure)
    - *oid4_09:* fourth powdery mildew infection note for 2009 (empty for missing information)
    - *oid5_10:* fifth powdery mildew infection note for 2010 (empty for missing information)
    - *oid5_11:* fifth powdery mildew infection note for 2011 (empty for missing information)
    - *oid4_12:* fourth powdery mildew infection note for 2012 (empty for missing information)
    - *gel_13:* individual affected by frost in 2013 ("0" = no, "1"= yes, empty = missing information)
    - *oid2_13:* second powdery mildew infection note for 2013 (empty for missing information)
    - *oid_16:* powdery mildew infection note for 2016 (empty for missing information)
    - *oid_17:* powdery mildew infection note for 2017 (empty for missing information)
    - *pgland:* weight of the acorn (grams)
    - *date_em:* date of acorn raising
    - *Hfin09:* height measured at the end of 2009 (cm)
    - *Hfin10:* height measured at the end of 2010 (cm)
    - *Hfin11:* height measured at the end of 2011 (cm)
    - *Hfin12:* height measured at the end of 2012 (cm)
    - *Hdeb14:* height measured at the beginning of 2014 (cm)
    - *Hdeb15:* height measured at the beginning of 2015 (cm)
    - *Hdeb16:* height measured at the beginning of 2016 (cm)
    - *Hdeb17:* height measured at the beginning of 2017 (cm)
    - *diam16:* diameter measured in 2016 (cm)
    - *an_mort:* year of the death of the seedlings (from "2009"" to "2017"; "vivant" if the seedling was still alive in 2017)
    - *H09v:* size at the end of the year 2009 in cm (measured from the highest living shoot of the individual)
    - *H10v:* size at the end of the year 2010 in cm (measured from the highest living shoot of the individual)
    - *H11v:* size at the end of the year 2011 in cm (measured from the highest living shoot of the individual)
    - *H12v:* size at the end of the year 2012 in cm (measured from the highest living shoot of the individual)
    - *H14v:* size at the end of the year 2014 in cm (measured from the highest living shoot of the individual)
    - *H15v:* size at the end of the year 2015 in cm (measured from the highest living shoot of the individual)
    - *H16v:* size at the end of the year 2016 in cm (measured from the highest living shoot of the individual)
    - *H17v:* size at the end of the year 2017 in cm (measured from the highest living shoot of the individual)
    - *exp:* final exposure category "exp" for natural powdery mildew exposure and "low" for plot with protected exposure by the use of fungicide
    - *oid_moy:* ean powdery mildew infection score between 2009 and 2013
    - *survie2017:* dead or alive status of the plant in 2017 coded as a binary variable (0=dead; 1=alive). Empty = missing information 

    
- **RT_pop_tot.csv:** the data set which includes only individuals with genotypic information. In addition to the variables already described for the *RT_pop.csv* data set, there are four variables that estimate the level of heterozygosity/homozygosity in each individual:
    - *PHt:* number of heterozygous loci / number of genotyped loci
    - *Hs_obs:* standardized heterozygosity based on the mean observed heterozygosity
    - *Hs_exp:* standardized heterozygosity based on the mean expected heterozygosity
    - *IR:* Internal relatedness
    - *HL:* Homozygosity by locus
    

## SAS scripts
In this section, you will find a brief description of the script used for analyses with the SAS software.

- **RealTime_pgm_analyses_tot.sas:** this file lists all the scripts that have been used to perform statistical analyses and model fitting in the SAS framework. Each script begins with a line describing the general purpose of the analysis. All scripts have to be run on the *RT_pop.csv* data set, except for the two in the "with heterozygosity" section which should be run on the *RT_pop_tot.csv* data set. 


## Supplementary Data sets (/data/dataSup)
In this section, you will find the list of the supplementary data sets used in this study. The data file can be found in the "dataSup" subfolder. There are 4 supplementary data sets linked to this study.

- **Data_S1.txt:** this file lists the flanking sequences of the 819 SNP used in the study. These sequences were used to map the SNP markers on the *Q. robur* genome. This is a fasta file format. 
- **Ind_Geno_Qual.txt:** dataset providing quality statistics on genotyping results for genotyped individuals. These data are used to produce Figure S4
- **SNP_Ind_coord_FDT.txt** and **SNP_Stat_Tab_T.txt:** two datasets that are the output of the clustering permformed using the genotyping module of the BeadStudio/GenomeStudio package (Illumina, San Diego, CA, USA). This datasets are used to produce Figures S5 and S6 and can be used to produce the polar plot for all SNPs markers. 


## Citation
You will be soon (hopefully) able to cite the related study as follow: 
+ Barrès B., Dutech C., Saint-Jean G., Bodénès C., Burban C., Fievet V., Lepoittevin C., Garnier-Géré P. and Desprez-Loustau M.-L. [Demographic and genetic impacts of powdery mildew in a young oak cohort. ***Submitted to PCI Forest and Wood Sciences***.](https://)

If you want to use (some of) the code found on this page or if you want to cite this repository: 
+ Benoit Barrès and Marie-Laure Desprez-Loustau. [Supporting data and code for: Demographic and genetic impacts of powdery mildew in a young oak cohort. Zenodo.](https://zenodo.org/badge/latestdoi/33980368)
