[![DOI](https://zenodo.org/badge/41293576.svg)](https://zenodo.org/badge/latestdoi/41293576) \# Significant demographic and genetic impacts of powdery mildew in a young oak cohort *This repository contains the R code and SAS scripts used for the data analyses and production of the figures of the related article*

![alt text](https://am3pap005files.storage.live.com/y4mXUw4rIo7I6I-UHrfoPc32YcaZYoOIM_R2-WK8ZLDnMrfurYfJ5FV2WjlTIh_idbaCJMyDFBGPAnA3lwmFOaY_6M5ra2vfzfOiobz1ENwdeA1QbGCTFYXnkZznKUDZXRVNRsXyB7KCZzHJLIFI1B8rqivCd0_12NtVbUs-X7a5FWMEHdFaMUfnwWvsiDjm8JU?width=1584&height=588&cropmode=none)

## Context

## Datasets

In this section, you will find the list of the data sets used in this study. The data files can be found in the "data" folder. For the data tables, the name of the different variables are listed and explained as well. There are 5 data sets used in this study.

-   **datatot.txt:** the main data set which contains all the . Each line correspond to one individuals and the following information for each individuals can be found in this table:
    -   *Sample_ID*:
    -   *family*:
    -   *family_simp*:
    -   *parent_id*:
    -   *bloc*:
    -   *PU*:
    -   *rg*:
    -   *n*:
    -   *coord_X*:
    -   *coord_Y*:
    -   *na.9micro*:
    -   *na.12micro*:
    -   *barcode*:
    -   *A11_all1* to *S19_all2* (24 columns):
    -   *nb_conta_12*:
    -   *nb_conta_09*:
    -   *SNPage*:
    -   *Quality_SNPage*:
    -   *pb_robot_SNPage*:
    -   *mother.snp*:
    -   *father.snp*:
    -   *live_bin*:
    -   *live_year*:
    -   *na.snp*:
    -   *REF\~CL371CT472_02-163* to *REF\~CL9715CT16278_02-801* (819 columns):
    -   
-   **lim.hmp.txt.txt and nat.hmp.txt:** two data sets formatted for the GWAS analyses, for the protected and natural treatment respectively
    -   *indiv_ID*:

## R scripts

In this section, you will find the list of the different scripts used in the article with a brief description of their purpose.

-   **RealTime_load.R:** the script to load the different data sets, functions and packages that are necessary for the data analyses and representation in the R environment.
-   **RealTime_evolpheno.R:** script to analyze and plot the evolutionary patterns of the main phenotypic traits over the years.
-   **RealTime_plot_SurvShannon.R:** script to compute the Shannon index evolution and to plot the Figure of the survival and Shannon index evolution.
-   **RealTime_Fstcomparison.R:** script to compute F-statistics and compare them across the different treatment and time. The output files are stored in a Genepop folder created within the output folder.
-   **RealTime_GENHET.R:** script to compute intra-individual heterozygosity indices. The code to produce the related Figure is also included.
-   **RealTime_GWAS.R:** script to perform the Genome-Wide Association analyses for both treatment. This script will produce two folders each containing results from the GWAS analyses. These output files are necessary for plotting the manhattan and GWAS plot.
-   **RealTime_plot_manhattan.R:** script to plot the manhattan plot Figure.
-   **RealTime_plot_GWAS.R:** script to plot the GWAS results for the SNP of interest.

## R session info

For reproducibility purpose, you will find all the information about the versions of R, Rstudio, OS etc., as well as the list and version number of the packages used at the time of publishing this script in the **session_info.txt** file.

## Citation

You will be soon (hopefully) able to cite the related study as follow: + Barrès B., Saint-Jean G., Lepoittevin C., Burban C., Garnier-Géré P., Dutech C. and Desprez-Loustau M.-L. [Significant demographic and genetic impacts of powdery mildew in a young oak cohort. *journal name*.](https://)

If you want to use (some of) the code found on this page or if you want to cite this repository: + Benoit Barrès and Marie-Laure Desprez-Loustau. [Supporting data and code for: Significant demographic and genetic impacts of powdery mildew in a young oak cohort. Zenodo.](https://zenodo.org/badge/latestdoi/sss)
