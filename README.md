[![DOI](https://zenodo.org/badge/41293576.svg)](https://zenodo.org/badge/latestdoi/41293576)
# Significant demographic and genetic impacts of powdery mildew in a young oak cohort
*This repository contains the R code and SAS scripts used for the data analyses and production of the figures of the related article*

![alt text](https://am3pap005files.storage.live.com/y4mEo2XF6NRE3Q7vYX-H1vCjsYjXvfngA-U2S8sHSvrAqdHnIyhOpRKxwWaxZOtRUsiAsc2fql-lFhE2QGTQ7E6ecA5g7yITvUmXp1Bi-prt1TCqOI4wBOXzJzvhSRhrUmlL04cHbQKf7MRwKonEBrn6j5-BaTGjuehNdwQzywuS3w5logfETAYnODOXoUe99fx?width=1584&height=588&cropmode=none)


## Context
 


## Datasets
In this section, you will find the list of the data sets used in this study. The data files can be found in the "data" folder. For the data tables, the name of the different variables are listed and explained as well. There are 5 data sets used in this study.  

+ **AgrAph5.dat:** the first data set contains the data for all the indivudals analyzed. Each line correspond to one individuals and the following information for each individuals can be found in this table: 
  + *indiv_ID*: individual's ID, this is a unique string of character
  + *data_batch*: the global dataset can be divided in several sub-projects (AgrAphid, RAtransect, Rpp and Zepeda)



## R scripts
In this section, you will find the list of the different scripts used in the article with a brief description of their purpose.

+ **RealTime_load.R:** the script to load the different data sets, functions and packages that are necessary for the data analyses and representation in the R environment. 
+ **RealTime_fun.R:** functions to perform the delta-K analysis. 
+ **RealTime_strplot_fun.R:** a function to plot beautiful STRUCTURE-like plot with several parameters to control the output. 


## R session info
For reproducibility purpose, you will find all the information about the versions of R, Rstudio, OS etc., as well as the list and version number of the packages used at the time of publishing this script in the **session_info.txt** file.

## Citation
You will be soon (hopefully) able to cite the related study as follow: 
+ Barrès B., Saint-Jean G.,Garnier-Géré P., Dutech C. and Desprez-Loustau M.-L.
[Significant demographic and genetic impacts of powdery mildew in a young oak cohort. *journal name*.](https://)

If you want to use (some of) the code found on this page or if you want to cite this repository:
+ Benoit Barrès and Marie-Laure Desprez-Loustau. [Supporting data and code for: Significant demographic and genetic impacts of powdery mildew in a young oak cohort. Zenodo.](https://zenodo.org/badge/latestdoi/sss)
