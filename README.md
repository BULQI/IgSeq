# IgSeq v1.0.0
Code for IgSeq data analysis

Reference: Characterizing adjuvants’ effects at the murine immunoglobulin repertoire level ([Preprint](https://doi.org/10.1101/2022.11.19.517218 ))

([Final publication](https://doi.org/10.1016/j.isci.2023.108749))

## ABSTRACT
High-throughput immunoglobulin sequencing (IgSeq) has been developed and applied to study the adaptive immune response extensively for more than a decade. However, generating large-scale, high-fidelity sequencing data is still challenging, and furthermore, not much has been done to characterize adjuvants' effects at the repertoire level. Thus, we developed an improved library prep protocol and standardized the data analysis pipeline for accurate repertoire profiling. In addition, two metrics were implemented to assess repertoire clone properties. We then studied systemically the effects of two adjuvants, CpG and Alum, on the Ig heavy chain repertoire using the ovalbumin (OVA) challenged mouse model. Ig repertoires of different tissues (spleen and bone marrow) and isotypes (IgG and IgM) were examined and compared in terms of sequence mutation frequency, IGHV gene usage, CDR3 length, rescaled Hill numbers for clonal diversity, and clone selection strength. As a result, Ig repertoires of different tissues or isotypes exhibited distinguishable profiles at the non-immunized steady state. Adjuvanted immunizations further resulted in statistically significant alterations in Ig repertoire compared with PBS or OVA alone immunized groups. Lastly, we applied unsupervised machine learning techniques — multiple factor analysis and clustering — to identify Ig repertoire signatures in different compartments and under varying immunizations. We found that the IGH repertoires of distinct tissue-isotype compartments or under varying immunizations differed in unique sets of properties. Notably, Alum and CpG effects on the Ig repertoire exhibited different tissue and isotype preferences. The former led to increased diversity of abundant clones of both isotypes in BM only, and the latter promoted the selection of IgG clones only but in both tissues. The patterns of Ig repertoire changes likely reflected possible action mechanisms of these two adjuvants.

Contact: Feng Feng (ffeng@bu.edu)

## Explanation of the project

There are two folders: DataPreprocessing_code and R_code. The former contains code to preprocessing the raw read data and the latter contains R code for generating figures and results in the manuscripts.

### DataPreprocessing_code


## Raw IgSeq data

The raw IgSeq data were submitted to the NIH GEO repository. It can be found at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228856 . The record number is GSE228856.

## Docker containers

- IgSeq R3 docker image can be accessed at https://hub.docker.com/repository/docker/ffeng23/igseqr3/general

This is an Rstudio container for running R code files for plotting the figures in the manuscripts.
You need to provide the data (see the ReadMe files in each figure data folder for details).

- Umi merge docker image can be accessed at https://hub.docker.com/repository/docker/ffeng23/umi_merge/general

This is python docker container for running umi merge processing. 

