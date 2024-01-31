Cloanalyst package for doing both VDJ annotation and clonotyping.

## Reference 
Cloanalyst reference:
 Kepler TB. Reconstructing a B-cell clonal lineage. I. Statistical inference of unobserved ancestors [version 1; peer review: 2 approved, 1 approved with reservations]. F1000Research 2013, 2:103 (https://doi.org/10.12688/f1000research.2-103.v1) 

## Code github repository (windows version): 
https://github.com/BULQI/Cloanalyst


## Code package (linux version).

We also in this project have revised the software into a linux version and use the linux version software to run
analysis for the manuscript.

The linux version of Cloanalyst can be accessed at : https://drive.google.com/file/d/1RO2fkR2w3TdWW2fVEprmbk8ckMgsVR58/view?usp=sharing


### Usage of Cloanalyst (linux version)

- 1) download the archive and unzip the folder.

- 2) dotnet sdk 3.1 is needed to run the program.

- 3) Test the program with the example data set contained in folder by running
    
    dotnet ./Cloanalyst.dll parse -s "Mus musculus" -c Heavy  --thread 1 --excl 5 -g "AR20170307 FPC-F" -r2 S1_R2_IgG3200.fasta S1_R1_IgG3200.fasta

## docker container 

The docker container for running Cloanalyst under the linux environment is available at ffeng23/cloanalyst (docker hub).
