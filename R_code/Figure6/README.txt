#README for reproducing the figure 6 in the manuscript

##Data

To prepare the data for reproducing the figure 6, do one of the following:
	
	- download the data from . and put the data in the ../Data/Figure5
	
	- run the data precessing scripts (../DataPrepare), which will generate the data and save at the correct location.

##Steps to reproduce the figure
	
	- Run the R script "IntraClonal_batch.R" to estimate two metrics: intra clonal diversity (pairwise clonal difference) and also mean clonal mutation frequency.

	- Run the R script "intraClonaDiversity_process.R" to process the data and estimate the clonal selection strength and save the data for plotting.

	- Run the R script "intraClonalDiversity_process_figure.r" to to run anova and plot the comparsion among different groups (Figure 6).


