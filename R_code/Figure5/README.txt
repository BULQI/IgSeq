#README for reproducing the figure 5 in the manuscript

##Data 
To prepare the data for reproducing the figure 5, do one of the following:
	
	- download the data and put the data in the ../Data/Figure5 (see the links contained in the file ../../Data/Figure5/ReadMe.txt)
	
	- run the data precessing scripts (../DataPrepare), which will generate the data and save at the correct location.

##Steps to reproduce the figure
	
	- Run the R script "Diversity_sampleDepth_subsample.R" to do subsampling of sequencing data and compare the sample size effect on diversity estimations

	- Run the R script "sizeDistribution.R" to prepare the data for plotting the clone size distribution.

	- Run the R script "Diversity_Subsample_RCMD_v2.0.r" to generate diversity estimation on subsampled sequence data.

	- Run the script "Disersity_figure_subSample_v2.1_arcsin.r" to run anova and plot the comparsion among different groups (Figure 5).

	
