# README for the r code in this folder.

The R code files in this folder are used to plot Figure 3 in the manuscript. 
Please run the code file in the order as shown below to generate the figures.

- 1) DataAnalysis_v1.0_MutationFreq_figure.R
this module reads the data in the ../../Data/Figure3 folder, processes them, generated intermediate files/output for the next R code file to plot the final figures. Run these r code in this folder first.
It assume you have Ig recombinationation files and sample infomation files in the ../../Data file. Make you have downloaded them or run the data preprocessing modules already. (see details in ../../Data/Figure3/ReadMe.txt)

- 2) figurePlotting_figure_v5.0.r
this module reads the input files generated in the previous module and simply plot the figures. Make sure you have run the first code file before this.
